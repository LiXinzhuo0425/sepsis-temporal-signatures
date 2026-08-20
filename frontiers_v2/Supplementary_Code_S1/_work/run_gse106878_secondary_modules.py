#!/usr/bin/env python3
"""Secondary clinical and repeated-measurement analyses for GSE106878.

This cohort remains a single independent replication cohort.  The analyses
below do not enter the six-cohort primary meta-analysis.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
import statsmodels.formula.api as smf


ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "_work" / "analysis" / "gse106878"
META = ROOT / "_sources" / "cohort_audit" / "GSE106878_verified_sample_mapping.csv"
OUT = ROOT / "_work" / "analysis"
B = 2000

PRIMARY = ["SIG001", "SIG002", "SIG003", "SIG004"]
ALL = PRIMARY + ["SIG022", "SIG023", "SIG033", "SIG034"]
NAMES = {
    "SIG001": "Sepsis MetaScore", "SIG002": "SeptiCyte LAB", "SIG003": "FAIM3:PLAC8", "SIG004": "sNIP",
    "SIG022": "Bacterial/Viral MetaScore", "SIG023": "Herberg two-gene DRS",
    "SIG033": "Lin seven-gene mortality score", "SIG034": "Severe-or-Mild score",
}


def seed_for(*parts: str) -> int:
    digest = hashlib.sha256("|".join(parts).encode()).digest()
    return int.from_bytes(digest[:8], "big") % (2**32 - 1)


def holm(series: pd.Series) -> pd.Series:
    output = pd.Series(np.nan, index=series.index, dtype=float)
    ordered = series.dropna().sort_values()
    running = 0.0
    for rank, (index, value) in enumerate(ordered.items()):
        running = max(running, min(1.0, (len(ordered) - rank) * float(value)))
        output.loc[index] = running
    return output


def bootstrap_adjusted_status(frame: pd.DataFrame, seed: int) -> tuple[float, float, float]:
    rng = np.random.default_rng(seed)
    favorable = frame[frame["adverse"] == 0]
    adverse = frame[frame["adverse"] == 1]
    draws = np.empty(B)
    for b in range(B):
        sample = pd.concat([
            favorable.iloc[rng.integers(0, len(favorable), len(favorable))],
            adverse.iloc[rng.integers(0, len(adverse), len(adverse))],
        ], ignore_index=True)
        design = np.column_stack([
            np.ones(len(sample)), sample["adverse"].to_numpy(float), sample["hydrocortisone"].to_numpy(float)
        ])
        coef, *_ = np.linalg.lstsq(design, sample["delta_z"].to_numpy(float), rcond=None)
        draws[b] = coef[1]
    return float(draws.std(ddof=1)), float(np.quantile(draws, 0.025)), float(np.quantile(draws, 0.975))


def rank_metrics(frame: pd.DataFrame, seed: int) -> dict:
    baseline = frame["baseline_score"].to_numpy(float)
    followup = frame["followup_score"].to_numpy(float)
    n = len(frame)

    def metrics(x: np.ndarray, y: np.ndarray) -> tuple[float, float]:
        rho = float(stats.spearmanr(x, y).statistic)
        rx = (stats.rankdata(x, method="average") - 1) / (len(x) - 1)
        ry = (stats.rankdata(y, method="average") - 1) / (len(y) - 1)
        return rho, float(np.median(np.abs(ry - rx)))

    rho, displacement = metrics(baseline, followup)
    rng = np.random.default_rng(seed)
    rho_draws, displacement_draws = [], []
    for _ in range(B):
        idx = rng.integers(0, n, n)
        r, d = metrics(baseline[idx], followup[idx])
        if np.isfinite(r) and np.isfinite(d):
            rho_draws.append(r)
            displacement_draws.append(d)
    if len(rho_draws) < 0.95 * B:
        raise RuntimeError("Fewer than 95% finite rank bootstrap replicates")
    return {
        "n_pairs": n,
        "spearman_rho": rho,
        "rho_ci95_lower": float(np.quantile(rho_draws, 0.025)),
        "rho_ci95_upper": float(np.quantile(rho_draws, 0.975)),
        "median_rank_displacement": displacement,
        "displacement_ci95_lower": float(np.quantile(displacement_draws, 0.025)),
        "displacement_ci95_upper": float(np.quantile(displacement_draws, 0.975)),
        "finite_bootstrap_replicates": len(rho_draws),
    }


def main() -> None:
    pairs = pd.read_csv(SRC / "GSE106878_patient_paired_changes.csv")
    metadata = pd.read_csv(META, dtype=str)
    baseline = metadata[metadata["timepoint"] == "Pre"][["patient_id", "survival_(28_days)", "treatment"]].copy()
    baseline["adverse"] = (baseline["survival_(28_days)"] == "died").astype(int)
    baseline["hydrocortisone"] = (baseline["treatment"] == "hydrocortisone").astype(int)
    if baseline["patient_id"].duplicated().any() or len(baseline) != 47:
        raise RuntimeError("GSE106878 baseline outcome mapping is not one-to-one")
    pairs = pairs.merge(baseline, on="patient_id", how="left", validate="many_to_one", suffixes=("", "_meta"))
    if pairs[["adverse", "hydrocortisone"]].isna().any().any():
        raise RuntimeError("Missing GSE106878 clinical mapping")

    clinical_rows, rank_rows = [], []
    for signature in ALL:
        frame = pairs[pairs["signature_id"] == signature].copy()
        fit = smf.ols("delta_z ~ adverse + hydrocortisone", frame).fit(cov_type="HC3")
        ancova = frame.copy()
        ancova["baseline_z"] = ancova["baseline_score"] / ancova["baseline_sd"]
        ancova["followup_z"] = ancova["followup_score"] / ancova["baseline_sd"]
        afit = smf.ols("followup_z ~ baseline_z + adverse + hydrocortisone", ancova).fit(cov_type="HC3")
        seed = seed_for("GSE106878", signature, "mortality")
        boot_se, boot_lo, boot_hi = bootstrap_adjusted_status(frame, seed)
        clinical_rows.append({
            "dataset": "GSE106878", "signature_id": signature, "signature_name": NAMES[signature],
            "endpoint": "28-day mortality", "analysis_role": "secondary treatment-adjusted clinical context",
            "survivor_n": int((frame["adverse"] == 0).sum()), "nonsurvivor_n": int((frame["adverse"] == 1).sum()),
            "placebo_n": int((frame["hydrocortisone"] == 0).sum()), "hydrocortisone_n": int((frame["hydrocortisone"] == 1).sum()),
            "beta_adverse_minus_favorable": float(fit.params["adverse"]),
            "se_hc3": float(fit.bse["adverse"]), "ci95_lower": float(fit.conf_int().loc["adverse", 0]),
            "ci95_upper": float(fit.conf_int().loc["adverse", 1]), "p_value": float(fit.pvalues["adverse"]),
            "bootstrap_se": boot_se, "bootstrap_ci95_lower": boot_lo, "bootstrap_ci95_upper": boot_hi,
            "ancova_beta": float(afit.params["adverse"]), "ancova_p_hc3": float(afit.pvalues["adverse"]),
            "bootstrap_replicates": B, "seed": seed,
        })
        rank = rank_metrics(frame, seed_for("GSE106878", signature, "rank"))
        rank_rows.append({"dataset": "GSE106878", "signature_id": signature, "signature_name": NAMES[signature], **rank})

    clinical = pd.DataFrame(clinical_rows)
    clinical["holm_adjusted_p_primary_four"] = np.nan
    ix = clinical[clinical["signature_id"].isin(PRIMARY)].index
    clinical.loc[ix, "holm_adjusted_p_primary_four"] = holm(clinical.loc[ix, "p_value"])
    ranks = pd.DataFrame(rank_rows)

    # Pairwise diagnostic change concordance in the same 47 patients.
    wide = pairs[pairs["signature_id"].isin(PRIMARY)].pivot(index="patient_id", columns="signature_id", values="delta_z")
    concordance = []
    for i, a in enumerate(PRIMARY):
        for b in PRIMARY[i + 1:]:
            rho = float(stats.spearmanr(wide[a], wide[b]).statistic)
            concordance.append({"dataset": "GSE106878", "signature_a": a, "signature_b": b, "n_pairs": len(wide), "spearman_rho": rho})
    concordance = pd.DataFrame(concordance)

    clinical.to_csv(OUT / "07_gse106878_treatment_adjusted_mortality.csv", index=False)
    ranks.to_csv(OUT / "07_gse106878_rank_persistence.csv", index=False)
    concordance.to_csv(OUT / "07_gse106878_diagnostic_concordance.csv", index=False)
    summary = {
        "patients": int(pairs["patient_id"].nunique()), "signatures": int(pairs["signature_id"].nunique()),
        "survivors": int((baseline["adverse"] == 0).sum()), "nonsurvivors": int((baseline["adverse"] == 1).sum()),
        "placebo": int((baseline["hydrocortisone"] == 0).sum()), "hydrocortisone": int((baseline["hydrocortisone"] == 1).sum()),
        "status": "PASS",
    }
    (OUT / "07_gse106878_secondary_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()

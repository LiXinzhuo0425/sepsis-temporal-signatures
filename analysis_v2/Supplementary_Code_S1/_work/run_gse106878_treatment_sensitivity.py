#!/usr/bin/env python3
"""Treatment-arm sensitivity for four diagnostic signatures in GSE106878.

This module tests whether directional recurrence depends on one randomized
arm. It is not an efficacy or causal treatment-effect analysis.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests


ROOT = Path(__file__).resolve().parents[1]
INPUT = ROOT / "_work" / "analysis" / "gse106878" / "GSE106878_patient_paired_changes.csv"
FROZEN_INPUT = ROOT / "_sources" / "revision_addition" / "GSE106878_patient_paired_changes.csv"
OUTPUT = ROOT / "_work" / "analysis" / "gse106878_treatment_sensitivity"
SIGNATURES = {
    "SIG001": "Sepsis MetaScore",
    "SIG002": "SeptiCyte LAB",
    "SIG003": "FAIM3:PLAC8",
    "SIG004": "sNIP",
}
ALL_SIGNATURE_ORDER = ["SIG001", "SIG002", "SIG003", "SIG004", "SIG022", "SIG023", "SIG033", "SIG034"]
SEED = 20260720
N_BOOT = 2000


def bootstrap_mean(values: np.ndarray, rng: np.random.Generator) -> tuple[float, float, float]:
    draws = np.empty(N_BOOT, dtype=float)
    for index in range(N_BOOT):
        draws[index] = rng.choice(values, size=len(values), replace=True).mean()
    return float(draws.std(ddof=1)), float(np.quantile(draws, 0.025)), float(np.quantile(draws, 0.975))


def arm_estimates(data: pd.DataFrame) -> pd.DataFrame:
    rng = np.random.default_rng(SEED)
    rows = []
    for signature_id in ALL_SIGNATURE_ORDER:
        frame = data.loc[data["signature_id"].eq(signature_id)].sort_values("baseline_sample_id", kind="stable")
        if len(frame) != 47:
            raise ValueError(f"{signature_id}: expected 47 paired patients")
        for stratum, subset in [("all", frame), *list(frame.groupby("treatment", sort=True))]:
            values = subset["delta_z"].to_numpy(float)
            se, low, high = bootstrap_mean(values, rng)
            rows.append({
                "signature_id": signature_id,
                "stratum": stratum,
                "n_pairs": len(values),
                "mean_delta_z": float(values.mean()),
                "bootstrap_se": se,
                "bootstrap_ci_low": low,
                "bootstrap_ci_high": high,
                "bootstrap_replicates": N_BOOT,
                "seed": SEED,
            })
    return pd.DataFrame(rows)


def main() -> None:
    data = pd.read_csv(INPUT) if INPUT.exists() else pd.read_csv(FROZEN_INPUT)
    required = {"patient_id", "signature_id", "treatment", "baseline_sample_id", "delta_z"}
    if required - set(data.columns):
        raise ValueError("Paired input is missing required columns")
    arm = arm_estimates(data)
    model_rows = []
    coefficient_rows = []
    for signature_id, signature in SIGNATURES.items():
        frame = data.loc[data["signature_id"].eq(signature_id), ["patient_id", "treatment", "delta_z"]].copy()
        frame["hydrocortisone"] = frame["treatment"].eq("hydrocortisone").astype(int)
        fit = smf.ols("delta_z ~ hydrocortisone", data=frame).fit(cov_type="HC3")
        intervals = fit.conf_int(alpha=0.05)
        for term in ["Intercept", "hydrocortisone"]:
            coefficient_rows.append({
                "signature_id": signature_id,
                "signature": signature,
                "term": term,
                "estimate": float(fit.params[term]),
                "robust_se_hc3": float(fit.bse[term]),
                "ci95_low": float(intervals.loc[term, 0]),
                "ci95_high": float(intervals.loc[term, 1]),
                "p_value": float(fit.pvalues[term]),
                "n": int(fit.nobs),
                "r_squared": float(fit.rsquared),
            })
        model_rows.append({
            "signature_id": signature_id,
            "signature": signature,
            "n_total": int(fit.nobs),
            "n_placebo": int(frame["hydrocortisone"].eq(0).sum()),
            "n_hydrocortisone": int(frame["hydrocortisone"].eq(1).sum()),
            "treatment_beta_hydrocortisone_minus_placebo": float(fit.params["hydrocortisone"]),
            "robust_se_hc3": float(fit.bse["hydrocortisone"]),
            "ci95_low": float(intervals.loc["hydrocortisone", 0]),
            "ci95_high": float(intervals.loc["hydrocortisone", 1]),
            "p_value_raw": float(fit.pvalues["hydrocortisone"]),
            "r_squared": float(fit.rsquared),
        })

    summary = pd.DataFrame(model_rows)
    summary["p_value_holm"] = multipletests(summary["p_value_raw"], method="holm")[1]
    selected_arm = arm.loc[arm["signature_id"].isin(SIGNATURES)]
    wide = selected_arm.pivot(index="signature_id", columns="stratum")
    for stratum in ["placebo", "hydrocortisone"]:
        summary[f"{stratum}_mean_delta_z"] = summary["signature_id"].map(wide["mean_delta_z"][stratum])
        summary[f"{stratum}_bootstrap_ci_low"] = summary["signature_id"].map(wide["bootstrap_ci_low"][stratum])
        summary[f"{stratum}_bootstrap_ci_high"] = summary["signature_id"].map(wide["bootstrap_ci_high"][stratum])
    summary["direction_consistent_between_arms"] = (
        np.sign(summary["placebo_mean_delta_z"]) == np.sign(summary["hydrocortisone_mean_delta_z"])
    )
    coefficients = pd.DataFrame(coefficient_rows).merge(
        summary[["signature_id", "p_value_holm"]], on="signature_id", how="left"
    )
    coefficients.loc[coefficients["term"].ne("hydrocortisone"), "p_value_holm"] = np.nan

    OUTPUT.mkdir(parents=True, exist_ok=True)
    summary.to_csv(OUTPUT / "GSE106878_treatment_arm_sensitivity.csv", index=False)
    coefficients.to_csv(OUTPUT / "GSE106878_treatment_arm_model_coefficients.csv", index=False)
    selected_arm.to_csv(OUTPUT / "GSE106878_treatment_arm_bootstrap_effects.csv", index=False)
    (OUTPUT / "GSE106878_treatment_arm_sensitivity.json").write_text(
        json.dumps({
            "analysis": "delta_z ~ hydrocortisone; HC3 covariance; Holm across four treatment coefficients",
            "bootstrap_replicates": N_BOOT,
            "seed": SEED,
            "scope": "directional-replication sensitivity, not efficacy",
        }, indent=2) + "\n",
        encoding="utf-8",
    )
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()

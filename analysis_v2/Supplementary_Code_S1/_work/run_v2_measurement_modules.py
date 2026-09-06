#!/usr/bin/env python3
"""Run locked Longitudinal V2 rank, concordance and contribution modules."""

from __future__ import annotations

import hashlib
import json
import math
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import optimize, stats


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "_sources" / "repo_subset_v1.3.1" / "data"
OUT = ROOT / "_work" / "analysis"
OUT.mkdir(parents=True, exist_ok=True)

WINDOW_LABEL = {"T1": "T24", "T2": "T48", "T3": "T72", "T4": "Day 5"}
SIGNATURE_NAMES = {
    "SIG001": "Sepsis MetaScore",
    "SIG002": "SeptiCyte LAB",
    "SIG003": "FAIM3:PLAC8",
    "SIG004": "sNIP",
    "SIG022": "Bacterial/Viral MetaScore",
    "SIG023": "Herberg two-gene DRS",
    "SIG033": "Lin seven-gene mortality score",
    "SIG034": "Severe-or-Mild score",
}
DIAGNOSTIC = ["SIG001", "SIG002", "SIG003", "SIG004"]
B = 2000


def stable_seed(*parts: str) -> int:
    digest = hashlib.sha256("|".join(map(str, parts)).encode("utf-8")).digest()
    return int.from_bytes(digest[:8], "big") % (2**32 - 1)


def holm(p_values: pd.Series) -> pd.Series:
    values = pd.to_numeric(p_values, errors="coerce")
    out = pd.Series(np.nan, index=values.index, dtype=float)
    valid = values.dropna().sort_values()
    m = len(valid)
    running = 0.0
    for rank, (idx, value) in enumerate(valid.items()):
        adjusted = min(1.0, (m - rank) * float(value))
        running = max(running, adjusted)
        out.loc[idx] = running
    return out


def rowwise_corr(x: np.ndarray, y: np.ndarray) -> np.ndarray:
    x = x - x.mean(axis=1, keepdims=True)
    y = y - y.mean(axis=1, keepdims=True)
    denominator = np.sqrt(np.sum(x * x, axis=1) * np.sum(y * y, axis=1))
    result = np.full(len(x), np.nan, dtype=float)
    valid = denominator > 0
    result[valid] = np.sum(x[valid] * y[valid], axis=1) / denominator[valid]
    return result


def bootstrap_spearman(x: np.ndarray, y: np.ndarray, seed: int) -> tuple[float, float, float, float, float, float, int, str]:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    rho = float(stats.spearmanr(x, y).statistic)
    rng = np.random.default_rng(seed)
    indices = rng.integers(0, len(x), size=(B, len(x)))
    # Re-rank within every resample so ties induced by resampling are handled as Spearman ties.
    rx = stats.rankdata(x[indices], axis=1, method="average")
    ry = stats.rankdata(y[indices], axis=1, method="average")
    boot = rowwise_corr(rx, ry)
    valid_mask = np.isfinite(boot)
    boot = boot[valid_mask]
    if len(boot) < int(0.95 * B):
        return rho, math.nan, math.nan, math.nan, math.nan, math.nan, len(boot), "INSUFFICIENT_BOOTSTRAP"
    lo, hi = np.quantile(boot, [0.025, 0.975])
    z = np.arctanh(np.clip(boot, -0.999999, 0.999999))
    se_z = float(z.std(ddof=1))
    se_source = "bootstrap"
    if not np.isfinite(se_z) or se_z <= 0:
        se_z = float(1 / math.sqrt(max(len(x) - 3, 1)))
        se_source = "analytic_fallback_1_over_sqrt_n_minus_3"
    n = len(x)
    px = (rx - 1) / max(n - 1, 1) * 100
    py = (ry - 1) / max(n - 1, 1) * 100
    boot_displacement = np.median(np.abs(py - px), axis=1)[valid_mask]
    disp_lo, disp_hi = np.quantile(boot_displacement, [0.025, 0.975])
    return rho, float(lo), float(hi), se_z, float(disp_lo), float(disp_hi), len(boot), se_source


def reml_hk(y: np.ndarray, se: np.ndarray) -> dict[str, float]:
    y = np.asarray(y, dtype=float)
    vi = np.maximum(np.asarray(se, dtype=float) ** 2, 1e-12)
    k = len(y)
    if k < 2:
        return {"k": k, "mu": math.nan, "se": math.nan, "lower": math.nan, "upper": math.nan, "tau2": math.nan, "I2": math.nan, "prediction_lower": math.nan, "prediction_upper": math.nan, "p_value": math.nan}

    def objective(tau2: float) -> float:
        weights = 1.0 / (vi + tau2)
        mu = float(np.sum(weights * y) / np.sum(weights))
        return 0.5 * (float(np.sum(np.log(vi + tau2))) + math.log(float(np.sum(weights))) + float(np.sum(weights * (y - mu) ** 2)))

    upper = max(10.0, float(np.var(y, ddof=1) * 20 + np.max(vi) * 20))
    opt = optimize.minimize_scalar(objective, bounds=(0.0, upper), method="bounded", options={"xatol": 1e-12})
    tau2 = max(0.0, float(opt.x))
    if objective(0.0) <= objective(tau2) + 1e-10:
        tau2 = 0.0
    weights = 1.0 / (vi + tau2)
    mu = float(np.sum(weights * y) / np.sum(weights))
    q_hk = float(np.sum(weights * (y - mu) ** 2) / (k - 1))
    hk_se = math.sqrt(max(q_hk, 1e-12) / float(np.sum(weights)))
    crit = float(stats.t.ppf(0.975, k - 1))
    lower, upper_ci = mu - crit * hk_se, mu + crit * hk_se
    p = float(2 * stats.t.sf(abs(mu / hk_se), k - 1)) if hk_se > 0 else math.nan
    fixed_weights = 1.0 / vi
    fixed = float(np.sum(fixed_weights * y) / np.sum(fixed_weights))
    q = float(np.sum(fixed_weights * (y - fixed) ** 2))
    i2 = max(0.0, (q - (k - 1)) / q * 100) if q > 0 else 0.0
    if k >= 3:
        pred_crit = float(stats.t.ppf(0.975, k - 2))
        pred_se = math.sqrt(tau2 + hk_se**2)
        pred_lo, pred_hi = mu - pred_crit * pred_se, mu + pred_crit * pred_se
    else:
        pred_lo, pred_hi = math.nan, math.nan
    modified_se = math.sqrt(max(q_hk, 1.0) / float(np.sum(weights)))
    modified_lower, modified_upper = mu - crit * modified_se, mu + crit * modified_se
    modified_p = float(2 * stats.t.sf(abs(mu / modified_se), k - 1)) if modified_se > 0 else math.nan
    if k >= 3:
        modified_pred_se = math.sqrt(tau2 + modified_se**2)
        modified_pred_lo, modified_pred_hi = mu - pred_crit * modified_pred_se, mu + pred_crit * modified_pred_se
    else:
        modified_pred_lo, modified_pred_hi = math.nan, math.nan
    return {
        "k": k, "mu": mu, "se": hk_se, "lower": lower, "upper": upper_ci, "tau2": tau2, "I2": i2,
        "prediction_lower": pred_lo, "prediction_upper": pred_hi, "p_value": p,
        "modified_se": modified_se, "modified_lower": modified_lower, "modified_upper": modified_upper,
        "modified_prediction_lower": modified_pred_lo, "modified_prediction_upper": modified_pred_hi,
        "modified_p_value": modified_p, "q_hk": q_hk,
    }


def primary_sets() -> dict[tuple[str, str], set[str]]:
    meta = pd.read_csv(DATA / "source_tables" / "stage3" / "03_05_signature_level_meta_analysis.csv")
    meta = meta[(meta["analysis_set"] == "PRIMARY_INDEPENDENT") & meta["time_window"].isin(["T1", "T2"])]
    return {(row.signature_id, row.time_window): set(str(row.cohorts).split(";")) for row in meta.itertuples()}


def rank_module(paired: pd.DataFrame, sets: dict[tuple[str, str], set[str]]) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    rows = []
    displacement_rows = []
    eligible = paired[paired["time_window"].isin(["T1", "T2"]) & paired["analysis_eligible"].astype("string").str.upper().eq("YES")]
    for (dataset, signature, window), group in eligible.groupby(["dataset", "signature_id", "time_window"], sort=True):
        if len(group) < 10:
            continue
        group = group.sort_values("patient_id").copy()
        seed = stable_seed("rank", dataset, signature, window)
        rho, lo, hi, se_z, disp_boot_lo, disp_boot_hi, valid_b, se_source = bootstrap_spearman(group["baseline_score"].to_numpy(float), group["followup_score"].to_numpy(float), seed)
        n = len(group)
        baseline_pct = (stats.rankdata(group["baseline_score"], method="average") - 1) / max(n - 1, 1) * 100
        followup_pct = (stats.rankdata(group["followup_score"], method="average") - 1) / max(n - 1, 1) * 100
        displacement = np.abs(followup_pct - baseline_pct)
        for patient, b_pct, f_pct, value in zip(group["patient_id"], baseline_pct, followup_pct, displacement):
            displacement_rows.append({
                "dataset": dataset, "patient_id": patient, "signature_id": signature,
                "signature_name": SIGNATURE_NAMES[signature], "time_window": WINDOW_LABEL[window],
                "baseline_percentile": b_pct, "followup_percentile": f_pct,
                "absolute_rank_displacement_pp": value,
            })
        rows.append({
            "dataset": dataset,
            "signature_id": signature,
            "signature_name": SIGNATURE_NAMES[signature],
            "time_window": WINDOW_LABEL[window],
            "paired_n": n,
            "spearman_rho": rho,
            "bootstrap_ci95_lower": lo,
            "bootstrap_ci95_upper": hi,
            "bootstrap_se_fisher_z": se_z,
            "fisher_z_se_source": se_source,
            "median_absolute_rank_displacement_pp": float(np.median(displacement)),
            "q1_absolute_rank_displacement_pp": float(np.quantile(displacement, 0.25)),
            "q3_absolute_rank_displacement_pp": float(np.quantile(displacement, 0.75)),
            "bootstrap_ci95_lower_median_rank_displacement_pp": disp_boot_lo,
            "bootstrap_ci95_upper_median_rank_displacement_pp": disp_boot_hi,
            "bootstrap_replicates": B,
            "valid_bootstrap_replicates": valid_b,
            "seed": seed,
            "primary_independent_eligible": dataset in sets[(signature, window)],
        })
    cohort = pd.DataFrame(rows)
    meta_rows = []
    for (signature, landmark), group in cohort[cohort["primary_independent_eligible"]].groupby(["signature_id", "time_window"], sort=True):
        group = group[np.isfinite(group["bootstrap_se_fisher_z"]) & (group["bootstrap_se_fisher_z"] > 0)]
        if len(group) < 3:
            continue
        values = np.arctanh(np.clip(group["spearman_rho"].to_numpy(float), -0.999999, 0.999999))
        result = reml_hk(values, group["bootstrap_se_fisher_z"].to_numpy(float))
        meta_rows.append({
            "signature_id": signature,
            "signature_name": SIGNATURE_NAMES[signature],
            "time_window": landmark,
            "cohort_n": int(result["k"]),
            "patient_n": int(group["paired_n"].sum()),
            "pooled_spearman_rho": float(np.tanh(result["mu"])),
            "ci95_lower": float(np.tanh(result["lower"])),
            "ci95_upper": float(np.tanh(result["upper"])),
            "modified_hk_ci95_lower": float(np.tanh(result["modified_lower"])),
            "modified_hk_ci95_upper": float(np.tanh(result["modified_upper"])),
            "tau2_fisher_z": result["tau2"],
            "I2_percent": result["I2"],
            "prediction_lower": float(np.tanh(result["prediction_lower"])),
            "prediction_upper": float(np.tanh(result["prediction_upper"])),
            "modified_hk_prediction_lower": float(np.tanh(result["modified_prediction_lower"])),
            "modified_hk_prediction_upper": float(np.tanh(result["modified_prediction_upper"])),
            "p_value": result["p_value"],
            "modified_hk_p_value": result["modified_p_value"],
            "q_hk": result["q_hk"],
            "cohorts": ";".join(group["dataset"].tolist()),
            "median_cohort_rank_displacement_pp": float(group["median_absolute_rank_displacement_pp"].median()),
        })
    meta = pd.DataFrame(meta_rows)
    if not meta.empty:
        meta["holm_adjusted_p_diagnostic_family"] = np.nan
        for landmark, idx in meta[meta["signature_id"].isin(DIAGNOSTIC)].groupby("time_window").groups.items():
            meta.loc[idx, "holm_adjusted_p_diagnostic_family"] = holm(meta.loc[idx, "modified_hk_p_value"])
    return cohort, meta, pd.DataFrame(displacement_rows)


def concordance_module(paired: pd.DataFrame, sets: dict[tuple[str, str], set[str]]) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    eligible = paired[paired["time_window"].isin(["T1", "T2"]) & paired["signature_id"].isin(DIAGNOSTIC)]
    cohort_rows = []
    patient_rows = []
    for (dataset, window), group in eligible.groupby(["dataset", "time_window"], sort=True):
        wide = group.pivot_table(index="patient_id", columns="signature_id", values="delta_z", aggfunc="first").dropna(subset=DIAGNOSTIC)
        if len(wide) < 10:
            continue
        signs = np.sign(wide[DIAGNOSTIC].to_numpy(float))
        all_same = ((signs > 0).all(axis=1) | (signs < 0).all(axis=1))
        for patient, same, row_signs in zip(wide.index, all_same, signs):
            patient_rows.append({
                "dataset": dataset,
                "patient_id": patient,
                "time_window": WINDOW_LABEL[window],
                "all_four_same_nonzero_direction": bool(same),
                **{f"{signature}_direction": int(value) for signature, value in zip(DIAGNOSTIC, row_signs)},
            })
        for i, a in enumerate(DIAGNOSTIC):
            for b in DIAGNOSTIC[i + 1:]:
                seed = stable_seed("concordance", dataset, window, a, b)
                rho, lo, hi, se_z, _, _, valid_b, se_source = bootstrap_spearman(wide[a].to_numpy(float), wide[b].to_numpy(float), seed)
                sign_a = np.sign(wide[a].to_numpy(float))
                sign_b = np.sign(wide[b].to_numpy(float))
                exact_zero = (sign_a == 0) | (sign_b == 0)
                cohort_rows.append({
                    "dataset": dataset,
                    "time_window": WINDOW_LABEL[window],
                    "signature_a": a,
                    "signature_b": b,
                    "signature_pair": f"{a} vs {b}",
                    "paired_n": len(wide),
                    "spearman_rho_delta_z": rho,
                    "bootstrap_ci95_lower": lo,
                    "bootstrap_ci95_upper": hi,
                    "bootstrap_se_fisher_z": se_z,
                    "fisher_z_se_source": se_source,
                    "sign_agreement_fraction": float(np.mean(sign_a == sign_b)),
                    "exact_zero_in_either_fraction": float(np.mean(exact_zero)),
                    "all_four_same_direction_fraction": float(np.mean(all_same)),
                    "bootstrap_replicates": B,
                    "valid_bootstrap_replicates": valid_b,
                    "seed": seed,
                    "primary_independent_intersection_eligible": dataset in (sets[(a, window)] & sets[(b, window)]),
                })
    cohort = pd.DataFrame(cohort_rows)
    meta_rows = []
    for (landmark, a, b), group in cohort[cohort["primary_independent_intersection_eligible"]].groupby(["time_window", "signature_a", "signature_b"], sort=True):
        group = group[np.isfinite(group["bootstrap_se_fisher_z"]) & (group["bootstrap_se_fisher_z"] > 0)]
        if len(group) < 3:
            continue
        z = np.arctanh(np.clip(group["spearman_rho_delta_z"].to_numpy(float), -0.999999, 0.999999))
        result = reml_hk(z, group["bootstrap_se_fisher_z"].to_numpy(float))
        meta_rows.append({
            "time_window": landmark,
            "signature_a": a,
            "signature_b": b,
            "signature_pair": f"{a} vs {b}",
            "cohort_n": int(result["k"]),
            "patient_n": int(group["paired_n"].sum()),
            "pooled_spearman_rho_delta_z": float(np.tanh(result["mu"])),
            "ci95_lower": float(np.tanh(result["lower"])),
            "ci95_upper": float(np.tanh(result["upper"])),
            "modified_hk_ci95_lower": float(np.tanh(result["modified_lower"])),
            "modified_hk_ci95_upper": float(np.tanh(result["modified_upper"])),
            "tau2_fisher_z": result["tau2"],
            "I2_percent": result["I2"],
            "prediction_lower": float(np.tanh(result["prediction_lower"])),
            "prediction_upper": float(np.tanh(result["prediction_upper"])),
            "modified_hk_prediction_lower": float(np.tanh(result["modified_prediction_lower"])),
            "modified_hk_prediction_upper": float(np.tanh(result["modified_prediction_upper"])),
            "p_value": result["p_value"],
            "modified_hk_p_value": result["modified_p_value"],
            "q_hk": result["q_hk"],
            "cohorts": ";".join(group["dataset"].tolist()),
            "median_sign_agreement_fraction": float(group["sign_agreement_fraction"].median()),
        })
    meta = pd.DataFrame(meta_rows)
    if not meta.empty:
        meta["holm_adjusted_p_six_pairs"] = np.nan
        for landmark, idx in meta.groupby("time_window").groups.items():
            meta.loc[idx, "holm_adjusted_p_six_pairs"] = holm(meta.loc[idx, "modified_hk_p_value"])
    return cohort, meta, pd.DataFrame(patient_rows)


def contribution_module(paired: pd.DataFrame, gene: pd.DataFrame, sets: dict[tuple[str, str], set[str]]) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    sums = gene.groupby(["dataset", "patient_id", "signature_id", "time_window"], as_index=False)["standardized_contribution"].sum().rename(columns={"standardized_contribution": "reconstructed_delta_z"})
    target = paired[["dataset", "patient_id", "signature_id", "time_window", "delta_z"]]
    qc = sums.merge(target, on=["dataset", "patient_id", "signature_id", "time_window"], how="outer", validate="one_to_one", indicator=True)
    qc["absolute_error"] = (qc["reconstructed_delta_z"] - qc["delta_z"]).abs()
    qc["status"] = np.where((qc["_merge"] == "both") & (qc["absolute_error"] < 1e-8), "PASS", "FAIL")

    patient_rows = []
    for keys, group in gene.groupby(["dataset", "patient_id", "signature_id", "time_window"], sort=False):
        values = group["standardized_contribution"].to_numpy(float)
        abs_values = np.abs(values)
        abs_sum = float(abs_values.sum())
        index = int(np.argmax(abs_values))
        patient_rows.append({
            "dataset": keys[0], "patient_id": keys[1], "signature_id": keys[2], "time_window": keys[3],
            "total_delta_z": float(values.sum()),
            "absolute_contribution_sum": abs_sum,
            "dominance_ratio": float(abs_values[index] / abs_sum) if abs_sum > 0 else 0.0,
            "cancellation_index": float(1 - abs(values.sum()) / abs_sum) if abs_sum > 0 else 0.0,
            "dominant_gene": str(group.iloc[index]["gene"]),
        })
    patient_metrics = pd.DataFrame(patient_rows)

    gene_means = (
        gene.groupby(["dataset", "signature_id", "time_window", "gene"], as_index=False)["standardized_contribution"]
        .mean()
        .rename(columns={"standardized_contribution": "cohort_mean_signed_contribution"})
    )
    leader_lookup: dict[tuple[str, str, str], dict[str, object]] = {}
    for keys, group in gene_means.groupby(["dataset", "signature_id", "time_window"], sort=False):
        absolute = group["cohort_mean_signed_contribution"].abs()
        maximum = float(absolute.max())
        leaders = sorted(group.loc[np.abs(absolute - maximum) <= 1e-12, "gene"].astype(str).tolist())
        leader_lookup[keys] = {
            "cohort_mean_leading_gene_set": ";".join(leaders),
            "cohort_mean_leading_gene_tie_n": len(leaders),
            "cohort_mean_leading_absolute_contribution": maximum,
        }

    cohort_rows = []
    for keys, group in patient_metrics.groupby(["dataset", "signature_id", "time_window"], sort=False):
        counter = Counter(group["dominant_gene"])
        dominant_gene, dominant_n = counter.most_common(1)[0]
        top_n = dominant_n
        tied = sorted([gene_name for gene_name, count in counter.items() if count == top_n])
        cohort_rows.append({
            "dataset": keys[0], "signature_id": keys[1], "signature_name": SIGNATURE_NAMES[keys[1]],
            "time_window": WINDOW_LABEL.get(keys[2], keys[2]), "paired_n": len(group),
            "median_dominance_ratio": float(group["dominance_ratio"].median()),
            "median_cancellation_index": float(group["cancellation_index"].median()),
            "median_absolute_contribution_sum": float(group["absolute_contribution_sum"].median()),
            "cohort_dominant_gene": dominant_gene,
            "cohort_dominant_gene_patient_n": dominant_n,
            "cohort_dominant_gene_patient_fraction": dominant_n / len(group),
            "cohort_dominant_gene_tie_n": len(tied),
            "cohort_dominant_gene_tied_candidates": ";".join(tied),
            **leader_lookup[keys],
            "primary_independent_eligible": keys[0] in sets.get((keys[1], keys[2]), set()),
        })
    cohort = pd.DataFrame(cohort_rows)
    recurrence_rows = []
    for (signature, landmark), group in cohort[(cohort["time_window"].isin(["T24", "T48"])) & cohort["primary_independent_eligible"]].groupby(["signature_id", "time_window"], sort=True):
        leader_sets = [set(str(value).split(";")) for value in group["cohort_mean_leading_gene_set"]]
        all_genes = sorted(set().union(*leader_sets))
        counts = {gene_name: sum(gene_name in leaders for leaders in leader_sets) for gene_name in all_genes}
        top_n = max(counts.values())
        top_candidates = sorted([gene_name for gene_name, count in counts.items() if count == top_n])
        recurrence_rows.append({
            "signature_id": signature,
            "signature_name": SIGNATURE_NAMES[signature],
            "time_window": landmark,
            "eligible_cohort_n": len(group),
            "recurring_cohort_mean_leading_gene_set": ";".join(top_candidates),
            "cohort_mean_leading_gene_recurrence": top_n / len(group),
            "cohort_count_with_recurring_leader": top_n,
            "recurrence_tie_n": len(top_candidates),
            "recurrence_tied_candidates": ";".join(top_candidates),
            "cohorts": ";".join(group["dataset"].tolist()),
        })
    return qc.drop(columns="_merge"), patient_metrics, cohort, pd.DataFrame(recurrence_rows)


def main() -> None:
    paired = pd.read_parquet(DATA / "derived_patient_level" / "paired_changes_all_windows.parquet")
    gene = pd.read_parquet(DATA / "derived_patient_level" / "gene_contributions.parquet")
    sets = primary_sets()

    rank_cohort, rank_meta, displacement = rank_module(paired, sets)
    concordance_cohort, concordance_meta, direction = concordance_module(paired, sets)
    qc, patient_metrics, contribution_cohort, recurrence = contribution_module(paired, gene, sets)

    outputs = {
        "02_rank_persistence_cohort.csv": rank_cohort,
        "02_rank_persistence_meta.csv": rank_meta,
        "02_rank_displacement_patient.csv": displacement,
        "03_cross_signature_concordance_cohort.csv": concordance_cohort,
        "03_cross_signature_concordance_meta.csv": concordance_meta,
        "03_cross_signature_patient_direction.csv": direction,
        "04_contribution_reconstruction_qc.csv": qc,
        "04_contribution_patient_metrics.csv": patient_metrics,
        "04_contribution_cohort_dominant.csv": contribution_cohort,
        "04_contribution_recurrence.csv": recurrence,
    }
    for filename, frame in outputs.items():
        frame.to_csv(OUT / filename, index=False)

    failed = int((qc["status"] != "PASS").sum())
    summary = {
        "rank_cohort_rows": len(rank_cohort),
        "rank_meta_rows": len(rank_meta),
        "concordance_cohort_rows": len(concordance_cohort),
        "concordance_meta_rows": len(concordance_meta),
        "contribution_reconstruction_rows": len(qc),
        "contribution_reconstruction_failures": failed,
        "maximum_contribution_reconstruction_error": float(qc["absolute_error"].max()),
        "bootstrap_replicates": B,
        "status": "PASS" if failed == 0 else "PAUSE",
    }
    (OUT / "02_04_measurement_modules_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()

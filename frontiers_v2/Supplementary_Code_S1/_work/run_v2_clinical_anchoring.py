#!/usr/bin/env python3
"""Run source-defined clinical-course anchoring under the frozen V2 SAP."""

from __future__ import annotations

import hashlib
import json
import math
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
import statsmodels.formula.api as smf


ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "_sources" / "repo_subset_v1.3.1" / "data" / "derived_patient_level"
META = ROOT / "_sources" / "cohort_audit"
OUT = ROOT / "_work" / "analysis"
OUT.mkdir(parents=True, exist_ok=True)

PRIMARY_SIGNATURES = ["SIG001", "SIG002", "SIG003", "SIG004"]
ALL_SIGNATURES = PRIMARY_SIGNATURES + ["SIG022", "SIG023", "SIG033", "SIG034"]
SIGNATURE_NAMES = {
    "SIG001": "Sepsis MetaScore", "SIG002": "SeptiCyte LAB", "SIG003": "FAIM3:PLAC8", "SIG004": "sNIP",
    "SIG022": "Bacterial/Viral MetaScore", "SIG023": "Herberg two-gene DRS",
    "SIG033": "Lin seven-gene mortality score", "SIG034": "Severe-or-Mild score",
}
B = 2000

ENDPOINTS = {
    "GSE54514": {
        "endpoint": "source-defined survival status",
        "definition": "Survivor or non-survivor label deposited with GSE54514; analyzed cohort-specifically because the public GEO record does not state a mortality horizon.",
        "favorable": "Survivor", "adverse": "Non-survivor", "role": "secondary source-defined survival context",
        "source": "GEO GSE54514; Parnell et al., Shock 2013, doi:10.1097/SHK.0b013e31829ee604",
    },
    "GSE95233": {
        "endpoint": "28-day survival",
        "definition": "Day-28 survivor or non-survivor label deposited with GSE95233.",
        "favorable": "Survivor", "adverse": "Non-survivor", "role": "primary clinical-course anchor",
        "source": "GEO GSE95233; doi:10.1016/j.jinf.2018.11.012",
    },
    "GSE110487": {
        "endpoint": "source-defined organ-function course group",
        "definition": "Non-response was defined by source investigators as SOFA(T1)-SOFA(T2)<5 together with SOFA(T2)>8; all other patients were responders.",
        "favorable": "R", "adverse": "NR", "role": "primary clinical-course anchor",
        "source": "GEO GSE110487; Barcella et al., Crit Care 2018, doi:10.1186/s13054-018-2242-3",
    },
    "GSE57065": {
        "endpoint": "baseline SAPS II severity category",
        "definition": "Source-deposited SAPSII-High or SAPSII-Low category at septic-shock onset; treated as severity-stratified trajectory context.",
        "favorable": "SAPSII-Low", "adverse": "SAPSII-High", "role": "secondary severity context",
        "source": "GEO GSE57065; Cazalis et al., Intensive Care Med Exp 2014, doi:10.1186/s40635-014-0020-3",
    },
}


def stable_seed(*parts: str) -> int:
    digest = hashlib.sha256("|".join(map(str, parts)).encode()).digest()
    return int.from_bytes(digest[:8], "big") % (2**32 - 1)


def holm(values: pd.Series) -> pd.Series:
    values = pd.to_numeric(values, errors="coerce")
    result = pd.Series(np.nan, index=values.index, dtype=float)
    valid = values.dropna().sort_values()
    running = 0.0
    m = len(valid)
    for rank, (index, value) in enumerate(valid.items()):
        running = max(running, min(1.0, (m - rank) * float(value)))
        result.loc[index] = running
    return result


def bootstrap_group_difference(frame: pd.DataFrame, seed: int, adjust_day: bool = False) -> tuple[float, float, float]:
    rng = np.random.default_rng(seed)
    favorable = frame[frame["adverse"] == 0].copy()
    adverse = frame[frame["adverse"] == 1].copy()
    draws = np.empty(B, dtype=float)
    if not adjust_day:
        favorable_values = favorable["delta_z"].to_numpy(float)
        adverse_values = adverse["delta_z"].to_numpy(float)
        favorable_draws = favorable_values[rng.integers(0, len(favorable_values), size=(B, len(favorable_values)))].mean(axis=1)
        adverse_draws = adverse_values[rng.integers(0, len(adverse_values), size=(B, len(adverse_values)))].mean(axis=1)
        draws[:] = adverse_draws - favorable_draws
    else:
        # Status-stratified patient bootstrap with a fast closed OLS solve.  This
        # is numerically identical to refitting delta_z ~ adverse + day3, but
        # avoids thousands of high-overhead statsmodels object constructions.
        for iteration in range(B):
            sample = pd.concat([
                favorable.iloc[rng.integers(0, len(favorable), len(favorable))],
                adverse.iloc[rng.integers(0, len(adverse), len(adverse))],
            ], ignore_index=True)
            design = np.column_stack([
                np.ones(len(sample)),
                sample["adverse"].to_numpy(float),
                sample["followup_day3"].to_numpy(float),
            ])
            beta, *_ = np.linalg.lstsq(design, sample["delta_z"].to_numpy(float), rcond=None)
            draws[iteration] = beta[1]
    return float(draws.std(ddof=1)), float(np.quantile(draws, 0.025)), float(np.quantile(draws, 0.975))


def fit_contrast(frame: pd.DataFrame, dataset: str, signature: str, contrast: str, endpoint: dict, independent: bool) -> dict:
    long = pd.concat([
        frame[["patient_id", "clinical_status", "adverse", "baseline_score", "baseline_sd"]].assign(followup=0).rename(columns={"baseline_score": "score"}),
        frame[["patient_id", "clinical_status", "adverse", "followup_score", "baseline_sd"]].assign(followup=1).rename(columns={"followup_score": "score"}),
    ], ignore_index=True)
    long["score_z"] = long["score"] / long["baseline_sd"]
    if "followup_day3" in frame.columns:
        day_lookup = frame.set_index("patient_id")["followup_day3"].to_dict()
        long["followup_day3"] = long["patient_id"].map(day_lookup).astype(int)
        mixed_formula = "score_z ~ followup * adverse + followup * followup_day3 + followup_day3"
    else:
        mixed_formula = "score_z ~ followup * adverse"
    warning_text = ""
    model = None
    method_used = ""
    for method in ["lbfgs", "bfgs", "powell"]:
        try:
            with warnings.catch_warnings(record=True) as caught:
                warnings.simplefilter("always")
                candidate = smf.mixedlm(mixed_formula, long, groups=long["patient_id"]).fit(reml=True, method=method, maxiter=1000, disp=False)
            if bool(getattr(candidate, "converged", False)):
                model = candidate
                method_used = method
                warning_text = " | ".join(str(item.message) for item in caught)
                break
        except Exception as exc:
            warning_text = str(exc)
    if model is None:
        raise RuntimeError(f"{dataset} {signature} {contrast}: MixedLM failed")

    interaction = "followup:adverse"
    estimate = float(model.params[interaction])
    se = float(model.bse[interaction])
    p_value = float(model.pvalues[interaction])
    zcrit = float(stats.norm.ppf(0.975))

    change = frame[["patient_id", "clinical_status", "adverse", "delta_z"]].copy()
    if "followup_day3" in frame.columns:
        change["followup_day3"] = frame["followup_day3"].to_numpy(int)
        change_formula = "delta_z ~ adverse + followup_day3"
        adjustment = "follow-up day (D3 vs D2)"
    else:
        change_formula = "delta_z ~ adverse"
        adjustment = "none"
    ols = smf.ols(change_formula, change).fit(cov_type="HC3")
    seed = stable_seed("clinical", dataset, signature, contrast)
    boot_se, boot_lo, boot_hi = bootstrap_group_difference(change, seed, adjust_day="followup_day3" in change.columns)

    ancova = frame[["baseline_score", "followup_score", "baseline_sd", "adverse"]].copy()
    ancova["baseline_z"] = ancova["baseline_score"] / ancova["baseline_sd"]
    ancova["followup_z"] = ancova["followup_score"] / ancova["baseline_sd"]
    if "followup_day3" in frame.columns:
        ancova["followup_day3"] = frame["followup_day3"].to_numpy(int)
        ancova_formula = "followup_z ~ baseline_z + adverse + followup_day3"
    else:
        ancova_formula = "followup_z ~ baseline_z + adverse"
    ancova_fit = smf.ols(ancova_formula, ancova).fit(cov_type="HC3")

    return {
        "dataset": dataset,
        "signature_id": signature,
        "signature_name": SIGNATURE_NAMES[signature],
        "contrast": contrast,
        "endpoint": endpoint["endpoint"],
        "endpoint_definition": endpoint["definition"],
        "analysis_role": endpoint["role"],
        "contrast_role": (
            "primary source-defined contrast" if dataset == "GSE54514" and contrast == "T48"
            else "secondary source-defined contrast" if dataset == "GSE54514" and contrast == "T24"
            else "primary within-endpoint contrast"
        ),
        "independence_status": "PRIMARY_INDEPENDENT" if independent else "OVERLAP_FLAGGED_EXPLORATORY",
        "favorable_group": endpoint["favorable"],
        "adverse_group": endpoint["adverse"],
        "favorable_n": int((change["adverse"] == 0).sum()),
        "adverse_n": int((change["adverse"] == 1).sum()),
        "patient_n": int(change["patient_id"].nunique()),
        "observation_n": len(long),
        "interaction_beta_adverse_minus_favorable": estimate,
        "interaction_se": se,
        "ci95_lower": estimate - zcrit * se,
        "ci95_upper": estimate + zcrit * se,
        "p_value": p_value,
        "mixedlm_converged": bool(model.converged),
        "mixedlm_method": method_used,
        "mixedlm_warning": warning_text,
        "primary_inference_method": "paired-change contrast with HC3 covariance",
        "primary_adjustment": adjustment,
        "paired_change_hc3_beta": float(ols.params["adverse"]),
        "paired_change_hc3_se": float(ols.bse["adverse"]),
        "paired_change_hc3_ci95_lower": float(ols.conf_int().loc["adverse", 0]),
        "paired_change_hc3_ci95_upper": float(ols.conf_int().loc["adverse", 1]),
        "paired_change_hc3_p": float(ols.pvalues["adverse"]),
        "primary_beta_adverse_minus_favorable": float(ols.params["adverse"]),
        "primary_se_hc3": float(ols.bse["adverse"]),
        "primary_ci95_lower": float(ols.conf_int().loc["adverse", 0]),
        "primary_ci95_upper": float(ols.conf_int().loc["adverse", 1]),
        "primary_p_value": float(ols.pvalues["adverse"]),
        "bootstrap_se": boot_se,
        "bootstrap_ci95_lower": boot_lo,
        "bootstrap_ci95_upper": boot_hi,
        "bootstrap_replicates": B,
        "seed": seed,
        "ancova_adverse_beta": float(ancova_fit.params["adverse"]),
        "ancova_adverse_se_hc3": float(ancova_fit.bse["adverse"]),
        "ancova_adverse_p_hc3": float(ancova_fit.pvalues["adverse"]),
        "source_provenance": endpoint["source"],
    }


def assemble_contrasts(paired: pd.DataFrame, mapping: dict[str, pd.DataFrame]) -> list[tuple[str, str, pd.DataFrame]]:
    contrasts: list[tuple[str, str, pd.DataFrame]] = []
    # GSE54514: retain strict 24-h and 48-h paired contrasts separately.
    for window, label in [("T1", "T24"), ("T2", "T48")]:
        frame = paired[(paired["dataset"] == "GSE54514") & (paired["time_window"] == window)].copy()
        frame["clinical_status"] = frame["outcome"]
        contrasts.append(("GSE54514", label, frame))

    # GSE95233: every patient contributes one source-defined D2 or D3 follow-up.
    frame = paired[(paired["dataset"] == "GSE95233") & paired["time_window"].isin(["T1", "T2"])].copy()
    frame["clinical_status"] = frame["outcome"]
    frame["followup_day3"] = (frame["time_window"] == "T2").astype(int)
    if frame.groupby(["patient_id", "signature_id"]).size().max() != 1:
        raise RuntimeError("GSE95233 has duplicated second measurements")
    contrasts.append(("GSE95233", "Admission to D2/D3", frame))

    frame = paired[(paired["dataset"] == "GSE110487") & (paired["time_window"] == "T2")].copy()
    frame["clinical_status"] = frame["outcome"]
    contrasts.append(("GSE110487", "T48", frame))

    # GSE57065 severity comes from verified GEO sample metadata.
    severity = mapping["GSE57065"][["sample_id", "clinical_status"]].drop_duplicates()
    for window, label in [("T1", "T24"), ("T2", "T48")]:
        frame = paired[(paired["dataset"] == "GSE57065") & (paired["time_window"] == window)].copy()
        frame = frame.merge(severity.rename(columns={"sample_id": "baseline_sample_id", "clinical_status": "source_clinical_status"}), on="baseline_sample_id", how="left", validate="many_to_one")
        frame["clinical_status"] = frame["source_clinical_status"]
        contrasts.append(("GSE57065", label, frame))
    return contrasts


def metadata_audit(score: pd.DataFrame, mapping: dict[str, pd.DataFrame]) -> pd.DataFrame:
    rows = []
    sig1 = score[score["signature_id"] == "SIG001"]
    for dataset, source in mapping.items():
        registry = sig1[sig1["dataset"] == dataset][["sample_id", "patient_id", "time_window", "relative_hours", "outcome", "disease_status", "main_longitudinal_case"]].drop_duplicates()
        merged = registry.merge(source, on="sample_id", how="left", suffixes=("_pipeline", "_source"), indicator=True, validate="one_to_one")
        for _, row in merged.iterrows():
            endpoint = ENDPOINTS[dataset]
            source_status = row.get("clinical_status", "")
            if dataset == "GSE57065":
                clinical_status = source_status
            else:
                clinical_status = row["outcome"]
            case = str(row["main_longitudinal_case"]).upper() == "YES"
            mapping_status = row.get("mapping_status", "unmatched") if row["_merge"] == "both" else "unmatched"
            eligible = case and bool(clinical_status) and str(clinical_status).lower() != "nan" and mapping_status != "unmatched"
            rows.append({
                "dataset": dataset,
                "sample_id": row["sample_id"],
                "pipeline_patient_id": row["patient_id_pipeline"],
                "source_patient_id": row.get("patient_id_source", ""),
                "pipeline_time_window": row["time_window"],
                "relative_hours": row["relative_hours_pipeline"],
                "clinical_endpoint": endpoint["endpoint"],
                "clinical_status": clinical_status,
                "mapping_status": mapping_status,
                "analysis_eligible": "YES" if eligible else "NO",
                "exclusion_reason": "" if eligible else ("not a longitudinal case" if not case else "missing clinical label or mapping"),
                "source_provenance": endpoint["source"],
            })
    return pd.DataFrame(rows)


def trajectory_summary(score: pd.DataFrame, mapping: dict[str, pd.DataFrame]) -> pd.DataFrame:
    rows = []
    for dataset in ENDPOINTS:
        endpoint = ENDPOINTS[dataset]
        data = score[(score["dataset"] == dataset) & score["signature_id"].isin(PRIMARY_SIGNATURES) & (score["main_longitudinal_case"] == "YES")].copy()
        if dataset == "GSE57065":
            status = mapping[dataset][["sample_id", "clinical_status"]].drop_duplicates()
            data = data.merge(status, on="sample_id", how="left", validate="many_to_one")
        else:
            data["clinical_status"] = data["outcome"]
        data = data[data["clinical_status"].isin([endpoint["favorable"], endpoint["adverse"]])]
        if dataset == "GSE95233":
            data["display_time"] = np.where(data["time_window"] == "T0", "Admission", "D2/D3")
            # One source-defined second measurement per patient.
        else:
            display = {"T0": "Baseline", "T1": "T24", "T2": "T48", "T3": "T72", "T4": "Day 5"}
            data["display_time"] = data["time_window"].map(display)
        for signature, sig_data in data.groupby("signature_id"):
            baseline = sig_data[sig_data["time_window"] == "T0"]["oriented_score"]
            mean, sd = float(baseline.mean()), float(baseline.std(ddof=1))
            sig_data = sig_data.assign(score_z=(sig_data["oriented_score"] - mean) / sd)
            for (time, status), group in sig_data.groupby(["display_time", "clinical_status"], sort=False):
                values = group["score_z"].to_numpy(float)
                seed = stable_seed("trajectory", dataset, signature, time, status)
                rng = np.random.default_rng(seed)
                draws = values[rng.integers(0, len(values), size=(B, len(values)))].mean(axis=1)
                rows.append({
                    "dataset": dataset, "signature_id": signature, "signature_name": SIGNATURE_NAMES[signature],
                    "display_time": time, "clinical_status": status, "n": len(values),
                    "mean_score_z": float(values.mean()), "se_score_z": float(values.std(ddof=1) / math.sqrt(len(values))) if len(values) > 1 else math.nan,
                    "bootstrap_ci95_lower": float(np.quantile(draws, 0.025)), "bootstrap_ci95_upper": float(np.quantile(draws, 0.975)),
                    "endpoint": endpoint["endpoint"], "analysis_role": endpoint["role"], "seed": seed,
                })
    return pd.DataFrame(rows)


def main() -> None:
    paired = pd.read_parquet(BASE / "paired_changes_all_windows.parquet")
    score = pd.read_parquet(BASE / "score_matrix.parquet")
    mapping = {dataset: pd.read_csv(META / f"{dataset}_verified_sample_mapping.csv", dtype=str) for dataset in ENDPOINTS}
    meta = pd.read_csv(ROOT / "_sources" / "repo_subset_v1.3.1" / "data" / "source_tables" / "stage3" / "03_05_signature_level_meta_analysis.csv")
    meta = meta[(meta["analysis_set"] == "PRIMARY_INDEPENDENT") & meta["time_window"].isin(["T1", "T2"])]
    independent_sets = {(row.signature_id, row.time_window): set(str(row.cohorts).split(";")) for row in meta.itertuples()}

    audit = metadata_audit(score, mapping)
    unresolved_eligible = audit[(audit["analysis_eligible"] == "YES") & audit["mapping_status"].isin(["unmatched", "unresolved"])]
    if len(unresolved_eligible):
        raise RuntimeError("Unresolved sample mappings entered the clinical analysis")

    rows = []
    for dataset, contrast, frame in assemble_contrasts(paired, mapping):
        endpoint = ENDPOINTS[dataset]
        frame = frame[frame["clinical_status"].isin([endpoint["favorable"], endpoint["adverse"]])].copy()
        frame["adverse"] = (frame["clinical_status"] == endpoint["adverse"]).astype(int)
        for signature in ALL_SIGNATURES:
            sig_frame = frame[frame["signature_id"] == signature].copy()
            if sig_frame["adverse"].nunique() != 2:
                raise RuntimeError(f"{dataset} {contrast} {signature}: clinical groups incomplete")
            if dataset == "GSE95233":
                independent = dataset in independent_sets[(signature, "T1")] and dataset in independent_sets[(signature, "T2")]
            else:
                source_window = "T1" if contrast == "T24" else "T2"
                independent = dataset in independent_sets[(signature, source_window)]
            rows.append(fit_contrast(sig_frame, dataset, signature, contrast, endpoint, independent))
    results = pd.DataFrame(rows)
    results["holm_adjusted_p_primary_four"] = np.nan
    for (dataset, contrast), index in results[results["signature_id"].isin(PRIMARY_SIGNATURES)].groupby(["dataset", "contrast"]).groups.items():
        results.loc[index, "holm_adjusted_p_primary_four"] = holm(results.loc[index, "paired_change_hc3_p"])

    # Prespecified D2- and D3-specific sensitivity estimates for GSE95233.
    day_rows = []
    endpoint = ENDPOINTS["GSE95233"]
    for window, label in [("T1", "Admission to D2"), ("T2", "Admission to D3")]:
        day_frame = paired[(paired["dataset"] == "GSE95233") & (paired["time_window"] == window)].copy()
        day_frame["clinical_status"] = day_frame["outcome"]
        day_frame = day_frame[day_frame["clinical_status"].isin([endpoint["favorable"], endpoint["adverse"]])].copy()
        day_frame["adverse"] = (day_frame["clinical_status"] == endpoint["adverse"]).astype(int)
        for signature in ALL_SIGNATURES:
            sig_frame = day_frame[day_frame["signature_id"] == signature].copy()
            independent = "GSE95233" in independent_sets[(signature, window)]
            fitted = fit_contrast(sig_frame, "GSE95233", signature, label, endpoint, independent)
            fitted["analysis_role"] = "follow-up-day-stratified sensitivity"
            fitted["contrast_role"] = "sensitivity"
            day_rows.append(fitted)
    day_results = pd.DataFrame(day_rows)
    day_results["holm_adjusted_p_primary_four"] = np.nan
    for contrast, index in day_results[day_results["signature_id"].isin(PRIMARY_SIGNATURES)].groupby("contrast").groups.items():
        day_results.loc[index, "holm_adjusted_p_primary_four"] = holm(day_results.loc[index, "primary_p_value"])

    trajectories = trajectory_summary(score, mapping)
    audit.to_csv(OUT / "05_clinical_metadata_audit.csv", index=False)
    results.to_csv(OUT / "05_clinical_course_interactions.csv", index=False)
    day_results.to_csv(OUT / "05_gse95233_day_stratified_sensitivity.csv", index=False)
    trajectories.to_csv(OUT / "05_clinical_course_trajectories.csv", index=False)
    summary = {
        "clinical_mapping_rows": len(audit),
        "unresolved_eligible_mappings": len(unresolved_eligible),
        "interaction_rows": len(results),
        "primary_interaction_rows": int(results["signature_id"].isin(PRIMARY_SIGNATURES).sum()),
        "all_mixed_models_converged": bool(results["mixedlm_converged"].all()),
        "independent_course_anchor_cohorts": 2,
        "secondary_severity_context_cohorts": 1,
        "secondary_source_defined_survival_cohorts": 1,
        "gse95233_day_stratified_rows": len(day_results),
        "status": "PASS" if results["mixedlm_converged"].all() and len(unresolved_eligible) == 0 else "PAUSE",
    }
    (OUT / "05_clinical_course_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()

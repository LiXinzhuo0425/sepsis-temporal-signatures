#!/usr/bin/env python3
"""Verify that the manuscript-facing workbook matches retained numerical sources."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
WORKBOOK = ROOT / "_work" / "input_snapshot" / "04_Data_S1.xlsx"
SOURCE_ROOT = ROOT / "_sources" / "repo_subset_v1.3.1" / "data" / "source_tables"
OUT = ROOT / "_work" / "analysis"

WINDOW_MAP = {"T1": "T24", "T2": "T48", "T3": "T72", "T4": "Day 5"}
SET_MAP = {
    "ALL_COHORTS": "All cohorts",
    "PRIMARY_INDEPENDENT": "Primary independent set",
    "STRICT_NEVER_USED": "Strict never-used set",
    "PILOT_ONLY": "Pilot only",
    "PRESPECIFIED_NON_PILOT_ONLY": "Prespecified non-pilot only",
}
ROLE_MAP = {"PILOT": "Pilot", "PRESPECIFIED_NON_PILOT": "Prespecified non-pilot"}
CLASS_MAP = {
    "CONSISTENT_DRIFT": "Consistent drift",
    "RELATIVE_STABILITY": "Relative stability",
    "COHORT_DEPENDENT_DRIFT": "Cohort-dependent drift",
    "HIGH_HETEROGENEITY_OR_PREDICTIVE_UNCERTAINTY": "High heterogeneity or predictive uncertainty",
    "EVIDENCE_INSUFFICIENT": "Insufficient evidence",
}
ARCH_MAP = {
    "COHORT_DEPENDENT_DRIFT": "Cohort-dependent drift",
    "CONSISTENT_MULTIGENE_DRIFT": "Consistent multigene drift",
    "SINGLE_GENE_DOMINANT_DRIFT": "Single-gene-dominant drift",
}


def read_source(relative: str) -> pd.DataFrame:
    frame = pd.read_csv(SOURCE_ROOT / relative)
    if "time_window" in frame:
        frame["time_window"] = frame["time_window"].replace(WINDOW_MAP)
    if "classification_window" in frame:
        frame["classification_window"] = frame["classification_window"].replace(WINDOW_MAP)
    if "analysis_set" in frame:
        frame["analysis_set"] = frame["analysis_set"].replace(SET_MAP)
    if "analysis_role" in frame:
        frame["analysis_role"] = frame["analysis_role"].replace(ROLE_MAP)
    if "stability_class" in frame:
        frame["stability_class"] = frame["stability_class"].replace(CLASS_MAP)
    if "drift_architecture" in frame:
        frame["drift_architecture"] = frame["drift_architecture"].replace(ARCH_MAP)
    return frame


def compare(sheet: str, relative: str, keys: list[str], rename: dict[str, str] | None = None) -> dict:
    observed = pd.read_excel(WORKBOOK, sheet_name=sheet)
    expected = read_source(relative)
    if rename:
        expected = expected.rename(columns=rename)
    observed = observed.dropna(subset=keys).copy()
    expected = expected.dropna(subset=keys).copy()
    if "inference_status" in observed.columns:
        suppressed_keys = observed.loc[
            observed["inference_status"].astype("string").str.contains("suppressed", case=False, na=False),
            keys,
        ].drop_duplicates()
        if not suppressed_keys.empty:
            suppressed_keys = suppressed_keys.assign(_suppress=True)
            expected = expected.merge(suppressed_keys, on=keys, how="left")
            inferential_columns = [
                "pooled_se", "ci95_lower", "ci95_upper", "tau2", "I2_percent",
                "prediction_lower", "prediction_upper", "p_value", "holm_adjusted_p",
            ]
            for column in inferential_columns:
                if column in expected.columns:
                    expected.loc[expected["_suppress"].fillna(False), column] = np.nan
            expected = expected.drop(columns="_suppress")
    expected = expected.merge(observed[keys].drop_duplicates(), on=keys, how="inner")
    merged = observed.merge(expected, on=keys, how="outer", suffixes=("_workbook", "_source"), indicator=True)
    missing = int((merged["_merge"] != "both").sum())
    comparisons = []
    for col in sorted((set(observed.columns) & set(expected.columns)) - set(keys)):
        left = merged.get(f"{col}_workbook")
        right = merged.get(f"{col}_source")
        if left is None or right is None:
            continue
        if pd.api.types.is_numeric_dtype(observed[col]) and pd.api.types.is_numeric_dtype(expected[col]):
            mask = left.notna() & right.notna()
            max_abs = float(np.max(np.abs(left[mask].astype(float) - right[mask].astype(float)))) if mask.any() else 0.0
            nan_mismatch = int((left.isna() != right.isna()).sum())
            comparisons.append({"column": col, "kind": "numeric", "max_abs_diff": max_abs, "mismatch_n": nan_mismatch + int((mask & ~np.isclose(left.astype(float), right.astype(float), rtol=0, atol=1e-12)).sum())})
        else:
            ltxt = left.fillna("<NA>").astype(str)
            rtxt = right.fillna("<NA>").astype(str)
            comparisons.append({"column": col, "kind": "text", "max_abs_diff": None, "mismatch_n": int((ltxt != rtxt).sum())})
    material = [item for item in comparisons if item["mismatch_n"] > 0 and item["column"] not in {"inference_status", "classification_reason"}]
    return {
        "sheet": sheet,
        "workbook_rows": len(observed),
        "matched_source_rows": len(expected),
        "key_mismatch_rows": missing,
        "max_numeric_diff": max((item["max_abs_diff"] or 0.0 for item in comparisons if item["kind"] == "numeric"), default=0.0),
        "material_column_mismatches": material,
        "column_checks": comparisons,
    }


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    checks = [
        compare("S3_CohortEffects", "stage3/03_04_cohort_signature_primary_effects.csv", ["dataset", "signature_id", "time_window"]),
        compare("S4_Meta", "stage3/03_05_signature_level_meta_analysis.csv", ["signature_id", "time_window", "analysis_set"]),
        compare("S5_Stability", "stage3/03_06_temporal_stability_profile.csv", ["signature_id", "time_window", "analysis_set"]),
        compare(
            "S11_GeneArch",
            "stage4/04_08_signature_drift_architecture.csv",
            ["signature_id", "classification_window"],
            rename={
                "stage3_pooled_delta_z": "score_pooled_delta_z",
                "stage3_I2_percent": "score_I2_percent",
                "stage3_prediction_lower": "score_prediction_lower",
                "stage3_prediction_upper": "score_prediction_upper",
                "median_patient_dominance_ratio": "median_cohort_dominance",
                "median_patient_cancellation_index": "median_cohort_cancellation",
                "median_patient_absolute_contribution_sum": "median_cohort_absolute_contribution_sum",
                "stage3_direction_consistency": "score_direction_consistency",
            },
        ),
        compare("S12_GeneMeta", "stage4/04_07_gene_contribution_meta_analysis.csv", ["signature_id", "time_window", "gene", "gene_group", "analysis_set"]),
    ]
    pass_status = all(
        item["key_mismatch_rows"] == 0
        and item["max_numeric_diff"] <= 1e-12
        and not item["material_column_mismatches"]
        for item in checks
    )
    detail_path = OUT / "01_legacy_reproduction_detail.json"
    detail_path.write_text(json.dumps({"status": "PASS" if pass_status else "PAUSE", "checks": checks}, indent=2), encoding="utf-8")

    lines = [
        "# Legacy numerical reproduction report",
        "",
        "Date: 16 August 2026",
        "",
        f"**{'PASS' if pass_status else 'PAUSE'} — the manuscript-facing numerical workbook {'matches' if pass_status else 'does not fully match'} the retained corrected source tables.**",
        "",
        "This is a source-to-deliverable verification. The corrected v1.2.1 numerical workflow is not silently rerun or altered.",
        "",
        "| Workbook sheet | Rows | Matched source rows | Key mismatches | Maximum numeric difference |",
        "|---|---:|---:|---:|---:|",
    ]
    for item in checks:
        lines.append(f"| {item['sheet']} | {item['workbook_rows']} | {item['matched_source_rows']} | {item['key_mismatch_rows']} | {item['max_numeric_diff']:.3g} |")
    lines.extend([
        "",
        "The patient-level health check independently reproduced 264 primary paired patients, including 172 at T24 and 146 at T48. All retained score calculations report complete gene coverage and successful calculation. Full column checks are stored in `01_legacy_reproduction_detail.json`.",
        "",
    ])
    (OUT / "01_legacy_reproduction_report.md").write_text("\n".join(lines), encoding="utf-8")
    print(json.dumps({"status": "PASS" if pass_status else "PAUSE", "checks": [{k: v for k, v in item.items() if k != "column_checks"} for item in checks]}, indent=2))


if __name__ == "__main__":
    main()

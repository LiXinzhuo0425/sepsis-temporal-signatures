#!/usr/bin/env python3
"""Descriptive common-cohort sensitivity for T24/T48 signature changes.

This script imports the frozen Stage 3 REML/Hartung-Knapp implementation,
verifies that it reconstructs the stored primary results, and then evaluates
two transparent cohort restrictions without altering any manuscript files.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd


RELEASE_ROOT = Path(__file__).resolve().parents[1]
STAGE3_ROOT = RELEASE_ROOT / "reproduction_project/longitudinal_stage3"
STAGE3_CODE = RELEASE_ROOT / "code/portable_analysis/03_00_09_environment_lock/run_stage3.py"
EFFECTS_CSV = STAGE3_ROOT / "source_data/03_04_cohort_signature_primary_effects.csv"
META_CSV = STAGE3_ROOT / "source_data/03_05_signature_level_meta_analysis.csv"
DEFAULT_OUTPUT_DIR = RELEASE_ROOT / "reproducibility_evidence/common_cohort_sensitivity"

SIGNATURE_NAMES = {
    "SIG001": "Sepsis MetaScore",
    "SIG002": "SeptiCyte LAB",
    "SIG003": "FAIM3:PLAC8 ratio",
    "SIG004": "sNIP",
    "SIG022": "Bacterial/Viral MetaScore",
    "SIG023": "Herberg Disease Risk Score",
    "SIG033": "Lin seven-gene mortality score",
    "SIG034": "Severe-or-Mild score",
}
SIGNATURES = tuple(SIGNATURE_NAMES)
WINDOW_LABELS = {"T1": "T24", "T2": "T48"}
LITERAL_COMMON = ("GSE54514", "GSE57065", "GSE95233")
DIAGNOSTIC = {"SIG001", "SIG002", "SIG003", "SIG004"}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def release_relative(path: Path) -> str:
    """Return a stable POSIX path relative to the extracted release root."""
    return path.resolve().relative_to(RELEASE_ROOT).as_posix()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
        help=(
            "Directory for the five evidence files. Relative paths are resolved "
            "from the extracted release root."
        ),
    )
    return parser.parse_args()


def resolve_output_dir(path: Path) -> Path:
    return path.resolve() if path.is_absolute() else (RELEASE_ROOT / path).resolve()


def sign_label(value: float, tol: float = 1e-12) -> str:
    if value > tol:
        return "increase"
    if value < -tol:
        return "decrease"
    return "zero"


def import_frozen_stage3():
    os.environ["SEPSIS_SIGNATURE_ANALYSIS_ROOT"] = str(STAGE3_ROOT)
    sys.path.insert(0, str(STAGE3_CODE.parent))
    spec = importlib.util.spec_from_file_location("frozen_stage3_common_cohort", STAGE3_CODE)
    if spec is None or spec.loader is None:
        raise RuntimeError("Unable to load frozen Stage 3 module")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def parse_cohorts(value: str) -> tuple[str, ...]:
    return tuple(sorted(item for item in str(value).split(";") if item))


def meta_row(stage3, subset: pd.DataFrame) -> dict[str, float | int | str | None]:
    subset = subset.sort_values("dataset").copy()
    if subset.empty:
        raise ValueError("Empty meta-analysis subset")
    result = stage3.reml_meta(
        subset["mean_delta_z"].to_numpy(dtype=float),
        subset["bootstrap_se"].to_numpy(dtype=float),
    )
    pooled_sign = np.sign(float(result["mu"]))
    directions = np.sign(subset["mean_delta_z"].to_numpy(dtype=float))
    direction_consistency = (
        float(np.mean(directions == pooled_sign)) if pooled_sign != 0 else float(np.mean(directions == 0))
    )
    formal = len(subset) >= 3
    return {
        "cohort_n": int(len(subset)),
        "patient_n": int(subset["paired_n"].sum()),
        "cohorts": ";".join(subset["dataset"].tolist()),
        "pooled_delta_z": float(result["mu"]),
        "pooled_se_hk": float(result["se"]) if formal else None,
        "ci95_lower_hk": float(result["lower"]) if formal else None,
        "ci95_upper_hk": float(result["upper"]) if formal else None,
        "tau2_reml": float(result["tau2"]),
        "I2_percent": float(result["I2"]),
        "prediction_lower": float(result["prediction_lower"]) if formal else None,
        "prediction_upper": float(result["prediction_upper"]) if formal else None,
        "direction_consistency": direction_consistency,
        "direction": sign_label(float(result["mu"])),
        "inference_status": (
            "Formal REML random-effects point estimate with standard HK interval"
            if formal
            else "Descriptive two-cohort REML point estimate; formal HK interval suppressed"
        ),
    }


def main(output_dir: Path) -> None:
    stage3 = import_frozen_stage3()
    effects = pd.read_csv(EFFECTS_CSV)
    stored = pd.read_csv(META_CSV)
    primary = stored[
        (stored["analysis_set"] == "PRIMARY_INDEPENDENT")
        & stored["signature_id"].isin(SIGNATURES)
        & stored["time_window"].isin(WINDOW_LABELS)
    ].copy()
    if len(primary) != 16:
        raise AssertionError(f"Expected 16 stored primary comparisons, found {len(primary)}")

    reconstruction = []
    primary_cohorts: dict[tuple[str, str], tuple[str, ...]] = {}
    for row in primary.itertuples(index=False):
        cohorts = parse_cohorts(row.cohorts)
        primary_cohorts[(row.signature_id, row.time_window)] = cohorts
        subset = effects[
            (effects["signature_id"] == row.signature_id)
            & (effects["time_window"] == row.time_window)
            & effects["dataset"].isin(cohorts)
        ]
        rebuilt = stage3.reml_meta(
            subset["mean_delta_z"].to_numpy(dtype=float),
            subset["bootstrap_se"].to_numpy(dtype=float),
        )
        reconstruction.append(
            {
                "signature_id": row.signature_id,
                "time_window": WINDOW_LABELS[row.time_window],
                "mu_abs_error": abs(float(rebuilt["mu"]) - float(row.pooled_delta_z)),
                "tau2_abs_error": abs(float(rebuilt["tau2"]) - float(row.tau2)),
                "hk_se_abs_error": abs(float(rebuilt["se"]) - float(row.pooled_se)),
            }
        )
    recon = pd.DataFrame(reconstruction)
    component_max_errors = {
        column: float(recon[column].max())
        for column in ("mu_abs_error", "tau2_abs_error", "hk_se_abs_error")
    }
    max_error = max(component_max_errors.values())
    if max_error > 1e-8:
        raise AssertionError(f"Frozen method failed primary-result reconstruction: max error={max_error}")

    result_rows: list[dict] = []
    definitions: dict[str, dict[str, tuple[str, ...]]] = {
        "FIXED_THREE_LANDMARK_COHORTS": {signature: LITERAL_COMMON for signature in SIGNATURES},
        "SHARED_PRIMARY_ELIGIBLE_COHORTS": {
            signature: tuple(
                sorted(
                    set(primary_cohorts[(signature, "T1")])
                    & set(primary_cohorts[(signature, "T2")])
                )
            )
            for signature in SIGNATURES
        },
    }

    for definition, sets_by_signature in definitions.items():
        for signature in SIGNATURES:
            selected = sets_by_signature[signature]
            for internal_window, display_window in WINDOW_LABELS.items():
                subset = effects[
                    (effects["signature_id"] == signature)
                    & (effects["time_window"] == internal_window)
                    & effects["dataset"].isin(selected)
                ].copy()
                found = tuple(sorted(subset["dataset"].tolist()))
                if found != tuple(sorted(selected)):
                    raise AssertionError(
                        f"Missing expected cohort effects for {definition}, {signature}, {display_window}: "
                        f"expected={selected}, found={found}"
                    )
                full_cohorts = set(primary_cohorts[(signature, internal_window)])
                excluded_from_primary = tuple(sorted(set(selected) - full_cohorts))
                item = {
                    "analysis_definition": definition,
                    "signature_id": signature,
                    "signature_name": SIGNATURE_NAMES[signature],
                    "signature_family": "diagnostic" if signature in DIAGNOSTIC else "secondary",
                    "time_window": display_window,
                    "internal_time_window": internal_window,
                    **meta_row(stage3, subset),
                    "cohorts_excluded_from_primary_for_this_signature_window": ";".join(excluded_from_primary),
                    "primary_eligibility_preserved": len(excluded_from_primary) == 0,
                }
                result_rows.append(item)

    results = pd.DataFrame(result_rows)
    primary_display = primary.copy()
    primary_display["time_window"] = primary_display["time_window"].map(WINDOW_LABELS)

    comparison_rows: list[dict] = []
    for definition in definitions:
        for signature in SIGNATURES:
            restricted = results[
                (results["analysis_definition"] == definition)
                & (results["signature_id"] == signature)
            ].set_index("time_window")
            full = primary_display[primary_display["signature_id"] == signature].set_index("time_window")
            full_t24 = float(full.loc["T24", "pooled_delta_z"])
            full_t48 = float(full.loc["T48", "pooled_delta_z"])
            res_t24 = float(restricted.loc["T24", "pooled_delta_z"])
            res_t48 = float(restricted.loc["T48", "pooled_delta_z"])
            full_pair = f"{sign_label(full_t24)}->{sign_label(full_t48)}"
            res_pair = f"{sign_label(res_t24)}->{sign_label(res_t48)}"
            full_abs_relation = "T48 larger" if abs(full_t48) > abs(full_t24) else "T24 larger_or_equal"
            res_abs_relation = "T48 larger" if abs(res_t48) > abs(res_t24) else "T24 larger_or_equal"
            comparison_rows.append(
                {
                    "analysis_definition": definition,
                    "signature_id": signature,
                    "signature_name": SIGNATURE_NAMES[signature],
                    "signature_family": "diagnostic" if signature in DIAGNOSTIC else "secondary",
                    "full_T24_delta_z": full_t24,
                    "restricted_T24_delta_z": res_t24,
                    "T24_direction_concordant": sign_label(full_t24) == sign_label(res_t24),
                    "full_T48_delta_z": full_t48,
                    "restricted_T48_delta_z": res_t48,
                    "T48_direction_concordant": sign_label(full_t48) == sign_label(res_t48),
                    "full_direction_pair": full_pair,
                    "restricted_direction_pair": res_pair,
                    "direction_pair_concordant": full_pair == res_pair,
                    "full_T48_minus_T24": full_t48 - full_t24,
                    "restricted_T48_minus_T24": res_t48 - res_t24,
                    "cross_window_difference_direction_concordant": sign_label(full_t48 - full_t24)
                    == sign_label(res_t48 - res_t24),
                    "full_absolute_magnitude_relation": full_abs_relation,
                    "restricted_absolute_magnitude_relation": res_abs_relation,
                    "absolute_magnitude_relation_concordant": full_abs_relation == res_abs_relation,
                    "primary_eligibility_preserved": bool(
                        restricted["primary_eligibility_preserved"].all()
                    ),
                    "cohorts": restricted.loc["T24", "cohorts"],
                    "T24_patient_n": int(restricted.loc["T24", "patient_n"]),
                    "T48_patient_n": int(restricted.loc["T48", "patient_n"]),
                }
            )
    comparisons = pd.DataFrame(comparison_rows)

    common_inputs = effects[
        effects["signature_id"].isin(SIGNATURES)
        & effects["time_window"].isin(WINDOW_LABELS)
        & effects["dataset"].isin(LITERAL_COMMON)
    ].copy()
    common_inputs["internal_time_window"] = common_inputs["time_window"]
    common_inputs["time_window"] = common_inputs["time_window"].map(WINDOW_LABELS)
    common_inputs["sampling_variance"] = common_inputs["bootstrap_se"].astype(float) ** 2
    common_inputs["primary_independent_eligible"] = common_inputs.apply(
        lambda row: row["dataset"]
        in primary_cohorts[(row["signature_id"], row["internal_time_window"])],
        axis=1,
    )
    common_inputs["primary_exclusion_reason"] = np.where(
        common_inputs["primary_independent_eligible"],
        "",
        "POSSIBLE_SAME_MARS_PROGRAM",
    )
    common_inputs = common_inputs[
        [
            "signature_id",
            "time_window",
            "internal_time_window",
            "dataset",
            "analysis_role",
            "paired_n",
            "mean_delta_z",
            "bootstrap_se",
            "sampling_variance",
            "ci95_lower",
            "ci95_upper",
            "primary_independent_eligible",
            "primary_exclusion_reason",
        ]
    ].sort_values(["signature_id", "time_window", "dataset"])
    if len(common_inputs) != 48 or common_inputs.duplicated(
        ["signature_id", "time_window", "dataset"]
    ).any():
        raise AssertionError("Common-cohort input table must contain 48 unique rows")

    summaries: dict[str, dict] = {}
    for definition in definitions:
        block = comparisons[comparisons["analysis_definition"] == definition]
        diag = block[block["signature_family"] == "diagnostic"]
        summaries[definition] = {
            "all_signatures_landmark_direction_concordance": {
                "numerator": int(block[["T24_direction_concordant", "T48_direction_concordant"]].to_numpy().sum()),
                "denominator": int(block.shape[0] * 2),
            },
            "all_signatures_direction_pair_concordance": {
                "numerator": int(block["direction_pair_concordant"].sum()),
                "denominator": int(block.shape[0]),
            },
            "all_signatures_cross_window_difference_direction_concordance": {
                "numerator": int(block["cross_window_difference_direction_concordant"].sum()),
                "denominator": int(block.shape[0]),
            },
            "all_signatures_absolute_magnitude_relation_concordance": {
                "numerator": int(block["absolute_magnitude_relation_concordant"].sum()),
                "denominator": int(block.shape[0]),
            },
            "diagnostic_landmark_direction_concordance": {
                "numerator": int(diag[["T24_direction_concordant", "T48_direction_concordant"]].to_numpy().sum()),
                "denominator": int(diag.shape[0] * 2),
            },
            "diagnostic_direction_pair_concordance": {
                "numerator": int(diag["direction_pair_concordant"].sum()),
                "denominator": int(diag.shape[0]),
            },
            "diagnostic_cross_window_difference_direction_concordance": {
                "numerator": int(diag["cross_window_difference_direction_concordant"].sum()),
                "denominator": int(diag.shape[0]),
            },
            "diagnostic_absolute_magnitude_relation_concordance": {
                "numerator": int(diag["absolute_magnitude_relation_concordant"].sum()),
                "denominator": int(diag.shape[0]),
            },
        }

    output_dir.mkdir(parents=True, exist_ok=True)
    results.to_csv(output_dir / "common_cohort_results.csv", index=False)
    comparisons.to_csv(output_dir / "common_cohort_comparison.csv", index=False)
    common_inputs.to_csv(output_dir / "common_cohort_inputs.csv", index=False)
    recon.to_csv(output_dir / "primary_reconstruction_check.csv", index=False)
    summary = {
        "analysis_date": "2026-08-08",
        "purpose": "Descriptive common-cohort sensitivity; no new multiplicity testing",
        "canonical_method": "Frozen Stage 3 REML random-effects point estimate; standard Hartung-Knapp interval only when k>=3",
        "source_files": {
            release_relative(EFFECTS_CSV): sha256(EFFECTS_CSV),
            release_relative(META_CSV): sha256(META_CSV),
            release_relative(STAGE3_CODE): sha256(STAGE3_CODE),
        },
        "primary_reconstruction": {
            "comparisons": int(len(recon)),
            "max_absolute_numeric_error": max_error,
            "component_max_absolute_errors": component_max_errors,
            "tolerance": 1e-8,
            "status": "PASS",
        },
        "definitions": {
            "FIXED_THREE_LANDMARK_COHORTS": {
                "cohorts": list(LITERAL_COMMON),
                "note": "Fixed three-cohort definition; includes GSE54514 for SIG002-SIG004 even though that cohort was excluded from their primary syntheses because of possible same-program reuse.",
            },
            "SHARED_PRIMARY_ELIGIBLE_COHORTS": {
                "cohorts_by_signature": {
                    signature: list(definitions["SHARED_PRIMARY_ELIGIBLE_COHORTS"][signature])
                    for signature in SIGNATURES
                },
                "note": "Intersection of the signature-specific PRIMARY_INDEPENDENT cohort sets at T24 and T48; preserves original eligibility rules.",
            },
        },
        "concordance": summaries,
    }
    (output_dir / "common_cohort_summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    arguments = parse_args()
    main(resolve_output_dir(arguments.output_dir))

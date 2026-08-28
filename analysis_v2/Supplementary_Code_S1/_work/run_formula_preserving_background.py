#!/usr/bin/env python3
"""Formula-preserving temporal background benchmark in GSE106878.

Gene replacement is matched only on baseline expression mean and variance.
Published algebra, coefficients, component count and score orientation remain
fixed.  No outcome, treatment or follow-up information enters matching.
"""

from __future__ import annotations

import csv
import gzip
import hashlib
import importlib.util
import json
import re
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
DL = ROOT / "_sources" / "gse106878_reconstruction" / "downloads"
PLATFORM = DL / "GPL10295_self_full.txt"
MATRIX = DL / "GSE106878_series_matrix.txt.gz"
SIG_PY = ROOT / "_sources" / "repo_subset_v1.3.1" / "code" / "portable_analysis" / "03_00_09_environment_lock" / "signatures.py"
PAIR = ROOT / "_work" / "analysis" / "gse106878" / "GSE106878_patient_paired_changes.csv"
OUT = ROOT / "_work" / "analysis"
ITERATIONS = 1000
ALIASES = {"KIAA1370": "FAM214A", "C9ORF103": "IDNK", "C9ORF95": "NMRK1", "NALP1": "NLRP1", "FCMR": "FAIM3"}


def stable_seed(*parts: str) -> int:
    digest = hashlib.sha256("|".join(parts).encode()).digest()
    return int.from_bytes(digest[:8], "big") % (2**32 - 1)


def load_signatures():
    spec = importlib.util.spec_from_file_location("frozen_signatures_bg", SIG_PY)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def canonical(text: str) -> str:
    symbol = str(text).strip().upper()
    return ALIASES.get(symbol, symbol)


def platform_universe() -> dict[str, list[str]]:
    mapping: dict[str, list[str]] = defaultdict(list)
    inside = False
    with PLATFORM.open("rt", errors="replace") as handle:
        for line in handle:
            if line.startswith("!platform_table_begin"):
                inside = True
                header = next(handle).rstrip("\r\n").split("\t")
                index = {name: pos for pos, name in enumerate(header)}
                continue
            if line.startswith("!platform_table_end"):
                break
            if not inside:
                continue
            row = line.rstrip("\r\n").split("\t")
            if len(row) <= max(index["ID"], index["Symbol"]):
                continue
            raw_symbols = [canonical(x) for x in re.split(r"\s*(?:///|;|\|)\s*", row[index["Symbol"]]) if x.strip()]
            symbols = sorted(set(s for s in raw_symbols if re.fullmatch(r"[A-Z0-9][A-Z0-9._-]*", s) and s not in {"NA", "N/A", "NULL"}))
            # Restrict the background to unambiguous platform annotations.
            if len(symbols) == 1:
                mapping[symbols[0]].append(row[index["ID"]].strip())
    return mapping


def read_matrix(probes: set[str]) -> tuple[pd.DataFrame, pd.DataFrame]:
    metadata: dict[str, list[str]] = {}
    header = None
    rows = []
    with gzip.open(MATRIX, "rt", errors="replace") as handle:
        for line in handle:
            if line.startswith("!Sample_title"):
                metadata["title"] = next(csv.reader([line], delimiter="\t"))[1:]
            elif line.startswith("!Sample_geo_accession"):
                metadata["sample_id"] = next(csv.reader([line], delimiter="\t"))[1:]
            elif line.startswith("!series_matrix_table_begin"):
                header = next(csv.reader([next(handle)], delimiter="\t"))
                for row in csv.reader(handle, delimiter="\t"):
                    if row and row[0].startswith("!series_matrix_table_end"):
                        break
                    if row and row[0] in probes:
                        rows.append(row)
                break
    if header is None:
        raise RuntimeError("Series matrix table not found")
    matrix = pd.DataFrame(rows, columns=header).set_index("ID_REF").astype(float)
    samples = pd.DataFrame(metadata)
    parsed = samples["title"].str.extract(r"^(?P<patient_id>C\d+)\s+(?P<timepoint>Pre|Post\(24h\))\s+(?P<treatment>placebo|hydrocortisone)$")
    samples = pd.concat([samples, parsed], axis=1)
    if samples[["patient_id", "timepoint", "treatment"]].isna().any().any():
        raise RuntimeError("Sample title parsing failed")
    return matrix, samples


def assign_bins(values: pd.Series) -> pd.Series:
    # Rank-first qcut produces ten balanced baseline-only bins even with ties.
    return pd.qcut(values.rank(method="average"), 10, labels=False).astype(int)


def main() -> None:
    signatures = load_signatures()
    platform = platform_universe()
    probe, samples = read_matrix({p for ps in platform.values() for p in ps})
    gene = pd.DataFrame(index=sorted(platform), columns=probe.columns, dtype=float)
    for symbol in gene.index:
        available = [p for p in platform[symbol] if p in probe.index]
        if available:
            gene.loc[symbol] = probe.loc[available].mean(axis=0)
    gene = gene.dropna(axis=0, how="any")
    gene = gene[(gene > 0).all(axis=1)]

    sample_index = samples.set_index("sample_id")
    baseline_samples = sample_index.index[sample_index["timepoint"] == "Pre"].tolist()
    followup_samples = sample_index.index[sample_index["timepoint"] == "Post(24h)"].tolist()
    if len(baseline_samples) != 47 or len(followup_samples) != 47:
        raise RuntimeError("Expected 47 baseline and 47 follow-up samples")
    paired_samples = samples.pivot(index="patient_id", columns="timepoint", values="sample_id")

    metrics = pd.DataFrame({
        "baseline_mean": gene[baseline_samples].mean(axis=1),
        "baseline_variance": gene[baseline_samples].var(axis=1, ddof=1),
    })
    metrics["mean_bin"] = assign_bins(metrics["baseline_mean"])
    metrics["variance_bin"] = assign_bins(metrics["baseline_variance"])
    all_required = {g for genes in signatures.STAGE3_REQUIRED_GENES.values() for g in genes}
    missing = sorted(all_required - set(metrics.index))
    if missing:
        raise RuntimeError("Required components absent from background universe: " + ", ".join(missing))
    candidate_metrics = metrics.drop(index=list(all_required), errors="ignore")

    observed = pd.read_csv(PAIR).groupby("signature_id")["delta_z"].mean().to_dict()
    null_rows, mapping_rows, summary_rows = [], [], []
    for signature_id, function in signatures.STAGE3_FUNCTIONS.items():
        components = list(signatures.STAGE3_REQUIRED_GENES[signature_id])
        direction = signatures.SCORE_DIRECTION[signature_id]
        rng = np.random.default_rng(stable_seed("formula-preserving", "GSE106878", signature_id))
        successes = 0
        for iteration in range(1, ITERATIONS + 1):
            replacements: dict[str, str] = {}
            used: set[str] = set()
            exact_matches = 0
            for component in components:
                target_mean = int(metrics.at[component, "mean_bin"])
                target_var = int(metrics.at[component, "variance_bin"])
                available = candidate_metrics.loc[~candidate_metrics.index.isin(used)].copy()
                distance = (available["mean_bin"] - target_mean).abs() + (available["variance_bin"] - target_var).abs()
                minimum = int(distance.min())
                pool = available.index[distance == minimum].to_numpy()
                replacement = str(rng.choice(pool))
                replacements[component] = replacement
                used.add(replacement)
                exact_matches += int(minimum == 0)
                mapping_rows.append({
                    "signature_id": signature_id, "iteration": iteration, "component_gene": component,
                    "replacement_gene": replacement, "mean_bin": target_mean, "variance_bin": target_var,
                    "bin_distance": minimum,
                })
            try:
                pseudo_scores = {}
                for sample_id in gene.columns:
                    expression = {component: float(gene.at[replacement, sample_id]) for component, replacement in replacements.items()}
                    pseudo_scores[sample_id] = direction * float(function(expression))
                baseline_values = np.array([pseudo_scores[x] for x in baseline_samples], dtype=float)
                baseline_sd = float(baseline_values.std(ddof=1))
                if not np.isfinite(baseline_sd) or baseline_sd <= 0:
                    raise ValueError("non-positive baseline SD")
                deltas = []
                for _, row in paired_samples.iterrows():
                    deltas.append((pseudo_scores[row["Post(24h)"]] - pseudo_scores[row["Pre"]]) / baseline_sd)
                effect = float(np.mean(deltas))
                if not np.isfinite(effect):
                    raise ValueError("non-finite pseudo effect")
                successes += 1
                null_rows.append({
                    "dataset": "GSE106878", "signature_id": signature_id, "iteration": iteration,
                    "pseudo_mean_delta_z": effect, "exact_bin_match_fraction": exact_matches / len(components),
                    "status": "SUCCESS",
                })
            except Exception as exc:
                null_rows.append({
                    "dataset": "GSE106878", "signature_id": signature_id, "iteration": iteration,
                    "pseudo_mean_delta_z": np.nan, "exact_bin_match_fraction": exact_matches / len(components),
                    "status": f"FAILED: {type(exc).__name__}",
                })

        subset = pd.DataFrame([r for r in null_rows if r["signature_id"] == signature_id and r["status"] == "SUCCESS"])
        obs = float(observed[signature_id])
        abs_tail = int((subset["pseudo_mean_delta_z"].abs() >= abs(obs)).sum())
        summary_rows.append({
            "dataset": "GSE106878", "signature_id": signature_id,
            "observed_mean_delta_z": obs, "successful_iterations": successes,
            "success_fraction": successes / ITERATIONS,
            "median_exact_bin_match_fraction": float(subset["exact_bin_match_fraction"].median()),
            "null_mean": float(subset["pseudo_mean_delta_z"].mean()),
            "null_sd": float(subset["pseudo_mean_delta_z"].std(ddof=1)),
            "observed_absolute_percentile": float((subset["pseudo_mean_delta_z"].abs() <= abs(obs)).mean()),
            "empirical_two_sided_p": (abs_tail + 1) / (successes + 1),
            "gate_status": "PASS" if successes >= 0.95 * ITERATIONS else "UNAVAILABLE",
        })

    null = pd.DataFrame(null_rows)
    maps = pd.DataFrame(mapping_rows)
    summary = pd.DataFrame(summary_rows)
    null.to_csv(OUT / "08_formula_preserving_null_distribution.csv", index=False)
    maps.to_csv(OUT / "08_formula_preserving_gene_matches.csv", index=False)
    summary.to_csv(OUT / "08_formula_preserving_background_summary.csv", index=False)
    provenance = {
        "dataset": "GSE106878", "background_gene_universe": int(len(gene)),
        "series_matrix_url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE106nnn/GSE106878/matrix/GSE106878_series_matrix.txt.gz",
        "platform_table_url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL10295&targ=self&form=text&view=full",
        "baseline_samples": len(baseline_samples), "followup_samples": len(followup_samples),
        "matching_variables": ["baseline_mean_decile", "baseline_variance_decile"],
        "followup_or_outcome_used_for_matching": False,
        "iterations_per_signature": ITERATIONS,
        "all_signature_gates_passed": bool((summary["gate_status"] == "PASS").all()),
    }
    (OUT / "08_formula_preserving_background_provenance.json").write_text(json.dumps(provenance, indent=2), encoding="utf-8")
    print(json.dumps(provenance, indent=2))


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Reconstruct GSE106878 fixed-signature scores from public GEO files."""

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
DOWNLOADS = ROOT / "_sources" / "gse106878_reconstruction" / "downloads"
MATRIX = DOWNLOADS / "GSE106878_series_matrix.txt.gz"
PLATFORM = DOWNLOADS / "GPL10295_self_full.txt"
SIGNATURES_PY = ROOT / "_sources" / "repo_subset_v1.3.1" / "code" / "portable_analysis" / "03_00_09_environment_lock" / "signatures.py"
OUT = ROOT / "_work" / "analysis" / "gse106878"
OUT.mkdir(parents=True, exist_ok=True)

ALIAS_TO_CURRENT = {"KIAA1370": "FAM214A", "C9ORF103": "IDNK", "C9ORF95": "NMRK1", "NALP1": "NLRP1", "FCMR": "FAIM3"}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def load_signatures():
    spec = importlib.util.spec_from_file_location("frozen_signatures", SIGNATURES_PY)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def parse_platform(required: set[str]) -> tuple[dict[str, list[str]], pd.DataFrame]:
    gene_to_probes: dict[str, list[str]] = defaultdict(list)
    evidence = []
    inside = False
    with PLATFORM.open("rt", errors="replace") as handle:
        for line in handle:
            if line.startswith("!platform_table_begin"):
                inside = True
                header = next(handle).rstrip("\n").split("\t")
                idx = {name: i for i, name in enumerate(header)}
                continue
            if line.startswith("!platform_table_end"):
                break
            if not inside:
                continue
            row = line.rstrip("\n").split("\t")
            if len(row) <= max(idx["ID"], idx["Symbol"]):
                continue
            probe = row[idx["ID"]].strip()
            for raw in re.split(r"\s*(?:///|;|\|)\s*", row[idx["Symbol"]].strip()):
                canonical = ALIAS_TO_CURRENT.get(raw.upper(), raw.upper()) if raw else ""
                if canonical in required:
                    gene_to_probes[canonical].append(probe)
                    evidence.append({"probe_id": probe, "platform_symbol": raw, "canonical_symbol": canonical, "mapping_basis": "historical alias" if canonical != raw.upper() else "exact symbol"})
    return gene_to_probes, pd.DataFrame(evidence)


def parse_matrix(selected_probes: set[str]) -> tuple[pd.DataFrame, pd.DataFrame]:
    metadata: dict[str, list[str]] = {}
    rows = []
    header = None
    with gzip.open(MATRIX, "rt", errors="replace") as handle:
        for line in handle:
            if line.startswith("!Sample_title"):
                metadata["title"] = next(csv.reader([line], delimiter="\t"))[1:]
            elif line.startswith("!Sample_geo_accession"):
                metadata["geo_accession"] = next(csv.reader([line], delimiter="\t"))[1:]
            elif line.startswith("!series_matrix_table_begin"):
                header = next(csv.reader([next(handle)], delimiter="\t"))
                for row in csv.reader(handle, delimiter="\t"):
                    if row and row[0].startswith("!series_matrix_table_end"):
                        break
                    if row and row[0] in selected_probes:
                        rows.append(row)
                break
    if header is None:
        raise RuntimeError("GEO matrix table was not found")
    matrix = pd.DataFrame(rows, columns=header).set_index("ID_REF").astype(float)
    sample = pd.DataFrame(metadata)
    parsed = sample["title"].str.extract(r"^(?P<patient>C\d+)\s+(?P<timepoint>Pre|Post\(24h\))\s+(?P<treatment>placebo|hydrocortisone)$")
    sample = pd.concat([sample, parsed], axis=1)
    if sample[["patient", "timepoint", "treatment"]].isna().any().any():
        raise RuntimeError("Sample-title parsing failed")
    return matrix, sample


def bootstrap_mean(values: np.ndarray, rng: np.random.Generator, b: int = 2000) -> tuple[float, float, float]:
    draws = values[rng.integers(0, len(values), size=(b, len(values)))].mean(axis=1)
    return float(draws.std(ddof=1)), float(np.quantile(draws, 0.025)), float(np.quantile(draws, 0.975))


def main() -> None:
    signatures = load_signatures()
    required = {gene for genes in signatures.STAGE3_REQUIRED_GENES.values() for gene in genes}
    gene_to_probes, mapping = parse_platform(required)
    missing_mapping = sorted(required - set(gene_to_probes))
    if missing_mapping:
        raise RuntimeError("Incomplete platform mapping: " + ", ".join(missing_mapping))
    probe_matrix, sample = parse_matrix({probe for probes in gene_to_probes.values() for probe in probes})
    missing_matrix = sorted(gene for gene, probes in gene_to_probes.items() if not any(probe in probe_matrix.index for probe in probes))
    if missing_matrix:
        raise RuntimeError("Mapped probes absent from matrix: " + ", ".join(missing_matrix))

    gene_matrix = pd.DataFrame(index=sorted(required), columns=probe_matrix.columns, dtype=float)
    for gene in gene_matrix.index:
        probes = [probe for probe in gene_to_probes[gene] if probe in probe_matrix.index]
        gene_matrix.loc[gene] = probe_matrix.loc[probes].mean(axis=0)
    if (gene_matrix <= 0).any().any():
        raise RuntimeError("Non-positive expression prevents fixed-formula calculation")

    score_rows = []
    score_lookup: dict[tuple[str, str], float] = {}
    for signature, function in signatures.STAGE3_FUNCTIONS.items():
        for gsm in gene_matrix.columns:
            expression = {gene: float(gene_matrix.at[gene, gsm]) for gene in signatures.STAGE3_REQUIRED_GENES[signature]}
            value = signatures.SCORE_DIRECTION[signature] * float(function(expression))
            score_lookup[(signature, gsm)] = value
            record = sample.loc[sample["geo_accession"].eq(gsm)].iloc[0]
            score_rows.append({"dataset": "GSE106878", "patient_id": record.patient, "sample_id": gsm, "timepoint": record.timepoint, "treatment": record.treatment, "signature_id": signature, "oriented_score": value})
    scores = pd.DataFrame(score_rows)

    pair_rows = []
    effect_rows = []
    rng = np.random.default_rng(20260720)
    for signature in signatures.STAGE3_FUNCTIONS:
        subset = scores[scores["signature_id"] == signature]
        baseline_sd = float(subset.loc[subset["timepoint"] == "Pre", "oriented_score"].std(ddof=1))
        # Preserve GEO sample order so the seeded bootstrap reproduces the
        # archived S26 intervals bit-for-bit.
        for patient, group in subset.groupby("patient_id", sort=False):
            pre = group[group["timepoint"] == "Pre"].iloc[0]
            post = group[group["timepoint"] == "Post(24h)"].iloc[0]
            pair_rows.append({"dataset": "GSE106878", "patient_id": patient, "signature_id": signature, "time_window": "T24", "treatment": pre.treatment, "baseline_sample_id": pre.sample_id, "followup_sample_id": post.sample_id, "baseline_score": pre.oriented_score, "followup_score": post.oriented_score, "baseline_sd": baseline_sd, "delta_z": (post.oriented_score - pre.oriented_score) / baseline_sd})
        pair_frame = pd.DataFrame([row for row in pair_rows if row["signature_id"] == signature])
        for stratum, frame in [("all", pair_frame), *[(name, group) for name, group in pair_frame.groupby("treatment", sort=True)]]:
            values = frame["delta_z"].to_numpy(float)
            se, lo, hi = bootstrap_mean(values, rng)
            effect_rows.append({"dataset": "GSE106878", "signature_id": signature, "time_window": "T24", "stratum": stratum, "n_pairs": len(values), "baseline_sd": baseline_sd, "mean_delta_z": float(values.mean()), "bootstrap_se": se, "bootstrap_ci_low": lo, "bootstrap_ci_high": hi, "bootstrap_replicates": 2000, "seed": 20260720})

    pairs = pd.DataFrame(pair_rows)
    effects = pd.DataFrame(effect_rows)
    mapping.sort_values(["canonical_symbol", "probe_id"]).to_csv(OUT / "GSE106878_platform_feature_map.csv", index=False)
    sample.to_csv(OUT / "GSE106878_sample_registry.csv", index=False)
    gene_matrix.to_csv(OUT / "GSE106878_signature_gene_expression.csv")
    scores.to_csv(OUT / "GSE106878_patient_signature_scores.csv", index=False)
    pairs.to_csv(OUT / "GSE106878_patient_paired_changes.csv", index=False)
    effects.to_csv(OUT / "GSE106878_directional_replication_effects.csv", index=False)
    provenance = {
        "series_matrix_url": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE106nnn/GSE106878/matrix/GSE106878_series_matrix.txt.gz",
        "platform_table_url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL10295&targ=self&form=text&view=full",
        "matrix_sha256": sha256(MATRIX),
        "platform_sha256": sha256(PLATFORM),
        "patients": int(pairs["patient_id"].nunique()),
        "paired_rows": len(pairs),
        "signatures": int(pairs["signature_id"].nunique()),
        "required_genes": len(required),
        "status": "PASS",
    }
    (OUT / "GSE106878_reconstruction_provenance.json").write_text(json.dumps(provenance, indent=2), encoding="utf-8")
    print(json.dumps(provenance, indent=2))


if __name__ == "__main__":
    main()

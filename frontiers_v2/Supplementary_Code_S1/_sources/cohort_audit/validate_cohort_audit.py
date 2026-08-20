#!/usr/bin/env python3
"""Executable consistency checks for the cohort-feasibility audit."""

from __future__ import annotations

import csv
from collections import Counter, defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parent


def rows(name: str, delimiter: str = ",") -> list[dict[str, str]]:
    with (ROOT / name).open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter=delimiter))


def by_patient(data: list[dict[str, str]]) -> dict[str, list[dict[str, str]]]:
    grouped: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in data:
        grouped[row["patient_id"]].append(row)
    return grouped


def main() -> None:
    messages: list[str] = []

    g54514 = rows("GSE54514_verified_sample_mapping.csv")
    baseline = [r for r in g54514 if r["relative_hours"] == "0" and r["clinical_status"] != "control"]
    assert len(baseline) == 35
    assert Counter(r["clinical_status"] for r in baseline) == Counter({"survivor": 26, "nonsurvivor": 9})
    orphan = [r for r in g54514 if r["mapping_status"] == "series_count_conflict_orphan_day4"]
    assert len(orphan) == 1 and orphan[0]["patient_id"] == "20" and orphan[0]["relative_hours"] == "72"
    messages.append("PASS GSE54514: 35 baseline patients (26/9); one Day-4-only orphan ID20 flagged")

    g95233 = [r for r in rows("GSE95233_verified_sample_mapping.csv") if r["clinical_status"] != "control"]
    p95233 = by_patient(g95233)
    assert len(p95233) == 51 and all(len(v) == 2 for v in p95233.values())
    assert Counter(v[0]["clinical_status"] for v in p95233.values()) == Counter({"survivor": 34, "non_survivor": 17})
    follow = [r for r in g95233 if r["relative_hours"] in {"24", "48"}]
    assert Counter(r["relative_hours"] for r in follow) == Counter({"24": 20, "48": 31})
    messages.append("PASS GSE95233: 51 complete pairs; 34/17 D28 outcome; 20 D2 and 31 D3 follow-ups")

    g110487 = rows("GSE110487_verified_sample_mapping.csv")
    p110487 = by_patient(g110487)
    assert len(p110487) == 31 and all(len(v) == 2 for v in p110487.values())
    assert all({r["relative_hours"] for r in v} == {"0", "48"} for v in p110487.values())
    assert Counter(v[0]["clinical_status"] for v in p110487.values()) == Counter({"R": 17, "NR": 14})
    assert Counter(r["mapping_status"] for r in g110487)["resolved_title_source_name_characteristic_conflict"] == 21
    assert all(len({r["platform_id"] for r in v}) == 1 for v in p110487.values())
    messages.append("PASS GSE110487: 31 complete pairs; 17/14 R/NR; 21 GEO characteristic conflicts isolated")

    g57065 = [r for r in rows("GSE57065_verified_sample_mapping.csv") if r["clinical_status"] != "control"]
    p57065 = by_patient(g57065)
    assert len(p57065) == 28
    assert Counter(v[0]["clinical_status"] for v in p57065.values()) == Counter({"SAPSII-High": 14, "SAPSII-Low": 14})
    assert sum("0" in {r["relative_hours"] for r in v} for v in p57065.values()) == 28
    assert sum("24" in {r["relative_hours"] for r in v} for v in p57065.values()) == 28
    assert sum("48" in {r["relative_hours"] for r in v} for v in p57065.values()) == 26
    messages.append("PASS GSE57065: 28 patient maps; 14/14 SAPSII class; H0/H24=28 and H48=26")

    g106878 = rows("GSE106878_verified_sample_mapping.csv")
    p106878 = by_patient(g106878)
    assert len(p106878) == 47 and all(len(v) == 2 for v in p106878.values())
    assert all({r["relative_hours"] for r in v} == {"0", "24"} for v in p106878.values())
    assert Counter(v[0]["treatment"] for v in p106878.values()) == Counter({"hydrocortisone": 24, "placebo": 23})
    assert Counter(v[0]["survival_(28_days)"] for v in p106878.values()) == Counter({"survived": 34, "died": 13})
    messages.append("PASS GSE106878: 47 complete pairs; 24/23 treatment; 34/13 D28 outcome")

    prj = rows("PRJEB111201_verified_sample_mapping.csv")
    patients = [r for r in prj if r["clinical_group"] != "healthy_control"]
    pprj = by_patient(patients)
    assert len(pprj) == 11
    assert sum({r["source_day"] for r in v} == {"1", "3", "7"} for v in pprj.values()) == 10
    assert {r["source_day"] for r in pprj["septic_shock_00"]} == {"1"}
    total_bytes = sum(int(r["fastq_total_bytes"]) for r in prj)
    assert total_bytes == 174_571_503_968
    analysis_lines = (ROOT / "PRJEB111201_ENA_analysis.tsv").read_text(encoding="utf-8").splitlines()
    assert len(analysis_lines) == 1
    messages.append("PASS PRJEB111201: 10 complete trajectories plus one D1-only patient; 162.58 GiB; no ENA analysis row")

    for accession in ("GSE54514", "GSE95233", "GSE110487", "GSE57065"):
        qc = rows(f"{accession}_existing_score_qc.csv")
        assert len(qc) == 8
        assert all(r["status"] == "PASS" and r["score_missing_n"] == "0" for r in qc)
    g106_qc = rows("GSE106878_existing_eligibility_summary.csv")
    assert g106_qc[0]["required_genes"] == g106_qc[0]["covered_genes"] == "74"
    messages.append("PASS score reconstruction: all 8 signatures complete in four primary cohorts; GSE106878 covers 74/74 genes")

    (ROOT / "AUDIT_QA.txt").write_text("\n".join(messages) + "\n", encoding="utf-8")
    print("\n".join(messages))


if __name__ == "__main__":
    main()

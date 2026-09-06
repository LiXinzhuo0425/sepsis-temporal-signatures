#!/usr/bin/env python3
"""Extract sample-level metadata needed for the Longitudinal V2 feasibility audit.

This script is deliberately read-only with respect to the established analysis
workspace. It writes compact, reviewable CSV extracts beside this file.
"""

from __future__ import annotations

import csv
import gzip
import hashlib
import re
from collections import Counter
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
OUT = Path(__file__).resolve().parent
EXPR = ROOT / "longitudinal_stage3/inputs/expression"
GSE106878_MATRIX = (
    ROOT
    / "stage10_submission/10_14_decision_gate/source_evidence/"
    "GSE106878_series_matrix.txt.gz"
)


def read_series_header(path: Path) -> list[dict[str, str]]:
    """Read aligned !Sample_* fields from a GEO series matrix header."""
    fields: dict[str, list[list[str]]] = {}
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("!series_matrix_table_begin"):
                break
            if not line.startswith("!Sample_"):
                continue
            row = next(csv.reader([line.rstrip("\n")], delimiter="\t"))
            fields.setdefault(row[0], []).append(row[1:])

    accessions = fields["!Sample_geo_accession"][0]
    records: list[dict[str, str]] = []
    for idx, accession in enumerate(accessions):
        rec = {"sample_id": accession}
        for key in ("!Sample_title", "!Sample_source_name_ch1", "!Sample_platform_id"):
            if key in fields:
                rec[key.removeprefix("!Sample_")] = fields[key][0][idx]
        for characteristic_row in fields.get("!Sample_characteristics_ch1", []):
            value = characteristic_row[idx].strip()
            if not value:
                continue
            if ":" in value:
                k, v = value.split(":", 1)
                rec[k.strip().lower().replace(" ", "_")] = v.strip()
            else:
                rec.setdefault("unparsed_characteristic", value)
        records.append(rec)
    return records


def parse_soft_samples(path: Path) -> list[dict[str, str]]:
    """Read sample metadata blocks from a GEO family SOFT file."""
    records: list[dict[str, str]] = []
    current: dict[str, str] | None = None
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            line = line.rstrip("\n")
            if line.startswith("^SAMPLE = "):
                if current:
                    records.append(current)
                current = {"sample_id": line.split("=", 1)[1].strip()}
                continue
            if current is None or not line.startswith("!Sample_"):
                continue
            key, value = line.split("=", 1)
            key = key.removeprefix("!Sample_").strip()
            value = value.strip()
            if key == "characteristics_ch1" and ":" in value:
                ck, cv = value.split(":", 1)
                current[ck.strip().lower().replace(" ", "_")] = cv.strip()
            elif key in {"title", "source_name_ch1", "platform_id", "instrument_model"}:
                current[key] = value
        if current:
            records.append(current)
    return records


def write_csv(name: str, rows: list[dict[str, str]]) -> None:
    keys: list[str] = []
    for row in rows:
        for key in row:
            if key not in keys:
                keys.append(key)
    with (OUT / name).open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=keys)
        writer.writeheader()
        writer.writerows(rows)


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_provenance() -> None:
    sources = {
        "GSE54514_GEO_quick.txt": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE54514&targ=self&form=text&view=quick",
        "GSE95233_GEO_quick.txt": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE95233&targ=self&form=text&view=quick",
        "GSE110487_GEO_quick.txt": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE110487&targ=self&form=text&view=quick",
        "GSE57065_GEO_quick.txt": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE57065&targ=self&form=text&view=quick",
        "GSE106878_GEO_quick.txt": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE106878&targ=self&form=text&view=quick",
        "GSE110487_family.soft.gz": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE110nnn/GSE110487/soft/GSE110487_family.soft.gz",
        "PMC6249814_BioC.xml": "https://www.ncbi.nlm.nih.gov/research/bionlp/RESTful/pmcoa.cgi/BioC_xml/PMC6249814/unicode",
        "PMC4512996_BioC.xml": "https://www.ncbi.nlm.nih.gov/research/bionlp/RESTful/pmcoa.cgi/BioC_xml/PMC4512996/unicode",
        "PMC13335145_BioC.xml": "https://www.ncbi.nlm.nih.gov/research/bionlp/RESTful/pmcoa.cgi/BioC_xml/PMC13335145/unicode",
        "PMID23807251_PubMed.xml": "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pubmed&id=23807251&retmode=xml",
        "PRJEB111201_ENA_read_run.tsv": "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJEB111201&result=read_run",
        "PRJEB111201_ENA_analysis.tsv": "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJEB111201&result=analysis",
        "ENA_read_run_return_fields.tsv": "https://www.ebi.ac.uk/ena/portal/api/returnFields?result=read_run",
        "ENA_analysis_return_fields.tsv": "https://www.ebi.ac.uk/ena/portal/api/returnFields?result=analysis",
        "existing_longitudinal_sample_manifest.csv": str(
            ROOT
            / "stage10_submission/10_04_10_05_cohort_search_missingness_audit/"
            "00_input_locked/longitudinal_sample_manifest.csv"
        ),
        "GSE106878_existing_sample_registry.csv": str(
            ROOT
            / "stage10_submission/10_15_final_freeze_v1_1_1/qa/"
            "gse106878_portable_recheck/GSE106878_sample_registry.csv"
        ),
        "GSE106878_existing_eligibility_summary.csv": str(
            ROOT
            / "stage10_submission/10_15_final_freeze_v1_1_1/qa/"
            "gse106878_portable_recheck/GSE106878_eligibility_summary.csv"
        ),
        "GSE106878_existing_directional_replication_effects.csv": str(
            ROOT
            / "stage10_submission/10_17_figure6a_redesign_v1_2_1/release_work/"
            "srep-sepsis-temporal-signatures-v1.2.1/reproducibility_evidence/"
            "S26_GSE106878_directional_replication_effects.csv"
        ),
    }
    for accession in ("GSE54514", "GSE95233", "GSE110487", "GSE57065"):
        sources[f"{accession}_existing_score_qc.csv"] = str(
            ROOT / f"longitudinal_stage3/intermediate/score_qc_{accession}.csv"
        )
    for accession in (
        "GSE54514",
        "GSE95233",
        "GSE110487",
        "GSE57065",
        "GSE106878",
        "PRJEB111201",
    ):
        sources[f"{accession}_verified_sample_mapping.csv"] = (
            "generated by extract_verified_metadata.py from the primary deposit and locked local input"
        )

    rows = []
    for name, source in sources.items():
        path = OUT / name
        if not path.exists():
            continue
        rows.append(
            {
                "local_file": name,
                "sha256": file_sha256(path),
                "source": source,
                "audit_date": "2026-08-16",
            }
        )
    write_csv("PROVENANCE.csv", rows)


def gse54514() -> list[dict[str, str]]:
    records = read_series_header(EXPR / "GSE54514_series_matrix.txt.gz")
    pattern = re.compile(
        r"PAXgene whole blood, (Control|sepsis_nonsurvivor|sepsis_survivor), "
        r"Day_(\d+), ID=(\d+)"
    )
    for rec in records:
        match = pattern.fullmatch(rec["title"])
        if not match:
            rec["mapping_status"] = "unresolved_title"
            continue
        group, day, pid = match.groups()
        rec["patient_id"] = pid
        rec["source_day"] = day
        rec["relative_hours"] = str((int(day) - 1) * 24)
        rec["clinical_status"] = {
            "sepsis_survivor": "survivor",
            "sepsis_nonsurvivor": "nonsurvivor",
            "Control": "control",
        }[group]
        rec["mapping_status"] = "resolved"

    # GEO states 26 survivors and 9 nonsurvivors, but title-level parsing yields
    # a tenth nonsurvivor ID (20) represented only at Day 4. Preserve and flag it.
    for rec in records:
        if rec.get("clinical_status") == "nonsurvivor" and rec.get("patient_id") == "20":
            rec["mapping_status"] = "series_count_conflict_orphan_day4"
    return records


def gse95233() -> list[dict[str, str]]:
    records = read_series_header(EXPR / "GSE95233_series_matrix.txt.gz")
    title_pattern = re.compile(r"Blood-(?P<patient>(?:SV|NS|CS|PC)_\d+)_D(?P<day>\d{2})")
    for rec in records:
        match = title_pattern.fullmatch(rec.get("title", ""))
        rec["patient_id"] = match.group("patient") if match else ""
        timepoint = rec.get("time_point", "")
        rec["source_day"] = timepoint
        rec["relative_hours"] = {"D01": "0", "D02": "24", "D03": "48"}.get(
            timepoint, ""
        )
        survival = rec.get("survival", "")
        rec["clinical_status"] = survival.lower().replace(" ", "_") if survival != "NA" else "control"
        rec["mapping_status"] = "resolved" if rec["patient_id"] and timepoint else "unresolved"
    return records


def gse57065() -> list[dict[str, str]]:
    records = read_series_header(EXPR / "GSE57065_series_matrix.txt.gz")
    patient_pattern = re.compile(r"Blood_(P\d+)_H(00|24|48)")
    control_pattern = re.compile(r"Blood_(HV\d+)")
    for rec in records:
        title = rec.get("title", "")
        if match := patient_pattern.fullmatch(title):
            rec["patient_id"], rec["source_time"] = match.groups()
            rec["relative_hours"] = str(int(rec["source_time"]))
            rec["clinical_status"] = rec.get("sapsii", "")
            rec["mapping_status"] = "resolved"
        elif match := control_pattern.fullmatch(title):
            rec["patient_id"] = match.group(1)
            rec["clinical_status"] = "control"
            rec["mapping_status"] = "control"
        else:
            rec["mapping_status"] = "unresolved_title"
    return records


def gse110487() -> list[dict[str, str]]:
    records = parse_soft_samples(OUT / "GSE110487_family.soft.gz")
    for rec in records:
        rec["patient_id"] = rec.get("patient", "")
        rec["deposited_timepoint_characteristic"] = rec.get("timepoint", "")
        title_match = re.search(r"(T[12])$", rec.get("title", ""))
        source_match = re.search(r"\b(T[12])\b", rec.get("source_name_ch1", ""))
        title_time = title_match.group(1) if title_match else ""
        source_time = source_match.group(1) if source_match else ""
        # The raw-count column labels and sample source names agree for all 62
        # records. The separate GEO timepoint characteristic is misassigned for
        # a subset, so title/source-name concordance is the recoverable mapping.
        rec["source_time"] = title_time if title_time == source_time else ""
        rec["relative_hours"] = {"T1": "0", "T2": "48"}.get(
            rec["source_time"], ""
        )
        rec["clinical_status"] = rec.get("clinical_classification", "")
        required = rec["patient_id"] and rec["source_time"] and rec["clinical_status"]
        if not required:
            rec["mapping_status"] = "unresolved"
        elif rec["deposited_timepoint_characteristic"] != rec["source_time"]:
            rec["mapping_status"] = "resolved_title_source_name_characteristic_conflict"
        else:
            rec["mapping_status"] = "resolved"
    return records


def gse106878() -> list[dict[str, str]]:
    records = read_series_header(GSE106878_MATRIX)
    pattern = re.compile(
        r"(?P<patient>[A-Za-z0-9]+) (?P<time>Pre|Post\(24h\)) "
        r"(?P<treatment>placebo|hydrocortisone)"
    )
    for rec in records:
        match = pattern.fullmatch(rec.get("title", ""))
        if not match:
            rec["mapping_status"] = "unresolved_title"
            continue
        rec.update(match.groupdict())
        rec["patient_id"] = rec.pop("patient")
        rec["source_time"] = rec.pop("time")
        rec["relative_hours"] = "0" if rec["source_time"] == "Pre" else "24"
        rec["mapping_status"] = "resolved"
    return records


def prjeb111201() -> list[dict[str, str]]:
    path = OUT / "PRJEB111201_ENA_read_run.tsv"
    with path.open(newline="", encoding="utf-8") as handle:
        records = list(csv.DictReader(handle, delimiter="\t"))
    patient_pattern = re.compile(
        r"(?P<group>Sepsis|Septic shock) patient (?P<number>\d+) "
        r"hospital day (?P<day>\d+) - whole blood RNA-seq"
    )
    control_pattern = re.compile(r"Healthy control (?P<number>\d+) - whole blood RNA-seq")
    for rec in records:
        title = rec.get("sample_title", "")
        if match := patient_pattern.fullmatch(title):
            group = match.group("group").lower().replace(" ", "_")
            rec["clinical_group"] = group
            # Keep 00 distinct from 0: the paper states that one of six shock
            # patients contributed only a day-1 sample, matching this record.
            rec["patient_id"] = f"{group}_{match.group('number')}"
            rec["source_day"] = match.group("day")
            rec["mapping_status"] = "resolved"
        elif match := control_pattern.fullmatch(title):
            rec["clinical_group"] = "healthy_control"
            rec["patient_id"] = f"healthy_control_{match.group('number')}"
            rec["source_day"] = ""
            rec["mapping_status"] = "resolved_control"
        else:
            rec["mapping_status"] = "unresolved_title"
        byte_parts = [int(value) for value in rec.get("fastq_bytes", "").split(";") if value]
        rec["fastq_total_bytes"] = str(sum(byte_parts))
        rec["fastq_total_gib"] = f"{sum(byte_parts) / 1024**3:.6f}"
    return records


def main() -> None:
    datasets = {
        "GSE54514": gse54514(),
        "GSE95233": gse95233(),
        "GSE110487": gse110487(),
        "GSE57065": gse57065(),
        "GSE106878": gse106878(),
        "PRJEB111201": prjeb111201(),
    }
    for accession, rows in datasets.items():
        write_csv(f"{accession}_verified_sample_mapping.csv", rows)
        print(accession, len(rows), Counter(r.get("mapping_status", "") for r in rows))
    write_provenance()


if __name__ == "__main__":
    main()

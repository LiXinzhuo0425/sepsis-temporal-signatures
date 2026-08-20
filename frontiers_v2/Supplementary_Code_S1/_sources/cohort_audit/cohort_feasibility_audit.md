# Frontiers V2 cohort feasibility audit

Audit date: 2026-08-16. This document separates what is directly recoverable from primary deposits from what remains conditional. No FASTQ files were downloaded.

## Decision summary

| Cohort | Decision | Defensible V2 role | Main condition |
|---|---|---|---|
| GSE95233 | **GREEN** | Primary 28-day mortality anchoring | Preserve D2 versus D3 follow-up timing |
| GSE110487 | **GREEN** | Primary SOFA-response anchoring | Use title/source-name time mapping, not the erroneous characteristic field |
| GSE57065 | **GREEN** | Primary binary SAPSII severity anchoring | Use High/Low class only; no patient-level mortality claim |
| GSE106878 | **GREEN** | Independent directional replication; secondary mortality anchoring | Account for hydrocortisone/placebo assignment |
| GSE54514 | **YELLOW** | Secondary source-defined survival trajectory | Endpoint horizon and one orphan patient ID remain unresolved |
| PRJEB111201 | **YELLOW** | Conditional new D1-to-D3/D7 replication | Public archive contains raw FASTQ only, about 162.58 GiB |

The three planned clinical-anchor families are therefore feasible without new wet-lab work. They are not interchangeable endpoints. GSE95233 contributes explicit 28-day mortality, GSE110487 contributes early organ-function response, and GSE57065 contributes baseline SAPSII severity class. Each should retain a separate multiplicity family and interpretation.

## GSE54514

GEO describes daily whole-blood sampling for up to five days in 26 survivors, 9 nonsurvivors, and 18 controls. The public series matrix contains 163 samples. Thirty-five sepsis samples occur at baseline, followed by 31, 28, 19, and 14 at 24, 48, 72, and 96 hours. Patient and day labels are encoded in every sample title, and the existing fixed-formula workflow generated all eight scores without missing values for 163/163 samples.

Two constraints prevent a green clinical-anchor decision. First, neither the GEO record nor the accessible PubMed abstract specifies the survival horizon. The endpoint must therefore be called **source-defined survival status**, rather than 28-day mortality. Second, title parsing yields a Day-4-only nonsurvivor labelled ID 20, although GEO states that there were only nine nonsurvivors and 35 sepsis patients. The nine baseline nonsurvivors are IDs 19 and 21 to 28. ID 20 should remain excluded from baseline-defined clinical models and logged as an unresolved orphan record.

Later follow-up is strongly depleted, especially among nonsurvivors. A Day-1-to-Day-5 mixed model can be reported as secondary, with the number at risk shown at every day. A paired early-landmark analysis is safer for the main text.

## GSE95233

This is the cleanest mortality anchor. GEO explicitly defines day-28 survival, and the public matrix contains 51 complete septic-shock pairs plus 22 controls. Among patients, 34 survived and 17 died by day 28. Twenty second samples were collected on D2 and 31 on D3. Patient IDs, exact day, and survival status are recoverable without ambiguity. All eight fixed signatures were already reconstructed for 124/124 samples.

The D2/D3 mixture should be handled analytically. A robust two-visit formulation is paired score change regressed on day-28 status and exact follow-up day, with D2- and D3-stratified estimates as sensitivity analyses. Calling every second measurement “T48” would overstate the timing precision.

## GSE110487

GEO and the primary full text report 31 septic-shock patients, each sampled twice. T1 occurred within 16 hours of ICU admission and T2 occurred 48 hours after study enrolment. The deposited classifications comprise 17 responders and 14 nonresponders. The exact rule is important: a patient was a nonresponder when **SOFA(T1) − SOFA(T2) was below 5 and SOFA(T2) was above 8**; all other patients were classified as responders.

The public `GSE110487_rawcounts.xlsx` contains 58,096 Ensembl rows plus a header and all 62 sample columns. The deposit describes GRCh38, GENCODE v25, and featureCounts. Fifty-six samples used HiSeq 2500 and six used Genome Analyzer IIx, with each patient's two samples remaining on the same platform. Existing processing reconstructed all eight fixed signatures in every sample.

A metadata defect must be locked before new analysis. For 21 of 62 GEO sample records, the standalone `timepoint` characteristic conflicts with both the sample title and source name. The raw-count column names also agree with the title. The defensible authority is therefore the concordant title, source name, and raw-count column label. The conflict list is preserved in `GSE110487_verified_sample_mapping.csv`.

## GSE57065

The dataset contains 28 septic-shock patients and 25 controls. All patients have H0 and H24 samples; 26 have H48 samples. GEO titles provide patient IDs and exact hours. The patient-level SAPSII class is also present in the series matrix for every septic-shock sample: 14 patients are SAPSII-Low and 14 are SAPSII-High. The primary paper defines the classes using the cohort median, with scores below 45 and above 45, respectively.

The public deposit does not expose exact patient-level SAPSII values. The paper reports five 28-day deaths in aggregate, one in the Low group and four in the High group, but GEO does not map those deaths to patients. This cohort is green for a binary severity-trajectory interaction and unsuitable for patient-level mortality anchoring. All eight signatures were already reconstructed for 107/107 samples.

## GSE106878

The transcriptomic subset comprises 47 CORTICUS septic-shock patients with complete pre-treatment and 24-hour samples. GEO metadata resolve patient, time, randomized arm, 28-day survival, ACTH response, and baseline IFNG/IL10 stratum. There are 24 hydrocortisone and 23 placebo patients; 34 survived and 13 died by day 28.

This cohort is already stronger than a merely hypothetical replication. The existing audit recovered all 74 component genes, reconstructed all eight formulas, and produced treatment-stratified effects for all 47 pairs. The four diagnostic signatures reproduced the prespecified T24 direction. It should stay outside the original meta-analysis as an independent directional replication.

Because treatment may alter the 24-hour transcriptome, any added survival-trajectory model should include treatment assignment or report treatment-stratified sensitivity. The cohort can support more than one analysis, but it remains one independent patient sample and must not be counted twice when describing validation breadth.

## PRJEB111201

The primary paper enrolled 11 patients, five with sepsis and six with septic shock, plus four healthy controls. ENA exposes 35 paired-end RNA-seq runs. Ten patients have D1, D3, and D7 samples. A distinct record named `Septic shock patient 00` has only a D1 sample, matching the paper's statement that one shock patient died before later collection. Thus, D1-to-D3 and D1-to-D7 replication would contain ten complete patients, not eleven.

The public run metadata are complete, but the analysis is not plug-and-play. The 70 compressed FASTQ files total 174,571,503,968 bytes, or 162.58 GiB. An ENA `analysis` query returned no records. The paper states that HISAT2 and StringTie were used to obtain raw counts and TPM, while processed data and scripts are available only from the corresponding author on request. No public count matrix was identified.

The recommended first action is an author request for the StringTie gene-count or TPM matrix and its sample key. Full raw-data processing remains technically feasible, but it requires a separate compute and storage gate. Fixed-score reconstruction cannot be marked complete until gene coverage, identifier mapping, normalization, and the existing cross-platform formula rules pass QC.

## Immediate analysis gates

1. Freeze the verified sample maps in this folder before fitting any new models.
2. Use GSE95233 for explicit 28-day mortality, GSE110487 for the deposited response definition, and GSE57065 for the deposited SAPSII class. Keep these endpoint families separate.
3. Keep GSE54514 secondary and label its outcome as source-defined survival.
4. Retain GSE106878 as the independent directional replication. A mortality analysis is secondary and treatment-aware.
5. Contact the PRJEB111201 authors for processed counts before authorizing a 162.58-GiB raw download.

The machine-readable decisions are in `cohort_feasibility_audit.csv`; sample-level reconciliations are in the `*_verified_sample_mapping.csv` files.

## Primary-source links

- GEO: [GSE54514](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE54514), [GSE95233](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE95233), [GSE110487](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE110487), [GSE57065](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE57065), and [GSE106878](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE106878).
- Primary papers/full text: [PMID 23807251](https://pubmed.ncbi.nlm.nih.gov/23807251/), [PMC6249814](https://pmc.ncbi.nlm.nih.gov/articles/PMC6249814/), [PMC4512996](https://pmc.ncbi.nlm.nih.gov/articles/PMC4512996/), and [PMC13335145](https://pmc.ncbi.nlm.nih.gov/articles/PMC13335145/).
- ENA study: [PRJEB111201](https://www.ebi.ac.uk/ena/browser/view/PRJEB111201).

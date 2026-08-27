# isoform_CWOW

**Long-read transcriptome analysis across the Alzheimer’s and Lewy body disease neuropathological spectrum**

This repository contains the workflows and downstream analyses used to characterize full-length transcript isoforms in human anterior cingulate cortex using PacBio Kinnex long-read RNA sequencing.

The study includes **276 human brain samples** spanning control, Alzheimer’s disease (AD), and Lewy body disease (LBD) neuropathological groups. The overall workflow processes PacBio Kinnex reads into sample-level transcript annotations, constructs a cohort-wide transcript catalog, performs transcript classification and quality control, and generates transcript-level measurements for downstream differential transcript expression (DTE) and differential transcript usage (DTU) analyses.

---

## Workflow overview

```text
PacBio Revio Kinnex sequencing
        │
        ▼
Demultiplexing and read generation
        │
        ▼
Full-length non-concatemer (FLNC) reads
        │
        ▼
Iso-Seq clustering and alignment
        │
        ▼
Sample-level transcript models
        │
        ▼
SQANTI3 filtering
        │
        ▼
276 sample-level filtered GTFs
        │
        ▼
IsoMatch merge
        │
        ▼
7,770,780 cohort transcript models
        │
        ▼
IsoMatch classify
        │
        ▼
Sample-support QC
        │
        ▼
Cohort confidence filtering
        │
        ▼
Final transcript catalog
        │
        ▼
Transcript × sample quantification
        │
        ├───────────────┐
        ▼               ▼
      DTE              DTU
        │               │
        └───────┬───────┘
                ▼
     Long-read / short-read
           comparison
```

**Current stage:** sample-support QC and cohort confidence filtering.

---

## Reference genome and annotation

PacBio provides a prebuilt human reference set for Iso-Seq/Pigeon workflows:

https://downloads.pacbcloud.com/public/dataset/MAS-Seq/REF-pigeon_ref_sets/Human_hg38_Gencode_v39/

For the CWOW analysis, alignment and transcript annotation use GRCh38 with GENCODE v38.

Sex-informed reference genomes are used during alignment:

```text
Male:   GRCh38_YPARsmasked_XY.fa
Female: GRCh38_Ymasked_XX.fa
```

Annotation:

```text
gencode.v38.annotation.gtf
```

This approach retains the appropriate sex chromosomes while minimizing ambiguous alignment within homologous X/Y regions.

---

# 1. PacBio read generation

Initial PacBio Revio output was processed by TGen's long-read sequencing workflow.

TGen workflow:

https://github.com/tgen/tgs_readgeneration?tab=readme-ov-file#pacbio-revio

The initial processing follows Steps 1–3 of the PacBio Iso-Seq workflow:

https://isoseq.how/clustering/cli-workflow.html

The resulting files provided for downstream analysis are **unmapped BAM (`u.bam`) files**.

Based on the BAM headers, preprocessing consisted of:

1. **`lima`** — demultiplexing at the barcoded-adapter level.
2. **`skera`** — deconcatenating Kinnex arrays into S-reads.
3. **`lima`** — primer removal and additional demultiplexing.
4. **`isoseq refine`** — poly(A)-tail removal and generation of full-length non-concatemer reads.

These processed reads serve as input to the cohort-level Iso-Seq workflow described below.

---

# 2. TGen Jetstream / GrandCanyon workflow

The primary long-read processing workflow is executed on TGen's HPC infrastructure using the Jetstream/GrandCanyon workflow.

## Step 1 — Create cohort configuration

Generate the cohort JSON containing sample information, sequencing inputs, and sex information.

```bash
python3 create_kinnex_json_with_sex.py \
  --base /tgen_labs/jfryer/cores/tgen \
  --cohort-csv CWOW_full_long_read_cohort_n276.csv \
  --name CWOW_Kinnex_LongRead_with_sex_info \
  --study CWOW \
  --account YOUR_SLURM_ACCOUNT \
  --email YOUR_EMAIL \
  --check-pbi
```

## Step 2 — Audit input BAMs

Confirm that expected BAM files are present and correctly associated with cohort samples.

```bash
python3 audit_kinnex_bams_pandas.py \
  --base /tgen_labs/jfryer/cores/tgen \
  --csv CWOW_sample_submission_with_replacement_samples.csv
```

## Step 3 — Confirm sex-informed configuration

Generate or verify the sex-informed cohort configuration.

```bash
python3 create_kinnex_json_with_sex.py \
  --base /tgen_labs/jfryer/cores/tgen \
  --cohort-csv CWOW_full_long_read_cohort_n276.csv \
  --name CWOW_Kinnex_LongRead_with_sex_info \
  --study CWOW \
  --account YOUR_SLURM_ACCOUNT \
  --email YOUR_EMAIL \
  --check-pbi
```

## Step 4 — Generate GrandCanyon submission commands

Generate sample-specific workflow JSON files and SLURM submission commands.

```bash
python3 create_all_cwow_jsons_npid_fixed_2bams.py \
  --cohort-csv CWOW_full_long_read_cohort_n276.csv \
  --fastq-to-json fastq_to_json.py \
  --template grandcanyon_template.json \
  --bam-base /tgen_labs/jfryer/cores/tgen \
  --output-dir output/jsons \
  --results-base grandcanyon/results \
  --male-reference /tgen_labs/jfryer/projects/references/human/GRCh38/GRCh38_YPARsmasked_XY.fa \
  --female-reference /tgen_labs/jfryer/projects/references/human/GRCh38/GRCh38_Ymasked_XX.fa \
  --gtf /tgen_labs/jfryer/projects/references/human/GRCh38/gencode.v38.annotation.gtf \
  --sj-root /tgen_labs/jfryer/kolney/LBD_CWOW/isoform_CWOW/STAR_SJ \
  --account YOUR_SLURM_ACCOUNT \
  --dry-run > submit_runs_second_half.sh
```

The workflow performs Iso-Seq processing, alignment, transcript collapsing, and SQANTI3-based transcript QC for each sample.

---

# 3. Cohort transcript catalog

Following sample-level processing, SQANTI3-filtered transcript annotations are combined to construct a unified CWOW transcript catalog.

## Collect SQANTI3-filtered GTF files

Define the GrandCanyon results directory:

```bash
RESULTS_BASE="/tgen_labs/jfryer/kolney/LBD_CWOW/isoform_CWOW/scripts/grandcanyon/results"
```

Identify all sample-level filtered transcript annotations:

```bash
find "$RESULTS_BASE" \
  -type f \
  -path "*/rna/isoforms/isoseq/*/*.filtered.gtf" \
  | sort \
  > CWOW_filtered_gtfs.txt
```

Verify the resulting file list:

```bash
wc -l CWOW_filtered_gtfs.txt
head CWOW_filtered_gtfs.txt
```

The expected input consists of **276 sample-level SQANTI3-filtered GTF files**.

---

# 4. IsoMatch cohort merge

[IsoMatch](https://github.com/zhengxinchang/isomatch) is used to merge sample-level transcript annotations into a unified cohort transcriptome.

Version used:

```text
IsoMatch v0.6.2
```

Download:

```bash
wget https://github.com/zhengxinchang/isomatch/releases/download/v0.6.2/isomatch-v0.6.2-linux-x86_64.tar.gz
```

## Merge sample transcriptomes

```bash
isomatch merge \
  --ref-fa "$REF_FA" \
  -o "$OUTDIR/CWOW_cohort" \
  $(cat "$GTF_LIST")
```

Example output directory:

```bash
OUTDIR="/tgen_labs/jfryer/kolney/LBD_CWOW/isoform_CWOW/isomatch_cohort"
```

Inspect merge metadata:

```bash
cat "$OUTDIR/CWOW_cohort.merged_info.json"
cat "$OUTDIR/CWOW_cohort.merged_params.json"
```

Inspect the merged transcript catalog:

```bash
zcat "$OUTDIR/CWOW_cohort.merged.gtf.gz" | head -20
```

Inspect transcript tracking and sample-support matrices:

```bash
zcat "$OUTDIR/CWOW_cohort.track.tsv.gz" | head
zcat "$OUTDIR/CWOW_cohort.present_absent.tsv.gz" | head
```

Count transcript models:

```bash
zcat "$OUTDIR/CWOW_cohort.merged.gtf.gz" \
  | awk '$3=="transcript"' \
  | wc -l
```

Initial cohort merge:

```text
7,770,783 transcript models
```

---

# 5. IsoMatch transcript classification

The merged cohort transcriptome is classified relative to the reference annotation.

Submit classification:

```bash
sbatch iso_class.sh
```

## Classification results

| Category | Transcripts | Percent |
|---|---:|---:|
| FSM | 1,453,870 | 18.71% |
| ISM | 2,916,921 | 37.54% |
| NIC | 822,683 | 10.59% |
| NNC | 886,488 | 11.41% |
| Genic intron | 815,752 | 10.50% |
| Genic | 300,676 | 3.87% |
| Antisense | 270,467 | 3.48% |
| Intergenic | 255,812 | 3.29% |
| Fusion | 47,916 | 0.62% |

Classification was successful for essentially the entire catalog.

```text
Classified models: 7,770,780
Canonical models:  7,767,421
Canonical rate:    ~99.96%
```

Approximately **1.74 million models (22.4%) are mono-exonic**.

The very high canonical splice-site rate indicates that noncanonical splice junctions are not the primary source of the large initial transcript catalog. Instead, the catalog requires additional filtering based on sample support, transcript structure, abundance, and reproducibility across the cohort.

---

# 6. Cohort confidence filtering

The raw IsoMatch merge intentionally represents a permissive union of transcript structures observed across all samples.

The next stage is therefore **sample-support and cohort-level confidence filtering**.

```text
7.77 million merged transcript models
                 │
                 ▼
       Sample-support filtering
                 │
                 ▼
       Structural confidence QC
                 │
                 ▼
       Abundance / read support
                 │
                 ▼
       Final CWOW transcriptome
```

The goal is to distinguish reproducibly observed transcript structures from low-support sample-specific models before downstream statistical testing.

Filtering criteria and thresholds will be documented here once finalized.

---

# 7. Transcript quantification and downstream analysis

Following cohort-level filtering, the final transcript catalog will be used to construct a:

```text
transcript × sample abundance matrix
```

Primary downstream analyses include:

### Differential transcript expression (DTE)

Tests whether the abundance of an individual transcript changes between neuropathological groups.

### Differential transcript usage (DTU)

Tests whether the **relative proportion of transcripts within a gene** changes between groups, capturing isoform switching and other changes in transcript composition.

### Differential gene expression (DGE)

Gene-level expression will be evaluated alongside transcript-level results to distinguish changes in total gene abundance from changes specifically affecting transcript composition.

Of particular interest are genes showing **DTU without corresponding DGE**, which may indicate disease-associated transcript remodeling that would not be detected using conventional gene-level RNA-seq analysis.

---

# 8. Long-read vs. short-read comparison

Long-read transcript results will be compared with accompanying short-read RNA-seq data from the same cohort.

This enables evaluation of:

- transcript structures supported by both technologies;
- isoform-level changes detectable only with full-length sequencing;
- agreement between gene-level expression estimates;
- disease-associated transcript usage;
- novel transcript and splice-junction support; and
- potential cryptic exon or disease-associated splicing events.

---

# Repository structure

```text
isoform_CWOW/
├── scripts/
│   ├── grandcanyon/
│   ├── isomatch/
│   ├── quantification/
│   └── downstream_analysis/
│
├── metadata/
│   └── cohort metadata and sample manifests
│
├── output/
│   └── generated workflow configuration files
│
├── results/
│   ├── transcript_catalog/
│   ├── DGE/
│   ├── DTE/
│   └── DTU/
│
└── README.md
```

> **Note:** Large sequencing files and intermediate HPC outputs are not stored directly in this Git repository. Paths shown in this README refer to TGen HPC storage and may not be accessible outside the TGen computing environment.

---

# Key software and resources

| Resource | Purpose |
|---|---|
| PacBio Iso-Seq | Full-length transcript processing |
| Kinnex / MAS-Seq | High-throughput long-read RNA sequencing |
| lima | Barcode demultiplexing and primer removal |
| skera | Kinnex array deconcatenation |
| Iso-Seq refine | Poly(A) trimming / FLNC generation |
| SQANTI3 | Transcript annotation and quality control |
| IsoMatch v0.6.2 | Cohort transcript merging and classification |
| GRCh38 | Human reference genome |
| GENCODE v38 | Primary transcript annotation |

---

# Analysis status

```text
[✓] PacBio Revio sequencing
[✓] Kinnex read generation
[✓] Sample demultiplexing
[✓] Iso-Seq processing
[✓] Alignment
[✓] Sample-level transcript collapsing
[✓] SQANTI3 filtering
[✓] Collection of 276 filtered transcriptomes
[✓] IsoMatch cohort merge
[✓] IsoMatch transcript classification
[ ] Sample-support QC
[ ] Cohort confidence filtering
[ ] Final transcript catalog
[ ] Cohort transcript quantification
[ ] Differential gene expression (DGE)
[ ] Differential transcript expression (DTE)
[ ] Differential transcript usage (DTU)
[ ] Long-read / short-read integration
```

---

## Project summary

**Cohort:** 276 human anterior cingulate cortex samples  
**Technology:** PacBio Revio Kinnex long-read RNA sequencing  
**Reference:** GRCh38 / GENCODE v38  
**Initial cohort catalog:** ~7.77 million transcript models  
**Primary objective:** identify disease-associated changes in full-length transcript expression and usage across the Alzheimer’s and Lewy body disease neuropathological spectrum.
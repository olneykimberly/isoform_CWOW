

# Human reference 
https://downloads.pacbcloud.com/public/dataset/MAS-Seq/REF-pigeon_ref_sets/Human_hg38_Gencode_v39/

# isoform_CWOW
Long and short read isoform differences among human Lewy body dementia brains

Output form TGen:
The results are in unmapped BAM file format (u.bam) and the workflow used to process the results from the sequencer is located here https://github.com/tgen/tgs_readgeneration?tab=readme-ov-file#pacbio-revio. Steps 1-3 were followed from PacBio's workflow https://isoseq.how/clustering/cli-workflow.html.

From the samtools header\
-lima for demultiplexing at the barcoded adapters level\
-skera for deconcatenating the Kinnex arrays into S-reads\
-lima for primer removal and demultiplexing\
-isoseq refine for removal of the poly-A tails\

# TGen's jetstream GrandCanyon workflow
### Step 1. create the .json
python3 create_kinnex_json_with_sex.py \
  --base /tgen_labs/jfryer/cores/tgen \
  --cohort-csv CWOW_full_long_read_cohort_n276.csv \
  --name CWOW_Kinnex_LongRead_with_sex_info \
  --study CWOW \
  --account YOUR_SLURM_ACCOUNT \
  --email kolney@tgen.org \
  --check-pbi

### Step 2. inspect
python3 audit_kinnex_bams_pandas.py \
  --base /tgen_labs/jfryer/cores/tgen \
  --csv CWOW_sample_submission_with_replacement_samples.csv

### Step 3. sex informed reference 
python3 create_kinnex_json_with_sex.py \
  --base /tgen_labs/jfryer/cores/tgen \
  --cohort-csv CWOW_full_long_read_cohort_n276.csv \
  --name CWOW_Kinnex_LongRead_with_sex_info \
  --study CWOW \
  --account YOUR_SLURM_ACCOUNT \
  --email kolney@tgen.org \
  --check-pbi


### Step 4. generate submit commands
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
  --account tgen-206000 \
  --dry-run > submit_runs_second_half.sh


# Down stream anaysis   
### Create a table of individual SQANTI3 outputs into a single text file 
RESULTS_BASE="/tgen_labs/jfryer/kolney/LBD_CWOW/isoform_CWOW/scripts/grandcanyon/results"

find "$RESULTS_BASE" \
    -type f \
    -path "*/rna/isoforms/isoseq/*/*.filtered.gtf" \
    | sort \
    > CWOW_filtered_gtfs.txt

wc -l CWOW_filtered_gtfs.txt
head CWOW_filtered_gtfs.txt

### IsoSeq Match to merge individual files to generate a cohort matrix expression 
wget https://github.com/zhengxinchang/isomatch/releases/download/v0.6.2/isomatch-v0.6.2-linux-x86_64.tar.gz 


isomatch merge \
    --ref-fa "$REF_FA" \
    -o "$OUTDIR/CWOW_cohort" \
    $(cat "$GTF_LIST")

OUTDIR="/tgen_labs/jfryer/kolney/LBD_CWOW/isoform_CWOW/isomatch_cohort"
cat "$OUTDIR/CWOW_cohort.merged_info.json"
cat "$OUTDIR/CWOW_cohort.merged_params.json"
zcat "$OUTDIR/CWOW_cohort.merged.gtf.gz" | head -20
zcat "$OUTDIR/CWOW_cohort.track.tsv.gz" | head
zcat "$OUTDIR/CWOW_cohort.present_absent.tsv.gz" | head
zcat "$OUTDIR/CWOW_cohort.merged.gtf.gz" \
    | awk '$3=="transcript"' \
    | wc -l
7,770,783

### IsoSeq classify 
sbatch iso_class.sh
#--- class complete
Category	N	%
FSM	1,453,870	18.71%
ISM	2,916,921	37.54%
NIC	822,683	10.59%
NNC	886,488	11.41%
Genic intron	815,752	10.50%
Genic	300,676	3.87%
Antisense	270,467	3.48%
Intergenic	255,812	3.29%
Fusion	47,916	0.62%

There are also 1.74 million mono-exon models (22.4%).
The good news is that essentially all models classified successfully, and 7,767,421 / 7,770,780 ≈ 99.96% are canonical, so splice-site canonicality itself is not the major issue.

276 sample-level SQANTI3-filtered GTFs
                 ↓
          IsoMatch merge
                 ✓
        7,770,780 models
                 ↓
        IsoMatch classify
                 ✓
                 ↓
       SAMPLE-SUPPORT QC       ← WE ARE HERE
                 ↓
     cohort confidence filtering
                 ↓
      final transcript catalog
                 ↓
      transcript × sample counts
                 ↓
         DTE / DTU analyses
                 ↓
 long-read vs short-read comparison


# Human reference 
https://downloads.pacbcloud.com/public/dataset/MAS-Seq/REF-pigeon_ref_sets/Human_hg38_Gencode_v39/

# isoform_CWOW
Long and short read isoform differences among human Lewy body dementia brains

Output form TGen:
The results are in unmapped BAM file format (u.bam) and the workflow used to process the results from the sequencer is located here https://github.com/tgen/tgs_readgeneration?tab=readme-ov-file#pacbio-revio. Steps 1-3 were followed from PacBio's workflow https://isoseq.how/clustering/cli-workflow.html.



# Samtools header 
-lima for demultiplexing at the barcoded adapters level
-skera for deconcatenating the Kinnex arrays into S-reads
-lima for primer removal and demultiplexing
-isoseq refine for removal of the poly-A tails

#-------
conda activate isoseq 


ls /tgen_labs/jfryer/cores/tgen/r84132_20250721_215156/1_D01/CWOW_NA07150/CWOW_NA07150_1_BR_LBD_C1_PKFLR_A38419_SMRT3337_ACAGTC_L004.hifi.u.bam /tgen_labs/jfryer/cores/tgen/r84132_20250726_004614/1_C01/CWOW_NA07150/CWOW_NA07150_1_BR_LBD_C1_PKFLR_A38419_SMRT3340_ACAGTC_L003.hifi.u.bam > CWOW_NA07150.fofn


#-------- TGen
1. lima for demultiplexing at the barcoded adapters level
2. skera for deconcatenating the Kinnex arrays into S-reads
3. lima for primer removal and demultiplexing
4. isoseq refine for removal of the poly-A tails

#-------- Me
merged the SMRT cells
isoseq cluster2
aligned via pbmm2.


Prepare reference files for pigeon


isoseq collapse --do-not-collapse-extra-5exons NA15-031_CONTROL_mapped.bam NA15-031_CONTROL_collapsed.gff
pigeon sort NA15-031_CONTROL_collapsed.gff
pigeon classify NA15-031_CONTROL_collapsed.gff.sorted ../GRCh38/gencode.v38.annotation.sorted.gtf ../GRCh38/GRCh38.primary_assembly.genome.fa




#-------- SQANTI3
wget https://github.com/ConesaLab/SQANTI3/releases/download/v5.5/SQANTI3_v5.5.zip
mkdir sqanti3
unzip SQANTI3_v5.5.zip -d sqanti3

The SQANTI3 tool is designed to enable quality control and filtering of long read-defined transcriptomes, which are often rich in artifacts and false-positive isoforms. Therefore, a good curation of the transcriptome is indispensable to proceed with FIT analysis and produce valid, biologically sound conclusions/hypothesis.

the files head NA18-285_LBD_collapsed_sorted_classification.filtered_lite_classification.txt contain the counts 
After running the pigeon filter pipeline step, the file that contains the counts for each isoform is named with the extension _filtered.classification.txt.

This file contains isoform annotation information, which includes columns for full-length (FL) read counts. Specifically, the columns FL.<Sample> and fl_assoc provide the full-length read counts for each isoform.


#------- FLAIR https://flair.readthedocs.io/en/latest/modules.html
conda activate flair

Run flair correct and collapse individually on each sample, then combine your transcriptomes using flair combine. This method will be faster and easier, but you may miss some low-expression transcripts.
Your other option is to merge your samples before running FLAIR. If using PacBio reads, be careful doing this, as the reads may not have unique IDs. You may need to label each read with its sample ID to keep the read IDs unique. You can either merge the FASTA/FASTQ files before running FLAIR (simplest, recommended).

Modules must be run in order (align, correct, collapse).

** If you want to compare multiple samples, there are two primary ways of doing this:
Run FLAIR align, correct, and collapse (or FLAIR transcriptome) on each sample separately (better for large sets of samples)

#------ August readme

# create the .json
python3 create_kinnex_json_with_sex.py \
  --base /tgen_labs/jfryer/cores/tgen \
  --cohort-csv CWOW_full_long_read_cohort_n276.csv \
  --name CWOW_Kinnex_LongRead_with_sex_info \
  --study CWOW \
  --account YOUR_SLURM_ACCOUNT \
  --email kolney@tgen.org \
  --check-pbi

# add the grandcanyon tasks to the .json file 
python3 create_kinnex_json_with_sex.py \
  --base /tgen_labs/jfryer/cores/tgen \
  --cohort-csv CWOW_full_long_read_cohort_n276.csv \
  --name CWOW_Kinnex_LongRead_with_sex_info \
  --study CWOW \
  --account YOUR_SLURM_ACCOUNT \
  --email kolney@tgen.org \
  --template kinnex_template.json \
  --check-pbi


# audit bams
python3 audit_kinnex_bams_pandas.py \
  --base /tgen_labs/jfryer/cores/tgen \
  --csv CWOW_sample_submission_with_replacement_samples.csv

# add sex 
python3 create_kinnex_json_with_sex.py \
  --base /tgen_labs/jfryer/cores/tgen \
  --cohort-csv CWOW_full_long_read_cohort_n276.csv \
  --name CWOW_Kinnex_LongRead_with_sex_info \
  --study CWOW \
  --account YOUR_SLURM_ACCOUNT \
  --email kolney@tgen.org \
  --check-pbi


#---- updated
python3 create_all_cwow_jsons.py \
  --cohort-csv CWOW_full_long_read_cohort_n276.csv \
  --fastq-to-json fastq_to_json.py \
  --template grandcanyon_template \
  --bam-base /tgen_labs/jfryer/cores/tgen \
  --output-dir output/jsons \
  --results-base grandcanyon/results \
  --male-reference /path/to/male/reference.fa \
  --female-reference /path/to/female/reference.fa \
  --sj-root /path/to/short_read/results \
  --sj-pattern "{sample}/SJ.out.tab" \
  --allow-missing-sj \
  --overwrite \
  --dry-run

python3 create_all_cwow_jsons.py \
  --cohort-csv CWOW_full_long_read_cohort_n276.csv \
  --fastq-to-json fastq_to_json.py \
  --template grandcanyon_template.json \
  --bam-base /tgen_labs/jfryer/cores/tgen \
  --output-dir output/jsons \
  --results-base grandcanyon/results \
  --male-reference /tgen_labs/jfryer/projects/references/human/GRCh38/GRCh38_YPARsmasked_XY.fa \
  --female-reference /tgen_labs/jfryer/projects/references/human/GRCh38/GRCh38_Ymasked_XX.fa \
  --sj-root /tgen_labs/jfryer/kolney/LBD_CWOW/bulkRNA/starAligned_SCC 

python3 create_all_cwow_jsons_npid_fixed.py \
  --cohort-csv CWOW_full_long_read_cohort_n276.csv \
  --fastq-to-json fastq_to_json.py \
  --template grandcanyon_template.json \
  --bam-base /tgen_labs/jfryer/cores/tgen \
  --output-dir output/jsons \
  --results-base grandcanyon/results \
  --male-reference /tgen_labs/jfryer/projects/references/human/GRCh38/GRCh38_YPARsmasked_XY.fa \
  --female-reference /tgen_labs/jfryer/projects/references/human/GRCh38/GRCh38_Ymasked_XX.fa \
  --sj-root /tgen_labs/jfryer/kolney/LBD_CWOW/isoform_CWOW/STAR_SJ \
  --account tgen-206000 \
  --dry-run

#----
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

    
#--- merge 
RESULTS_BASE="/tgen_labs/jfryer/kolney/LBD_CWOW/isoform_CWOW/scripts/grandcanyon/results"

find "$RESULTS_BASE" \
    -type f \
    -path "*/rna/isoforms/isoseq/*/*.filtered.gtf" \
    | sort \
    > CWOW_filtered_gtfs.txt

wc -l CWOW_filtered_gtfs.txt
head CWOW_filtered_gtfs.txt

wget https://github.com/zhengxinchang/isomatch/releases/download/v0.6.2/isomatch-v0.6.2-linux-x86_64.tar.gz 

# see script
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
7770783

# classify 
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

 #-- 
# We need to know the exact names, but based on IsoMatch/SQANTI-style output, we can make this robustly in R.

#I recommend extracting the GTF attributes to a manageable TSV first:

zcat "$OUTDIR/CWOW_cohort.merged.gtf.gz" \
| awk -F'\t' '$3=="transcript" {
    match($9, /transcript_id "([^"]+)"/, a)
    match($9, /gene_id "([^"]+)"/, b)
    match($9, /ISOM_EXONS "([^"]+)"/, c)
    match($9, /ISOM_COUNT "([^"]+)"/, d)
    match($9, /ISOM_SAMPLE_CNT "([^"]+)"/, e)


    print a[1] "\t" b[1] "\t" c[1] "\t" d[1] "\t" e[1]
}' OFS='\t' \
> "$OUTDIR/CWOW_cohort_support.tsv"
sed -i '1i transcript_id\tisom_gene_id\tn_exons\tn_source_transcripts\tn_samples' \
    "$OUTDIR/CWOW_cohort_support.tsv"

head "$OUTDIR/CWOW_cohort_support.tsv"

wc -l "$OUTDIR/CWOW_cohort_support.tsv"
# should have 7,770,781
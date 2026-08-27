#!/bin/sh
#SBATCH --signal=USR2
#SBATCH --ntasks=1
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8
#SBATCH --time=04:00:00
#SBATCH --output=/tgen_labs/jfryer/kolney/LBD_CWOW/isoform_CWOW/scripts/iso_class.job.%j
#SBATCH --job-name=iso_merge # Job name
#SBATCH --mail-type=NONE     # Mail events (NONE, BEGIN, END, FAIL, ALL)


set -euo pipefail
cd /tgen_labs/jfryer/kolney/LBD_CWOW/isoform_CWOW/scripts/
OUTDIR="/tgen_labs/jfryer/kolney/LBD_CWOW/isoform_CWOW/isomatch_cohort"
REF_FA="/tgen_labs/jfryer/projects/references/human/GRCh38/GRCh38_YPARsmasked_XY.fa"
REF_GTF="/tgen_labs/jfryer/projects/references/human/GRCh38/gencode.v38.annotation.gtf"
ISOMATCH="/tgen_labs/jfryer/kolney/tools/isomatch"

echo "Started: $(date)"

"$ISOMATCH" classify \
    --ref-fa "$REF_FA" \
    --ref-gtf "$REF_GTF" \
    --out "$OUTDIR/CWOW_cohort_gencodev38" \
    "$OUTDIR/CWOW_cohort.merged.gtf.gz"

echo "Finished: $(date)"

echo
echo "Output files:"
ls -lh "$OUTDIR"/CWOW_cohort_gencodev38*
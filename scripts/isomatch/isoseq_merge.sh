#!/bin/sh
#SBATCH --signal=USR2
#SBATCH --ntasks=1
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8
#SBATCH --time=24:00:00
#SBATCH --output=/tgen_labs/jfryer/kolney/LBD_CWOW/QTL/LBD_CWOW_QTL/isoform_CWOW/scripts/iso_merge.job.%j
#SBATCH --job-name=iso_merge # Job name
#SBATCH --mail-type=NONE     # Mail events (NONE, BEGIN, END, FAIL, ALL)

cd /tgen_labs/jfryer/kolney/LBD_CWOW/isoform_CWOW/scripts/
GTF_LIST="CWOW_filtered_gtfs.txt"
REF_FA="/tgen_labs/jfryer/projects/references/human/GRCh38/GRCh38_YPARsmasked_XY.fa"
OUTDIR="/tgen_labs/jfryer/kolney/LBD_CWOW/isoform_CWOW/isomatch_cohort"
mkdir -p "$OUTDIR"

echo "========================================"
echo "CWOW IsoMatch cohort merge"
echo "Started: $(date)"
echo "========================================"

echo "IsoMatch version:"
/tgen_labs/jfryer/kolney/tools/isomatch --version

echo "Starting merge..."

/tgen_labs/jfryer/kolney/tools/isomatch merge \
    --ref-fa "$REF_FA" \
    -o "$OUTDIR/CWOW_cohort" \
    $(cat "$GTF_LIST")


echo "========================================"
echo "Merge completed: $(date)"
echo "========================================"

echo
echo "Output files:"
ls -lh "$OUTDIR"
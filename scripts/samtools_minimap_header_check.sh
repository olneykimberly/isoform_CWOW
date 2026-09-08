#!/bin/bash

BASE="/tgen_labs/jfryer/kolney/LBD_CWOW/isoform_CWOW/scripts/grandcanyon/results"

OUT="${BASE}/smrt_cell_header_check.tsv"

echo -e "sample\tn_read_groups\tread_groups\tstatus" > "$OUT"

find "$BASE" \
    -type f \
    -path "*/rna/alignment/minimap2/*.cram" \
    | sort \
    | while read -r cram; do

        relative="${cram#${BASE}/}"
        sample="${relative%%/*}"

        echo "Checking: $sample"

        rg_lines=$(samtools view -H "$cram" | grep '^@RG')

        n_rg=$(echo "$rg_lines" | grep -c '^@RG')

        rg_ids=$(echo "$rg_lines" \
            | tr '\t' '\n' \
            | grep '^ID:' \
            | sed 's/^ID://' \
            | paste -sd ',' -)

        if [ "$n_rg" -eq 2 ]; then
            status="PASS"
        else
            status="CHECK"
        fi

        echo -e "${sample}\t${n_rg}\t${rg_ids}\t${status}" >> "$OUT"

    done

echo
echo "Finished."
echo "Results:"
echo "$OUT"

echo
echo "Samples requiring inspection:"
awk -F'\t' 'NR==1 || $4 != "PASS"' "$OUT"


#!/usr/bin/env bash

output="out/mapping_percentages.txt"
echo -e "sample\ttotal_reads\tmapped\tpercent_mapped" > "$output"

for f in out/trimmed/*_cutadapt_report.txt; do
    # strip path + suffix
    fname=$(basename "$f")
    sample="${fname%_cutadapt_report.txt}"

    # --- 1) READS WRITTEN (line 13 of cutadapt report) ---
    line="$(sed -n '13p' "$f")"

    reads=$(
        echo "$line" |
        awk '{
            for (i = 1; i <= NF; i++) {
                gsub(/,/, "", $i)
                if ($i ~ /^[0-9]+$/) { print $i; break }
            }
        }'
    )

    # --- 2) MAPPED (sum over featureCounts summary) ---
    summary="out/counts/mito/longFeatures/${sample}_mito_featureCounts.txt.summary"

    if [[ -f "$summary" ]]; then
        mapped=$(awk 'NR>1 {s+=$2} END{print s+0}' "$summary")
    else
        mapped="NA"
    fi

    # --- 3) percent mapped ---
    if [[ "$mapped" != "NA" ]] && [[ "$reads" != "NA" ]]; then
        percent=$(awk -v m="$mapped" -v r="$reads" 'BEGIN {printf "%.2f", (m/r)*100}')
    else
        percent="NA"
    fi

    echo -e "${sample}\t${reads}\t${mapped}\t${percent}" >> "$output"
done

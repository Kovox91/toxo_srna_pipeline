#!/bin/bash

bedtools bamtobed -split -i out/filtered/par1_mito_filtered.bam \
| awk 'BEGIN{OFS="\t"}
{
  if ($6 == "+") {
    s = $2 + 1; e = $3
  } else {
    s = $3;     e = $2 + 1
  }
  start[$1 FS s FS $6]++
  end[$1 FS e FS $6]++
}
END {
  for (k in start) {
    split(k, a, FS)
    print a[1], a[2], start[k], a[3], "start"
  }
  for (k in end) {
    split(k, a, FS)
    print a[1], a[2], end[k], a[3], "end"
  }
}' > strand_start_end_counts.tsv

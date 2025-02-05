#!/bin/bash

# Set the directory containing BAM files
bam_dir="/Genomics/argo/users/ed7982/crvi-analysis/acute/dna/results/dmel-all-chromosome-r6.61.fasta/bams"

# Loop through each BAM file in the directory
for bam in "$bam_dir"/*.bam; do
    if [[ -f "$bam" ]]; then
        avg_depth=$(samtools depth -a "$bam" | awk '{sum+=$3; count++} END {if (count>0) print sum/count; else print "0"}')
        echo -e "$(basename "$bam")\t$avg_depth" >> avg_depths.txt
    fi
done

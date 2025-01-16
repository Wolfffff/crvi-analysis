module load samtools

dir_of_interest="/Genomics/argo/users/ed7982/crvi-analysis/acute/dna/results/dmel-all-chromosome-r6.61.fasta/bams"
output_dir="/Genomics/argo/users/ed7982/crvi-analysis/acute/dna/manual_samstats/stats"

# Create the output directory if it doesn't exist

# Iterate over BAM files
for file in $dir_of_interest/*_final.bam; do
    # Extract the base name of the file (without path)
    base_name=$(basename $file .bam)
    stats_file=$output_dir/${base_name}.stats

    # Submit a SLURM job
    sbatch --mem=8G --wrap="echo $base_name >> $stats_file && \
        samtools stats $file \
        | grep ^SN \
        | cut -f 2- \
        | sed 's/^[^:]*://; s/#.*$//; s/[[:space:]\t]//g' >> $stats_file"
done

#sq | grep "mem=8G" | awk '{print $1}' | xargs scancel
#find . -type f -empty -delete

echo "sample name" > header.txt
samtools stats /scratch/tmp/ed7982/shradda/illumina_snpcalling/06_recalibrated_bam/Plt3A3_Illumina.clean.sort.dedup.recal.bam \
    | grep ^SN \
    | cut -f 2- \
    | cut -d':' -f1 >> header.txt

cd stats

paste * > ../all_stats.txt
cd ..
paste header.txt all_stats.txt > all_stats_with_header.txt
mv all_stats_with_header.txt illumina_stats_with_header.txt

dir_of_interest="/scratch/tmp/ed7982/shradda/tn5_readcounting/06_recalibrated_bam"
output_dir="/scratch/tmp/ed7982/shradda/samtoolsstats/tn5_stats"

# Create the output directory if it doesn't exist

# Iterate over BAM files
for file in $dir_of_interest/*.bam; do
    # Extract the base name of the file (without path)
    base_name=$(basename $file .bam)
    stats_file=$output_dir/${base_name}.stats

    # Submit a SLURM job
    sbatch --mem=8G --wrap="echo $base_name >> $stats_file && \
        samtools stats $file \
        | grep ^SN \
        | cut -f 2- \
        | sed 's/^[^:]*://; s/#.*$//; s/[[:space:]\t]//g' >> $stats_file"
done

cd stats

paste * > ../tn5_all_stats.txt
paste header.txt tn5_all_stats.txt > tn5_stats_with_header.txt

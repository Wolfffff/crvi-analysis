module load R
for start in $(seq 1 5000 1119387); do
  sbatch --job-name=bb_$start \
         --output=logs/bb_$start.out \
         --error=logs/bb_$start.err \
         --cpus-per-task=1 \
         --mem=32G \
         --time=1-0 \
         --wrap="Rscript /Genomics/argo/users/ed7982/crvi-analysis/longevity_bp_resolution/bbin/beta_binomial.R $start"
done

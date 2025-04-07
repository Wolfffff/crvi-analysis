module load R
for start in $(seq 1 5000 1119387); do
  sbatch --job-name=bbp_$start \
         --output=perm_logs/bbp_$start.out \
         --error=perm_logs/bbp_$start.err \
         --cpus-per-task=1 \
         --mem=32G \
         --time=1-0 \
         --wrap="Rscript /Genomics/argo/users/ed7982/crvi-analysis/longevity_bp_resolution/bbin/beta_binomial_permuted.R $start"
done

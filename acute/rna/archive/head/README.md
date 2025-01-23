export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/mambaforge/23.3.1-1/lib/ 
snakemake --cores all --use-conda --conda-frontend conda -s ../rna-seq-star-deseq2/workflow/Snakefile -d ./

biomaRt::biomartCacheClear()
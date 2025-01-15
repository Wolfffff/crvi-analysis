# crvi-analysis

.
├── snpArcher
├── project_1/
│   ├── config/
│   │   ├── config.yaml
│   │   └── samples.csv
│   ├── data
│   └── results
└── project_2/
    ├── config/
    │   ├── config.yaml
    │   └── samples.csv
    └── data

#must delete slurm partition and account from slurm config file. also manually added slurm_logs

mamba activate snparcher
snakemake -s ./snpArcher/workflow/Snakefile -d ./longevity --cores all --use-conda --workflow-profile ./longevity/config/profiles/slurm

snakemake -s ./snpArcher/workflow/Snakefile -d ./acute/dna --cores all --use-conda --workflow-profile ./acute/dna/config/profiles/slurm


RNA
mamba create -c conda-forge -c bioconda --name snakemake snakemake snakedeploy
conda activate snakemake

cd  /Genomics/argo/users/ed7982/crvi-analysis/acute/rna/body
snakedeploy deploy-workflow https://github.com/snakemake-workflows/rna-seq-star-deseq2 . --tag v2.1.2

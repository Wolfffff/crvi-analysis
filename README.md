# crvi-analysis


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

snakemake -s ./snpArcher/workflow/Snakefile -d ./acute/dna --cores all --use-conda --workflow-profile ./acute/dna/config/profiles/slurm --rerun-incomplete --dry-run


--dry-run


RNA
mamba create -c conda-forge -c bioconda --name snakemake snakemake snakedeploy
conda activate snakemake

cd  /Genomics/argo/users/ed7982/crvi-analysis/acute/rna/body
snakedeploy deploy-workflow https://github.com/snakemake-workflows/rna-seq-star-deseq2 . --tag v2.1.2


snakemake -s ./acute/rna/body/workflow/Snakefile -d ./acute/rna/body --cores all --use-conda --rerun-incomplete 
snakemake --cores all --use-conda --rerun-incomplete 


cd /Genomics/argo/users/ed7982/crvi-analysis/acute/rna

head -n 3 rna_sample_sheet_head.tsv > ./head/rna_sample_sheet_head_top3.tsv
head -n 3 rna_unit_sheet_head.tsv > ./head/rna_unit_sheet_head_top3.tsv 



snakemake -s ./acute/rna/body/rna-seq-star-deseq2/workflow/Snakefile -d ./acute/rna/head --cores all --use-conda --workflow-profile ./acute/rna/head/config/profiles/slurm -F


<!-- conda install r-base r-ashr r-stringr rseqc r-tidyverse r-dbplyr -->
<!-- + biomart + DESeq2 -->

mamba create -c conda-forge -c bioconda --name stardeseq snakemake snakedeploy conda r-base r-ashr r-stringr rseqc r-tidyverse r-dbplyr bioconductor-deseq2 bioconductor-biomart=2.56


snakemake -s ./snpArcher/workflow/Snakefile -d ./longevity --cores all --use-conda --workflow-profile ./longevity/config/profiles/slurm --rerun-incomplete

bcftools +setGT acute_dna_raw_toedit.vcf.gz -- -t q -n . -e 'FMT/DP>=4' > acute_dna_raw_toedit_replaced.vcf


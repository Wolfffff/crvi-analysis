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

mamba activate snparcher
snakemake -s ./snpArcher/workflow/Snakefile -d ./longevity --cores all --use-conda --workflow-profile ./longevity/config/profiles/slurm --report report.html --dry-run

snakemake -s ./snpArcher/workflow/Snakefile -d ./acute/dna --cores all --use-conda --workflow-profile ./acute/dna/config/profiles/slurm --report report.html --dry-run



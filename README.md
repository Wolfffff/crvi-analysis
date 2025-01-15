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


snakemake -s ./snpArcher/workflow/Snakefile -d ./project_2 <other CLI options>


# PRIMED: Joint Loss Estimation with Selective Penalization for Multi-Reader Imaging Data with Binary Outcomes

## Overview

PRIMED (Penalized Regularization for Integrated Multi-path Estimation with Dispersion) is a joint loss framework for multi-reader imaging data that integrates three information channels: consensus, disagreement, and the structural relationship between them.

## Repository Structure

```
PRIMED_submission/
├── paper/                    Manuscript (LaTeX)
│   ├── primed_paper.tex        Main manuscript
│   └── primed_paper.pdf        Compiled PDF
├── code/                     Analysis code
│   ├── simulation/             Simulation studies (Julia + R)
│   │   ├── sim_modified_primed_core.jl       PRIMED optimizer (K <= 50)
│   │   ├── sim_modified_primed_core_highdim.jl  PRIMED optimizer (K = 100)
│   │   ├── sim_full_single_job.R             Single simulation job (R wrapper)
│   │   ├── run_full_parallel.sh              Parallel execution script
│   │   └── sim_all_with_sd.R                 Result aggregation
│   ├── rectal/                 Rectal cancer MRI application
│   │   ├── main.R                CD Group LASSO analysis
│   │   ├── apply_primed_real_data.R          PRIMED (linear dispersion)
│   │   ├── apply_quadratic_primed_real_data.R  PRIMED (quadratic dispersion)
│   │   ├── make_scatter_consensus_disagreement.R  Scatter plot figure
│   │   └── 01_rectal_bootstrap_stability.R   Bootstrap stability analysis
│   └── lidc/                   LIDC-IDRI lung nodule application
│       ├── 03_extract_radiomics_parallel.py  PyRadiomics feature extraction
│       ├── 04_apply_primed.R                 PRIMED application
│       ├── 06_cv10_evaluation_v2.R           10-fold nested CV evaluation
│       └── 02_lidc_bootstrap_stability.R     Bootstrap stability analysis
├── data/                     Input data
│   ├── data.xlsx               Rectal cancer multi-reader MRI data
│   └── lidc_annotations.csv    LIDC-IDRI lung nodule annotations
├── results/                  Analysis results
│   ├── simulation/             Simulation results (CSV)
│   └── lidc/                   LIDC-IDRI results
└── figures/                  Figures for manuscript
```

## Requirements

- **R** (>= 4.0) with packages: `tidyverse`, `readxl`, `grpreg`, `pROC`, `glmnet`, `parallel`
- **Julia** (>= 1.8) with standard libraries
- **Python** (>= 3.8) with `pyradiomics`, `SimpleITK`, `numpy`, `pandas` (for LIDC feature extraction only)

## Running the Code

### Simulation

```bash
cd code/simulation
bash run_full_parallel.sh        # 40 jobs, 8 parallel
Rscript sim_all_with_sd.R        # Aggregate results
```

### Rectal Cancer Application

```bash
cd code/rectal
Rscript apply_primed_real_data.R                  # Linear PRIMED
Rscript apply_quadratic_primed_real_data.R        # Linear vs Quadratic comparison
Rscript main.R                                     # CD Group LASSO baseline
```

### LIDC-IDRI Application

```bash
cd code/lidc
Rscript 06_cv10_evaluation_v2.R   # 10-fold nested CV (requires radiomics_cd.csv)
```

## Data Availability

- **LIDC-IDRI**: Publicly available at [The Cancer Imaging Archive (TCIA)](https://www.cancerimagingarchive.net/collection/lidc-idri/)
- **Rectal cancer MRI**: Available from the corresponding author on reasonable request

## License

This project is licensed under the MIT License.

## Citation

If you use this code, please cite:

> Kwon Y, Han K. PRIMED: Joint Loss Estimation with Selective Penalization for Multi-Reader Imaging Data with Binary Outcomes. *BMC Medical Research Methodology*. (submitted)

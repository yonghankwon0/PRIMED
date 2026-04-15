# PRIMED: Joint Loss Estimation with Selective Penalization for Multi-Reader Imaging Data with Binary Outcomes

## Overview

PRIMED (Penalized Regularization for Integrated Multi-path Estimation with Dispersion) is a joint loss framework for multi-reader imaging data that integrates three information channels: consensus, disagreement, and the structural relationship between them.

## Repository Structure

```
PRIMED/
├── paper/                        Manuscript (LaTeX)
│   ├── primed_paper.tex            Main manuscript
│   └── primed_paper.pdf            Compiled PDF
├── code/
│   ├── simulation/                 Simulation studies
│   │   ├── 01_primed_core.jl        PRIMED optimizer (K <= 50)
│   │   ├── 01_primed_core_highdim.jl  PRIMED optimizer (K = 100)
│   │   ├── 02_run_single_job.R      Single simulation job (R wrapper)
│   │   ├── 03_run_parallel.sh       Parallel execution (40 jobs)
│   │   └── 04_aggregate_results.R   Result aggregation
│   ├── rectal/                     Rectal cancer MRI application
│   │   ├── 01_cd_group_lasso.R      CD Group LASSO baseline
│   │   ├── 02_primed_linear.R       PRIMED (linear dispersion)
│   │   ├── 03_primed_quadratic.R    PRIMED (quadratic dispersion)
│   │   ├── 04_bootstrap_stability.R Bootstrap stability analysis
│   │   └── 05_scatter_plot.R        Consensus-disagreement scatter plot
│   └── lidc/                       LIDC-IDRI lung nodule application
│       ├── 01_extract_radiomics.py  PyRadiomics feature extraction
│       ├── 02_apply_primed.R        PRIMED application
│       ├── 03_cv10_evaluation.R     10-fold nested CV evaluation
│       └── 04_bootstrap_stability.R Bootstrap stability analysis
├── data/
│   ├── data.xlsx                   Rectal cancer multi-reader MRI data
│   └── lidc_annotations.csv        LIDC-IDRI lung nodule annotations
├── results/                      Analysis results
│   ├── simulation/
│   └── lidc/
└── figures/
```

## Requirements

- **R** (>= 4.0) with packages: `tidyverse`, `readxl`, `grpreg`, `pROC`, `glmnet`, `parallel`
- **Julia** (>= 1.8) with standard libraries
- **Python** (>= 3.8) with `pyradiomics`, `SimpleITK`, `numpy`, `pandas` (for LIDC feature extraction only)

## Running the Code

### Simulation

```bash
cd code/simulation
bash 03_run_parallel.sh          # 40 jobs, 8 parallel
Rscript 04_aggregate_results.R   # Aggregate results
```

### Rectal Cancer Application

```bash
cd code/rectal
Rscript 01_cd_group_lasso.R      # CD Group LASSO baseline
Rscript 02_primed_linear.R       # PRIMED (linear dispersion)
Rscript 03_primed_quadratic.R    # PRIMED (quadratic dispersion)
```

### LIDC-IDRI Application

```bash
cd code/lidc
Rscript 03_cv10_evaluation.R     # 10-fold nested CV (requires radiomics_cd.csv)
```

## Data Availability

- **LIDC-IDRI**: Publicly available at [The Cancer Imaging Archive (TCIA)](https://www.cancerimagingarchive.net/collection/lidc-idri/)
- **Rectal cancer MRI**: Available from the corresponding author on reasonable request

## License

This project is licensed under the MIT License.

## Citation

If you use this code, please cite:

> Kwon Y, Han K. PRIMED: Joint Loss Estimation with Selective Penalization for Multi-Reader Imaging Data with Binary Outcomes. *BMC Medical Research Methodology*. (submitted)

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
│   ├── simulation/                 Simulation studies (Tables 1-9)
│   │   ├── 01_primed_core.jl        PRIMED optimizer (K <= 50)
│   │   ├── 01_primed_core_highdim.jl  PRIMED optimizer (K = 100)
│   │   ├── 02_run_single_job.R      Single simulation job
│   │   ├── 03_run_parallel.sh       Parallel execution (40 jobs)
│   │   └── 04_aggregate_results.R   Result aggregation
│   ├── rectal/                     Rectal cancer application (Tables 4-7)
│   │   ├── 01_primed_linear.R       PRIMED with linear dispersion
│   │   ├── 02_primed_quadratic.R    Linear vs quadratic comparison
│   │   └── 03_scatter_plot.R        Consensus-disagreement scatter plot
│   └── lidc/                       LIDC-IDRI application (Tables 6, 8)
│       ├── 01_extract_radiomics.py  PyRadiomics feature extraction
│       └── 02_cv10_evaluation.R     10-fold nested CV evaluation
├── data/
│   ├── data.xlsx                   Rectal cancer multi-reader MRI data
│   └── lidc_annotations.csv        LIDC-IDRI lung nodule annotations
├── results/
└── figures/
```

## Requirements

- **R** (>= 4.0) with packages: `tidyverse`, `readxl`, `grpreg`, `pROC`, `glmnet`, `parallel`
- **Julia** (>= 1.8) with standard libraries
- **Python** (>= 3.8) with `pyradiomics`, `SimpleITK` (for LIDC feature extraction only)

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
Rscript 01_primed_linear.R       # PRIMED + baselines (Tables 4-5)
Rscript 02_primed_quadratic.R    # Quadratic comparison (Tables 6-7)
```

### LIDC-IDRI Application

```bash
cd code/lidc
Rscript 02_cv10_evaluation.R     # 10-fold nested CV (Tables 6, 8)
```

## Data Availability

- **LIDC-IDRI**: Publicly available at [The Cancer Imaging Archive (TCIA)](https://www.cancerimagingarchive.net/collection/lidc-idri/)
- **Rectal cancer MRI**: Available from the corresponding author on reasonable request

## License

This project is licensed under the MIT License.

## Citation

If you use this code, please cite:

> Kwon Y, Han K. PRIMED: Joint Loss Estimation with Selective Penalization for Multi-Reader Imaging Data with Binary Outcomes. *BMC Medical Research Methodology*. (submitted)

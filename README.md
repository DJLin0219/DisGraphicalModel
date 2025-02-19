# Distributional Graphical Models for Functional Data

This repository contains the implementation of the paper  
**Graphical Models for Random Distributions of Functional Data with Application to fMRI**  
by Dianjun Lin, Bing Li and Lingzhou Xue. The paper introduces a novel framework for estimating graphical models where each node represents a distribution of functional data rather than a single averaged time series. By leveraging kernel mean embeddings and reproducing kernel Hilbert spaces (RKHS), the method captures complex dependencies among regions of interest (ROIs) in fMRI data without losing higher‐order information.

## Features
- Implements the Distributional Additive Semigraphoid Graphical Model (DASG) framework.
- Models each ROI’s signal as a distribution of random functions via kernel mean embeddings.
- Provides several sparse graph estimation methods including:
  - **MAPO**: Mean Additive Precision Operator estimator.
  - **MARO**: Mean Additive Regression Operator via neighborhood selection.
  - **MGLASSO**: Mean Functional Graphical Lasso.
- Supports both balanced and unbalanced sampling schemes.
- Evaluates performance through simulation studies and real rsfMRI (ADHD) data analysis.

## Dependencies
The following R packages are required:
- `Matrix`
- `parallel`
- `fda`

Install them from CRAN:
```r
install.packages(c("Matrix", "parallel", "fda"))





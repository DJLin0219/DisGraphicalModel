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

## Installation
1. Clone the repository:
   ```sh
   git clone https://github.com/your_username/DistributionalGraphicalModels.git
   cd DistributionalGraphicalModels
   ```

2. Install required R packages:
   ```r
   install.packages(c("Matrix", "parallel", "fda"))
   ```

3. Compile the R source files:
   ```r
   file_path <- "/path/to/repository"
   source(paste(file_path, "MAPO.cpp", sep="/"))
   ```
   
## Results
- Simulation Studies: Demonstrate the effectiveness of the proposed methods in recovering true graph structures under various dependency models.
- rsfMRI Data Analysis: Applied to ADHD data, the methods reveal interpretable brain connectivity patterns that differentiate ADHD subjects from controls.

## Citation
If you use this code in your research, please cite (The paper will be posted to arxiv soon).

## Contact

For any questions or issues, please open an issue in this repository or contact dzl5618@psu.edu.

## License
This project is licensed under the MIT License.



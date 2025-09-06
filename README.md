# Gower-kPrototypes Predictive Mean Matching (GKP-PMM)

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![R](https://img.shields.io/badge/R-4.0+-blue.svg)](https://www.r-project.org/)

A custom imputation method for the MICE package that integrates Gower's distance for refined donor selection within a Predictive Mean Matching (PMM) framework. This method supports mixed data types and includes optional pre-clustering for improved performance.

## Features

- **Mixed Data Types**: Handles numeric, categorical, and ordinal variables seamlessly
- **Gower's Distance**: Uses appropriate distance metrics for different data types
- **Optional Pre-clustering**: k-prototypes clustering to improve donor selection
- **MICE Integration**: Fully compatible with the Multiple Imputation by Chained Equations framework
- **Robust Validation**: Automatic model selection and error handling

## Installation

### Prerequisites
```r
install.packages(c("mice", "cluster", "clustMixType", "MASS", "nnet",
                   "glmnet", "survey", "dplyr", "ggplot2", "tidyr"))
```

### Loading the Method
```r
source("mice.impute.gkp_pmm.r")
```

## Usage

### Basic Example
```r
library(mice)
data(nhanes)

# Register the custom method
mice.impute.gkp_pmm <- function(...) { ... }  # Load the function

# Perform imputation
imp <- mice(nhanes, method = "gkp_pmm", donors = 5, k_pre_clusters = 3)
```

### Advanced Usage
```r
# With custom parameters
imp <- mice(nhanes,
            method = "gkp_pmm",
            donors = 7,
            k_pre_clusters = 5,
            predictive_model = "auto",
            pmm_pool_factor = 10)
```

## Parameters

- `donors`: Number of nearest neighbors for final imputation (default: 5)
- `k_pre_clusters`: Number of clusters for pre-clustering (default: 0, disabled)
- `predictive_model`: Predictive model type ("auto", "lm", "logit", "polr", "multinom")
- `pmm_pool_factor`: Factor for initial PMM candidate pool size (default: 5)

## Project Structure

```
├── mice.impute.gkp_pmm.r    # Main imputation function
├── source.r                 # Two-phase sampling methodology
├── thesees code.r           # Setup and utility scripts
├── Imputer paper.ipynb      # Example notebook
├── enhance simulation.ipynb # Enhanced simulation framework
├── Thesis coded.ipynb       # Thesis implementation
├── README.md                # This file
└── .gitignore              # Git ignore rules
```

## Methodology

1. **Data Preparation**: Separate observed and missing data
2. **Model Training**: Train predictive model on observed data
3. **Optional Pre-clustering**: Use k-prototypes for donor grouping
4. **Donor Selection**: Find nearest neighbors using Gower's distance
5. **Imputation**: Sample from selected donors

## Related Work

This implementation is part of a PhD thesis on optimal selection and variance estimation with multiple auxiliary variables in two-phase sampling.

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Citation

If you use this method in your research, please cite:

```
[Your Name]. (2025). A Methodology for Optimal Selection and Variance Estimation
with Multiple Auxiliary Variables in Two-Phase Sampling [PhD Thesis].
```

## Contributing

1. Fork the repository
2. Create a feature branch
3. Make your changes
4. Submit a pull request

## Contact

[Your Name] - [tertimothy@gmail.com]

Project Link: [https://github.com/timmoyaj/gkp-pmm](https://github.com/Timoyaj/gkp-pmm)

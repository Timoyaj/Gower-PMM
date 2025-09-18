# Gower-PMM Benchmark Framework

A comprehensive benchmarking framework for evaluating the Gower-PMM imputation method against alternative approaches for handling missing data in mixed-type datasets. Designed specifically for thesis research and methodological evaluation.

## Overview

This benchmark framework provides a **two-part evaluation structure** for comprehensive methodological assessment:

### Part 1: Distance Engine Comparison
- **Engines Compared**: Gower-PMM, FD::gowdis, cluster::daisy, StatMatch::gower.dist, ade4::dist.ktab, Euclidean, Manhattan
- **Focus**: Distance/dissimilarity measure quality for mixed-type data
- **Metrics**: Rank correlation balance, nearest neighbor preservation, computational efficiency
- **Data**: Complete datasets (no missing values)

### Part 2: Imputation Method Comparison
- **Methods Compared**: Gower-PMM variants, MICE PMM/CART/RF, VIM k-NN/IRM, mean/mode, random, hot deck, regression
- **Focus**: Imputation quality and performance under different missing data patterns
- **Metrics**: RMSE, MAE, bias, distribution preservation, correlation preservation, computational efficiency
- **Data**: Datasets with various missing data patterns (MCAR, MAR, MNAR)

### Framework Features
- **Multiple Evaluation Metrics**: Comprehensive assessment of quality and efficiency
- **Diverse Simulation Scenarios**: Various dataset sizes, variable types, correlation structures, and noise levels
- **Thesis-Ready Reports**: Automated generation of publication-quality tables, figures, and statistical analysis
- **Parallel Processing**: Efficient execution using multiple CPU cores
- **Reproducibility**: Full control over random seeds and comprehensive logging

## Quick Start

### Complete Two-Part Benchmark (Recommended for Thesis)

```r
# Load the package
library(gowerpmm)

# Run complete two-part benchmark (distance engines + imputation methods)
results <- run_complete_two_part_benchmark(
  n_replications_distance = 50,
  n_replications_imputation = 50,
  parallel = TRUE
)

# Results are automatically saved and reports generated
# Check the output directories for comprehensive analysis
```

### Individual Part Benchmarks

```r
# Part 1: Distance Engine Comparison
distance_config <- create_distance_benchmark_config(n_replications = 100)
distance_results <- run_distance_comparison_benchmark(distance_config, parallel = TRUE)
generate_distance_comparison_report(distance_results, "distance_report")

# Part 2: Imputation Method Comparison
imputation_config <- create_imputation_benchmark_config(n_replications = 100)
imputation_results <- run_imputation_comparison_benchmark(imputation_config, parallel = TRUE)
generate_imputation_comparison_report(imputation_results, "imputation_report")
```

## Framework Components

### Core Framework (`R/benchmark_framework.R`)

- **BenchmarkConfig Class**: S4 class for storing benchmark configurations
- **Simulation Scenarios**: Pre-defined and custom data generation scenarios
- **Missing Data Patterns**: Various mechanisms for introducing missing values
- **Method Specifications**: Configuration for different methods

### Distance Engine Comparison (`R/benchmark_distance_comparison.R`)

- **Distance Engines**: Gower-PMM, FD::gowdis, cluster::daisy, StatMatch::gower.dist, ade4::dist.ktab, Euclidean, Manhattan
- **Quality Metrics**: Rank correlation balance, nearest neighbor preservation
- **Performance Analysis**: Computational efficiency and scalability
- **Specialized Reporting**: Distance-specific analysis and visualization

### Imputation Method Comparison (`R/benchmark_imputation_comparison.R`)

- **Imputation Methods**: Gower-PMM variants, MICE PMM/CART/RF, VIM k-NN/IRM, mean/mode, random, hot deck, regression
- **Quality Metrics**: RMSE, MAE, bias, distribution preservation, correlation preservation
- **Missing Pattern Analysis**: MCAR, MAR, MNAR performance evaluation
- **Specialized Reporting**: Imputation-specific analysis and visualization

### Evaluation Metrics (`R/benchmark_evaluation.R`)

- **Distance Quality Metrics**: Rank correlation balance, nearest neighbor preservation
- **Imputation Quality Metrics**: RMSE, MAE, bias, distribution preservation, correlation preservation
- **Computational Performance**: Execution time, memory usage, convergence analysis

### Execution Engine (`R/benchmark_execution.R`)

- **Single Run Execution**: `run_single_benchmark()` for individual replications
- **Batch Execution**: `run_benchmark_suite()` for complete benchmark suites
- **Parallel Processing**: Multi-core execution support
- **Error Handling**: Robust error handling and logging

### Analysis and Visualization (`R/benchmark_analysis.R`)

- **Statistical Analysis**: T-tests, ranking analysis, sensitivity analysis
- **Publication-Quality Plots**: Box plots, efficiency frontiers, heatmaps
- **Thesis-Ready Tables**: Formatted summary tables for academic papers
- **Automated Reporting**: R Markdown report generation

### Usage Examples (`R/benchmark_examples.R`)

- **Two-Part Benchmark**: Complete evaluation combining distance and imputation analysis
- **Individual Benchmarks**: Separate distance engine or imputation method evaluation
- **Thesis Benchmark**: Comprehensive evaluation for dissertation work
- **Specialized Benchmarks**: Distance comparison, efficiency analysis, validation
- **Sensitivity Analysis**: Parameter influence assessment

## Key Features

### Simulation Scenarios

- **Small datasets** (n=100): Quick testing and validation
- **Medium datasets** (n=500): Balanced evaluation scenarios
- **Large datasets** (n=1000+): Performance and scalability testing
- **High-dimensional data**: Many variables with mixed types
- **Correlated structures**: Various correlation patterns from low to very high
- **Noisy data**: Different noise levels for robustness testing

### Imputation Methods

- **Gower-PMM Variants**:
  - Auto-optimized weights (GA)
  - Equal weights
  - Range vs IQR scaling
- **Traditional Methods**:
  - FD::gowdis + PMM
  - MICE PMM, CART, RF
  - VIM k-NN, IRMI
- **Extensible Framework**: Easy to add new methods

### Evaluation Metrics

#### Imputation Quality
- **Numeric**: RMSE, MAE, bias, NRMSE, R²
- **Categorical**: Accuracy, mean level difference
- **Distribution**: KS test, mean/SD preservation
- **Correlation**: Preservation of variable relationships

#### Distance Quality
- **Balance**: Rank correlation across variable types
- **Preservation**: Nearest neighbor structure maintenance
- **Efficiency**: Computational performance metrics

## Usage Examples

### Complete Two-Part Benchmark (Recommended)

```r
# Run both distance engine and imputation method comparisons
results <- run_complete_two_part_benchmark(
  n_replications_distance = 50,
  n_replications_imputation = 50,
  parallel = TRUE
)

# Automatically generates:
# - Distance comparison report (two_part_distance_report/)
# - Imputation comparison report (two_part_imputation_report/)
# - Combined summary report (two_part_combined_report/)
```

### Individual Part Benchmarks

#### Distance Engine Comparison
```r
# Configure distance engine comparison
distance_config <- create_distance_benchmark_config(
  n_replications = 100,
  engines = c("gowerpmm", "fd_gowdis", "daisy", "euclidean")
)

# Run distance comparison
distance_results <- run_distance_comparison_benchmark(distance_config, parallel = TRUE)

# Generate specialized distance report
generate_distance_comparison_report(distance_results, "distance_analysis")
```

#### Imputation Method Comparison
```r
# Configure imputation method comparison
imputation_config <- create_imputation_benchmark_config(
  n_replications = 100,
  methods = c("gowerpmm_auto", "mice_pmm", "vim_knn", "mean_mode")
)

# Run imputation comparison
imputation_results <- run_imputation_comparison_benchmark(imputation_config, parallel = TRUE)

# Generate specialized imputation report
generate_imputation_comparison_report(imputation_results, "imputation_analysis")
```

### Custom Benchmark Configuration

```r
# Create custom distance engine configuration
custom_distance_config <- create_distance_benchmark_config(
  n_replications = 30,
  engines = c("gowerpmm", "fd_gowdis"),  # Only compare these two
  output_dir = "custom_distance_results"
)

# Create custom imputation configuration
custom_imputation_config <- create_imputation_benchmark_config(
  n_replications = 30,
  methods = c("gowerpmm_auto", "gowerpmm_equal"),  # Only Gower-PMM variants
  output_dir = "custom_imputation_results"
)
```

## Output and Results

### Two-Part Benchmark File Structure

```
two_part_distance_results/
├── distance_engine_analysis.rds     # Distance analysis results
├── distance_balance_comparison.png  # Balance score plots
├── distance_efficiency_frontier.png # Efficiency frontier
├── distance_engine_rankings.png     # Ranking plots
└── distance_comparison_report.html  # Distance report

two_part_imputation_results/
├── imputation_method_analysis.rds   # Imputation analysis results
├── imputation_rmse_comparison.png   # RMSE comparison plots
├── imputation_efficiency_frontier.png # Efficiency frontier
├── imputation_method_rankings.png   # Ranking plots
└── imputation_comparison_report.html # Imputation report

two_part_combined_report/
├── combined_performance.png         # Combined analysis plot
└── combined_benchmark_report.html   # Combined summary report
```

### Individual Benchmark File Structure

#### Distance Engine Results
```
distance_benchmark_results/
├── distance_engine_analysis.rds     # Analysis results
├── distance_balance_comparison.png  # Balance score comparison
├── distance_efficiency_frontier.png # Quality vs speed
├── distance_engine_rankings.png     # Performance rankings
└── distance_comparison_report.html  # Specialized report
```

#### Imputation Method Results
```
imputation_benchmark_results/
├── imputation_method_analysis.rds   # Analysis results
├── imputation_rmse_comparison.png   # RMSE comparison
├── imputation_efficiency_frontier.png # Quality vs speed
├── imputation_method_rankings.png   # Performance rankings
└── imputation_comparison_report.html # Specialized report
```

### Key Output Files

- **RDS Files**: Complete R objects for further analysis
- **PNG Files**: Publication-quality figures for thesis/dissertation
- **HTML Reports**: Comprehensive analysis reports with tables and figures
- **Analysis Objects**: Specialized analysis results for each benchmark type

## Statistical Analysis

### Comparative Tests

- **Pairwise t-tests**: Compare method performance differences
- **Ranking Analysis**: Overall method rankings across scenarios
- **Sensitivity Analysis**: Parameter influence assessment

### Validation

- **Reproducibility Testing**: Multiple runs for consistency
- **Convergence Analysis**: Ranking stability with more replications
- **Error Analysis**: Failure mode identification

## Performance Considerations

### Computational Requirements

- **Small benchmarks** (10-50 replications): Minutes on standard hardware
- **Thesis benchmarks** (100+ replications): Hours with parallel processing
- **Memory usage**: Scales with dataset size and number of methods

### Parallel Processing

```r
# Use all available cores
results <- run_benchmark_suite(config, parallel = TRUE)

# Specify number of cores
results <- run_benchmark_suite(config, parallel = TRUE, n_cores = 4)
```

### Memory Management

- Large datasets may require significant RAM
- Consider reducing replications for memory-constrained systems
- Monitor memory usage with `gc()` and system tools

## Extending the Framework

### Adding New Methods

```r
# Define new method in config
new_method <- list(
  name = "My New Method",
  type = "custom",
  params = list(
    parameter1 = value1,
    parameter2 = value2
  )
)

# Implement imputation function
impute_my_method <- function(data, params) {
  # Your imputation logic here
  return(imputed_data)
}

# Add to execution engine
run_custom_imputation <- function(data, params) {
  # Call your method
  imputed_data <- impute_my_method(data, params)
  list(imputed_data = imputed_data, convergence_info = list(...))
}
```

### Adding New Metrics

```r
# Define custom evaluation function
evaluate_custom_metric <- function(original, imputed, missing_mask) {
  # Your metric calculation
  return(metric_value)
}

# Add to evaluation pipeline
custom_results <- evaluate_custom_metric(original_data, imputed_data, missing_mask)
```

## Dependencies

### Required Packages

- **Core**: `mice`, `FD`, `GA`, `MASS`, `moments`
- **Analysis**: `dplyr`, `tidyr`, `ggplot2`, `kableExtra`
- **Parallel**: `parallel`, `doParallel` (optional)
- **Reporting**: `rmarkdown` (optional)

### Installation

```r
# Install required packages
install.packages(c("mice", "FD", "GA", "MASS", "moments",
                   "dplyr", "tidyr", "ggplot2", "kableExtra"))

# Optional for parallel processing
install.packages(c("parallel", "doParallel"))

# Optional for reporting
install.packages("rmarkdown")
```

## Troubleshooting

### Common Issues

1. **Memory errors**: Reduce dataset sizes or replications
2. **Package not found**: Install missing dependencies
3. **Parallel failures**: Try sequential execution (`parallel = FALSE`)
4. **Long runtimes**: Use smaller configurations for testing

### Debugging

```r
# Run single replication for testing
single_result <- run_single_benchmark(scenario, method, pattern, 1, config)

# Check for errors
if (!single_result$success) {
  print(single_result$error_message)
}
```

## Citation

If you use this benchmark framework in your research, please cite:

```
T. Yaji. (2024). Gower-PMM Benchmark Framework: Comprehensive Evaluation
of Imputation Methods for Mixed-Type Data. R package version 0.1.0.
```

## Support

For issues, questions, or contributions:

1. Check the documentation and examples
2. Review existing issues on the project repository
3. Create a new issue with reproducible examples
4. Contact the maintainer: tertimothy@gmail.com

## License

This benchmark framework is part of the Gower-PMM package and follows the same GPL-3 license.
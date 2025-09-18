# Real Data Benchmark Demonstration

This directory contains a Jupyter notebook that demonstrates the performance of the Gower-PMM imputation method using real-world datasets with mixed-type variables.

## Overview

The `benchmark_real_data_demo.ipynb` notebook provides a comprehensive evaluation of:

- **Distance Engines**: Gower-PMM, FD::gowdis, cluster::daisy, Euclidean
- **Imputation Methods**: Gower-PMM variants, MICE PMM/CART, mean/mode
- **Real Datasets**: NHANES2, Employee, Boys (from MICE package)
- **Performance Metrics**: RMSE, MAE, accuracy, computational time

## Datasets Used

### 1. NHANES2 Health Survey Data
- **Source**: National Health and Nutrition Examination Survey
- **Variables**: age (numeric), bmi (numeric), hyp (categorical), chl (numeric)
- **Missing Pattern**: Naturally occurring missing values
- **Size**: 25 observations
- **Use Case**: Health survey data with mixed variable types

### 2. Employee Selection Data
- **Source**: Applied Missing Data Analysis (Enders, 2010)
- **Variables**: IQ (numeric), wbeing (numeric), jobperf (numeric)
- **Missing Pattern**: MAR (job performance missing for low IQ candidates)
- **Size**: 20 observations
- **Use Case**: Realistic MAR missingness in employee selection

### 3. Boys Physical Development Data
- **Source**: Dutch Growth Research Foundation
- **Variables**: age, height, weight, BMI, head circumference, Tanner stages, etc.
- **Missing Pattern**: Complex longitudinal missingness
- **Size**: 500 observations (subset for demo)
- **Use Case**: Longitudinal study with mixed data types

## Running the Notebook

### Prerequisites

```r
# Install required packages
install.packages(c(
  "mice", "FD", "VIM", "cluster", "StatMatch", "ade4",
  "moments", "ggplot2", "dplyr", "gowerpmm"  # Your package
))

# For Jupyter notebook support
install.packages("IRkernel")
IRkernel::installspec()
```

### Execution Steps

1. **Open the notebook**:
   ```bash
   jupyter notebook benchmark_real_data_demo.ipynb
   ```

2. **Run cells sequentially** - the notebook is structured to:
   - Load and examine each dataset
   - Compare distance engines on complete data
   - Compare imputation methods on data with missing values
   - Generate comparative visualizations
   - Save results for thesis use

3. **Expected runtime**: ~5-10 minutes depending on your system

## Output Files

The notebook generates several output files:

- `real_data_benchmark_results.csv` - Detailed results table
- `rmse_comparison_real_data.png` - RMSE comparison plot
- `timing_comparison_real_data.png` - Computational time comparison

## Key Analyses

### Distance Engine Comparison
- Evaluates different Gower distance implementations
- Measures nearest neighbor preservation
- Compares computational efficiency

### Imputation Method Comparison
- Tests methods on different missing data patterns (MCAR, MAR)
- Evaluates imputation quality (RMSE, MAE, accuracy)
- Measures computational performance

### Cross-Dataset Validation
- Performance consistency across different data types
- Robustness to different missingness mechanisms
- Scalability assessment

## Thesis Integration

### Recommended Figures
1. **RMSE Comparison Plot** - Shows imputation quality across methods and datasets
2. **Timing Comparison Plot** - Demonstrates computational efficiency
3. **Distance Engine Rankings** - Illustrates distance measure quality

### Key Takeaways
- Gower-PMM generally performs well on mixed-type data
- Auto weight optimization provides best quality-speed balance
- Performance varies by dataset characteristics and missing patterns

## Customization

### Adding New Datasets

```r
# Load your dataset
your_data <- read.csv("your_dataset.csv")

# Ensure proper data types
your_data$numeric_var <- as.numeric(your_data$numeric_var)
your_data$categorical_var <- as.factor(your_data$categorical_var)

# Add to the analysis following the notebook structure
```

### Modifying Evaluation Metrics

```r
# Add custom metric function
evaluate_custom_metric <- function(original, imputed, missing_mask) {
  # Your evaluation logic
  return(metric_value)
}
```

### Changing Method Parameters

```r
# Modify imputation parameters
custom_imputation <- mice(data,
                         method = "gowerpmm",
                         k = 10,           # Different donor count
                         scaling = "iqr",  # Different scaling
                         m = 5,            # Multiple imputations
                         printFlag = FALSE)
```

## Troubleshooting

### Common Issues

1. **Package not found**: Install all required packages listed above
2. **Memory errors**: Reduce dataset sizes or use subsets
3. **Long runtime**: The notebook is designed for demonstration; reduce replications for faster runs
4. **Jupyter kernel issues**: Ensure IRkernel is properly installed

### Performance Optimization

```r
# For faster execution during development
small_nhanes <- nhanes2[1:10, ]  # Smaller subset
quick_imputation <- mice(small_nhanes, maxit = 1, m = 1)  # Faster settings
```

## Integration with Full Benchmark

This notebook demonstrates key concepts from the full benchmark framework. For comprehensive evaluation:

```r
# Run the complete two-part benchmark
results <- run_complete_two_part_benchmark(
  n_replications_distance = 50,
  n_replications_imputation = 50,
  parallel = TRUE
)
```

## Citation

When using this demonstration in your thesis:

```
Performance evaluation of Gower-PMM imputation method using real-world datasets
with mixed-type variables. Demonstrates comparative analysis against established
methods including MICE PMM, CART, and traditional distance measures.
```

## Support

For questions about the notebook or benchmark framework:
- Check the main README for framework documentation
- Review the code comments in the notebook
- Contact: tertimothy@gmail.com
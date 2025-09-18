#' Benchmark Framework Usage Examples
#'
#' This file provides comprehensive examples for using the Gower-PMM benchmark framework
#' in thesis research and methodological evaluation.

# --- Two-Part Benchmark Structure ---

#' Run Complete Two-Part Benchmark
#'
#' Executes both distance engine comparison and imputation method comparison
#' as a comprehensive evaluation for thesis research
#'
#' @param n_replications_distance Number of replications for distance comparison
#' @param n_replications_imputation Number of replications for imputation comparison
#' @param parallel Whether to use parallel processing
#' @return List with both benchmark results
#' @export
run_complete_two_part_benchmark <- function(n_replications_distance = 50,
                                          n_replications_imputation = 50,
                                          parallel = TRUE) {

  message("Starting Complete Two-Part Benchmark Analysis")
  message("==============================================")

  results <- list()

  # Part 1: Distance Engine Comparison
  message("\nPart 1: Distance Engine Comparison")
  message("-----------------------------------")

  distance_config <- create_distance_benchmark_config(
    n_replications = n_replications_distance,
    output_dir = "two_part_distance_results"
  )

  results$distance <- run_distance_comparison_benchmark(distance_config, parallel = parallel)

  # Generate distance report
  generate_distance_comparison_report(results$distance, "two_part_distance_report")

  # Part 2: Imputation Method Comparison
  message("\nPart 2: Imputation Method Comparison")
  message("-------------------------------------")

  imputation_config <- create_imputation_benchmark_config(
    n_replications = n_replications_imputation,
    output_dir = "two_part_imputation_results"
  )

  results$imputation <- run_imputation_comparison_benchmark(imputation_config, parallel = parallel)

  # Generate imputation report
  generate_imputation_comparison_report(results$imputation, "two_part_imputation_report")

  # Generate combined summary
  generate_combined_summary_report(results, "two_part_combined_report")

  message("\nBenchmark Complete!")
  message("===================")
  message("Results saved in:")
  message("- Distance comparison: two_part_distance_results/")
  message("- Imputation comparison: two_part_imputation_results/")
  message("- Combined report: two_part_combined_report/")

  results
}

#' Generate Combined Summary Report
#' @keywords internal
generate_combined_summary_report <- function(results, output_dir = "combined_report") {

  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  # Create combined plots
  combined_plot <- create_combined_performance_plot(results)

  # Save combined plot
  ggplot2::ggsave(file.path(output_dir, "combined_performance.png"), combined_plot,
                  width = 12, height = 8, dpi = 300)

  # Generate combined report
  rmd_content <- create_combined_rmd_report(results)

  rmd_file <- file.path(output_dir, "combined_benchmark_report.Rmd")
  writeLines(rmd_content, rmd_file)

  # Try to render
  tryCatch({
    rmarkdown::render(rmd_file, output_dir = output_dir)
    message("Combined report generated successfully")
  }, error = function(e) {
    warning("Could not render combined report: ", e$message)
  })
}

#' Create Combined Performance Plot
#' @keywords internal
create_combined_performance_plot <- function(results) {

  # Extract key metrics from both benchmarks
  distance_data <- results$distance@summary_stats$by_engine %>%
    dplyr::select(method, mean_balance_score, mean_nn_preservation, mean_time) %>%
    dplyr::mutate(benchmark = "Distance Engines")

  imputation_data <- results$imputation@summary_stats$by_method %>%
    dplyr::select(method, mean_rmse, mean_mae, mean_time) %>%
    dplyr::mutate(benchmark = "Imputation Methods")

  # Create a combined visualization
  p <- ggplot2::ggplot() +
    # Distance engines
    ggplot2::geom_point(data = distance_data,
                       ggplot2::aes(x = mean_time, y = mean_balance_score,
                                   color = "Distance Engines", shape = "Distance Engines"),
                       size = 3, alpha = 0.7) +
    ggplot2::geom_text(data = distance_data,
                      ggplot2::aes(x = mean_time, y = mean_balance_score, label = method),
                      vjust = -1, size = 3, color = "blue") +
    # Imputation methods
    ggplot2::geom_point(data = imputation_data,
                       ggplot2::aes(x = mean_time, y = mean_rmse,
                                   color = "Imputation Methods", shape = "Imputation Methods"),
                       size = 3, alpha = 0.7) +
    ggplot2::geom_text(data = imputation_data,
                      ggplot2::aes(x = mean_time, y = mean_rmse, label = method),
                      vjust = 1, size = 3, color = "red") +
    ggplot2::labs(
      title = "Combined Benchmark Performance Overview",
      subtitle = "Distance Engines (Balance Score) vs Imputation Methods (RMSE)",
      x = "Mean Execution Time (seconds)",
      y = "Performance Metric",
      color = "Benchmark Type",
      shape = "Benchmark Type"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::scale_x_log10()

  p
}

#' Create Combined RMD Report
#' @keywords internal
create_combined_rmd_report <- function(results) {

  rmd <- sprintf('---
title: "Combined Benchmark Analysis Report"
author: "Gower-PMM Two-Part Evaluation"
date: "`r Sys.Date()`"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = FALSE, warning = FALSE, message = FALSE)
library(ggplot2)
library(dplyr)
library(kableExtra)
```

# Combined Benchmark Analysis: Distance Engines & Imputation Methods

This report provides a comprehensive evaluation combining both parts of the benchmark framework.

## Overview

### Distance Engine Comparison
- **Engines Compared**: %s
- **Scenarios Tested**: %d
- **Total Runs**: %d
- **Successful Runs**: %d (%.1f%%)

### Imputation Method Comparison
- **Methods Compared**: %s
- **Scenarios Tested**: %d
- **Missing Patterns**: %d
- **Total Runs**: %d
- **Successful Runs**: %d (%.1f%%)

## Key Findings

### Distance Engine Rankings
```{r distance-rankings}
distance_summary <- results$distance@summary_stats$by_engine
kable(distance_summary, caption = "Distance Engine Performance Summary") %>%
  kable_styling(bootstrap_options = c("striped", "hover"))
```

### Imputation Method Rankings
```{r imputation-rankings}
imputation_summary <- results$imputation@summary_stats$by_method
kable(imputation_summary, caption = "Imputation Method Performance Summary") %>%
  kable_styling(bootstrap_options = c("striped", "hover"))
```

## Combined Performance Overview

### Performance Comparison
![Combined Performance](combined_performance.png)

## Conclusions

### Distance Engines
- **Top Performer**: %s (Balance Score: %.3f)
- **Fastest**: %s (Time: %.3f sec)
- **Most Efficient**: %s (Efficiency: %.3f)

### Imputation Methods
- **Best Quality**: %s (RMSE: %.3f)
- **Fastest**: %s (Time: %.3f sec)
- **Most Efficient**: %s (Efficiency: %.3f)

## Technical Details

- **R Version**: `r R.version.string`
- **Analysis Date**: `r format(results$distance@timestamp, "%%Y-%%m-%%d %%H:%%M")`
- **Total Distance Benchmark Time**: `r round(results$distance@execution_time, 2)` seconds
- **Total Imputation Benchmark Time**: `r round(results$imputation@execution_time, 2)` seconds
',
                 # Distance engine stats
                 paste(names(results$distance@config@methods), collapse = ", "),
                 length(results$distance@config@scenarios),
                 results$distance@summary_stats$overall$total_runs,
                 results$distance@summary_stats$overall$successful_runs,
                 results$distance@summary_stats$overall$success_rate * 100,

                 # Imputation method stats
                 paste(names(results$imputation@config@methods), collapse = ", "),
                 length(results$imputation@config@scenarios),
                 length(results$imputation@config@missing_patterns),
                 results$imputation@summary_stats$overall$total_runs,
                 results$imputation@summary_stats$overall$successful_runs,
                 results$imputation@summary_stats$overall$success_rate * 100,

                 # Top performers
                 distance_summary$method[which.min(distance_summary$mean_balance_score)],
                 min(distance_summary$mean_balance_score, na.rm = TRUE),
                 distance_summary$method[which.min(distance_summary$mean_time)],
                 min(distance_summary$mean_time, na.rm = TRUE),
                 distance_summary$method[which.max(distance_summary$mean_efficiency)],
                 max(distance_summary$mean_efficiency, na.rm = TRUE),

                 imputation_summary$method[which.min(imputation_summary$mean_rmse)],
                 min(imputation_summary$mean_rmse, na.rm = TRUE),
                 imputation_summary$method[which.min(imputation_summary$mean_time)],
                 min(imputation_summary$mean_time, na.rm = TRUE),
                 imputation_summary$method[which.max(imputation_summary$mean_efficiency)],
                 max(imputation_summary$mean_efficiency, na.rm = TRUE))

  rmd
}

# --- Quick Start Example ---

#' Quick Benchmark Example
#'
#' Demonstrates a minimal benchmark setup for rapid testing
#'
#' @return BenchmarkResults object
#' @export
run_quick_benchmark <- function() {

  # Create a minimal configuration
  config <- create_default_benchmark_config(
    n_replications = 10,  # Small number for quick testing
    output_dir = "quick_benchmark_results"
  )

  # Select only a few methods and scenarios for quick testing
  config@methods <- config@methods[c("gowerpmm_auto", "gowerpmm_equal", "mice_pmm")]
  config@scenarios <- config@scenarios[c("small_balanced", "small_mixed")]
  config@missing_patterns <- config@missing_patterns[c("mcar_20", "mar_20")]

  message("Running quick benchmark with:")
  message("- Methods: ", paste(names(config@methods), collapse = ", "))
  message("- Scenarios: ", paste(names(config@scenarios), collapse = ", "))
  message("- Missing patterns: ", paste(names(config@missing_patterns), collapse = ", "))
  message("- Replications: ", config@n_replications)

  # Run the benchmark
  results <- run_benchmark_suite(config, parallel = FALSE)

  # Save results
  save_benchmark_results(results)

  # Generate basic report
  generate_thesis_report(results, output_dir = "quick_benchmark_report",
                        report_title = "Quick Benchmark Test Report")

  results
}

# --- Comprehensive Thesis Benchmark ---

#' Thesis-Quality Benchmark Setup
#'
#' Creates a comprehensive benchmark configuration suitable for thesis research
#'
#' @param n_replications Number of Monte Carlo replications (default: 100)
#' @param include_additional_methods Whether to include additional comparison methods
#' @return BenchmarkConfig object
#' @export
create_thesis_benchmark_config <- function(n_replications = 100,
                                         include_additional_methods = TRUE) {

  config <- create_default_benchmark_config(
    n_replications = n_replications,
    output_dir = "thesis_benchmark_results"
  )

  if (include_additional_methods) {
    # Add more sophisticated comparison methods
    additional_methods <- list(
      # Advanced MICE methods
      mice_norm = list(
        name = "MICE Norm",
        type = "mice_standard",
        params = list(method = "norm")
      ),
      mice_norm_predict = list(
        name = "MICE Norm Predict",
        type = "mice_standard",
        params = list(method = "norm.predict")
      ),

      # More VIM methods
      vim_hotdeck = list(
        name = "VIM Hot Deck",
        type = "vim",
        params = list(method = "hotdeck")
      ),

      # Custom distance-based methods
      gowerpmm_robust = list(
        name = "Gower-PMM (Robust)",
        type = "gowerpmm",
        params = list(weights = "auto", k = 3, scaling = "iqr")
      )
    )

    config@methods <- c(config@methods, additional_methods)
  }

  # Add more challenging scenarios
  additional_scenarios <- list(
    # Very large datasets
    large_sparse = list(n = 2000, p_num = 8, p_cat = 4, p_ord = 2,
                       correlation = "low", noise_level = "medium"),

    # High-dimensional scenarios
    high_dim_mixed = list(n = 300, p_num = 15, p_cat = 8, p_ord = 5,
                         correlation = "moderate", noise_level = "high"),

    # Extreme correlation scenarios
    perfect_corr = list(n = 400, p_num = 4, p_cat = 2, p_ord = 1,
                       correlation = "very_high", noise_level = "low"),

    # Complex missing patterns
    irregular_missing = list(n = 500, p_num = 6, p_cat = 3, p_ord = 2,
                           correlation = "mixed", noise_level = "medium")
  )

  config@scenarios <- c(config@scenarios, additional_scenarios)

  # Add more missing data patterns
  additional_patterns <- list(
    # Different missing rates
    mcar_5 = list(type = "MCAR", rate = 0.05, description = "5% MCAR"),
    mcar_40 = list(type = "MCAR", rate = 0.40, description = "40% MCAR"),

    # More complex MAR patterns
    mar_5 = list(type = "MAR", rate = 0.05, description = "5% MAR"),
    mar_40 = list(type = "MAR", rate = 0.40, description = "40% MAR"),

    # MNAR patterns
    mnar_5 = list(type = "MNAR", rate = 0.05, description = "5% MNAR"),
    mnar_40 = list(type = "MNAR", rate = 0.40, description = "40% MNAR")
  )

  config@missing_patterns <- c(config@missing_patterns, additional_patterns)

  config
}

#' Run Comprehensive Thesis Benchmark
#'
#' Executes a full thesis-quality benchmark analysis
#'
#' @param n_replications Number of replications
#' @param parallel Whether to use parallel processing
#' @param n_cores Number of cores (NULL for auto-detection)
#' @return BenchmarkResults object
#' @export
run_thesis_benchmark <- function(n_replications = 100, parallel = TRUE, n_cores = NULL) {

  message("Setting up comprehensive thesis benchmark...")
  config <- create_thesis_benchmark_config(n_replications)

  message(sprintf("Benchmark configuration:"))
  message(sprintf("- Methods: %d", length(config@methods)))
  message(sprintf("- Scenarios: %d", length(config@scenarios)))
  message(sprintf("- Missing patterns: %d", length(config@missing_patterns)))
  message(sprintf("- Total runs: %d",
                  length(config@methods) * length(config@scenarios) *
                  length(config@missing_patterns) * config@n_replications))

  # Run the benchmark
  results <- run_benchmark_suite(config, parallel = parallel, n_cores = n_cores)

  # Save results
  save_benchmark_results(results)

  # Generate comprehensive report
  generate_thesis_report(results, output_dir = "thesis_full_report",
                        report_title = "Comprehensive Gower-PMM Benchmark Analysis")

  results
}

# --- Custom Benchmark Examples ---

#' Distance Measure Comparison Benchmark
#'
#' Focuses specifically on comparing different distance measures for imputation
#'
#' @param n_replications Number of replications
#' @return BenchmarkResults object
#' @export
run_distance_comparison_benchmark <- function(n_replications = 50) {

  # Create focused configuration for distance measure comparison
  config <- create_default_benchmark_config(n_replications, "distance_comparison_results")

  # Focus on distance-based methods
  config@methods <- config@methods[c("gowerpmm_auto", "gowerpmm_equal",
                                   "gowerpmm_range", "gowerpmm_iqr",
                                   "fd_gowdis")]

  # Use a diverse set of scenarios
  config@scenarios <- config@scenarios[c("small_mixed", "medium_complex",
                                       "large_complex", "highly_correlated")]

  # Focus on moderate missing rates
  config@missing_patterns <- config@missing_patterns[c("mcar_20", "mar_20", "mnar_20")]

  message("Running distance measure comparison benchmark...")

  results <- run_benchmark_suite(config, parallel = TRUE)

  # Additional distance-specific analysis
  distance_analysis <- analyze_distance_measures(results)

  # Save extended results
  saveRDS(distance_analysis, file.path(config@output_dir, "distance_analysis.rds"))

  results
}

#' Analyze Distance Measure Performance
#' @keywords internal
analyze_distance_measures <- function(results) {
  df <- results@results
  df <- df[df$success, ]

  # Extract distance quality metrics
  distance_metrics <- df %>%
    dplyr::filter(method %in% c("gowerpmm_auto", "gowerpmm_equal",
                               "gowerpmm_range", "gowerpmm_iqr", "fd_gowdis")) %>%
    dplyr::group_by(method, scenario) %>%
    dplyr::summarise(
      mean_balance_score = mean(balance_score, na.rm = TRUE),
      mean_correlation = mean(mean_correlation, na.rm = TRUE),
      mean_nn_preservation = mean(nn_preservation_5, na.rm = TRUE),
      .groups = "drop"
    )

  distance_metrics
}

#' Computational Efficiency Benchmark
#'
#' Compares methods specifically on computational performance
#'
#' @param n_replications Number of replications
#' @return BenchmarkResults object
#' @export
run_efficiency_benchmark <- function(n_replications = 30) {

  config <- create_default_benchmark_config(n_replications, "efficiency_results")

  # Include all methods for efficiency comparison
  # Use larger datasets to see performance differences
  config@scenarios <- list(
    large_dataset = list(n = 1000, p_num = 10, p_cat = 5, p_ord = 3,
                        correlation = "moderate", noise_level = "medium"),
    very_large = list(n = 2000, p_num = 8, p_cat = 4, p_ord = 2,
                     correlation = "low", noise_level = "low")
  )

  config@missing_patterns <- config@missing_patterns["mcar_20"]

  message("Running computational efficiency benchmark...")

  results <- run_benchmark_suite(config, parallel = TRUE)

  # Generate efficiency-focused report
  efficiency_report <- analyze_efficiency(results)

  # Save efficiency analysis
  saveRDS(efficiency_report, file.path(config@output_dir, "efficiency_analysis.rds"))

  results
}

#' Analyze Computational Efficiency
#' @keywords internal
analyze_efficiency <- function(results) {
  df <- results@results
  df <- df[df$success, ]

  efficiency_summary <- df %>%
    dplyr::group_by(method, scenario) %>%
    dplyr::summarise(
      mean_time = mean(execution_time, na.rm = TRUE),
      sd_time = sd(execution_time, na.rm = TRUE),
      median_time = median(execution_time, na.rm = TRUE),
      min_time = min(execution_time, na.rm = TRUE),
      max_time = max(execution_time, na.rm = TRUE),
      mean_efficiency = mean(efficiency_score, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::arrange(scenario, mean_time)

  efficiency_summary
}

# --- Validation and Sensitivity Analysis ---

#' Validation Benchmark
#'
#' Tests benchmark framework validation and reproducibility
#'
#' @param n_replications Number of replications for validation
#' @return List with validation results
#' @export
run_validation_benchmark <- function(n_replications = 20) {

  message("Running benchmark validation...")

  # Run the same benchmark multiple times to check reproducibility
  validation_runs <- list()

  for (i in 1:3) {
    message(sprintf("Validation run %d/3...", i))

    config <- create_default_benchmark_config(n_replications, "validation_results")
    # Use minimal setup for validation
    config@methods <- config@methods[c("gowerpmm_auto", "mice_pmm")]
    config@scenarios <- config@scenarios["small_balanced"]
    config@missing_patterns <- config@missing_patterns["mcar_20"]

    results <- run_benchmark_suite(config, parallel = FALSE)
    validation_runs[[i]] <- results@results
  }

  # Analyze reproducibility
  validation_analysis <- analyze_validation(validation_runs)

  saveRDS(validation_analysis, "validation_results/validation_analysis.rds")

  validation_analysis
}

#' Analyze Validation Results
#' @keywords internal
analyze_validation <- function(validation_runs) {

  # Combine all runs
  all_results <- do.call(rbind, validation_runs)
  all_results$run_id <- rep(1:length(validation_runs),
                           each = nrow(validation_runs[[1]]))

  # Calculate reproducibility metrics
  reproducibility <- all_results %>%
    dplyr::group_by(method, scenario, missing_pattern) %>%
    dplyr::summarise(
      mean_rmse = mean(rmse, na.rm = TRUE),
      sd_rmse = sd(rmse, na.rm = TRUE),
      cv_rmse = sd_rmse / mean_rmse,  # Coefficient of variation
      mean_time = mean(execution_time, na.rm = TRUE),
      sd_time = sd(execution_time, na.rm = TRUE),
      cv_time = sd_time / mean_time,
      .groups = "drop"
    )

  list(
    reproducibility_metrics = reproducibility,
    all_results = all_results,
    summary = list(
      mean_cv_rmse = mean(reproducibility$cv_rmse, na.rm = TRUE),
      mean_cv_time = mean(reproducibility$cv_time, na.rm = TRUE),
      validation_runs = length(validation_runs)
    )
  )
}

#' Sensitivity Analysis Benchmark
#'
#' Tests how results change with different benchmark parameters
#'
#' @return List with sensitivity analysis results
#' @export
run_sensitivity_analysis <- function() {

  message("Running sensitivity analysis...")

  base_config <- create_default_benchmark_config(20, "sensitivity_results")
  base_config@methods <- base_config@methods[c("gowerpmm_auto", "mice_pmm")]
  base_config@scenarios <- base_config@scenarios["medium_balanced"]
  base_config@missing_patterns <- base_config@missing_patterns["mcar_20"]

  sensitivity_results <- list()

  # Test different replication numbers
  for (n_rep in c(10, 20, 50)) {
    message(sprintf("Testing with %d replications...", n_rep))
    config <- base_config
    config@n_replications <- n_rep

    results <- run_benchmark_suite(config, parallel = FALSE)
    sensitivity_results[[sprintf("replications_%d", n_rep)]] <- list(
      config = config,
      results = results@results,
      summary = results@summary_stats
    )
  }

  # Test different missing rates
  for (rate in c(0.1, 0.2, 0.3)) {
    message(sprintf("Testing with %.1f missing rate...", rate))
    config <- base_config
    config@missing_patterns <- list(
      test_rate = list(type = "MCAR", rate = rate,
                      description = sprintf("%.0f%% MCAR", rate * 100))
    )

    results <- run_benchmark_suite(config, parallel = FALSE)
    sensitivity_results[[sprintf("missing_rate_%.1f", rate)]] <- list(
      config = config,
      results = results@results,
      summary = results@summary_stats
    )
  }

  # Analyze sensitivity
  sensitivity_analysis <- analyze_sensitivity(sensitivity_results)

  saveRDS(sensitivity_analysis, "sensitivity_results/sensitivity_analysis.rds")

  sensitivity_analysis
}

#' Analyze Sensitivity Results
#' @keywords internal
analyze_sensitivity <- function(sensitivity_results) {

  # Extract key metrics for comparison
  sensitivity_summary <- lapply(names(sensitivity_results), function(param_set) {
    result <- sensitivity_results[[param_set]]

    summary_df <- result$results %>%
      dplyr::filter(success) %>%
      dplyr::group_by(method) %>%
      dplyr::summarise(
        mean_rmse = mean(rmse, na.rm = TRUE),
        sd_rmse = sd(rmse, na.rm = TRUE),
        mean_time = mean(execution_time, na.rm = TRUE),
        parameter_set = param_set,
        .groups = "drop"
      )

    summary_df
  })

  sensitivity_df <- do.call(rbind, sensitivity_summary)

  # Calculate sensitivity metrics
  sensitivity_metrics <- sensitivity_df %>%
    dplyr::group_by(method) %>%
    dplyr::summarise(
      rmse_range = max(mean_rmse) - min(mean_rmse),
      rmse_cv = sd(mean_rmse) / mean(mean_rmse),
      time_range = max(mean_time) - min(mean_time),
      time_cv = sd(mean_time) / mean(mean_time),
      .groups = "drop"
    )

  list(
    detailed_results = sensitivity_df,
    sensitivity_metrics = sensitivity_metrics,
    parameter_sets = names(sensitivity_results)
  )
}

# --- Utility Functions ---

#' Create Custom Benchmark Configuration
#'
#' Helper function to create custom benchmark configurations
#'
#' @param methods List of methods to include
#' @param scenarios List of scenarios to test
#' @param missing_patterns List of missing patterns
#' @param n_replications Number of replications
#' @param output_dir Output directory
#' @return BenchmarkConfig object
#' @export
create_custom_benchmark <- function(methods = NULL, scenarios = NULL,
                                   missing_patterns = NULL, n_replications = 50,
                                   output_dir = "custom_benchmark_results") {

  config <- create_default_benchmark_config(n_replications, output_dir)

  if (!is.null(methods)) {
    config@methods <- methods
  }

  if (!is.null(scenarios)) {
    config@scenarios <- scenarios
  }

  if (!is.null(missing_patterns)) {
    config@missing_patterns <- missing_patterns
  }

  config
}

#' Print Benchmark Summary
#'
#' Prints a human-readable summary of benchmark results
#'
#' @param results BenchmarkResults object
#' @export
print_benchmark_summary <- function(results) {

  cat("=== Benchmark Results Summary ===\n\n")

  cat("Configuration:\n")
  cat(sprintf("- Methods: %d (%s)\n", length(results@config@methods),
              paste(names(results@config@methods), collapse = ", ")))
  cat(sprintf("- Scenarios: %d (%s)\n", length(results@config@scenarios),
              paste(names(results@config@scenarios), collapse = ", ")))
  cat(sprintf("- Missing Patterns: %d (%s)\n", length(results@config@missing_patterns),
              paste(names(results@config@missing_patterns), collapse = ", ")))
  cat(sprintf("- Replications: %d\n", results@config@n_replications))
  cat(sprintf("- Total Runs: %d\n", results@summary_stats$overall$total_runs))

  cat("\nExecution Summary:\n")
  cat(sprintf("- Successful Runs: %d (%.1f%%)\n",
              results@summary_stats$overall$successful_runs,
              results@summary_stats$overall$success_rate * 100))
  cat(sprintf("- Total Execution Time: %.2f seconds\n",
              results@execution_time))
  cat(sprintf("- Average Time per Run: %.3f seconds\n",
              results@summary_stats$overall$mean_execution_time))

  cat("\nTop Performing Methods (by RMSE):\n")
  method_summary <- results@summary_stats$by_method %>%
    dplyr::arrange(mean_rmse) %>%
    head(5)

  for (i in 1:nrow(method_summary)) {
    cat(sprintf("%d. %s: RMSE = %.4f ± %.4f, Time = %.3f ± %.3f sec\n",
                i, method_summary$method[i], method_summary$mean_rmse[i],
                method_summary$sd_rmse[i], method_summary$mean_time[i],
                method_summary$sd_time[i]))
  }

  cat("\nResults saved to:", results@config@output_dir, "\n")
}
#' Benchmark Results Analysis and Visualization
#'
#' This module provides comprehensive analysis and visualization functions
#' for benchmark results, designed for thesis-quality reporting.

# --- Results Analysis Functions ---

#' Create Thesis-Ready Summary Tables
#'
#' Generates publication-quality summary tables comparing methods across scenarios
#'
#' @param results BenchmarkResults object
#' @param metric Primary metric for comparison ("rmse", "mae", "efficiency_score")
#' @param by_scenario Whether to create separate tables by scenario
#' @return List of formatted tables
#' @export
create_summary_tables <- function(results, metric = "rmse", by_scenario = TRUE) {

  df <- results@results

  # Filter to successful runs only
  df <- df[df$success, ]

  if (nrow(df) == 0) {
    warning("No successful runs found in results")
    return(NULL)
  }

  tables <- list()

  if (by_scenario) {
    # Create separate table for each scenario
    scenarios <- unique(df$scenario)

    for (scenario in scenarios) {
      scenario_data <- df[df$scenario == scenario, ]

      table_data <- scenario_data %>%
        dplyr::group_by(method, missing_pattern) %>%
        dplyr::summarise(
          mean_metric = mean(.data[[metric]], na.rm = TRUE),
          sd_metric = sd(.data[[metric]], na.rm = TRUE),
          n_replications = dplyr::n(),
          mean_time = mean(execution_time, na.rm = TRUE),
          .groups = "drop"
        ) %>%
        dplyr::mutate(
          metric_formatted = sprintf("%.4f (%.4f)", mean_metric, sd_metric),
          time_formatted = sprintf("%.3f", mean_time)
        ) %>%
        dplyr::select(method, missing_pattern, metric_formatted, time_formatted, n_replications)

      # Reshape to wide format
      wide_table <- table_data %>%
        tidyr::pivot_wider(
          names_from = missing_pattern,
          values_from = c(metric_formatted, time_formatted),
          names_glue = "{missing_pattern}_{.value}"
        )

      tables[[scenario]] <- wide_table
    }
  } else {
    # Create overall summary table
    overall_table <- df %>%
      dplyr::group_by(method) %>%
      dplyr::summarise(
        mean_metric = mean(.data[[metric]], na.rm = TRUE),
        sd_metric = sd(.data[[metric]], na.rm = TRUE),
        mean_time = mean(execution_time, na.rm = TRUE),
        sd_time = sd(execution_time, na.rm = TRUE),
        n_total = dplyr::n(),
        success_rate = mean(success),
        .groups = "drop"
      ) %>%
      dplyr::arrange(mean_metric)

    tables$overall <- overall_table
  }

  tables
}

#' Perform Statistical Comparisons
#'
#' Conducts statistical tests to compare method performance
#'
#' @param results BenchmarkResults object
#' @param metric Metric to compare
#' @param group_by Variables to group comparisons by
#' @return List of statistical test results
#' @export
perform_statistical_tests <- function(results, metric = "rmse",
                                    group_by = c("scenario", "missing_pattern")) {

  df <- results@results
  df <- df[df$success, ]

  if (nrow(df) == 0) {
    warning("No successful runs found in results")
    return(NULL)
  }

  # Group data for comparisons
  grouped_data <- df %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(group_by)))

  test_results <- list()

  # Perform pairwise t-tests for each group
  group_combinations <- grouped_data %>% dplyr::group_keys()

  for (i in 1:nrow(group_combinations)) {
    group_key <- group_combinations[i, ]
    group_data <- grouped_data %>% dplyr::group_split() %>% .[[i]]

    if (length(unique(group_data$method)) < 2) next

    # Pairwise comparisons
    methods <- unique(group_data$method)
    comparisons <- combn(methods, 2, simplify = FALSE)

    group_tests <- lapply(comparisons, function(comp) {
      method1_data <- group_data[group_data$method == comp[1], metric]
      method2_data <- group_data[group_data$method == comp[2], metric]

      if (length(method1_data) < 3 || length(method2_data) < 3) {
        return(list(
          method1 = comp[1],
          method2 = comp[2],
          test = "insufficient_data",
          p_value = NA,
          mean_diff = NA
        ))
      }

      # Perform t-test
      t_test <- tryCatch({
        t.test(method1_data, method2_data, alternative = "two.sided")
      }, error = function(e) {
        list(p.value = NA, estimate = c(NA, NA))
      })

      list(
        method1 = comp[1],
        method2 = comp[2],
        test = "t_test",
        p_value = t_test$p.value,
        mean_diff = diff(t_test$estimate),
        method1_mean = t_test$estimate[1],
        method2_mean = t_test$estimate[2]
      )
    })

    group_name <- paste(as.character(group_key), collapse = "_")
    test_results[[group_name]] <- list(
      group = group_key,
      comparisons = group_tests
    )
  }

  test_results
}

#' Create Performance Rankings
#'
#' Ranks methods by performance across different metrics and scenarios
#'
#' @param results BenchmarkResults object
#' @param metrics Vector of metrics to rank by
#' @param weights Optional weights for combining multiple metrics
#' @return Ranked performance table
#' @export
create_performance_rankings <- function(results, metrics = c("rmse", "mae", "execution_time"),
                                       weights = NULL) {

  df <- results@results
  df <- df[df$success, ]

  if (is.null(weights)) {
    weights <- rep(1/length(metrics), length(metrics))
  }

  # Normalize each metric (lower is better for all metrics)
  df_normalized <- df
  for (metric in metrics) {
    if (metric %in% colnames(df)) {
      metric_values <- df[[metric]]
      # Min-max normalization
      min_val <- min(metric_values, na.rm = TRUE)
      max_val <- max(metric_values, na.rm = TRUE)
      if (max_val > min_val) {
        df_normalized[[paste0(metric, "_norm")]] <- (metric_values - min_val) / (max_val - min_val)
      } else {
        df_normalized[[paste0(metric, "_norm")]] <- 0
      }
    }
  }

  # Create composite score
  norm_cols <- paste0(metrics, "_norm")
  df_normalized$composite_score <- rowSums(
    sapply(norm_cols, function(col) df_normalized[[col]] * weights[match(gsub("_norm$", "", col), metrics)])
  )

  # Rank by composite score within each scenario and missing pattern
  rankings <- df_normalized %>%
    dplyr::group_by(scenario, missing_pattern) %>%
    dplyr::mutate(
      rank = rank(composite_score),
      best_method = method[which.min(composite_score)]
    ) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(method, scenario, missing_pattern) %>%
    dplyr::summarise(
      mean_composite = mean(composite_score, na.rm = TRUE),
      mean_rank = mean(rank, na.rm = TRUE),
      n_wins = sum(rank == 1),
      .groups = "drop"
    ) %>%
    dplyr::arrange(scenario, missing_pattern, mean_rank)

  rankings
}

# --- Visualization Functions ---

#' Create Performance Comparison Plot
#'
#' Generates publication-quality plots comparing method performance
#'
#' @param results BenchmarkResults object
#' @param metric Metric to plot
#' @param facet_by Variable to facet by ("scenario" or "missing_pattern")
#' @param plot_type Type of plot ("boxplot", "violin", "bar")
#' @return ggplot object
#' @export
create_performance_plot <- function(results, metric = "rmse",
                                   facet_by = "scenario", plot_type = "boxplot") {

  df <- results@results
  df <- df[df$success, ]

  if (!(metric %in% colnames(df))) {
    stop("Metric '", metric, "' not found in results")
  }

  # Create the base plot
  p <- ggplot2::ggplot(df, ggplot2::aes_string(x = "method", y = metric, fill = "method")) +
    ggplot2::labs(
      title = paste("Method Comparison:", toupper(metric)),
      x = "Method",
      y = toupper(metric)
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )

  # Add plot type
  if (plot_type == "boxplot") {
    p <- p + ggplot2::geom_boxplot(alpha = 0.7, outlier.size = 0.5)
  } else if (plot_type == "violin") {
    p <- p + ggplot2::geom_violin(alpha = 0.7) +
      ggplot2::geom_boxplot(width = 0.1, fill = "white", alpha = 0.5)
  } else if (plot_type == "bar") {
    # Calculate means for bar plot
    means_df <- df %>%
      dplyr::group_by(method, !!rlang::sym(facet_by)) %>%
      dplyr::summarise(mean_val = mean(.data[[metric]], na.rm = TRUE),
                      se_val = sd(.data[[metric]], na.rm = TRUE) / sqrt(dplyr::n()),
                      .groups = "drop")

    p <- ggplot2::ggplot(means_df, ggplot2::aes_string(x = "method", y = "mean_val", fill = "method")) +
      ggplot2::geom_bar(stat = "identity", alpha = 0.7) +
      ggplot2::geom_errorbar(ggplot2::aes(ymin = mean_val - se_val, ymax = mean_val + se_val),
                            width = 0.2)
  }

  # Add faceting
  if (facet_by == "scenario") {
    p <- p + ggplot2::facet_wrap(~ scenario, scales = "free_y")
  } else if (facet_by == "missing_pattern") {
    p <- p + ggplot2::facet_wrap(~ missing_pattern, scales = "free_y")
  }

  p
}

#' Create Efficiency Frontier Plot
#'
#' Plots methods on an efficiency frontier (quality vs. speed)
#'
#' @param results BenchmarkResults object
#' @param quality_metric Quality metric (default: "rmse")
#' @param time_metric Time metric (default: "execution_time")
#' @return ggplot object
#' @export
create_efficiency_frontier <- function(results, quality_metric = "rmse",
                                      time_metric = "execution_time") {

  df <- results@results
  df <- df[df$success, ]

  # Calculate mean performance for each method
  efficiency_data <- df %>%
    dplyr::group_by(method) %>%
    dplyr::summarise(
      mean_quality = mean(.data[[quality_metric]], na.rm = TRUE),
      mean_time = mean(.data[[time_metric]], na.rm = TRUE),
      sd_quality = sd(.data[[quality_metric]], na.rm = TRUE),
      sd_time = sd(.data[[time_metric]], na.rm = TRUE),
      .groups = "drop"
    )

  # Create efficiency frontier
  p <- ggplot2::ggplot(efficiency_data,
                      ggplot2::aes(x = mean_time, y = mean_quality, label = method)) +
    ggplot2::geom_point(size = 3, alpha = 0.7) +
    ggplot2::geom_text(vjust = -1, size = 3) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = mean_quality - sd_quality,
                                       ymax = mean_quality + sd_quality),
                          alpha = 0.3) +
    ggplot2::geom_errorbarh(ggplot2::aes(xmin = mean_time - sd_time,
                                        xmax = mean_time + sd_time),
                           alpha = 0.3) +
    ggplot2::labs(
      title = "Efficiency Frontier: Quality vs. Speed",
      x = paste("Mean", toupper(time_metric), "(seconds)"),
      y = paste("Mean", toupper(quality_metric))
    ) +
    ggplot2::theme_minimal() +
    ggplot2::scale_x_log10() +
    ggplot2::scale_y_log10()

  # Add Pareto frontier line (simplified)
  # Methods on the frontier have no other method that is better in both dimensions
  efficiency_data <- efficiency_data[order(efficiency_data$mean_time), ]
  frontier <- efficiency_data[1, ]  # Start with fastest method

  for (i in 2:nrow(efficiency_data)) {
    if (efficiency_data$mean_quality[i] < min(frontier$mean_quality)) {
      frontier <- rbind(frontier, efficiency_data[i, ])
    }
  }

  p <- p + ggplot2::geom_line(data = frontier,
                             ggplot2::aes(x = mean_time, y = mean_quality),
                             linetype = "dashed", color = "red", alpha = 0.7)

  p
}

#' Create Scenario Comparison Heatmap
#'
#' Creates a heatmap showing method performance across scenarios and missing patterns
#'
#' @param results BenchmarkResults object
#' @param metric Metric to visualize
#' @return ggplot object
#' @export
create_scenario_heatmap <- function(results, metric = "rmse") {

  df <- results@results
  df <- df[df$success, ]

  # Calculate mean performance for each combination
  heatmap_data <- df %>%
    dplyr::group_by(method, scenario, missing_pattern) %>%
    dplyr::summarise(
      mean_metric = mean(.data[[metric]], na.rm = TRUE),
      .groups = "drop"
    )

  # Create heatmap
  p <- ggplot2::ggplot(heatmap_data,
                      ggplot2::aes(x = scenario, y = method, fill = mean_metric)) +
    ggplot2::geom_tile(color = "white") +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.3f", mean_metric)),
                      color = "black", size = 3) +
    ggplot2::facet_wrap(~ missing_pattern) +
    ggplot2::scale_fill_gradient(low = "green", high = "red",
                                name = toupper(metric)) +
    ggplot2::labs(
      title = paste("Method Performance Heatmap:", toupper(metric)),
      x = "Scenario",
      y = "Method"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      axis.text.y = ggplot2::element_text(size = 8)
    )

  p
}

#' Create Convergence Analysis Plot
#'
#' Shows how method rankings change with increasing replications
#'
#' @param results BenchmarkResults object
#' @param metric Metric for ranking
#' @param max_replications Maximum replications to analyze
#' @return ggplot object
#' @export
create_convergence_plot <- function(results, metric = "rmse", max_replications = NULL) {

  df <- results@results
  df <- df[df$success, ]

  if (is.null(max_replications)) {
    max_replications <- max(df$replication_id)
  }

  convergence_data <- list()

  # Calculate rankings for increasing numbers of replications
  for (n_rep in seq(10, max_replications, by = 10)) {
    subset_data <- df[df$replication_id <= n_rep, ]

    rankings <- subset_data %>%
      dplyr::group_by(method) %>%
      dplyr::summarise(
        mean_metric = mean(.data[[metric]], na.rm = TRUE),
        n_replications = n_rep,
        .groups = "drop"
      ) %>%
      dplyr::arrange(mean_metric) %>%
      dplyr::mutate(rank = 1:dplyr::n())

    convergence_data[[as.character(n_rep)]] <- rankings
  }

  convergence_df <- do.call(rbind, convergence_data)

  # Create plot showing rank stability
  p <- ggplot2::ggplot(convergence_df,
                      ggplot2::aes(x = n_replications, y = rank,
                                  color = method, group = method)) +
    ggplot2::geom_line(alpha = 0.7) +
    ggplot2::geom_point(alpha = 0.7) +
    ggplot2::labs(
      title = "Ranking Convergence Analysis",
      subtitle = paste("Metric:", toupper(metric)),
      x = "Number of Replications",
      y = "Method Rank (1 = best)"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::scale_y_reverse()  # Lower rank is better

  p
}

# --- Report Generation ---

#' Generate Thesis-Ready Report
#'
#' Creates a comprehensive report with tables, figures, and statistical analysis
#'
#' @param results BenchmarkResults object
#' @param output_dir Directory to save report
#' @param report_title Title for the report
#' @return Path to generated report
#' @export
generate_thesis_report <- function(results, output_dir = "thesis_report",
                                  report_title = "Benchmark Analysis Report") {

  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  # Generate key tables and plots
  summary_tables <- create_summary_tables(results)
  performance_rankings <- create_performance_rankings(results)
  statistical_tests <- perform_statistical_tests(results)

  # Create plots
  rmse_plot <- create_performance_plot(results, "rmse")
  efficiency_plot <- create_efficiency_frontier(results)
  heatmap_plot <- create_scenario_heatmap(results)

  # Save plots
  ggplot2::ggsave(file.path(output_dir, "rmse_comparison.png"), rmse_plot,
                  width = 10, height = 6, dpi = 300)
  ggplot2::ggsave(file.path(output_dir, "efficiency_frontier.png"), efficiency_plot,
                  width = 8, height = 6, dpi = 300)
  ggplot2::ggsave(file.path(output_dir, "performance_heatmap.png"), heatmap_plot,
                  width = 12, height = 8, dpi = 300)

  # Generate R Markdown report
  rmd_content <- create_rmd_report(results, summary_tables, performance_rankings,
                                  statistical_tests, report_title)

  rmd_file <- file.path(output_dir, "benchmark_report.Rmd")
  writeLines(rmd_content, rmd_file)

  # Try to render the report
  tryCatch({
    rmarkdown::render(rmd_file, output_dir = output_dir)
    message("Report generated successfully: ", file.path(output_dir, "benchmark_report.html"))
  }, error = function(e) {
    warning("Could not render R Markdown report: ", e$message)
    message("R Markdown file saved to: ", rmd_file)
  })

  file.path(output_dir, "benchmark_report.html")
}

#' Create R Markdown Report Content
#' @keywords internal
create_rmd_report <- function(results, summary_tables, performance_rankings,
                             statistical_tests, title) {

  rmd <- sprintf('---
title: "%s"
author: "Gower-PMM Benchmark Analysis"
date: "`r Sys.Date()`"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = FALSE, warning = FALSE, message = FALSE)
library(ggplot2)
library(dplyr)
library(tidyr)
library(kableExtra)
```

# Executive Summary

This report presents a comprehensive benchmark analysis comparing the Gower-PMM imputation method against several alternative approaches for handling missing data in mixed-type datasets.

## Benchmark Overview

- **Total Runs**: %d
- **Successful Runs**: %d (%.1f%%)
- **Methods Compared**: %s
- **Scenarios Tested**: %s
- **Missing Patterns**: %s

## Key Findings

### Performance Rankings

```{r rankings}
kable(performance_rankings, caption = "Method Performance Rankings") %>%
  kable_styling(bootstrap_options = c("striped", "hover"))
```

### Summary Tables

#### Overall Performance
```{r overall-table}
# Overall summary table would go here
```

### Visual Comparisons

#### RMSE Comparison Across Scenarios
![RMSE Comparison](rmse_comparison.png)

#### Efficiency Frontier
![Efficiency Frontier](efficiency_frontier.png)

#### Performance Heatmap
![Performance Heatmap](performance_heatmap.png)

## Detailed Analysis

### Statistical Comparisons

```{r stats}
# Statistical test results would be summarized here
```

### Method Characteristics

- **Gower-PMM (Auto)**: Optimized weights using genetic algorithm
- **Gower-PMM (Equal)**: Equal weights for all variables
- **FD::gowdis + PMM**: Traditional Gower distance with PMM
- **MICE PMM**: Standard MICE predictive mean matching
- **MICE CART**: Classification and regression trees
- **MICE RF**: Random forest imputation
- **VIM k-NN**: k-nearest neighbors imputation
- **VIM IRMI**: Iterative robust model-based imputation

## Conclusions

[Conclusions would be written here based on the analysis]

## Technical Details

- **R Version**: `r R.version.string`
- **Benchmark Run Date**: `r format(results@timestamp, "%%Y-%%m-%%d %%H:%%M")`
- **Total Execution Time**: `r round(results@execution_time, 2)` seconds
',
                 title,
                 results@summary_stats$overall$total_runs,
                 results@summary_stats$overall$successful_runs,
                 results@summary_stats$overall$success_rate * 100,
                 paste(names(results@config@methods), collapse = ", "),
                 paste(names(results@config@scenarios), collapse = ", "),
                 paste(names(results@config@missing_patterns), collapse = ", "))

  rmd
}
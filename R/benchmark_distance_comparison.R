#' Distance Engine Comparison Benchmark
#'
#' Specialized benchmark for comparing different dissimilarity/distance engines
#' for mixed-type data: Gower, DAISY, Gowdis, etc.

# --- Distance Engine Implementations ---

#' Compute Gower Distance using cluster::daisy
#'
#' @param data Data frame with mixed-type variables
#' @param weights Optional variable weights
#' @return Distance matrix
#' @export
compute_daisy_distance <- function(data, weights = NULL) {
  if (!requireNamespace("cluster", quietly = TRUE)) {
    stop("cluster package is required for DAISY distance")
  }

  # DAISY uses Gower's distance by default for mixed data
  daisy_result <- cluster::daisy(data, metric = "gower", stand = FALSE)

  # Convert to matrix if needed
  as.matrix(daisy_result)
}

#' Compute Gower Distance using StatMatch::gower.dist
#'
#' @param data Data frame with mixed-type variables
#' @param weights Optional variable weights
#' @return Distance matrix
#' @export
compute_statmatch_gower <- function(data, weights = NULL) {
  if (!requireNamespace("StatMatch", quietly = TRUE)) {
    stop("StatMatch package is required for StatMatch Gower distance")
  }

  # Convert to matrix for StatMatch
  data_matrix <- as.matrix(data)

  # StatMatch::gower.dist expects a matrix
  distances <- StatMatch::gower.dist(data_matrix)

  distances
}

#' Compute Distance using ade4::dist.ktab
#'
#' @param data Data frame with mixed-type variables
#' @param weights Optional variable weights
#' @return Distance matrix
#' @export
compute_ade4_distance <- function(data, weights = NULL) {
  if (!requireNamespace("ade4", quietly = TRUE)) {
    stop("ade4 package is required for ktab distance")
  }

  # Convert to ktab object
  ktab_data <- ade4::ktab.data.frame(data)
  distances <- ade4::dist.ktab(ktab_data, type = "F")

  as.matrix(distances)
}

#' Compute Euclidean Distance (for comparison)
#'
#' @param data Data frame (numeric variables only)
#' @param weights Optional variable weights
#' @return Distance matrix
#' @export
compute_euclidean_distance <- function(data, weights = NULL) {
  # Only use numeric variables for Euclidean distance
  num_cols <- sapply(data, is.numeric)
  if (!any(num_cols)) {
    warning("No numeric variables found for Euclidean distance")
    return(matrix(0, nrow(data), nrow(data)))
  }

  num_data <- data[, num_cols, drop = FALSE]
  as.matrix(stats::dist(num_data, method = "euclidean"))
}

#' Compute Manhattan Distance (for comparison)
#'
#' @param data Data frame (numeric variables only)
#' @param weights Optional variable weights
#' @return Distance matrix
#' @export
compute_manhattan_distance <- function(data, weights = NULL) {
  # Only use numeric variables for Manhattan distance
  num_cols <- sapply(data, is.numeric)
  if (!any(num_cols)) {
    warning("No numeric variables found for Manhattan distance")
    return(matrix(0, nrow(data), nrow(data)))
  }

  num_data <- data[, num_cols, drop = FALSE]
  as.matrix(stats::dist(num_data, method = "manhattan"))
}

# --- Distance Engine Comparison Framework ---

#' Create Distance Engine Comparison Configuration
#'
#' @param n_replications Number of replications
#' @param engines Engines to compare (default: all available)
#' @param output_dir Output directory
#' @return BenchmarkConfig object
#' @export
create_distance_benchmark_config <- function(n_replications = 100,
                                           engines = c("gowerpmm", "fd_gowdis", "daisy", "statmatch", "euclidean", "manhattan"),
                                           output_dir = "distance_benchmark_results") {

  # Define distance engines
  available_engines <- list(
    gowerpmm = list(
      name = "Gower-PMM Engine",
      type = "distance_engine",
      compute_function = "gower_dist_engine",
      description = "Optimized Gower distance with C++ backend"
    ),
    fd_gowdis = list(
      name = "FD::gowdis",
      type = "distance_engine",
      compute_function = "FD::gowdis",
      description = "Traditional FD package Gower distance"
    ),
    daisy = list(
      name = "cluster::daisy",
      type = "distance_engine",
      compute_function = "compute_daisy_distance",
      description = "DAISY clustering distance (Gower-based)"
    ),
    statmatch = list(
      name = "StatMatch::gower.dist",
      type = "distance_engine",
      compute_function = "compute_statmatch_gower",
      description = "StatMatch Gower distance implementation"
    ),
    ade4 = list(
      name = "ade4::dist.ktab",
      type = "distance_engine",
      compute_function = "compute_ade4_distance",
      description = "ADE4 K-table distance"
    ),
    euclidean = list(
      name = "Euclidean",
      type = "distance_engine",
      compute_function = "compute_euclidean_distance",
      description = "Euclidean distance (numeric only)"
    ),
    manhattan = list(
      name = "Manhattan",
      type = "distance_engine",
      compute_function = "compute_manhattan_distance",
      description = "Manhattan distance (numeric only)"
    )
  )

  # Filter to requested engines
  selected_engines <- available_engines[engines]

  # Create methods list (treating distance engines as "methods" for the framework)
  methods <- lapply(names(selected_engines), function(engine_name) {
    engine <- selected_engines[[engine_name]]
    list(
      name = engine$name,
      type = "distance_engine",
      engine_name = engine_name,
      compute_function = engine$compute_function,
      description = engine$description
    )
  })
  names(methods) <- names(selected_engines)

  # Use diverse scenarios for distance comparison
  scenarios <- list(
    # Balanced mixed-type datasets
    balanced_mixed = list(n = 200, p_num = 4, p_cat = 3, p_ord = 2,
                         correlation = "moderate", noise_level = "low"),
    high_cat = list(n = 200, p_num = 2, p_cat = 5, p_ord = 3,
                   correlation = "moderate", noise_level = "low"),
    high_num = list(n = 200, p_num = 6, p_cat = 2, p_ord = 1,
                   correlation = "moderate", noise_level = "low"),

    # Different correlation structures
    uncorrelated = list(n = 200, p_num = 4, p_cat = 3, p_ord = 2,
                       correlation = "low", noise_level = "low"),
    correlated = list(n = 200, p_num = 4, p_cat = 3, p_ord = 2,
                     correlation = "high", noise_level = "low"),

    # Different sizes
    small = list(n = 100, p_num = 3, p_cat = 2, p_ord = 1,
                correlation = "moderate", noise_level = "low"),
    large = list(n = 500, p_num = 5, p_cat = 4, p_ord = 3,
                correlation = "moderate", noise_level = "low")
  )

  # No missing data patterns for distance comparison (we evaluate on complete data)
  missing_patterns <- list(
    complete = list(type = "NONE", rate = 0, description = "Complete data")
  )

  new("BenchmarkConfig",
      scenarios = scenarios,
      methods = methods,
      metrics = list(),  # Will be set by distance-specific metrics
      n_replications = n_replications,
      missing_patterns = missing_patterns,
      output_dir = output_dir)
}

#' Run Distance Engine Comparison Benchmark
#'
#' @param config BenchmarkConfig object
#' @param parallel Whether to use parallel processing
#' @param n_cores Number of cores
#' @return BenchmarkResults object
#' @export
run_distance_comparison_benchmark <- function(config, parallel = TRUE, n_cores = NULL) {

  message("Running Distance Engine Comparison Benchmark...")
  message("Comparing engines: ", paste(names(config@methods), collapse = ", "))

  # Modify the execution to focus on distance computation
  results <- run_distance_engine_suite(config, parallel = parallel, n_cores = n_cores)

  # Save results
  save_benchmark_results(results)

  # Generate distance-specific analysis
  distance_analysis <- analyze_distance_engine_performance(results)

  # Save extended analysis
  saveRDS(distance_analysis, file.path(config@output_dir, "distance_engine_analysis.rds"))

  results
}

#' Run Distance Engine Suite
#' @keywords internal
run_distance_engine_suite <- function(config, parallel = TRUE, n_cores = NULL) {

  start_time <- Sys.time()

  # Prepare combinations
  combinations <- expand.grid(
    scenario = names(config@scenarios),
    method = names(config@methods),
    replication = 1:config@n_replications,
    stringsAsFactors = FALSE
  )

  total_runs <- nrow(combinations)
  message(sprintf("Starting distance comparison with %d total runs...", total_runs))

  results_list <- vector("list", total_runs)

  if (parallel && requireNamespace("parallel", quietly = TRUE)) {
    if (is.null(n_cores)) n_cores <- max(1, parallel::detectCores() - 1)

    message(sprintf("Running in parallel with %d cores...", n_cores))

    cl <- parallel::makeCluster(n_cores)
    on.exit(parallel::stopCluster(cl))

    # Export functions
    parallel::clusterExport(cl, c(
      "run_distance_engine_single", "generate_simulation_data",
      "evaluate_distance_quality", "compute_daisy_distance",
      "compute_statmatch_gower", "compute_ade4_distance",
      "compute_euclidean_distance", "compute_manhattan_distance",
      "compute_rank_correlation_balance", "compute_distance_preservation",
      "compute_variable_type_balance", "compute_efficiency_score"
    ), envir = environment())

    # Load packages
    parallel::clusterEvalQ(cl, {
      library(gowerpmm)
      if (requireNamespace("FD", quietly = TRUE)) library(FD)
      if (requireNamespace("cluster", quietly = TRUE)) library(cluster)
      if (requireNamespace("StatMatch", quietly = TRUE)) library(StatMatch)
      if (requireNamespace("ade4", quietly = TRUE)) library(ade4)
      if (requireNamespace("moments", quietly = TRUE)) library(moments)
    })

    results_list <- parallel::parLapply(cl, 1:total_runs, function(i) {
      combo <- combinations[i, ]

      scenario <- config@scenarios[[combo$scenario]]
      scenario$name <- combo$scenario

      method <- config@methods[[combo$method]]
      method$name <- combo$method

      run_distance_engine_single(scenario, method, combo$replication)
    })

  } else {
    # Sequential execution
    for (i in 1:total_runs) {
      combo <- combinations[i, ]

      scenario <- config@scenarios[[combo$scenario]]
      scenario$name <- combo$scenario

      method <- config@methods[[combo$method]]
      method$name <- combo$method

      if (i %% 10 == 0) {
        message(sprintf("Completed %d/%d runs (%.1f%%)...", i, total_runs, (i/total_runs)*100))
      }

      results_list[[i]] <- run_distance_engine_single(scenario, method, combo$replication)
    }
  }

  # Compile results
  results_df <- compile_distance_results(results_list)

  # Summary stats
  summary_stats <- compute_distance_summary(results_df, config)

  end_time <- Sys.time()
  total_execution_time <- as.numeric(difftime(end_time, start_time, units = "secs"))

  benchmark_results <- new("BenchmarkResults",
                          config = config,
                          results = results_df,
                          summary_stats = summary_stats,
                          execution_time = total_execution_time,
                          timestamp = start_time)

  message(sprintf("Distance comparison completed in %.2f seconds", total_execution_time))

  benchmark_results
}

#' Run Single Distance Engine Evaluation
#' @keywords internal
run_distance_engine_single <- function(scenario, method, replication_id) {

  set.seed(replication_id * 1000)
  start_time <- Sys.time()

  tryCatch({
    # Generate complete data (no missing values for distance comparison)
    data <- generate_simulation_data(scenario, seed = replication_id)

    # Compute distance matrix using specified engine
    distance_matrix <- compute_distance_with_engine(data, method)

    # Evaluate distance quality
    quality_metrics <- evaluate_distance_quality(data, distance_matrix, method$name)

    # Performance metrics
    end_time <- Sys.time()
    execution_time <- as.numeric(difftime(end_time, start_time, units = "secs"))

    performance_metrics <- evaluate_computational_performance(
      method_name = method$name,
      execution_time = execution_time,
      convergence_info = list(success = TRUE)
    )

    # Return results
    list(
      replication_id = replication_id,
      scenario = scenario$name,
      method = method$name,
      engine = method$engine_name,
      distance_quality = quality_metrics,
      computational_performance = performance_metrics,
      execution_time = execution_time,
      success = TRUE,
      error_message = NULL
    )

  }, error = function(e) {
    end_time <- Sys.time()
    execution_time <- as.numeric(difftime(end_time, start_time, units = "secs"))

    list(
      replication_id = replication_id,
      scenario = scenario$name,
      method = method$name,
      engine = method$engine_name,
      distance_quality = NULL,
      computational_performance = evaluate_computational_performance(
        method_name = method$name,
        execution_time = execution_time,
        convergence_info = list(error = as.character(e))
      ),
      execution_time = execution_time,
      success = FALSE,
      error_message = as.character(e)
    )
  })
}

#' Compute Distance with Specified Engine
#' @keywords internal
compute_distance_with_engine <- function(data, method) {

  compute_func <- method$compute_function

  if (compute_func == "gower_dist_engine") {
    # Our optimized engine
    distance_matrix <- gower_dist_engine(data)

  } else if (compute_func == "FD::gowdis") {
    # FD package
    if (!requireNamespace("FD", quietly = TRUE)) {
      stop("FD package required for FD::gowdis")
    }
    distance_matrix <- FD::gowdis(data)

  } else if (compute_func == "compute_daisy_distance") {
    distance_matrix <- compute_daisy_distance(data)

  } else if (compute_func == "compute_statmatch_gower") {
    distance_matrix <- compute_statmatch_gower(data)

  } else if (compute_func == "compute_ade4_distance") {
    distance_matrix <- compute_ade4_distance(data)

  } else if (compute_func == "compute_euclidean_distance") {
    distance_matrix <- compute_euclidean_distance(data)

  } else if (compute_func == "compute_manhattan_distance") {
    distance_matrix <- compute_manhattan_distance(data)

  } else {
    stop("Unknown distance computation function: ", compute_func)
  }

  # Ensure we return a matrix
  if (inherits(distance_matrix, "dist")) {
    as.matrix(distance_matrix)
  } else {
    distance_matrix
  }
}

#' Compile Distance Results
#' @keywords internal
compile_distance_results <- function(results_list) {

  df_rows <- lapply(results_list, function(result) {
    if (!result$success) {
      return(data.frame(
        replication_id = result$replication_id,
        scenario = result$scenario,
        method = result$method,
        engine = result$engine,
        success = FALSE,
        error_message = result$error_message,
        execution_time = result$execution_time,
        balance_score = NA, mean_correlation = NA, correlation_variance = NA,
        nn_preservation_1 = NA, nn_preservation_5 = NA,
        mean_distance = NA, sd_distance = NA, distance_skewness = NA, distance_kurtosis = NA,
        efficiency_score = NA,
        stringsAsFactors = FALSE
      ))
    }

    # Extract distance quality metrics
    dq <- result$distance_quality
    rank_corr <- dq$rank_correlation_balance
    dist_preserve <- dq$distance_preservation

    data.frame(
      replication_id = result$replication_id,
      scenario = result$scenario,
      method = result$method,
      engine = result$engine,
      success = TRUE,
      error_message = NA,
      execution_time = result$execution_time,

      # Distance quality metrics
      balance_score = rank_corr$balance_score %||% NA,
      mean_correlation = rank_corr$mean_correlation %||% NA,
      correlation_variance = rank_corr$correlation_variance %||% NA,
      nn_preservation_1 = dist_preserve$nearest_neighbor_preservation[1] %||% NA,
      nn_preservation_5 = dist_preserve$nearest_neighbor_preservation[5] %||% NA,
      mean_distance = dist_preserve$mean_distance %||% NA,
      sd_distance = dist_preserve$sd_distance %||% NA,
      distance_skewness = dist_preserve$distance_skewness %||% NA,
      distance_kurtosis = dist_preserve$distance_kurtosis %||% NA,

      # Performance
      efficiency_score = result$computational_performance$efficiency_score %||% NA,

      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, df_rows)
}

#' Compute Distance Summary Statistics
#' @keywords internal
compute_distance_summary <- function(results_df, config) {

  df <- results_df[results_df$success, ]

  summary_stats <- list()

  # Overall summary
  summary_stats$overall <- list(
    total_runs = nrow(results_df),
    successful_runs = sum(results_df$success),
    success_rate = mean(results_df$success),
    total_execution_time = sum(results_df$execution_time, na.rm = TRUE),
    mean_execution_time = mean(results_df$execution_time, na.rm = TRUE)
  )

  # By engine summary
  engine_summary <- df %>%
    dplyr::group_by(method, engine) %>%
    dplyr::summarise(
      n_runs = dplyr::n(),
      mean_balance_score = mean(balance_score, na.rm = TRUE),
      sd_balance_score = sd(balance_score, na.rm = TRUE),
      mean_nn_preservation = mean(nn_preservation_5, na.rm = TRUE),
      sd_nn_preservation = sd(nn_preservation_5, na.rm = TRUE),
      mean_time = mean(execution_time, na.rm = TRUE),
      sd_time = sd(execution_time, na.rm = TRUE),
      mean_efficiency = mean(efficiency_score, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::arrange(desc(mean_balance_score))

  summary_stats$by_engine <- engine_summary

  # By scenario summary
  scenario_summary <- df %>%
    dplyr::group_by(scenario) %>%
    dplyr::summarise(
      n_runs = dplyr::n(),
      mean_balance_score = mean(balance_score, na.rm = TRUE),
      mean_nn_preservation = mean(nn_preservation_5, na.rm = TRUE),
      .groups = "drop"
    )

  summary_stats$by_scenario <- scenario_summary

  summary_stats
}

#' Analyze Distance Engine Performance
#' @keywords internal
analyze_distance_engine_performance <- function(results) {

  df <- results@results
  df <- df[df$success, ]

  # Ranking by different metrics
  rankings <- list()

  # By balance score (higher is better)
  rankings$balance <- df %>%
    dplyr::group_by(scenario) %>%
    dplyr::mutate(rank = rank(-balance_score)) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(method) %>%
    dplyr::summarise(
      mean_rank = mean(rank),
      mean_balance = mean(balance_score, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::arrange(mean_rank)

  # By nearest neighbor preservation (higher is better)
  rankings$nn_preservation <- df %>%
    dplyr::group_by(scenario) %>%
    dplyr::mutate(rank = rank(-nn_preservation_5)) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(method) %>%
    dplyr::summarise(
      mean_rank = mean(rank),
      mean_nn_preserve = mean(nn_preservation_5, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::arrange(mean_rank)

  # By efficiency (higher is better)
  rankings$efficiency <- df %>%
    dplyr::group_by(scenario) %>%
    dplyr::mutate(rank = rank(-efficiency_score)) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(method) %>%
    dplyr::summarise(
      mean_rank = mean(rank),
      mean_efficiency = mean(efficiency_score, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::arrange(mean_rank)

  # Comparative analysis
  comparison <- df %>%
    dplyr::group_by(method) %>%
    dplyr::summarise(
      balance_score_mean = mean(balance_score, na.rm = TRUE),
      balance_score_sd = sd(balance_score, na.rm = TRUE),
      nn_preservation_mean = mean(nn_preservation_5, na.rm = TRUE),
      nn_preservation_sd = sd(nn_preservation_5, na.rm = TRUE),
      time_mean = mean(execution_time, na.rm = TRUE),
      time_sd = sd(execution_time, na.rm = TRUE),
      n_scenarios = length(unique(scenario)),
      .groups = "drop"
    )

  list(
    rankings = rankings,
    comparison = comparison,
    detailed_results = df
  )
}

#' Generate Distance Engine Comparison Report
#' @export
generate_distance_comparison_report <- function(results, output_dir = "distance_comparison_report") {

  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  # Load analysis
  analysis_file <- file.path(results@config@output_dir, "distance_engine_analysis.rds")
  if (file.exists(analysis_file)) {
    analysis <- readRDS(analysis_file)
  } else {
    analysis <- analyze_distance_engine_performance(results)
  }

  # Create plots
  balance_plot <- create_distance_balance_plot(results)
  efficiency_plot <- create_distance_efficiency_plot(results)
  ranking_plot <- create_distance_ranking_plot(analysis$rankings)

  # Save plots
  ggplot2::ggsave(file.path(output_dir, "distance_balance_comparison.png"), balance_plot,
                  width = 10, height = 6, dpi = 300)
  ggplot2::ggsave(file.path(output_dir, "distance_efficiency_frontier.png"), efficiency_plot,
                  width = 8, height = 6, dpi = 300)
  ggplot2::ggsave(file.path(output_dir, "distance_engine_rankings.png"), ranking_plot,
                  width = 12, height = 6, dpi = 300)

  # Generate report
  rmd_content <- create_distance_rmd_report(results, analysis)

  rmd_file <- file.path(output_dir, "distance_comparison_report.Rmd")
  writeLines(rmd_content, rmd_file)

  # Try to render
  tryCatch({
    rmarkdown::render(rmd_file, output_dir = output_dir)
    message("Distance comparison report generated successfully")
  }, error = function(e) {
    warning("Could not render R Markdown report: ", e$message)
    message("R Markdown file saved to: ", rmd_file)
  })

  file.path(output_dir, "distance_comparison_report.html")
}

#' Create Distance Balance Plot
#' @keywords internal
create_distance_balance_plot <- function(results) {

  df <- results@results
  df <- df[df$success, ]

  p <- ggplot2::ggplot(df, ggplot2::aes(x = method, y = balance_score, fill = method)) +
    ggplot2::geom_boxplot(alpha = 0.7) +
    ggplot2::facet_wrap(~ scenario) +
    ggplot2::labs(
      title = "Distance Engine Balance Score Comparison",
      subtitle = "Higher scores indicate better balance across variable types",
      x = "Distance Engine",
      y = "Balance Score"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )

  p
}

#' Create Distance Efficiency Plot
#' @keywords internal
create_distance_efficiency_plot <- function(results) {

  df <- results@results
  df <- df[df$success, ]

  efficiency_data <- df %>%
    dplyr::group_by(method) %>%
    dplyr::summarise(
      mean_balance = mean(balance_score, na.rm = TRUE),
      mean_time = mean(execution_time, na.rm = TRUE),
      .groups = "drop"
    )

  p <- ggplot2::ggplot(efficiency_data,
                      ggplot2::aes(x = mean_time, y = mean_balance, label = method)) +
    ggplot2::geom_point(size = 3, alpha = 0.7) +
    ggplot2::geom_text(vjust = -1, size = 3) +
    ggplot2::labs(
      title = "Distance Engine Efficiency Frontier",
      x = "Mean Execution Time (seconds)",
      y = "Mean Balance Score"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::scale_x_log10()

  p
}

#' Create Distance Ranking Plot
#' @keywords internal
create_distance_ranking_plot <- function(rankings) {

  # Combine rankings
  balance_ranks <- rankings$balance %>%
    dplyr::mutate(metric = "Balance Score")
  nn_ranks <- rankings$nn_preservation %>%
    dplyr::mutate(metric = "NN Preservation")

  combined_ranks <- rbind(balance_ranks, nn_ranks)

  p <- ggplot2::ggplot(combined_ranks,
                      ggplot2::aes(x = reorder(method, -mean_rank), y = mean_rank, fill = metric)) +
    ggplot2::geom_bar(stat = "identity", position = "dodge") +
    ggplot2::labs(
      title = "Distance Engine Rankings by Performance Metric",
      subtitle = "Lower rank numbers indicate better performance",
      x = "Distance Engine",
      y = "Mean Rank"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
    ggplot2::scale_y_reverse()  # Lower ranks are better

  p
}

#' Create Distance RMD Report
#' @keywords internal
create_distance_rmd_report <- function(results, analysis) {

  rmd <- sprintf('---
title: "Distance Engine Comparison Report"
author: "Gower-PMM Benchmark Analysis"
date: "`r Sys.Date()`"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = FALSE, warning = FALSE, message = FALSE)
library(ggplot2)
library(dplyr)
library(kableExtra)
```

# Distance Engine Comparison Analysis

This report compares different distance/dissimilarity engines for mixed-type data.

## Overview

- **Engines Compared**: %s
- **Scenarios Tested**: %s
- **Total Runs**: %d
- **Successful Runs**: %d (%.1f%%)

## Performance Rankings

### By Balance Score
```{r balance-rankings}
kable(analysis$rankings$balance, caption = "Ranking by Balance Score") %>%
  kable_styling(bootstrap_options = c("striped", "hover"))
```

### By Nearest Neighbor Preservation
```{r nn-rankings}
kable(analysis$rankings$nn_preservation, caption = "Ranking by NN Preservation") %>%
  kable_styling(bootstrap_options = c("striped", "hover"))
```

## Visual Comparisons

### Balance Score Comparison
![Balance Score Comparison](distance_balance_comparison.png)

### Efficiency Frontier
![Efficiency Frontier](distance_efficiency_frontier.png)

### Engine Rankings
![Engine Rankings](distance_engine_rankings.png)

## Detailed Comparison

```{r comparison-table}
kable(analysis$comparison, caption = "Detailed Engine Comparison") %>%
  kable_styling(bootstrap_options = c("striped", "hover"))
```

## Conclusions

[Conclusions would be written here based on the analysis]

## Technical Details

- **R Version**: `r R.version.string`
- **Analysis Date**: `r format(results@timestamp, "%%Y-%%m-%%d %%H:%%M")`
- **Total Execution Time**: `r round(results@execution_time, 2)` seconds
',
                 paste(names(results@config@methods), collapse = ", "),
                 paste(names(results@config@scenarios), collapse = ", "),
                 results@summary_stats$overall$total_runs,
                 results@summary_stats$overall$successful_runs,
                 results@summary_stats$overall$success_rate * 100)

  rmd
}
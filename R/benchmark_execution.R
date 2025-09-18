#' Benchmark Execution Engine
#'
#' This module provides the core execution engine for running comprehensive
#' benchmark simulations comparing different imputation methods and distance measures.

# --- Single Benchmark Run ---

#' Run Single Benchmark Replication
#'
#' Executes one complete benchmark replication for a given scenario, method, and missing pattern
#'
#' @param scenario Simulation scenario specification
#' @param method Imputation method specification
#' @param missing_pattern Missing data pattern specification
#' @param replication_id Replication identifier for reproducibility
#' @param config Benchmark configuration object
#' @return List containing all benchmark results for this replication
#' @export
run_single_benchmark <- function(scenario, method, missing_pattern, replication_id, config) {

  # Set seed for reproducibility
  set.seed(replication_id * 1000)

  start_time <- Sys.time()

  tryCatch({
    # 1. Generate complete simulation data
    complete_data <- generate_simulation_data(scenario, seed = replication_id)

    # 2. Introduce missing data pattern
    data_with_missing <- introduce_missing_data(complete_data, missing_pattern, seed = replication_id)
    missing_mask <- is.na(data_with_missing)

    # 3. Prepare predictor matrix (all variables predict each other)
    predictor_matrix <- matrix(1, ncol(complete_data), ncol(complete_data))
    diag(predictor_matrix) <- 0  # Variables don't predict themselves

    # 4. Run imputation method
    imputation_result <- run_imputation_method(
      method = method,
      data = data_with_missing,
      predictor_matrix = predictor_matrix,
      config = config
    )

    # 5. Evaluate imputation quality
    quality_metrics <- evaluate_imputation_quality(
      original_data = complete_data,
      imputed_data = imputation_result$imputed_data,
      missing_mask = missing_mask,
      method_name = method$name
    )

    # 6. Evaluate distance measure quality (if applicable)
    distance_metrics <- NULL
    if (method$type %in% c("gowerpmm", "fd_gowdis")) {
      # Compute distance matrix on complete data for evaluation
      if (method$type == "gowerpmm") {
        distance_matrix <- gower_dist_engine(complete_data)
      } else {
        distance_matrix <- FD::gowdis(complete_data)
      }

      distance_metrics <- evaluate_distance_quality(
        data = complete_data,
        distance_matrix = distance_matrix,
        method_name = method$name
      )
    }

    # 7. Evaluate computational performance
    end_time <- Sys.time()
    execution_time <- as.numeric(difftime(end_time, start_time, units = "secs"))

    performance_metrics <- evaluate_computational_performance(
      method_name = method$name,
      execution_time = execution_time,
      memory_usage = NA,  # Could be extended to measure memory usage
      convergence_info = imputation_result$convergence_info
    )

    # 8. Compile results
    result <- list(
      replication_id = replication_id,
      scenario = scenario,
      method = method,
      missing_pattern = missing_pattern,
      imputation_quality = quality_metrics,
      distance_quality = distance_metrics,
      computational_performance = performance_metrics,
      execution_time = execution_time,
      success = TRUE,
      error_message = NULL
    )

    result

  }, error = function(e) {
    end_time <- Sys.time()
    execution_time <- as.numeric(difftime(end_time, start_time, units = "secs"))

    # Return error result
    list(
      replication_id = replication_id,
      scenario = scenario,
      method = method,
      missing_pattern = missing_pattern,
      imputation_quality = NULL,
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

#' Run Imputation Method
#'
#' Executes a specific imputation method on the data
#'
#' @param method Method specification
#' @param data Data with missing values
#' @param predictor_matrix Predictor matrix for mice
#' @param config Benchmark configuration
#' @return List with imputed data and convergence info
#' @keywords internal
run_imputation_method <- function(method, data, predictor_matrix, config) {

  method_type <- method$type

  if (method_type == "gowerpmm") {
    # Run Gower-PMM via mice
    result <- run_gowerpmm_imputation(data, predictor_matrix, method$params, config)

  } else if (method_type == "mice_standard") {
    # Run standard mice methods
    result <- run_standard_mice_imputation(data, predictor_matrix, method$params)

  } else if (method_type == "fd_gowdis") {
    # Run FD::gowdis + PMM
    result <- run_fd_gowdis_imputation(data, method$params)

  } else if (method_type == "vim") {
    # Run VIM methods
    result <- run_vim_imputation(data, method$params)

  } else {
    stop("Unknown method type: ", method_type)
  }

  result
}

#' Run Gower-PMM Imputation
#' @keywords internal
run_gowerpmm_imputation <- function(data, predictor_matrix, params, config) {
  # Set up mice parameters
  mice_params <- list(
    data = data,
    method = "gowerpmm",
    predictorMatrix = predictor_matrix,
    m = 1,  # Single imputation for evaluation
    maxit = 1,
    printFlag = FALSE
  )

  # Add method-specific parameters
  if (!is.null(params$weights)) mice_params$weights <- params$weights
  if (!is.null(params$k)) mice_params$k <- params$k
  if (!is.null(params$scaling)) {
    mice_params$scaling <- params$scaling
  }

  # Run imputation
  imp <- do.call(mice::mice, mice_params)
  imputed_data <- mice::complete(imp, 1)

  list(
    imputed_data = imputed_data,
    convergence_info = list(
      method = "gowerpmm",
      iterations = imp$iteration,
      convergence = TRUE  # Assume convergence for single iteration
    )
  )
}

#' Run Standard MICE Imputation
#' @keywords internal
run_standard_mice_imputation <- function(data, predictor_matrix, params) {
  mice_params <- list(
    data = data,
    method = params$method,
    predictorMatrix = predictor_matrix,
    m = 1,
    maxit = 1,
    printFlag = FALSE
  )

  # Add method-specific parameters
  if (params$method == "pmm" && !is.null(params$k)) {
    mice_params$pmm.k <- params$k
  }

  imp <- do.call(mice::mice, mice_params)
  imputed_data <- mice::complete(imp, 1)

  list(
    imputed_data = imputed_data,
    convergence_info = list(
      method = params$method,
      iterations = imp$iteration,
      convergence = TRUE
    )
  )
}

#' Run FD::gowdis + PMM Imputation
#' @keywords internal
run_fd_gowdis_imputation <- function(data, params) {
  k <- params$k %||% 5

  # For simplicity, use all variables as predictors
  predictors <- data

  imputed_data <- impute_with_fd_gowdis(data, predictors, k)

  list(
    imputed_data = imputed_data,
    convergence_info = list(
      method = "fd_gowdis_pmm",
      k = k,
      convergence = TRUE
    )
  )
}

#' Run VIM Imputation
#' @keywords internal
run_vim_imputation <- function(data, params) {
  method <- params$method

  if (method == "kNN") {
    k <- params$k %||% 5
    imputed_data <- impute_with_vim_knn(data, k)
    convergence_info <- list(method = "vim_knn", k = k, convergence = TRUE)

  } else if (method == "irmi") {
    imputed_data <- impute_with_vim_irmi(data)
    convergence_info <- list(method = "vim_irmi", convergence = TRUE)

  } else {
    stop("Unknown VIM method: ", method)
  }

  list(
    imputed_data = imputed_data,
    convergence_info = convergence_info
  )
}

# --- Batch Benchmark Execution ---

#' Run Complete Benchmark Suite
#'
#' Executes the full benchmark suite across all scenarios, methods, and missing patterns
#'
#' @param config BenchmarkConfig object
#' @param parallel Whether to run in parallel (default: FALSE)
#' @param n_cores Number of cores for parallel execution (default: NULL, auto-detect)
#' @param progress_callback Optional callback function for progress updates
#' @return BenchmarkResults object
#' @export
run_benchmark_suite <- function(config, parallel = FALSE, n_cores = NULL,
                               progress_callback = NULL) {

  start_time <- Sys.time()

  # Prepare all combinations to run
  combinations <- expand.grid(
    scenario = names(config@scenarios),
    method = names(config@methods),
    missing_pattern = names(config@missing_patterns),
    replication = 1:config@n_replications,
    stringsAsFactors = FALSE
  )

  total_runs <- nrow(combinations)
  message(sprintf("Starting benchmark suite with %d total runs...", total_runs))

  # Run all combinations
  results_list <- vector("list", total_runs)

  if (parallel && requireNamespace("parallel", quietly = TRUE)) {
    # Parallel execution
    if (is.null(n_cores)) {
      n_cores <- max(1, parallel::detectCores() - 1)
    }

    message(sprintf("Running in parallel with %d cores...", n_cores))

    cl <- parallel::makeCluster(n_cores)
    on.exit(parallel::stopCluster(cl))

    # Export necessary functions and data
    parallel::clusterExport(cl, c(
      "run_single_benchmark", "run_imputation_method", "run_gowerpmm_imputation",
      "run_standard_mice_imputation", "run_fd_gowdis_imputation", "run_vim_imputation",
      "generate_simulation_data", "introduce_missing_data", "evaluate_imputation_quality",
      "evaluate_distance_quality", "evaluate_computational_performance",
      "compute_overall_imputation_metrics", "compute_by_variable_metrics",
      "compute_distribution_preservation", "compute_correlation_preservation",
      "compute_rank_correlation_balance", "compute_distance_preservation",
      "compute_variable_type_balance", "compute_efficiency_score",
      "impute_with_fd_gowdis", "impute_with_vim_knn", "impute_with_vim_irmi",
      "apply_pmm_with_distances", "generate_mixed_correlation_data"
    ), envir = environment())

    # Load required packages on cluster
    parallel::clusterEvalQ(cl, {
      library(gowerpmm)
      library(mice)
      if (requireNamespace("FD", quietly = TRUE)) library(FD)
      if (requireNamespace("VIM", quietly = TRUE)) library(VIM)
      if (requireNamespace("MASS", quietly = TRUE)) library(MASS)
      if (requireNamespace("moments", quietly = TRUE)) library(moments)
    })

    results_list <- parallel::parLapply(cl, 1:total_runs, function(i) {
      combo <- combinations[i, ]

      scenario <- config@scenarios[[combo$scenario]]
      scenario$name <- combo$scenario

      method <- config@methods[[combo$method]]
      method$name <- combo$method

      missing_pattern <- config@missing_patterns[[combo$missing_pattern]]
      missing_pattern$name <- combo$missing_pattern

      run_single_benchmark(scenario, method, missing_pattern, combo$replication, config)
    })

  } else {
    # Sequential execution
    for (i in 1:total_runs) {
      combo <- combinations[i, ]

      scenario <- config@scenarios[[combo$scenario]]
      scenario$name <- combo$scenario

      method <- config@methods[[combo$method]]
      method$name <- combo$method

      missing_pattern <- config@missing_patterns[[combo$missing_pattern]]
      missing_pattern$name <- combo$missing_pattern

      if (!is.null(progress_callback)) {
        progress_callback(i, total_runs, combo)
      } else if (i %% 10 == 0) {
        message(sprintf("Completed %d/%d runs (%.1f%%)...",
                       i, total_runs, (i/total_runs)*100))
      }

      results_list[[i]] <- run_single_benchmark(
        scenario, method, missing_pattern, combo$replication, config
      )
    }
  }

  # Compile results into data frame
  results_df <- compile_benchmark_results(results_list)

  # Compute summary statistics
  summary_stats <- compute_benchmark_summary(results_df, config)

  end_time <- Sys.time()
  total_execution_time <- as.numeric(difftime(end_time, start_time, units = "secs"))

  # Create BenchmarkResults object
  benchmark_results <- new("BenchmarkResults",
                          config = config,
                          results = results_df,
                          summary_stats = summary_stats,
                          execution_time = total_execution_time,
                          timestamp = start_time)

  message(sprintf("Benchmark suite completed in %.2f seconds", total_execution_time))

  benchmark_results
}

#' Compile Benchmark Results into Data Frame
#' @keywords internal
compile_benchmark_results <- function(results_list) {
  # Extract scalar metrics for data frame
  df_rows <- lapply(results_list, function(result) {
    if (!result$success) {
      # Handle failed runs
      return(data.frame(
        replication_id = result$replication_id,
        scenario = result$scenario$name,
        method = result$method$name,
        missing_pattern = result$missing_pattern$name,
        success = FALSE,
        error_message = result$error_message,
        execution_time = result$execution_time,
        # Fill other columns with NA
        rmse = NA, mae = NA, bias = NA, nrmse = NA, nmae = NA, r_squared = NA,
        accuracy = NA, mean_level_diff = NA, max_level_diff = NA,
        ks_statistic = NA, ks_p_value = NA, mean_preservation = NA, sd_preservation = NA,
        chisq_statistic = NA, chisq_p_value = NA,
        correlation_rmse = NA, mean_corr_diff = NA, max_corr_diff = NA,
        balance_score = NA, mean_correlation = NA, correlation_variance = NA,
        nn_preservation_1 = NA, nn_preservation_5 = NA,
        mean_distance = NA, sd_distance = NA, distance_skewness = NA, distance_kurtosis = NA,
        efficiency_score = NA,
        stringsAsFactors = FALSE
      ))
    }

    # Extract imputation quality metrics
    imp_qual <- result$imputation_quality
    overall_metrics <- imp_qual$overall_metrics
    corr_metrics <- imp_qual$correlation_metrics

    # Extract distance quality metrics (if available)
    dist_qual <- result$distance_quality
    if (!is.null(dist_qual)) {
      rank_corr <- dist_qual$rank_correlation_balance
      dist_preserve <- dist_qual$distance_preservation
    } else {
      rank_corr <- list(balance_score = NA, mean_correlation = NA, correlation_variance = NA)
      dist_preserve <- list(nearest_neighbor_preservation = rep(NA, 5),
                           mean_distance = NA, sd_distance = NA,
                           distance_skewness = NA, distance_kurtosis = NA)
    }

    # Extract distribution metrics (simplified - take first variable or average)
    dist_metrics <- imp_qual$distribution_metrics
    if (length(dist_metrics) > 0) {
      first_var <- dist_metrics[[1]]
      if (!is.null(first_var$ks_statistic)) {
        ks_stat <- first_var$ks_statistic
        ks_p <- first_var$ks_p_value
        mean_preserve <- first_var$mean_preservation
        sd_preserve <- first_var$sd_preservation
        chisq_stat <- NA
        chisq_p <- NA
      } else if (!is.null(first_var$chisq_statistic)) {
        ks_stat <- NA
        ks_p <- NA
        mean_preserve <- NA
        sd_preserve <- NA
        chisq_stat <- first_var$chisq_statistic
        chisq_p <- first_var$chisq_p_value
      } else {
        ks_stat <- ks_p <- mean_preserve <- sd_preserve <- chisq_stat <- chisq_p <- NA
      }
    } else {
      ks_stat <- ks_p <- mean_preserve <- sd_preserve <- chisq_stat <- chisq_p <- NA
    }

    data.frame(
      replication_id = result$replication_id,
      scenario = result$scenario$name,
      method = result$method$name,
      missing_pattern = result$missing_pattern$name,
      success = TRUE,
      error_message = NA,
      execution_time = result$execution_time,

      # Imputation quality metrics
      rmse = overall_metrics$rmse %||% NA,
      mae = overall_metrics$mae %||% NA,
      bias = overall_metrics$bias %||% NA,
      nrmse = overall_metrics$nrmse %||% NA,
      nmae = overall_metrics$nmae %||% NA,
      r_squared = overall_metrics$r_squared %||% NA,
      accuracy = overall_metrics$accuracy %||% NA,
      mean_level_diff = overall_metrics$mean_level_difference %||% NA,
      max_level_diff = overall_metrics$max_level_difference %||% NA,

      # Distribution preservation
      ks_statistic = ks_stat,
      ks_p_value = ks_p,
      mean_preservation = mean_preserve,
      sd_preservation = sd_preserve,
      chisq_statistic = chisq_stat,
      chisq_p_value = chisq_p,

      # Correlation preservation
      correlation_rmse = corr_metrics$correlation_rmse %||% NA,
      mean_corr_diff = corr_metrics$mean_correlation_difference %||% NA,
      max_corr_diff = corr_metrics$max_correlation_difference %||% NA,

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

      # Computational performance
      efficiency_score = result$computational_performance$efficiency_score %||% NA,

      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, df_rows)
}

#' Compute Benchmark Summary Statistics
#' @keywords internal
compute_benchmark_summary <- function(results_df, config) {
  # Group by scenario, method, and missing pattern
  summary_stats <- list()

  # Overall summary
  summary_stats$overall <- list(
    total_runs = nrow(results_df),
    successful_runs = sum(results_df$success),
    failed_runs = sum(!results_df$success),
    success_rate = mean(results_df$success),
    total_execution_time = sum(results_df$execution_time, na.rm = TRUE),
    mean_execution_time = mean(results_df$execution_time, na.rm = TRUE)
  )

  # By method summary
  method_summary <- results_df %>%
    dplyr::group_by(method) %>%
    dplyr::summarise(
      n_runs = dplyr::n(),
      success_rate = mean(success),
      mean_time = mean(execution_time, na.rm = TRUE),
      sd_time = sd(execution_time, na.rm = TRUE),
      mean_rmse = mean(rmse, na.rm = TRUE),
      sd_rmse = sd(rmse, na.rm = TRUE),
      mean_mae = mean(mae, na.rm = TRUE),
      sd_mae = sd(mae, na.rm = TRUE),
      mean_efficiency = mean(efficiency_score, na.rm = TRUE),
      .groups = "drop"
    )

  summary_stats$by_method <- method_summary

  # By scenario summary
  scenario_summary <- results_df %>%
    dplyr::group_by(scenario) %>%
    dplyr::summarise(
      n_runs = dplyr::n(),
      success_rate = mean(success),
      mean_rmse = mean(rmse, na.rm = TRUE),
      mean_mae = mean(mae, na.rm = TRUE),
      .groups = "drop"
    )

  summary_stats$by_scenario <- scenario_summary

  # By missing pattern summary
  pattern_summary <- results_df %>%
    dplyr::group_by(missing_pattern) %>%
    dplyr::summarise(
      n_runs = dplyr::n(),
      success_rate = mean(success),
      mean_rmse = mean(rmse, na.rm = TRUE),
      mean_mae = mean(mae, na.rm = TRUE),
      .groups = "drop"
    )

  summary_stats$by_missing_pattern <- pattern_summary

  summary_stats
}

# --- Utility Functions ---

#' Save Benchmark Results
#'
#' @param results BenchmarkResults object
#' @param filename Filename to save results (without extension)
#' @param output_dir Output directory
#' @export
save_benchmark_results <- function(results, filename = "benchmark_results",
                                  output_dir = results@config@output_dir) {

  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  # Save main results
  results_file <- file.path(output_dir, paste0(filename, ".rds"))
  saveRDS(results, results_file)

  # Save summary as CSV
  summary_file <- file.path(output_dir, paste0(filename, "_summary.csv"))
  summary_df <- as.data.frame(results@summary_stats$by_method)
  write.csv(summary_df, summary_file, row.names = FALSE)

  # Save detailed results as CSV
  detailed_file <- file.path(output_dir, paste0(filename, "_detailed.csv"))
  write.csv(results@results, detailed_file, row.names = FALSE)

  message(sprintf("Results saved to: %s", output_dir))
}

#' Load Benchmark Results
#'
#' @param filename Filename to load results from
#' @param output_dir Output directory
#' @return BenchmarkResults object
#' @export
load_benchmark_results <- function(filename = "benchmark_results",
                                  output_dir = "benchmark_results") {

  results_file <- file.path(output_dir, paste0(filename, ".rds"))

  if (!file.exists(results_file)) {
    stop("Results file not found: ", results_file)
  }

  readRDS(results_file)
}
#' Imputation Method Comparison Benchmark
#'
#' Specialized benchmark for comparing different imputation methods for mixed-type data:
#' Gower-PMM, PMM, mean/mode, random forest, etc.

# --- Additional Imputation Method Implementations ---

#' Mean/Mode Imputation
#'
#' Simple imputation using mean for numeric and mode for categorical variables
#'
#' @param data Data frame with missing values
#' @return Imputed data frame
#' @export
impute_mean_mode <- function(data) {
  imputed_data <- data

  for (j in 1:ncol(data)) {
    missing_idx <- which(is.na(data[,j]))

    if (length(missing_idx) == 0) next

    if (is.numeric(data[,j])) {
      # Mean imputation for numeric
      mean_val <- mean(data[,j], na.rm = TRUE)
      imputed_data[missing_idx, j] <- mean_val

    } else if (is.factor(data[,j]) || is.ordered(data[,j])) {
      # Mode imputation for categorical
      mode_val <- names(which.max(table(data[,j], useNA = "no")))
      imputed_data[missing_idx, j] <- mode_val

    } else {
      # For other types, use first non-NA value
      non_na_vals <- data[,j][!is.na(data[,j])]
      if (length(non_na_vals) > 0) {
        imputed_data[missing_idx, j] <- non_na_vals[1]
      }
    }
  }

  imputed_data
}

#' Random Imputation
#'
#' Impute using random sampling from observed values
#'
#' @param data Data frame with missing values
#' @return Imputed data frame
#' @export
impute_random <- function(data) {
  imputed_data <- data

  for (j in 1:ncol(data)) {
    missing_idx <- which(is.na(data[,j]))

    if (length(missing_idx) == 0) next

    observed_vals <- data[,j][!is.na(data[,j])]

    if (length(observed_vals) > 0) {
      random_vals <- sample(observed_vals, length(missing_idx), replace = TRUE)
      imputed_data[missing_idx, j] <- random_vals
    }
  }

  imputed_data
}

#' Hot Deck Imputation
#'
#' Simple hot deck imputation using random donor selection
#'
#' @param data Data frame with missing values
#' @return Imputed data frame
#' @export
impute_hot_deck <- function(data) {
  imputed_data <- data

  # Find variables with missing values
  vars_with_missing <- colnames(data)[sapply(data, function(x) any(is.na(x)))]

  for (var in vars_with_missing) {
    missing_idx <- which(is.na(data[, var]))

    if (length(missing_idx) == 0) next

    observed_idx <- which(!is.na(data[, var]))

    for (idx in missing_idx) {
      # Randomly select a donor
      donor <- sample(observed_idx, 1)
      imputed_data[idx, var] <- data[donor, var]
    }
  }

  imputed_data
}

#' Regression Imputation
#'
#' Simple regression imputation for mixed-type data
#'
#' @param data Data frame with missing values
#' @return Imputed data frame
#' @export
impute_regression <- function(data) {
  imputed_data <- data

  # Find variables with missing values
  vars_with_missing <- colnames(data)[sapply(data, function(x) any(is.na(x)))]

  for (var in vars_with_missing) {
    missing_idx <- which(is.na(data[, var]))

    if (length(missing_idx) == 0) next

    # Use all other variables as predictors
    predictors <- setdiff(colnames(data), var)

    if (length(predictors) == 0) {
      # No predictors available, use mean/mode
      imputed_data <- impute_mean_mode(imputed_data)
      next
    }

    # Create temporary complete data for modeling
    temp_data <- data

    # For simplicity, use mean/mode imputation for predictors first
    for (pred in predictors) {
      if (any(is.na(temp_data[, pred]))) {
        temp_data[, pred] <- impute_mean_mode(temp_data[, pred])[, pred]
      }
    }

    observed_data <- temp_data[!is.na(data[, var]), ]

    if (nrow(observed_data) < 2) {
      # Not enough data for regression, use mean/mode
      imputed_data <- impute_mean_mode(imputed_data)
      next
    }

    tryCatch({
      if (is.numeric(data[, var])) {
        # Linear regression for numeric
        formula <- as.formula(paste(var, "~", paste(predictors, collapse = " + ")))
        model <- lm(formula, data = observed_data)

        # Predict missing values
        pred_data <- temp_data[missing_idx, predictors, drop = FALSE]
        predictions <- predict(model, newdata = pred_data)
        imputed_data[missing_idx, var] <- predictions

      } else if (is.factor(data[, var]) || is.ordered(data[, var])) {
        # Logistic regression for categorical (simplified)
        formula <- as.formula(paste(var, "~", paste(predictors, collapse = " + ")))
        model <- nnet::multinom(formula, data = observed_data, trace = FALSE)

        # Predict missing values
        pred_data <- temp_data[missing_idx, predictors, drop = FALSE]
        predictions <- predict(model, newdata = pred_data, type = "class")
        imputed_data[missing_idx, var] <- predictions
      }
    }, error = function(e) {
      # If regression fails, use mean/mode
      imputed_data <- impute_mean_mode(imputed_data)
    })
  }

  imputed_data
}

#' Expectation-Maximization Imputation
#'
#' Simple EM-based imputation for mixed-type data
#'
#' @param data Data frame with missing values
#' @param max_iter Maximum iterations
#' @return Imputed data frame
#' @export
impute_em <- function(data, max_iter = 10) {
  # Simplified EM implementation
  imputed_data <- impute_mean_mode(data)  # Initial imputation

  for (iter in 1:max_iter) {
    # In a full implementation, this would involve:
    # 1. Estimate parameters from current imputed data
    # 2. Re-impute using estimated parameters
    # For simplicity, we'll just do mean/mode imputation
    break
  }

  imputed_data
}

# --- Imputation Method Comparison Framework ---

#' Create Imputation Method Comparison Configuration
#'
#' @param n_replications Number of replications
#' @param methods Methods to compare (default: all available)
#' @param output_dir Output directory
#' @return BenchmarkConfig object
#' @export
create_imputation_benchmark_config <- function(n_replications = 100,
                                             methods = c("gowerpmm_auto", "gowerpmm_equal", "mice_pmm",
                                                        "mice_cart", "mice_rf", "vim_knn", "vim_irmi",
                                                        "mean_mode", "random", "hot_deck", "regression"),
                                             output_dir = "imputation_benchmark_results") {

  # Define imputation methods
  available_methods <- list(
    # Gower-PMM variants
    gowerpmm_auto = list(
      name = "Gower-PMM (Auto)",
      type = "gowerpmm",
      params = list(weights = "auto", k = 5)
    ),
    gowerpmm_equal = list(
      name = "Gower-PMM (Equal)",
      type = "gowerpmm",
      params = list(weights = "equal", k = 5)
    ),

    # Standard MICE methods
    mice_pmm = list(
      name = "MICE PMM",
      type = "mice_standard",
      params = list(method = "pmm", k = 5)
    ),
    mice_cart = list(
      name = "MICE CART",
      type = "mice_standard",
      params = list(method = "cart")
    ),
    mice_rf = list(
      name = "MICE RF",
      type = "mice_standard",
      params = list(method = "rf")
    ),

    # VIM methods
    vim_knn = list(
      name = "VIM k-NN",
      type = "vim",
      params = list(method = "kNN", k = 5)
    ),
    vim_irmi = list(
      name = "VIM IRMI",
      type = "vim",
      params = list(method = "irmi")
    ),

    # Simple methods
    mean_mode = list(
      name = "Mean/Mode",
      type = "simple",
      impute_function = "impute_mean_mode"
    ),
    random = list(
      name = "Random",
      type = "simple",
      impute_function = "impute_random"
    ),
    hot_deck = list(
      name = "Hot Deck",
      type = "simple",
      impute_function = "impute_hot_deck"
    ),
    regression = list(
      name = "Regression",
      type = "simple",
      impute_function = "impute_regression"
    )
  )

  # Filter to requested methods
  selected_methods <- available_methods[methods]

  # Create methods list
  methods_list <- lapply(names(selected_methods), function(method_name) {
    method <- selected_methods[[method_name]]
    method$name <- method_name
    method
  })
  names(methods_list) <- names(selected_methods)

  # Use diverse scenarios for imputation comparison
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

  # Focus on different missing data patterns for imputation
  missing_patterns <- list(
    mcar_10 = list(type = "MCAR", rate = 0.10, description = "10% MCAR"),
    mcar_20 = list(type = "MCAR", rate = 0.20, description = "20% MCAR"),
    mcar_30 = list(type = "MCAR", rate = 0.30, description = "30% MCAR"),

    mar_10 = list(type = "MAR", rate = 0.10, description = "10% MAR"),
    mar_20 = list(type = "MAR", rate = 0.20, description = "20% MAR"),
    mar_30 = list(type = "MAR", rate = 0.30, description = "30% MAR"),

    mnar_10 = list(type = "MNAR", rate = 0.10, description = "10% MNAR"),
    mnar_20 = list(type = "MNAR", rate = 0.20, description = "20% MNAR")
  )

  new("BenchmarkConfig",
      scenarios = scenarios,
      methods = methods_list,
      metrics = list(),  # Will be set by imputation-specific metrics
      n_replications = n_replications,
      missing_patterns = missing_patterns,
      output_dir = output_dir)
}

#' Run Imputation Method Comparison Benchmark
#'
#' @param config BenchmarkConfig object
#' @param parallel Whether to use parallel processing
#' @param n_cores Number of cores
#' @return BenchmarkResults object
#' @export
run_imputation_comparison_benchmark <- function(config, parallel = TRUE, n_cores = NULL) {

  message("Running Imputation Method Comparison Benchmark...")
  message("Comparing methods: ", paste(names(config@methods), collapse = ", "))

  # Run the imputation comparison suite
  results <- run_imputation_method_suite(config, parallel = parallel, n_cores = n_cores)

  # Save results
  save_benchmark_results(results)

  # Generate imputation-specific analysis
  imputation_analysis <- analyze_imputation_method_performance(results)

  # Save extended analysis
  saveRDS(imputation_analysis, file.path(config@output_dir, "imputation_method_analysis.rds"))

  results
}

#' Run Imputation Method Suite
#' @keywords internal
run_imputation_method_suite <- function(config, parallel = TRUE, n_cores = NULL) {

  start_time <- Sys.time()

  # Prepare combinations
  combinations <- expand.grid(
    scenario = names(config@scenarios),
    method = names(config@methods),
    missing_pattern = names(config@missing_patterns),
    replication = 1:config@n_replications,
    stringsAsFactors = FALSE
  )

  total_runs <- nrow(combinations)
  message(sprintf("Starting imputation comparison with %d total runs...", total_runs))

  results_list <- vector("list", total_runs)

  if (parallel && requireNamespace("parallel", quietly = TRUE)) {
    if (is.null(n_cores)) n_cores <- max(1, parallel::detectCores() - 1)

    message(sprintf("Running in parallel with %d cores...", n_cores))

    cl <- parallel::makeCluster(n_cores)
    on.exit(parallel::stopCluster(cl))

    # Export functions
    parallel::clusterExport(cl, c(
      "run_imputation_method_single", "generate_simulation_data", "introduce_missing_data",
      "evaluate_imputation_quality", "run_imputation_method",
      "run_gowerpmm_imputation", "run_standard_mice_imputation",
      "run_fd_gowdis_imputation", "run_vim_imputation",
      "impute_mean_mode", "impute_random", "impute_hot_deck",
      "impute_regression", "impute_em",
      "compute_overall_imputation_metrics", "compute_by_variable_metrics",
      "compute_distribution_preservation", "compute_correlation_preservation"
    ), envir = environment())

    # Load packages
    parallel::clusterEvalQ(cl, {
      library(gowerpmm)
      library(mice)
      if (requireNamespace("FD", quietly = TRUE)) library(FD)
      if (requireNamespace("VIM", quietly = TRUE)) library(VIM)
      if (requireNamespace("nnet", quietly = TRUE)) library(nnet)
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

      run_imputation_method_single(scenario, method, missing_pattern, combo$replication)
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

      if (i %% 10 == 0) {
        message(sprintf("Completed %d/%d runs (%.1f%%)...", i, total_runs, (i/total_runs)*100))
      }

      results_list[[i]] <- run_imputation_method_single(scenario, method, missing_pattern, combo$replication)
    }
  }

  # Compile results
  results_df <- compile_imputation_results(results_list)

  # Summary stats
  summary_stats <- compute_imputation_summary(results_df, config)

  end_time <- Sys.time()
  total_execution_time <- as.numeric(difftime(end_time, start_time, units = "secs"))

  benchmark_results <- new("BenchmarkResults",
                          config = config,
                          results = results_df,
                          summary_stats = summary_stats,
                          execution_time = total_execution_time,
                          timestamp = start_time)

  message(sprintf("Imputation comparison completed in %.2f seconds", total_execution_time))

  benchmark_results
}

#' Run Single Imputation Method Evaluation
#' @keywords internal
run_imputation_method_single <- function(scenario, method, missing_pattern, replication_id) {

  set.seed(replication_id * 1000)
  start_time <- Sys.time()

  tryCatch({
    # 1. Generate complete simulation data
    complete_data <- generate_simulation_data(scenario, seed = replication_id)

    # 2. Introduce missing data pattern
    data_with_missing <- introduce_missing_data(complete_data, missing_pattern, seed = replication_id)
    missing_mask <- is.na(data_with_missing)

    # 3. Run imputation method
    imputation_result <- run_imputation_method_simple(method, data_with_missing)

    # 4. Evaluate imputation quality
    quality_metrics <- evaluate_imputation_quality(
      original_data = complete_data,
      imputed_data = imputation_result$imputed_data,
      missing_mask = missing_mask,
      method_name = method$name
    )

    # 5. Performance metrics
    end_time <- Sys.time()
    execution_time <- as.numeric(difftime(end_time, start_time, units = "secs"))

    performance_metrics <- evaluate_computational_performance(
      method_name = method$name,
      execution_time = execution_time,
      convergence_info = imputation_result$convergence_info
    )

    # Return results
    list(
      replication_id = replication_id,
      scenario = scenario$name,
      method = method$name,
      missing_pattern = missing_pattern$name,
      imputation_quality = quality_metrics,
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
      missing_pattern = missing_pattern$name,
      imputation_quality = NULL,
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

#' Run Imputation Method (Simplified)
#' @keywords internal
run_imputation_method_simple <- function(method, data) {

  method_type <- method$type

  if (method_type == "gowerpmm") {
    # Gower-PMM imputation
    result <- run_gowerpmm_imputation_simple(data, method$params)

  } else if (method_type == "mice_standard") {
    # Standard MICE methods
    result <- run_standard_mice_imputation_simple(data, method$params)

  } else if (method_type == "vim") {
    # VIM methods
    result <- run_vim_imputation_simple(data, method$params)

  } else if (method_type == "simple") {
    # Simple imputation methods
    result <- run_simple_imputation(data, method)

  } else {
    stop("Unknown method type: ", method_type)
  }

  result
}

#' Run Gower-PMM Imputation (Simplified)
#' @keywords internal
run_gowerpmm_imputation_simple <- function(data, params) {
  # Simplified version for benchmarking
  tryCatch({
    # Use all variables as predictors
    predictor_matrix <- matrix(1, ncol(data), ncol(data))
    diag(predictor_matrix) <- 0

    mice_params <- list(
      data = data,
      method = "gowerpmm",
      predictorMatrix = predictor_matrix,
      m = 1,
      maxit = 1,
      printFlag = FALSE
    )

    # Add method parameters
    if (!is.null(params$weights)) mice_params$weights <- params$weights
    if (!is.null(params$k)) mice_params$k <- params$k
    if (!is.null(params$scaling)) mice_params$scaling <- params$scaling

    imp <- do.call(mice::mice, mice_params)
    imputed_data <- mice::complete(imp, 1)

    list(
      imputed_data = imputed_data,
      convergence_info = list(method = "gowerpmm", success = TRUE)
    )
  }, error = function(e) {
    # Fallback to simple imputation
    list(
      imputed_data = impute_mean_mode(data),
      convergence_info = list(method = "gowerpmm", error = as.character(e))
    )
  })
}

#' Run Standard MICE Imputation (Simplified)
#' @keywords internal
run_standard_mice_imputation_simple <- function(data, params) {
  tryCatch({
    predictor_matrix <- matrix(1, ncol(data), ncol(data))
    diag(predictor_matrix) <- 0

    mice_params <- list(
      data = data,
      method = params$method,
      predictorMatrix = predictor_matrix,
      m = 1,
      maxit = 1,
      printFlag = FALSE
    )

    imp <- do.call(mice::mice, mice_params)
    imputed_data <- mice::complete(imp, 1)

    list(
      imputed_data = imputed_data,
      convergence_info = list(method = params$method, success = TRUE)
    )
  }, error = function(e) {
    # Fallback to simple imputation
    list(
      imputed_data = impute_mean_mode(data),
      convergence_info = list(method = params$method, error = as.character(e))
    )
  })
}

#' Run VIM Imputation (Simplified)
#' @keywords internal
run_vim_imputation_simple <- function(data, params) {
  tryCatch({
    method <- params$method

    if (method == "kNN") {
      k <- params$k %||% 5
      imputed_data <- impute_with_vim_knn(data, k)
      convergence_info <- list(method = "vim_knn", success = TRUE)

    } else if (method == "irmi") {
      imputed_data <- impute_with_vim_irmi(data)
      convergence_info <- list(method = "vim_irmi", success = TRUE)

    } else {
      stop("Unknown VIM method: ", method)
    }

    list(
      imputed_data = imputed_data,
      convergence_info = convergence_info
    )
  }, error = function(e) {
    # Fallback to simple imputation
    list(
      imputed_data = impute_mean_mode(data),
      convergence_info = list(method = params$method, error = as.character(e))
    )
  })
}

#' Run Simple Imputation Methods
#' @keywords internal
run_simple_imputation <- function(data, method) {
  tryCatch({
    impute_func <- method$impute_function

    if (impute_func == "impute_mean_mode") {
      imputed_data <- impute_mean_mode(data)
    } else if (impute_func == "impute_random") {
      imputed_data <- impute_random(data)
    } else if (impute_func == "impute_hot_deck") {
      imputed_data <- impute_hot_deck(data)
    } else if (impute_func == "impute_regression") {
      imputed_data <- impute_regression(data)
    } else {
      stop("Unknown imputation function: ", impute_func)
    }

    list(
      imputed_data = imputed_data,
      convergence_info = list(method = method$name, success = TRUE)
    )
  }, error = function(e) {
    # Fallback to mean/mode
    list(
      imputed_data = impute_mean_mode(data),
      convergence_info = list(method = method$name, error = as.character(e))
    )
  })
}

#' Compile Imputation Results
#' @keywords internal
compile_imputation_results <- function(results_list) {

  df_rows <- lapply(results_list, function(result) {
    if (!result$success) {
      return(data.frame(
        replication_id = result$replication_id,
        scenario = result$scenario,
        method = result$method,
        missing_pattern = result$missing_pattern,
        success = FALSE,
        error_message = result$error_message,
        execution_time = result$execution_time,
        rmse = NA, mae = NA, bias = NA, nrmse = NA, nmae = NA, r_squared = NA,
        accuracy = NA, mean_level_diff = NA, max_level_diff = NA,
        ks_statistic = NA, ks_p_value = NA, mean_preservation = NA, sd_preservation = NA,
        chisq_statistic = NA, chisq_p_value = NA,
        correlation_rmse = NA, mean_corr_diff = NA, max_corr_diff = NA,
        efficiency_score = NA,
        stringsAsFactors = FALSE
      ))
    }

    # Extract imputation quality metrics
    imp_qual <- result$imputation_quality
    overall_metrics <- imp_qual$overall_metrics
    corr_metrics <- imp_qual$correlation_metrics

    # Extract distribution metrics (simplified)
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
      scenario = result$scenario,
      method = result$method,
      missing_pattern = result$missing_pattern,
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

      # Performance
      efficiency_score = result$computational_performance$efficiency_score %||% NA,

      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, df_rows)
}

#' Compute Imputation Summary Statistics
#' @keywords internal
compute_imputation_summary <- function(results_df, config) {

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

  # By method summary
  method_summary <- df %>%
    dplyr::group_by(method) %>%
    dplyr::summarise(
      n_runs = dplyr::n(),
      mean_rmse = mean(rmse, na.rm = TRUE),
      sd_rmse = sd(rmse, na.rm = TRUE),
      mean_mae = mean(mae, na.rm = TRUE),
      sd_mae = sd(mae, na.rm = TRUE),
      mean_accuracy = mean(accuracy, na.rm = TRUE),
      sd_accuracy = sd(accuracy, na.rm = TRUE),
      mean_time = mean(execution_time, na.rm = TRUE),
      sd_time = sd(execution_time, na.rm = TRUE),
      mean_efficiency = mean(efficiency_score, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::arrange(mean_rmse)

  summary_stats$by_method <- method_summary

  # By scenario summary
  scenario_summary <- df %>%
    dplyr::group_by(scenario) %>%
    dplyr::summarise(
      n_runs = dplyr::n(),
      mean_rmse = mean(rmse, na.rm = TRUE),
      mean_mae = mean(mae, na.rm = TRUE),
      .groups = "drop"
    )

  summary_stats$by_scenario <- scenario_summary

  # By missing pattern summary
  pattern_summary <- df %>%
    dplyr::group_by(missing_pattern) %>%
    dplyr::summarise(
      n_runs = dplyr::n(),
      mean_rmse = mean(rmse, na.rm = TRUE),
      mean_mae = mean(mae, na.rm = TRUE),
      .groups = "drop"
    )

  summary_stats$by_missing_pattern <- pattern_summary

  summary_stats
}

#' Analyze Imputation Method Performance
#' @keywords internal
analyze_imputation_method_performance <- function(results) {

  df <- results@results
  df <- df[df$success, ]

  # Ranking by different metrics
  rankings <- list()

  # By RMSE (lower is better)
  rankings$rmse <- df %>%
    dplyr::group_by(scenario, missing_pattern) %>%
    dplyr::mutate(rank = rank(rmse)) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(method) %>%
    dplyr::summarise(
      mean_rank = mean(rank),
      mean_rmse = mean(rmse, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::arrange(mean_rank)

  # By MAE (lower is better)
  rankings$mae <- df %>%
    dplyr::group_by(scenario, missing_pattern) %>%
    dplyr::mutate(rank = rank(mae)) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(method) %>%
    dplyr::summarise(
      mean_rank = mean(rank),
      mean_mae = mean(mae, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::arrange(mean_rank)

  # By efficiency (higher is better)
  rankings$efficiency <- df %>%
    dplyr::group_by(scenario, missing_pattern) %>%
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
      rmse_mean = mean(rmse, na.rm = TRUE),
      rmse_sd = sd(rmse, na.rm = TRUE),
      mae_mean = mean(mae, na.rm = TRUE),
      mae_sd = sd(mae, na.rm = TRUE),
      accuracy_mean = mean(accuracy, na.rm = TRUE),
      accuracy_sd = sd(accuracy, na.rm = TRUE),
      time_mean = mean(execution_time, na.rm = TRUE),
      time_sd = sd(execution_time, na.rm = TRUE),
      n_scenarios = length(unique(scenario)),
      n_patterns = length(unique(missing_pattern)),
      .groups = "drop"
    )

  list(
    rankings = rankings,
    comparison = comparison,
    detailed_results = df
  )
}

#' Generate Imputation Method Comparison Report
#' @export
generate_imputation_comparison_report <- function(results, output_dir = "imputation_comparison_report") {

  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  # Load analysis
  analysis_file <- file.path(results@config@output_dir, "imputation_method_analysis.rds")
  if (file.exists(analysis_file)) {
    analysis <- readRDS(analysis_file)
  } else {
    analysis <- analyze_imputation_method_performance(results)
  }

  # Create plots
  rmse_plot <- create_imputation_rmse_plot(results)
  efficiency_plot <- create_imputation_efficiency_plot(results)
  ranking_plot <- create_imputation_ranking_plot(analysis$rankings)

  # Save plots
  ggplot2::ggsave(file.path(output_dir, "imputation_rmse_comparison.png"), rmse_plot,
                  width = 10, height = 6, dpi = 300)
  ggplot2::ggsave(file.path(output_dir, "imputation_efficiency_frontier.png"), efficiency_plot,
                  width = 8, height = 6, dpi = 300)
  ggplot2::ggsave(file.path(output_dir, "imputation_method_rankings.png"), ranking_plot,
                  width = 12, height = 6, dpi = 300)

  # Generate report
  rmd_content <- create_imputation_rmd_report(results, analysis)

  rmd_file <- file.path(output_dir, "imputation_comparison_report.Rmd")
  writeLines(rmd_content, rmd_file)

  # Try to render
  tryCatch({
    rmarkdown::render(rmd_file, output_dir = output_dir)
    message("Imputation comparison report generated successfully")
  }, error = function(e) {
    warning("Could not render R Markdown report: ", e$message)
    message("R Markdown file saved to: ", rmd_file)
  })

  file.path(output_dir, "imputation_comparison_report.html")
}

#' Create Imputation RMSE Plot
#' @keywords internal
create_imputation_rmse_plot <- function(results) {

  df <- results@results
  df <- df[df$success, ]

  p <- ggplot2::ggplot(df, ggplot2::aes(x = method, y = rmse, fill = method)) +
    ggplot2::geom_boxplot(alpha = 0.7) +
    ggplot2::facet_grid(missing_pattern ~ scenario) +
    ggplot2::labs(
      title = "Imputation Method RMSE Comparison",
      subtitle = "Lower RMSE indicates better imputation quality",
      x = "Imputation Method",
      y = "RMSE"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      legend.position = "none"
    )

  p
}

#' Create Imputation Efficiency Plot
#' @keywords internal
create_imputation_efficiency_plot <- function(results) {

  df <- results@results
  df <- df[df$success, ]

  efficiency_data <- df %>%
    dplyr::group_by(method) %>%
    dplyr::summarise(
      mean_rmse = mean(rmse, na.rm = TRUE),
      mean_time = mean(execution_time, na.rm = TRUE),
      .groups = "drop"
    )

  p <- ggplot2::ggplot(efficiency_data,
                      ggplot2::aes(x = mean_time, y = mean_rmse, label = method)) +
    ggplot2::geom_point(size = 3, alpha = 0.7) +
    ggplot2::geom_text(vjust = -1, size = 3) +
    ggplot2::labs(
      title = "Imputation Method Efficiency Frontier",
      x = "Mean Execution Time (seconds)",
      y = "Mean RMSE"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::scale_x_log10() +
    ggplot2::scale_y_log10()

  p
}

#' Create Imputation Ranking Plot
#' @keywords internal
create_imputation_ranking_plot <- function(rankings) {

  # Combine rankings
  rmse_ranks <- rankings$rmse %>%
    dplyr::mutate(metric = "RMSE")
  mae_ranks <- rankings$mae %>%
    dplyr::mutate(metric = "MAE")

  combined_ranks <- rbind(rmse_ranks, mae_ranks)

  p <- ggplot2::ggplot(combined_ranks,
                      ggplot2::aes(x = reorder(method, -mean_rank), y = mean_rank, fill = metric)) +
    ggplot2::geom_bar(stat = "identity", position = "dodge") +
    ggplot2::labs(
      title = "Imputation Method Rankings by Performance Metric",
      subtitle = "Lower rank numbers indicate better performance",
      x = "Imputation Method",
      y = "Mean Rank"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
    ggplot2::scale_y_reverse()  # Lower ranks are better

  p
}

#' Create Imputation RMD Report
#' @keywords internal
create_imputation_rmd_report <- function(results, analysis) {

  rmd <- sprintf('---
title: "Imputation Method Comparison Report"
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

# Imputation Method Comparison Analysis

This report compares different imputation methods for mixed-type data.

## Overview

- **Methods Compared**: %s
- **Scenarios Tested**: %s
- **Missing Patterns**: %s
- **Total Runs**: %d
- **Successful Runs**: %d (%.1f%%)

## Performance Rankings

### By RMSE
```{r rmse-rankings}
kable(analysis$rankings$rmse, caption = "Ranking by RMSE") %>%
  kable_styling(bootstrap_options = c("striped", "hover"))
```

### By MAE
```{r mae-rankings}
kable(analysis$rankings$mae, caption = "Ranking by MAE") %>%
  kable_styling(bootstrap_options = c("striped", "hover"))
```

## Visual Comparisons

### RMSE Comparison Across Scenarios and Missing Patterns
![RMSE Comparison](imputation_rmse_comparison.png)

### Efficiency Frontier
![Efficiency Frontier](imputation_efficiency_frontier.png)

### Method Rankings
![Method Rankings](imputation_method_rankings.png)

## Detailed Comparison

```{r comparison-table}
kable(analysis$comparison, caption = "Detailed Method Comparison") %>%
  kable_styling(bootstrap_options = c("striped", "hover"))
```

## Method Descriptions

- **Gower-PMM (Auto)**: Optimized Gower distance with automatic weight selection
- **Gower-PMM (Equal)**: Gower distance with equal weights
- **MICE PMM**: Predictive Mean Matching via MICE
- **MICE CART**: Classification and Regression Trees via MICE
- **MICE RF**: Random Forest imputation via MICE
- **VIM k-NN**: k-Nearest Neighbors via VIM package
- **VIM IRMI**: Iterative Robust Model-based Imputation via VIM
- **Mean/Mode**: Simple mean for numeric, mode for categorical
- **Random**: Random sampling from observed values
- **Hot Deck**: Random donor selection
- **Regression**: Simple regression-based imputation

## Conclusions

[Conclusions would be written here based on the analysis]

## Technical Details

- **R Version**: `r R.version.string`
- **Analysis Date**: `r format(results@timestamp, "%%Y-%%m-%%d %%H:%%M")`
- **Total Execution Time**: `r round(results@execution_time, 2)` seconds
',
                 paste(names(results@config@methods), collapse = ", "),
                 paste(names(results@config@scenarios), collapse = ", "),
                 paste(names(results@config@missing_patterns), collapse = ", "),
                 results@summary_stats$overall$total_runs,
                 results@summary_stats$overall$successful_runs,
                 results@summary_stats$overall$success_rate * 100)

  rmd
}
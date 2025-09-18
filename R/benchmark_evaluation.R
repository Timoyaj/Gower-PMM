#' Benchmark Evaluation Metrics and Comparison Methods
#'
#' This module provides comprehensive evaluation metrics for assessing the quality
#' of distance measures and imputation methods in mixed-type data scenarios.

# --- Distance Measure Evaluation ---

#' Evaluate Distance Measure Quality
#'
#' Computes various metrics to assess the quality of distance measures
#' for mixed-type data, focusing on rank correlation balance and preservation.
#'
#' @param data Complete dataset
#' @param distance_matrix Computed distance matrix
#' @param method_name Name of the distance method
#' @return List of distance quality metrics
#' @export
evaluate_distance_quality <- function(data, distance_matrix, method_name) {
  n <- nrow(data)
  p <- ncol(data)

  # Convert distance matrix to dist object if needed
  if (!inherits(distance_matrix, "dist")) {
    if (is.matrix(distance_matrix)) {
      distance_matrix <- as.dist(distance_matrix)
    }
  }

  # 1. Rank Correlation Balance
  # Assess how well the distance measure balances different variable types
  rank_corr_balance <- compute_rank_correlation_balance(data, distance_matrix)

  # 2. Distance Preservation
  # How well does the distance preserve the original data structure?
  dist_preservation <- compute_distance_preservation(data, distance_matrix)

  # 3. Variable Type Balance
  # Assess contribution of different variable types
  var_type_balance <- compute_variable_type_balance(data, distance_matrix)

  # 4. Computational Efficiency
  comp_efficiency <- list(
    method = method_name,
    n_observations = n,
    n_variables = p,
    distance_matrix_size = length(distance_matrix)
  )

  list(
    rank_correlation_balance = rank_corr_balance,
    distance_preservation = dist_preservation,
    variable_type_balance = var_type_balance,
    computational_efficiency = comp_efficiency
  )
}

#' Compute Rank Correlation Balance
#' @keywords internal
compute_rank_correlation_balance <- function(data, distance_matrix) {
  # For each variable type, compute how well distances correlate with variable differences
  var_types <- sapply(data, function(col) {
    if (is.ordered(col)) "ordinal"
    else if (is.factor(col)) "categorical"
    else if (is.numeric(col)) "numeric"
    else "other"
  })

  correlations <- list()

  # Numeric variables: correlation between distances and absolute differences
  if (any(var_types == "numeric")) {
    num_vars <- which(var_types == "numeric")
    for (i in num_vars) {
      if (length(unique(data[,i])) > 2) {  # Need variation
        var_distances <- as.matrix(stats::dist(data[,i]))
        overall_distances <- as.matrix(distance_matrix)

        # Compute correlation between variable-specific and overall distances
        corr <- tryCatch({
          cor(as.vector(var_distances), as.vector(overall_distances),
              method = "spearman", use = "pairwise.complete.obs")
        }, error = function(e) NA)

        correlations[[paste0("num_", i)]] <- corr
      }
    }
  }

  # Categorical variables: correlation between distances and matching
  if (any(var_types == "categorical")) {
    cat_vars <- which(var_types == "categorical")
    for (i in cat_vars) {
      cat_matrix <- outer(data[,i], data[,i], "!=") * 1.0
      overall_distances <- as.matrix(distance_matrix)

      corr <- tryCatch({
        cor(as.vector(cat_matrix), as.vector(overall_distances),
            method = "spearman", use = "pairwise.complete.obs")
      }, error = function(e) NA)

      correlations[[paste0("cat_", i)]] <- corr
    }
  }

  # Overall balance score: variance in correlations (lower is better balance)
  valid_corrs <- unlist(correlations)
  valid_corrs <- valid_corrs[!is.na(valid_corrs)]

  if (length(valid_corrs) > 1) {
    balance_score <- 1 / (1 + var(valid_corrs))  # Higher score = better balance
  } else {
    balance_score <- 0.5  # Neutral score
  }

  list(
    individual_correlations = correlations,
    balance_score = balance_score,
    mean_correlation = mean(valid_corrs, na.rm = TRUE),
    correlation_variance = var(valid_corrs, na.rm = TRUE)
  )
}

#' Compute Distance Preservation
#' @keywords internal
compute_distance_preservation <- function(data, distance_matrix) {
  n <- nrow(data)

  # 1. Nearest neighbor preservation
  # How often are the k nearest neighbors in original space also close in distance space?
  k <- min(5, floor(n/10))

  # Compute Euclidean distance for numeric variables only
  num_cols <- sapply(data, is.numeric)
  if (sum(num_cols) > 1) {
    euclidean_dist <- stats::dist(data[, num_cols, drop = FALSE])
    euclidean_matrix <- as.matrix(euclidean_dist)
  } else {
    euclidean_matrix <- matrix(0, n, n)
  }

  distance_matrix <- as.matrix(distance_matrix)

  preservation_scores <- numeric(k)
  for (i in 1:k) {
    # For each point, find k nearest neighbors in both spaces
    euclidean_neighbors <- apply(euclidean_matrix, 1, function(x) order(x)[2:(i+1)])
    distance_neighbors <- apply(distance_matrix, 1, function(x) order(x)[2:(i+1)])

    # Compute overlap
    overlaps <- sapply(1:n, function(j) {
      length(intersect(euclidean_neighbors[,j], distance_neighbors[,j]))
    })

    preservation_scores[i] <- mean(overlaps) / i
  }

  # 2. Distance distribution characteristics
  distance_values <- as.vector(distance_matrix[upper.tri(distance_matrix)])

  list(
    nearest_neighbor_preservation = preservation_scores,
    mean_distance = mean(distance_values),
    sd_distance = sd(distance_values),
    distance_range = range(distance_values),
    distance_skewness = moments::skewness(distance_values),
    distance_kurtosis = moments::kurtosis(distance_values)
  )
}

#' Compute Variable Type Balance
#' @keywords internal
compute_variable_type_balance <- function(data, distance_matrix) {
  var_types <- sapply(data, function(col) {
    if (is.ordered(col)) "ordinal"
    else if (is.factor(col)) "categorical"
    else if (is.numeric(col)) "numeric"
    else "other"
  })

  type_counts <- table(var_types)
  type_proportions <- type_counts / sum(type_counts)

  # Expected contribution if perfectly balanced
  expected_contribution <- 1 / length(type_counts)

  # Actual balance (inverse of variance in proportions)
  balance_score <- 1 / (1 + var(type_proportions))

  list(
    variable_type_counts = type_counts,
    variable_type_proportions = type_proportions,
    balance_score = balance_score,
    most_common_type = names(which.max(type_counts)),
    type_diversity = length(type_counts)
  )
}

# --- Imputation Quality Evaluation ---

#' Evaluate Imputation Quality
#'
#' Comprehensive assessment of imputation method performance
#'
#' @param original_data Complete original dataset
#' @param imputed_data Dataset with imputed values
#' @param missing_mask Logical matrix indicating missing values
#' @param method_name Name of imputation method
#' @return List of imputation quality metrics
#' @export
evaluate_imputation_quality <- function(original_data, imputed_data, missing_mask, method_name) {
  if (!all(dim(original_data) == dim(imputed_data))) {
    stop("Original and imputed data must have the same dimensions")
  }
  if (!all(dim(original_data) == dim(missing_mask))) {
    stop("Data and missing mask must have the same dimensions")
  }

  # Extract observed and imputed values
  observed_values <- original_data[missing_mask]
  imputed_values <- imputed_data[missing_mask]
  variable_indices <- rep(1:ncol(original_data), each = nrow(original_data))[missing_mask]

  results <- list(
    method = method_name,
    overall_metrics = compute_overall_imputation_metrics(observed_values, imputed_values),
    by_variable_metrics = compute_by_variable_metrics(original_data, imputed_data, missing_mask),
    distribution_metrics = compute_distribution_preservation(original_data, imputed_data, missing_mask),
    correlation_metrics = compute_correlation_preservation(original_data, imputed_data)
  )

  results
}

#' Compute Overall Imputation Metrics
#' @keywords internal
compute_overall_imputation_metrics <- function(observed, imputed) {
  # Handle different data types
  if (is.numeric(observed) && is.numeric(imputed)) {
    # Numeric metrics
    rmse <- sqrt(mean((observed - imputed)^2, na.rm = TRUE))
    mae <- mean(abs(observed - imputed), na.rm = TRUE)
    bias <- mean(imputed - observed, na.rm = TRUE)

    # Normalized metrics
    observed_sd <- sd(observed, na.rm = TRUE)
    nrmse <- rmse / observed_sd
    nmae <- mae / observed_sd

    # R-squared (coefficient of determination)
    ss_res <- sum((observed - imputed)^2, na.rm = TRUE)
    ss_tot <- sum((observed - mean(observed, na.rm = TRUE))^2, na.rm = TRUE)
    r_squared <- 1 - (ss_res / ss_tot)

    list(
      rmse = rmse, mae = mae, bias = bias,
      nrmse = nrmse, nmae = nmae, r_squared = r_squared
    )
  } else {
    # Categorical metrics
    accuracy <- mean(observed == imputed, na.rm = TRUE)
    # For ordered factors, compute mean absolute difference in levels
    if (is.ordered(observed) && is.ordered(imputed)) {
      level_diff <- abs(as.numeric(observed) - as.numeric(imputed))
      mean_level_diff <- mean(level_diff, na.rm = TRUE)
      max_level_diff <- max(level_diff, na.rm = TRUE)
    } else {
      mean_level_diff <- NA
      max_level_diff <- NA
    }

    list(
      accuracy = accuracy,
      mean_level_difference = mean_level_diff,
      max_level_difference = max_level_diff
    )
  }
}

#' Compute By-Variable Metrics
#' @keywords internal
compute_by_variable_metrics <- function(original_data, imputed_data, missing_mask) {
  p <- ncol(original_data)
  var_metrics <- vector("list", p)
  names(var_metrics) <- colnames(original_data)

  for (j in 1:p) {
    var_missing <- missing_mask[,j]
    if (any(var_missing)) {
      observed <- original_data[var_missing, j]
      imputed <- imputed_data[var_missing, j]

      var_metrics[[j]] <- compute_overall_imputation_metrics(observed, imputed)
      var_metrics[[j]]$n_missing <- sum(var_missing)
      var_metrics[[j]]$missing_rate <- mean(var_missing)
    } else {
      var_metrics[[j]] <- list(n_missing = 0, missing_rate = 0)
    }
  }

  var_metrics
}

#' Compute Distribution Preservation
#' @keywords internal
compute_distribution_preservation <- function(original_data, imputed_data, missing_mask) {
  p <- ncol(original_data)
  dist_metrics <- vector("list", p)
  names(dist_metrics) <- colnames(original_data)

  for (j in 1:p) {
    var_missing <- missing_mask[,j]
    observed_complete <- original_data[!var_missing, j]
    imputed_values <- imputed_data[var_missing, j]

    if (is.numeric(original_data[,j])) {
      # Kolmogorov-Smirnov test for numeric variables
      if (length(imputed_values) > 1 && length(observed_complete) > 1) {
        ks_test <- tryCatch({
          stats::ks.test(observed_complete, imputed_values)
        }, error = function(e) list(p.value = NA, statistic = NA))

        # Additional distribution metrics
        mean_preservation <- abs(mean(observed_complete) - mean(imputed_values)) / sd(observed_complete)
        sd_preservation <- abs(sd(observed_complete) - sd(imputed_values)) / sd(observed_complete)

        dist_metrics[[j]] <- list(
          ks_statistic = ks_test$statistic,
          ks_p_value = ks_test$p.value,
          mean_preservation = mean_preservation,
          sd_preservation = sd_preservation
        )
      }
    } else if (is.factor(original_data[,j]) || is.ordered(original_data[,j])) {
      # Chi-square test for categorical variables
      observed_table <- table(observed_complete)
      imputed_table <- table(imputed_values)

      # Combine levels that exist in both
      all_levels <- union(names(observed_table), names(imputed_table))
      obs_counts <- observed_table[all_levels]
      imp_counts <- imputed_table[all_levels]
      obs_counts[is.na(obs_counts)] <- 0
      imp_counts[is.na(imp_counts)] <- 0

      if (sum(obs_counts) > 0 && sum(imp_counts) > 0) {
        chisq_test <- tryCatch({
          stats::chisq.test(rbind(obs_counts, imp_counts))
        }, error = function(e) list(p.value = NA, statistic = NA))

        dist_metrics[[j]] <- list(
          chisq_statistic = chisq_test$statistic,
          chisq_p_value = chisq_test$p.value
        )
      }
    }
  }

  dist_metrics
}

#' Compute Correlation Preservation
#' @keywords internal
compute_correlation_preservation <- function(original_data, imputed_data) {
  # Compute correlations for numeric variables only
  num_cols <- sapply(original_data, is.numeric)

  if (sum(num_cols) < 2) {
    return(list(
      correlation_preservation = NA,
      correlation_difference = NA
    ))
  }

  original_corr <- tryCatch({
    cor(original_data[, num_cols], use = "pairwise.complete.obs")
  }, error = function(e) matrix(NA, sum(num_cols), sum(num_cols)))

  imputed_corr <- tryCatch({
    cor(imputed_data[, num_cols], use = "pairwise.complete.obs")
  }, error = function(e) matrix(NA, sum(num_cols), sum(num_cols)))

  # Compute correlation preservation metrics
  corr_diff <- abs(original_corr - imputed_corr)
  corr_diff <- corr_diff[upper.tri(corr_diff)]  # Upper triangle only

  list(
    mean_correlation_difference = mean(corr_diff, na.rm = TRUE),
    max_correlation_difference = max(corr_diff, na.rm = TRUE),
    correlation_rmse = sqrt(mean(corr_diff^2, na.rm = TRUE))
  )
}

# --- Comparison Method Implementations ---

#' Run FD::gowdis + PMM Imputation
#'
#' @param data Dataset with missing values
#' @param predictors Predictor variables
#' @param k Number of donors for PMM
#' @return Imputed values
#' @export
impute_with_fd_gowdis <- function(data, predictors, k = 5) {
  if (!requireNamespace("FD", quietly = TRUE)) {
    stop("FD package is required for this method")
  }

  # Compute Gower distances using FD package
  gower_dist <- FD::gowdis(predictors)

  # Apply PMM using the computed distances
  apply_pmm_with_distances(data, predictors, gower_dist, k)
}

#' Run VIM k-NN Imputation
#'
#' @param data Dataset with missing values
#' @param k Number of neighbors
#' @return Imputed dataset
#' @export
impute_with_vim_knn <- function(data, k = 5) {
  if (!requireNamespace("VIM", quietly = TRUE)) {
    stop("VIM package is required for this method")
  }

  VIM::kNN(data, k = k, imp_var = FALSE)
}

#' Run VIM IRMI Imputation
#'
#' @param data Dataset with missing values
#' @return Imputed dataset
#' @export
impute_with_vim_irmi <- function(data) {
  if (!requireNamespace("VIM", quietly = TRUE)) {
    stop("VIM package is required for this method")
  }

  VIM::irmi(data)
}

#' Apply PMM with Pre-computed Distances
#' @keywords internal
apply_pmm_with_distances <- function(data, predictors, distances, k = 5) {
  n <- nrow(data)
  imputed_data <- data

  # Find variables with missing values
  vars_with_missing <- colnames(data)[sapply(data, function(x) any(is.na(x)))]

  for (var in vars_with_missing) {
    missing_idx <- which(is.na(data[, var]))

    if (length(missing_idx) == 0) next

    observed_idx <- which(!is.na(data[, var]))

    for (idx in missing_idx) {
      # Find k nearest neighbors based on distance matrix
      distances_to_idx <- distances[idx, observed_idx]
      nearest_neighbors <- observed_idx[order(distances_to_idx)[1:k]]

      # Sample one donor randomly
      donor <- sample(nearest_neighbors, 1)
      imputed_data[idx, var] <- data[donor, var]
    }
  }

  imputed_data
}

# --- Computational Performance Evaluation ---

#' Evaluate Computational Performance
#'
#' @param method_name Name of the method
#' @param execution_time Time taken in seconds
#' @param memory_usage Memory used in MB
#' @param convergence_info Convergence information
#' @return List of performance metrics
#' @export
evaluate_computational_performance <- function(method_name, execution_time,
                                             memory_usage = NA, convergence_info = NULL) {
  list(
    method = method_name,
    execution_time_seconds = execution_time,
    memory_usage_mb = memory_usage,
    convergence_achieved = !is.null(convergence_info),
    convergence_info = convergence_info,
    efficiency_score = compute_efficiency_score(execution_time, memory_usage)
  )
}

#' Compute Efficiency Score
#' @keywords internal
compute_efficiency_score <- function(time, memory) {
  if (is.na(memory)) {
    # Time-based efficiency only
    1 / (1 + log(time + 1))
  } else {
    # Combined time and memory efficiency
    time_score <- 1 / (1 + log(time + 1))
    memory_score <- 1 / (1 + log(memory + 1))
    (time_score + memory_score) / 2
  }
}
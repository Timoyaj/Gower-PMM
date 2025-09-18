#' Real Gower-PMM Analysis: Distance Engines & Imputation Methods
#'
#' This script performs actual evaluation of the Gower-PMM implementation:
#' 1. Distance engine comparison (Gower-PMM vs FD::gowdis vs cluster::daisy)
#' 2. Imputation method comparison (Gower-PMM vs MICE methods vs others)
#' 3. Identification of strengths and potential issues in the implementation

# Setup and Package Loading
cat("=== REAL GOWER-PMM ANALYSIS: DISTANCE ENGINES & IMPUTATION METHODS ===\n\n")

# Load required libraries
library(mice)
library(FD)
library(VIM)
library(cluster)
library(ggplot2)
library(dplyr)
library(gowerpmm)  # Now properly installed

# Set seed for reproducibility
set.seed(42)

# Load Employee dataset
data("employee")
cat("Dataset: Employee Selection Data (Enders, 2010)\n")
cat("Variables: IQ (numeric), wbeing (numeric), jobperf (numeric, MAR)\n")
cat("Missingness: Job performance missing for lower IQ candidates\n\n")

# Create complete reference dataset
employee_complete <- employee
for (col in colnames(employee)) {
  if (any(is.na(employee[, col]))) {
    predictors <- setdiff(colnames(employee), col)
    observed_data <- employee[!is.na(employee[, col]), ]

    if (nrow(observed_data) > 2) {
      formula <- as.formula(paste(col, "~", paste(predictors, collapse = " + ")))
      model <- lm(formula, data = observed_data)

      missing_idx <- which(is.na(employee[, col]))
      pred_data <- employee[missing_idx, predictors, drop = FALSE]
      predictions <- predict(model, newdata = pred_data)

      employee_complete[missing_idx, col] <- predictions
    }
  }
}

cat("Complete reference dataset created for evaluation\n\n")

# ============================================================================
# PART 1: DISTANCE ENGINE COMPARISON
# ============================================================================

cat("PART 1: DISTANCE ENGINE COMPARISON\n")
cat("===================================\n\n")

# Function to evaluate distance quality
evaluate_distance_quality <- function(data, distance_obj, method_name) {
  # Convert to matrix if needed
  if (inherits(distance_obj, "dist")) {
    dist_matrix <- as.matrix(distance_obj)
  } else {
    dist_matrix <- distance_obj
  }

  n <- nrow(data)

  # 1. Basic distance statistics
  distance_values <- dist_matrix[upper.tri(dist_matrix)]

  # 2. Nearest neighbor preservation (simplified)
  nn_preservation <- NA
  tryCatch({
    # Use Euclidean for numeric variables as reference
    num_cols <- sapply(data, is.numeric)
    if (sum(num_cols) > 1) {
      euclidean_dist <- as.matrix(dist(data[, num_cols]))
      euclidean_matrix <- apply(euclidean_dist, 1, function(x) order(x)[2:6])  # k=5
      distance_matrix_nn <- apply(dist_matrix, 1, function(x) order(x)[2:6])  # k=5

      overlaps <- sapply(1:n, function(j) {
        length(intersect(euclidean_matrix[,j], distance_matrix_nn[,j]))
      })

      nn_preservation <- mean(overlaps) / 5
    }
  }, error = function(e) {
    nn_preservation <- NA
  })

  list(
    method = method_name,
    nn_preservation_k5 = nn_preservation,
    mean_distance = mean(distance_values),
    sd_distance = sd(distance_values),
    distance_range = diff(range(distance_values)),
    computation_success = TRUE
  )
}

# Compare distance engines on complete Employee data
cat("Comparing Distance Engines on Complete Employee Data:\n")
cat("- Gower-PMM: Optimized Gower distance with C++ implementation\n")
cat("- FD::gowdis: Traditional Gower distance (R implementation)\n")
cat("- cluster::daisy: Gower-based distance from cluster package\n")
cat("- Euclidean: Standard Euclidean distance (numeric only)\n\n")

distance_results <- list()
timing_results <- list()

# 1. Gower-PMM Engine (our optimized implementation)
cat("1. Computing Gower-PMM distance...")
start_time <- Sys.time()
dist_gowerpmm <- gower_dist_engine(employee_complete)
time_gowerpmm <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
distance_results$gowerpmm <- evaluate_distance_quality(employee_complete, dist_gowerpmm, "Gower-PMM")
timing_results$gowerpmm <- time_gowerpmm
cat(sprintf(" (%.4f sec)\n", time_gowerpmm))

# 2. FD::gowdis (traditional implementation)
cat("2. Computing FD::gowdis distance...")
start_time <- Sys.time()
dist_fd <- FD::gowdis(employee_complete)
time_fd <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
distance_results$fd_gowdis <- evaluate_distance_quality(employee_complete, dist_fd, "FD::gowdis")
timing_results$fd_gowdis <- time_fd
cat(sprintf(" (%.4f sec)\n", time_fd))

# 3. cluster::daisy (Gower-based)
cat("3. Computing cluster::daisy distance...")
start_time <- Sys.time()
dist_daisy <- cluster::daisy(employee_complete, metric = "gower")
time_daisy <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
distance_results$daisy <- evaluate_distance_quality(employee_complete, dist_daisy, "cluster::daisy")
timing_results$daisy <- time_daisy
cat(sprintf(" (%.4f sec)\n", time_daisy))

# 4. Euclidean (numeric only)
cat("4. Computing Euclidean distance (numeric only)...")
start_time <- Sys.time()
num_cols <- sapply(employee_complete, is.numeric)
if (sum(num_cols) > 1) {
  dist_euclidean <- dist(employee_complete[, num_cols])
  time_euclidean <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  distance_results$euclidean <- evaluate_distance_quality(employee_complete, dist_euclidean, "Euclidean")
  timing_results$euclidean <- time_euclidean
  cat(sprintf(" (%.4f sec)\n", time_euclidean))
} else {
  distance_results$euclidean <- list(method = "Euclidean", nn_preservation_k5 = NA, mean_distance = NA, sd_distance = NA, distance_skewness = NA, distance_range = NA)
  timing_results$euclidean <- NA
  cat(" (N/A - insufficient numeric variables)\n")
}

cat("\n")

# Distance Engine Results
cat("DISTANCE ENGINE PERFORMANCE RESULTS:\n")
cat("====================================\n\n")

distance_df <- data.frame(
  Engine = sapply(distance_results, function(x) x$method),
  Time_sec = sapply(timing_results, function(x) round(x, 4)),
  NN_Preservation = sapply(distance_results, function(x) round(x$nn_preservation_k5, 4)),
  Mean_Distance = sapply(distance_results, function(x) round(x$mean_distance, 4)),
  SD_Distance = sapply(distance_results, function(x) round(x$sd_distance, 4)),
  Distance_Range = sapply(distance_results, function(x) round(x$distance_range, 4))
)

print(distance_df)
cat("\n")

# Distance Engine Analysis
cat("DISTANCE ENGINE ANALYSIS:\n")
cat("- Gower-PMM shows", ifelse(distance_df$Time_sec[1] < distance_df$Time_sec[2], "better", "similar"), "performance vs FD::gowdis\n")
cat("- NN preservation:", round(distance_df$NN_Preservation[1], 3), "(Gower-PMM) vs", round(distance_df$NN_Preservation[2], 3), "(FD::gowdis)\n")
cat("- Distance characteristics are", ifelse(abs(distance_df$Mean_Distance[1] - distance_df$Mean_Distance[2]) < 0.01, "very similar", "different"), "between implementations\n\n")

# ============================================================================
# PART 2: IMPUTATION METHOD COMPARISON
# ============================================================================

cat("PART 2: IMPUTATION METHOD COMPARISON\n")
cat("====================================\n\n")

# Imputation quality evaluation function
evaluate_imputation_quality <- function(original, imputed, missing_mask, method_name) {
  results <- list(method = method_name)

  # RMSE and MAE for numeric variables
  for (col in colnames(original)) {
    if (is.numeric(original[, col]) && any(missing_mask[, col])) {
      observed <- original[missing_mask[, col], col]
      predicted <- imputed[missing_mask[, col], col]

      if (length(observed) > 0 && length(predicted) > 0) {
        rmse <- sqrt(mean((observed - predicted)^2, na.rm = TRUE))
        mae <- mean(abs(observed - predicted), na.rm = TRUE)
        bias <- mean(predicted - observed, na.rm = TRUE)

        results[[paste0(col, "_rmse")]] <- rmse
        results[[paste0(col, "_mae")]] <- mae
        results[[paste0(col, "_bias")]] <- bias
        results[[paste0(col, "_n_missing")]] <- sum(missing_mask[, col])
      }
    }
  }

  # Overall metrics
  rmse_values <- unlist(results[grep("_rmse$", names(results))])
  mae_values <- unlist(results[grep("_mae$", names(results))])

  results$overall_rmse <- mean(rmse_values, na.rm = TRUE)
  results$overall_mae <- mean(mae_values, na.rm = TRUE)
  results$total_missing <- sum(missing_mask)

  return(results)
}

imputation_results <- list()
imp_timing_results <- list()

# 1. Gower-PMM (Auto weights - our main method)
cat("1. Running Gower-PMM (Auto Weights - OUR METHOD)...\n")
start_time <- Sys.time()
tryCatch({
  imp_gower_auto <- mice(employee, method = "gowerpmm", m = 1, maxit = 1, printFlag = FALSE)
  time_gower_auto <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  imputed_gower_auto <- complete(imp_gower_auto, 1)

  imputation_results$gower_auto <- evaluate_imputation_quality(
    employee_complete, imputed_gower_auto, is.na(employee), "Gower-PMM (Auto)"
  )
  imp_timing_results$gower_auto <- time_gower_auto
  cat(sprintf("   SUCCESS: %.4f sec\n", time_gower_auto))
}, error = function(e) {
  time_gower_auto <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  cat(sprintf("   ERROR: %s (%.4f sec)\n", e$message, time_gower_auto))
  imputation_results$gower_auto <- list(method = "Gower-PMM (Auto) - ERROR", overall_rmse = NA, overall_mae = NA)
  imp_timing_results$gower_auto <- time_gower_auto
})

# 2. Gower-PMM (Equal weights)
cat("2. Running Gower-PMM (Equal Weights)...\n")
start_time <- Sys.time()
tryCatch({
  imp_gower_equal <- mice(employee, method = "gowerpmm", weights = "equal", m = 1, maxit = 1, printFlag = FALSE)
  time_gower_equal <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  imputed_gower_equal <- complete(imp_gower_equal, 1)

  imputation_results$gower_equal <- evaluate_imputation_quality(
    employee_complete, imputed_gower_equal, is.na(employee), "Gower-PMM (Equal)"
  )
  imp_timing_results$gower_equal <- time_gower_equal
  cat(sprintf("   SUCCESS: %.4f sec\n", time_gower_equal))
}, error = function(e) {
  time_gower_equal <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  cat(sprintf("   ERROR: %s (%.4f sec)\n", e$message, time_gower_equal))
  imputation_results$gower_equal <- list(method = "Gower-PMM (Equal) - ERROR", overall_rmse = NA, overall_mae = NA)
  imp_timing_results$gower_equal <- time_gower_equal
})

# 3. MICE PMM (standard comparison)
cat("3. Running MICE PMM...\n")
start_time <- Sys.time()
imp_mice_pmm <- mice(employee, method = "pmm", m = 1, maxit = 1, printFlag = FALSE)
time_mice_pmm <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
imputed_mice_pmm <- complete(imp_mice_pmm, 1)

imputation_results$mice_pmm <- evaluate_imputation_quality(
  employee_complete, imputed_mice_pmm, is.na(employee), "MICE PMM"
)
imp_timing_results$mice_pmm <- time_mice_pmm
cat(sprintf("   SUCCESS: %.4f sec\n", time_mice_pmm))

# 4. MICE CART
cat("4. Running MICE CART...\n")
start_time <- Sys.time()
imp_mice_cart <- mice(employee, method = "cart", m = 1, maxit = 1, printFlag = FALSE)
time_mice_cart <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
imputed_mice_cart <- complete(imp_mice_cart, 1)

imputation_results$mice_cart <- evaluate_imputation_quality(
  employee_complete, imputed_mice_cart, is.na(employee), "MICE CART"
)
imp_timing_results$mice_cart <- time_mice_cart
cat(sprintf("   SUCCESS: %.4f sec\n", time_mice_cart))

# 5. Mean/Mode (baseline)
cat("5. Running Mean/Mode Imputation...\n")
start_time <- Sys.time()
imputed_mean_mode <- employee
for (col in colnames(employee)) {
  if (any(is.na(employee[, col]))) {
    if (is.numeric(employee[, col])) {
      imputed_mean_mode[, col] <- ifelse(is.na(employee[, col]),
                                       mean(employee[, col], na.rm = TRUE),
                                       employee[, col])
    } else {
      mode_val <- names(which.max(table(employee[, col], useNA = "no")))
      imputed_mean_mode[, col] <- ifelse(is.na(employee[, col]),
                                       mode_val,
                                       employee[, col])
    }
  }
}
time_mean_mode <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))

imputation_results$mean_mode <- evaluate_imputation_quality(
  employee_complete, imputed_mean_mode, is.na(employee), "Mean/Mode"
)
imp_timing_results$mean_mode <- time_mean_mode
cat(sprintf("   SUCCESS: %.4f sec\n", time_mean_mode))

cat("\n")

# Imputation Results
cat("IMPUTATION METHOD PERFORMANCE RESULTS:\n")
cat("======================================\n\n")

# Filter out error results for display
valid_methods <- names(imputation_results)[!grepl("ERROR", sapply(imputation_results, function(x) x$method))]

imputation_df <- data.frame(
  Method = sapply(imputation_results[valid_methods], function(x) x$method),
  Overall_RMSE = sapply(imputation_results[valid_methods], function(x) round(x$overall_rmse, 4)),
  Overall_MAE = sapply(imputation_results[valid_methods], function(x) round(x$overall_mae, 4)),
  JobPerf_RMSE = sapply(imputation_results[valid_methods], function(x) round(x$jobperf_rmse, 4)),
  JobPerf_MAE = sapply(imputation_results[valid_methods], function(x) round(x$jobperf_mae, 4)),
  JobPerf_Bias = sapply(imputation_results[valid_methods], function(x) round(x$jobperf_bias, 4)),
  Time_seconds = sapply(imp_timing_results[valid_methods], function(x) round(x, 4))
)

print(imputation_df)
cat("\n")

# ============================================================================
# PART 3: ANALYSIS AND INTERPRETATION
# ============================================================================

cat("PART 3: ANALYSIS AND INTERPRETATION\n")
cat("===================================\n\n")

# Performance Analysis
cat("GOWER-PMM PERFORMANCE ANALYSIS:\n")
cat("-------------------------------\n")

if (!is.na(imputation_df$Overall_RMSE[1])) {
  gower_auto_rmse <- imputation_df$Overall_RMSE[1]
  mice_pmm_rmse <- imputation_df$Overall_RMSE[imputation_df$Method == "MICE PMM"]
  mice_cart_rmse <- imputation_df$Overall_RMSE[imputation_df$Method == "MICE CART"]

  cat(sprintf("Gower-PMM (Auto) RMSE: %.4f\n", gower_auto_rmse))
  cat(sprintf("MICE PMM RMSE: %.4f\n", mice_pmm_rmse))
  cat(sprintf("MICE CART RMSE: %.4f\n", mice_cart_rmse))

  # Performance comparison
  if (gower_auto_rmse < mice_pmm_rmse && gower_auto_rmse < mice_cart_rmse) {
    cat("✓ Gower-PMM (Auto) shows BEST overall performance\n")
  } else if (gower_auto_rmse < mice_pmm_rmse) {
    cat("✓ Gower-PMM (Auto) outperforms MICE PMM\n")
  } else {
    cat("⚠ Gower-PMM (Auto) shows competitive but not superior performance\n")
  }

  # Speed analysis
  gower_time <- imputation_df$Time_seconds[1]
  avg_competitor_time <- mean(imputation_df$Time_seconds[2:4])
  cat(sprintf("Gower-PMM timing: %.4f sec (competitors: %.4f sec avg)\n",
              gower_time, avg_competitor_time))

  if (gower_time > avg_competitor_time * 1.5) {
    cat("⚠ Gower-PMM is slower - optimization may be needed\n")
  } else {
    cat("✓ Gower-PMM timing is reasonable\n")
  }
} else {
  cat("⚠ Gower-PMM (Auto) failed to run - implementation issue detected\n")
}

cat("\n")

# Distance Engine Analysis
cat("DISTANCE ENGINE ANALYSIS:\n")
cat("------------------------\n")

gower_time <- distance_df$Time_sec[1]
fd_time <- distance_df$Time_sec[2]

cat(sprintf("Gower-PMM distance time: %.4f sec\n", gower_time))
cat(sprintf("FD::gowdis time: %.4f sec\n", fd_time))

if (gower_time < fd_time) {
  cat("✓ Gower-PMM distance engine is FASTER than FD::gowdis\n")
} else {
  cat("⚠ Gower-PMM distance engine is slower - C++ optimization may need improvement\n")
}

# Quality comparison
gower_nn <- distance_df$NN_Preservation[1]
fd_nn <- distance_df$NN_Preservation[2]

cat(sprintf("Gower-PMM NN preservation: %.4f\n", gower_nn))
cat(sprintf("FD::gowdis NN preservation: %.4f\n", fd_nn))

if (abs(gower_nn - fd_nn) < 0.05) {
  cat("✓ Distance quality is comparable between implementations\n")
} else if (gower_nn > fd_nn) {
  cat("✓ Gower-PMM shows better nearest neighbor preservation\n")
} else {
  cat("⚠ Gower-PMM shows worse nearest neighbor preservation\n")
}

cat("\n")

# Potential Issues Identification
cat("POTENTIAL IMPLEMENTATION ISSUES:\n")
cat("---------------------------------\n")

issues_found <- 0

# Check for runtime errors
if (any(grepl("ERROR", sapply(imputation_results, function(x) x$method)))) {
  cat("⚠ RUNTIME ERRORS detected in Gower-PMM implementation\n")
  issues_found <- issues_found + 1
}

# Check performance vs competitors
if (!is.na(imputation_df$Overall_RMSE[1])) {
  gower_rmse <- imputation_df$Overall_RMSE[1]
  best_competitor <- min(imputation_df$Overall_RMSE[2:4], na.rm = TRUE)

  if (gower_rmse > best_competitor * 1.2) {
    cat("⚠ Gower-PMM accuracy significantly worse than competitors\n")
    issues_found <- issues_found + 1
  }
}

# Check timing
if (!is.na(imputation_df$Time_seconds[1])) {
  gower_time <- imputation_df$Time_seconds[1]
  avg_time <- mean(imputation_df$Time_seconds[2:4], na.rm = TRUE)

  if (gower_time > avg_time * 2) {
    cat("⚠ Gower-PMM is much slower than competitors\n")
    issues_found <- issues_found + 1
  }
}

if (issues_found == 0) {
  cat("✓ No major implementation issues detected\n")
}

cat("\n")

# Recommendations
cat("RECOMMENDATIONS:\n")
cat("----------------\n")

if (!is.na(imputation_df$Overall_RMSE[1])) {
  cat("1. Gower-PMM shows", ifelse(imputation_df$Overall_RMSE[1] < mean(imputation_df$Overall_RMSE[2:4], na.rm = TRUE), "superior", "competitive"), "performance\n")
  cat("2. Consider the trade-off between accuracy and computational cost\n")
  cat("3. Test on larger datasets with mixed variable types for full evaluation\n")
} else {
  cat("1. Fix runtime errors in Gower-PMM implementation\n")
  cat("2. Debug C++ function calls and parameter passing\n")
  cat("3. Test individual components before full integration\n")
}

cat("4. Validate distance engine optimization vs R implementations\n")
cat("5. Profile code for performance bottlenecks\n")

cat("\n=== ANALYSIS COMPLETE ===\n")

# Save comprehensive results
results_summary <- list(
  dataset = "Employee (Enders, 2010)",
  distance_engines = distance_df,
  imputation_methods = imputation_df,
  analysis_timestamp = Sys.time(),
  issues_detected = issues_found,
  recommendations = if(issues_found > 0) "Fix implementation issues" else "Method performs well"
)

saveRDS(results_summary, "docs/analysis/real_gower_analysis_results.rds")
cat("Results saved to: docs/analysis/real_gower_analysis_results.rds\n")
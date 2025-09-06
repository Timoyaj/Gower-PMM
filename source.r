# =============================================================================
# Two-Phase Sampling Methodology with Gower's Pre-Selection and LASSO
# =============================================================================
#
# This script implements a methodology for optimal selection and variance
# estimation with multiple auxiliary variables in two-phase sampling. It uses
# Gower's distance for pre-selection of candidate variables followed by
# survey-weighted LASSO for final variable selection.
#
# Author: [Your Name]
# Date: September 6, 2025
# License: MIT
# =============================================================================

# Required libraries
library(glmnet)  # For LASSO regression
library(cluster) # For Gower's distance
library(survey)  # For survey-weighted analysis
library(dplyr)   # For data manipulation

# Set seed for reproducibility
set.seed(2025)


# 2. HELPER FUNCTION: GOWER'S PRE-SELECTION (STAGE 1)
# ----------------------------------------------------
# This function takes a dataframe of potential auxiliary variables, determines
# the optimal number of clusters using the silhouette width method, performs
# Partitioning Around Medoids (PAM), and returns the names of the medoids
# (representative variables) for each cluster.

gower_preselect <- function(data, max_clusters = 10) {
  # Ensure there are enough unique data points to cluster
  if (nrow(unique(data)) < 2) {
    warning("Not enough unique data points for clustering. Returning all variables.")
    return(colnames(data))
  }
  
  # Calculate Gower's distance matrix for mixed data types
  gower_dist <- daisy(data, metric = "gower")
  
  # Find the optimal number of clusters (k) using average silhouette width
  sil_width <- c(NA)
  # Iterate from 2 to a reasonable maximum number of clusters
  for (k in 2:min(max_clusters, nrow(data)-1)) {
    pam_fit <- pam(gower_dist, diss = TRUE, k = k)
    sil_width[k] <- pam_fit$silinfo$avg.width
  }
  
  # Select the k with the highest average silhouette width
  optimal_k <- which.max(sil_width)
  if (length(optimal_k) == 0 || is.na(optimal_k) || optimal_k < 2) {
    warning("Could not determine optimal clusters. Returning all variables.")
    return(colnames(data))
  }

  # Perform PAM clustering with the optimal k
  pam_final <- pam(gower_dist, diss = TRUE, k = optimal_k)
  
  # Get the names of the medoid variables
  selected_vars <- colnames(data)[pam_final$medoids]
  
  cat("Gower Pre-selection: \n")
  cat(" - Optimal number of clusters (k):", optimal_k, "\n")
  cat(" - Selected representative variables:", paste(selected_vars, collapse=", "), "\n\n")
  
  return(selected_vars)
}


# 3. MAIN FUNCTION: TWO-PHASE SELECTION & ESTIMATION
# ---------------------------------------------------
# This is the main wrapper function that executes the entire proposed methodology.
TwoPhaseSelect <- function(y_var, 
                           aux_vars, 
                           phase1_data, 
                           phase2_data, 
                           id_var_p1,
                           id_var_p2,
                           p1_wt_var,
                           p2_wt_var,
                           B = 500) {

  # --- PREPARATION ---
  # Combine weights to get the final full weight for phase 2 units
  phase2_data$final_weight <- phase2_data[[p1_wt_var]] * phase2_data[[p2_wt_var]]
  
  # --- STAGE 1: GOWER'S PRE-SELECTION ---
  # Perform on the first phase data to find candidate variables
  candidate_vars <- gower_preselect(phase1_data[aux_vars])
  
  # --- STAGE 2: SURVEY-WEIGHTED LASSO ---
  y <- phase2_data[[y_var]]
  X <- model.matrix(as.formula(paste("~", paste(candidate_vars, collapse = " + "))), data = phase2_data)[, -1]
  weights <- phase2_data$final_weight
  
  # Fit cross-validated LASSO using the final weights
  cv_lasso_fit <- cv.glmnet(X, y, weights = weights, alpha = 1)
  
  # Extract coefficients at lambda.1se for a more parsimonious model
  lasso_coefs <- coef(cv_lasso_fit, s = "lambda.1se")
  selected_vars_final <- rownames(lasso_coefs)[lasso_coefs[, 1] != 0]
  selected_vars_final <- selected_vars_final[selected_vars_final != "(Intercept)"]
  
  cat("Survey-Weighted LASSO Selection: \n")
  if (length(selected_vars_final) > 0) {
     cat(" - Final selected variables:", paste(selected_vars_final, collapse=", "), "\n\n")
  } else {
     cat(" - No variables selected by LASSO. Estimating population mean.\n\n")
  }

  # --- POINT ESTIMATION ---
  # Estimate the population mean using the Hajek-type difference estimator
  
  # Population size (estimated from sum of phase 1 weights)
  N_hat <- sum(phase1_data[[p1_wt_var]])
  
  # Initialize point estimate with the weighted mean from phase 2
  # This is the base of the estimator
  point_estimate <- sum(weights * y) / sum(weights)
  
  if (length(selected_vars_final) > 0) {
      # Extract coefficients for selected variables
      final_betas <- lasso_coefs[selected_vars_final, 1]
      
      # Calculate the weighted sum of X for phase 1 and phase 2 samples
      p1_X_sum <- colSums(phase1_data[selected_vars_final] * phase1_data[[p1_wt_var]])
      p2_X_sum <- colSums(phase2_data[selected_vars_final] * weights)
      
      # Add the adjustment term to the point estimate
      adjustment <- sum((p1_X_sum - p2_X_sum) * final_betas)
      point_estimate <- (sum(weights * y) + adjustment) / N_hat
  }


  # --- VARIANCE ESTIMATION (BOOTSTRAP) ---
  cat("Starting Bootstrap for Variance Estimation (B =", B, "replicates)...\n")
  bootstrap_estimates <- numeric(B)
  
  # Get unique first-phase IDs
  p1_ids <- unique(phase1_data[[id_var_p1]])

  for (b in 1:B) {
    # 1. Resample phase 1 units with replacement
    boot_p1_ids <- sample(p1_ids, size = length(p1_ids), replace = TRUE)
    
    # Create the bootstrap phase 1 sample by selecting rows
    boot_p1_sample <- phase1_data[match(boot_p1_ids, phase1_data[[id_var_p1]]), ]
    
    # 2. Resample phase 2 units within the resampled phase 1 units
    boot_p2_sample <- phase2_data %>%
      filter(.data[[id_var_p2]] %in% boot_p1_ids) %>%
      group_by(.data[[id_var_p2]]) %>%
      sample_n(size = n(), replace = TRUE) %>%
      ungroup()

    # --- REPEAT THE ENTIRE METHODOLOGY ON THE BOOTSTRAP SAMPLE ---
    tryCatch({
      # Stage 1: Gower
      boot_candidate_vars <- gower_preselect(boot_p1_sample[aux_vars])
      
      # Stage 2: LASSO
      boot_y <- boot_p2_sample[[y_var]]
      boot_X <- model.matrix(as.formula(paste("~", paste(boot_candidate_vars, collapse = " + "))), data = boot_p2_sample)[, -1]
      boot_weights <- boot_p2_sample$final_weight
      
      boot_cv_lasso <- cv.glmnet(boot_X, boot_y, weights = boot_weights, alpha = 1)
      boot_coefs <- coef(boot_cv_lasso, s = "lambda.1se")
      boot_selected_vars <- rownames(boot_coefs)[boot_coefs[, 1] != 0]
      boot_selected_vars <- boot_selected_vars[boot_selected_vars != "(Intercept)"]
      
      # Estimation
      boot_N_hat <- sum(boot_p1_sample[[p1_wt_var]])
      boot_est <- sum(boot_weights * boot_y) / sum(boot_weights)
      
      if (length(boot_selected_vars) > 0) {
        boot_betas <- boot_coefs[boot_selected_vars, 1]
        boot_p1_X_sum <- colSums(boot_p1_sample[boot_selected_vars] * boot_p1_sample[[p1_wt_var]])
        boot_p2_X_sum <- colSums(boot_p2_sample[boot_selected_vars] * boot_weights)
        boot_adjustment <- sum((boot_p1_X_sum - boot_p2_X_sum) * boot_betas)
        boot_est <- (sum(boot_weights * boot_y) + boot_adjustment) / boot_N_hat
      }
      
      bootstrap_estimates[b] <- boot_est
    }, error = function(e) {
      # In case of errors in a bootstrap replicate (e.g., no variables selected)
      # assign NA and handle later
      bootstrap_estimates[b] <- NA
    })
    
    # Progress indicator
    if (b %% 50 == 0) cat("  ... completed", b, "replicates\n")
  }
  
  # Calculate variance from the bootstrap estimates, ignoring any NAs
  variance_estimate <- var(bootstrap_estimates, na.rm = TRUE)
  
  # --- RESULTS ---
  se <- sqrt(variance_estimate)
  ci_lower <- point_estimate - 1.96 * se
  ci_upper <- point_estimate + 1.96 * se
  
  results <- list(
    point_estimate = point_estimate,
    variance_estimate = variance_estimate,
    standard_error = se,
    confidence_interval_95 = c(lower = ci_lower, upper = ci_upper),
    selected_variables = selected_vars_final
  )
  
  return(results)
}


# 4. EXAMPLE USAGE: California API Dataset
# -----------------------------------------
# The `apiclus2` dataset in the `survey` package comes from a two-phase
# cluster sample of California schools.
# Phase 1: 15 schools (clusters) were sampled.
# Phase 2: Students were sampled from within those 15 schools.

cat("\n\n--- Running Example on California API Dataset ---\n\n")

# Load the data
data(api)

# --- DATA PREPARATION ---
# Phase 2 sample (students from selected schools)
phase2_student_data <- apiclus2

# Phase 1 sample (the 15 unique schools selected)
phase1_school_data <- apisrs %>% 
  filter(dnum %in% unique(phase2_student_data$dnum))

# Define variables
# Response variable: Academic Performance Index in 2000
response_variable <- "api00"

# Potential auxiliary variables (school-level)
# Note: 'stype' is categorical, others are numeric. Our method handles this.
auxiliary_variables <- c("stype", "api99", "meals", "ell", "avg.ed", "col.grad", "full", "emer")

# Unique identifiers for schools
id_p1 <- "dnum" # School ID in phase 1 data
id_p2 <- "dnum" # School ID in phase 2 data

# --- WEIGHT CALCULATION ---
# Phase 1 weights (1 / probability of selecting a school)
# Total schools in population = 4421 (from apipop)
# Number of schools sampled = 15
N_schools <- 4421
n1_schools <- 15
phase1_school_data$p1_wt <- N_schools / n1_schools

# Phase 2 weights (1 / probability of selecting student WITHIN a school)
# This is provided in the dataset as `pw`
phase2_student_data$p2_wt <- phase2_student_data$pw

# Add phase 1 weights to the phase 2 student data by merging
phase2_student_data <- merge(phase2_student_data, phase1_school_data[, c("dnum", "p1_wt")], by = "dnum")


# --- RUN THE ANALYSIS ---
# For this example, let's use a smaller number of bootstrap reps for speed.
# For a real analysis, B should be 500 or higher.
api_results <- TwoPhaseSelect(
  y_var = response_variable,
  aux_vars = auxiliary_variables,
  phase1_data = phase1_school_data,
  phase2_data = phase2_student_data,
  id_var_p1 = id_p1,
  id_var_p2 = id_p2,
  p1_wt_var = "p1_wt",
  p2_wt_var = "p2_wt",
  B = 100 # Reduced for quick demonstration
)


# --- DISPLAY RESULTS ---
cat("\n--- FINAL RESULTS ---\n")
cat("Final Selected Auxiliary Variables:\n")
if (length(api_results$selected_variables) > 0) {
  print(api_results$selected_variables)
} else {
  cat("None\n")
}
cat("\n")
cat("Estimated Mean Population", toupper(response_variable), ":", round(api_results$point_estimate, 2), "\n")
cat("Bootstrap Standard Error:", round(api_results$standard_error, 2), "\n")
cat("95% Confidence Interval: [", round(api_results$confidence_interval_95[1], 2), ",", 
    round(api_results$confidence_interval_95[2], 2), "]\n")
cat("----------------------\n")

# =============================================================================
# Gower-kPrototypes Predictive Mean Matching (GKP-PMM) Imputation Method
# =============================================================================
#
# This script implements a custom imputation method for the MICE package that
# integrates Gower's distance for refined donor selection within a Predictive
# Mean Matching (PMM) framework. It supports mixed data types and includes
# optional pre-clustering for improved performance.
#
# Author: [Your Name]
# Date: September 6, 2025
# License: MIT
# =============================================================================

# Required libraries
library(mice)
library(cluster)      # For Gower's distance
library(clustMixType) # For k-prototypes clustering
library(MASS)         # For ordered logistic regression
library(nnet)         # For multinomial regression

# =============================================================================
# Main Imputation Function
# =============================================================================

#' Custom Imputation Function: Enhanced Gower-kPrototypes Predictive Mean Matching
#'
#' This function implements a custom imputation method that integrates Gower's distance
#' for refined donor selection within a PMM framework. This final, robust version validates
#' model compatibility to prevent crashes from incorrect user specifications.
#'
#' @param y The vector of the target variable (with NAs).
#' @param ry A logical vector indicating observed (TRUE) or missing (FALSE) values in `y`.
#' @param x A matrix of predictor variables for `y`.
#' @param donors The number of nearest neighbors (donors) for final imputation.
#' @param k_pre_clusters The number of clusters for optional k-prototypes pre-clustering.
#'                       Set to 0 to disable.
#' @param predictive_model A string specifying the predictive model for the initial PMM step.
#'                         Options: "auto", "lm", "logit", "polr", "multinom".
#'                         Default is "auto".
#' @param pmm_pool_factor A numeric factor to determine the size of the initial PMM candidate pool.
#'                        Default is 5.
#' @param ... Additional arguments.
#'
#' @return A vector of imputed values for the missing entries in `y`.
#'
#' @examples
#' # Example usage with MICE
#' library(mice)
#' data(nhanes)
#' imp <- mice(nhanes, method = "gkp_pmm", donors = 5, k_pre_clusters = 3)
#'
mice.impute.gkp_pmm <- function(y, ry, x, donors = 5, k_pre_clusters = 0,
                                predictive_model = "auto", pmm_pool_factor = 5, ...) {
  # ... [Paste the full, corrected gkp_pmm function code here] ...
    # --- 1. Prepare data for imputation ---
  x_donors_df <- as.data.frame(x[ry, , drop = FALSE])
  y_donors <- y[ry]
  x_recipients_df <- as.data.frame(x[!ry, , drop = FALSE])
  
  if (nrow(x_recipients_df) == 0) {
    return(y[!ry])
  }

  # --- 2. Train a predictive model ---
  model_type <- predictive_model
  if (model_type == "auto") {
    if (is.numeric(y)) model_type <- "lm"
    else if (is.factor(y) && nlevels(y) == 2) model_type <- "logit"
    else if (is.ordered(y)) model_type <- "polr"
    else if (is.factor(y) && !is.ordered(y) && nlevels(y) > 2) model_type <- "multinom"
    else model_type <- "lm" # Fallback for other data types
  }

  # *** NEW: Validate user-specified model against y's data type ***
  # If the model is incompatible, warn the user and fall back to "lm".
  if (model_type == "polr" && !is.ordered(y)) {
    warning(paste0("'polr' model requested for non-ordered variable. Falling back to 'lm'."), call. = FALSE)
    model_type <- "lm"
  } else if (model_type == "logit" && (!is.factor(y) || nlevels(y) != 2)) {
    warning(paste0("'logit' model requested for a non-binary variable. Falling back to 'lm'."), call. = FALSE)
    model_type <- "lm"
  } else if (model_type == "multinom" && (!is.factor(y) || nlevels(y) <= 2)) {
    warning(paste0("'multinom' model requested for a variable that is not a factor with >2 levels. Falling back to 'lm'."), call. = FALSE)
    model_type <- "lm"
  } else if (model_type %in% c("polr", "logit", "multinom") && !is.factor(y)) {
      warning(paste0("'", model_type, "' model requested for a non-factor variable. Falling back to 'lm'."), call. = FALSE)
      model_type <- "lm"
  }


  fit_data <- data.frame(y_target = y_donors)
  fit_data <- cbind(fit_data, x_donors_df)
  
  if (length(unique(y_donors)) < 2) {
      warning("Target 'y' has fewer than 2 unique levels. Falling back to random sampling from observed.", call. = FALSE)
      return(sample(y_donors, nrow(x_recipients_df), replace = TRUE))
  }

  switch(model_type,
    "lm" = {
      fit_data$y_target <- as.numeric(fit_data$y_target)
      if (var(fit_data$y_target, na.rm = TRUE) == 0) {
        return(sample(y_donors, nrow(x_recipients_df), replace = TRUE))
      }
      fit <- lm(y_target ~ ., data = fit_data)
      y_pred_donors <- predict(fit, newdata = x_donors_df)
      y_pred_recipients <- predict(fit, newdata = x_recipients_df)
    },
    "logit" = {
      fit <- glm(y_target ~ ., data = fit_data, family = "binomial")
      y_pred_donors <- predict(fit, newdata = x_donors_df, type = "response")
      y_pred_recipients <- predict(fit, newdata = x_recipients_df, type = "response")
    },
    "polr" = {
      fit <- MASS::polr(y_target ~ ., data = fit_data)
      ranks <- 1:nlevels(y_donors)
      pred_probs_donors <- predict(fit, newdata = x_donors_df, type = "probs")
      pred_probs_recipients <- predict(fit, newdata = x_recipients_df, type = "probs")
      y_pred_donors <- as.numeric(pred_probs_donors %*% ranks)
      y_pred_recipients <- as.numeric(pred_probs_recipients %*% ranks)
    },
    "multinom" = {
      fit <- nnet::multinom(y_target ~ ., data = fit_data, trace = FALSE)
      ranks <- 1:nlevels(y_donors)
      pred_probs_donors <- predict(fit, newdata = x_donors_df, type = "probs")
      pred_probs_recipients <- predict(fit, newdata = x_recipients_df, type = "probs")
      if (is.vector(pred_probs_donors)) pred_probs_donors <- matrix(pred_probs_donors, nrow = 1)
      if (is.vector(pred_probs_recipients)) pred_probs_recipients <- matrix(pred_probs_recipients, nrow = 1)
      y_pred_donors <- as.numeric(pred_probs_donors %*% ranks)
      y_pred_recipients <- as.numeric(pred_probs_recipients %*% ranks)
    },
    {
      warning(paste("Invalid predictive_model '", model_type, "'. Defaulting to 'lm'.", sep=""), call. = FALSE)
      fit_data$y_target <- as.numeric(fit_data$y_target)
      if (var(fit_data$y_target, na.rm = TRUE) == 0) {
        return(sample(y_donors, nrow(x_recipients_df), replace = TRUE))
      }
      fit <- lm(y_target ~ ., data = fit_data)
      y_pred_donors <- predict(fit, newdata = x_donors_df)
      y_pred_recipients <- predict(fit, newdata = x_recipients_df)
    }
  )

  # --- 3. Optional Pre-Clustering ---
  if (k_pre_clusters > 0 && nrow(x_donors_df) >= k_pre_clusters) {
    data_for_clustering <- rbind(x_donors_df, x_recipients_df)
    for(col_name in colnames(data_for_clustering)) {
      if (is.character(data_for_clustering[[col_name]])) {
        data_for_clustering[[col_name]] <- as.factor(data_for_clustering[[col_name]])
      }
    }
    kproto_fit <- clustMixType::kproto(data_for_clustering, k = k_pre_clusters, nstart = 1, type = "gower", verbose = FALSE)
    donor_clusters <- kproto_fit$cluster[1:nrow(x_donors_df)]
    recipient_clusters <- kproto_fit$cluster[(nrow(x_donors_df) + 1):nrow(data_for_clustering)]
    filtered_x_donors_list <- vector("list", nrow(x_recipients_df)); filtered_y_donors_list <- vector("list", nrow(x_recipients_df)); filtered_y_pred_donors_list <- vector("list", nrow(x_recipients_df))
    for (i in 1:nrow(x_recipients_df)) {
      current_recipient_cluster <- recipient_clusters[i]
      cluster_donors_idx <- which(donor_clusters == current_recipient_cluster)
      if (length(cluster_donors_idx) > 0) {
        filtered_x_donors_list[[i]] <- x_donors_df[cluster_donors_idx, , drop = FALSE]
        filtered_y_donors_list[[i]] <- y_donors[cluster_donors_idx]
        filtered_y_pred_donors_list[[i]] <- y_pred_donors[cluster_donors_idx]
      } else {
        filtered_x_donors_list[[i]] <- x_donors_df; filtered_y_donors_list[[i]] <- y_donors; filtered_y_pred_donors_list[[i]] <- y_pred_donors
      }
    }
  } else {
    filtered_x_donors_list <- rep(list(x_donors_df), nrow(x_recipients_df)); filtered_y_donors_list <- rep(list(y_donors), nrow(x_recipients_df)); filtered_y_pred_donors_list <- rep(list(y_pred_donors), nrow(x_recipients_df))
  }

  # --- 4. Enhanced Donor Identification and Imputation ---
  imputed_values <- vector("list", length = nrow(x_recipients_df))

  for (i in 1:nrow(x_recipients_df)) {
    current_recipient_x <- x_recipients_df[i, , drop = FALSE]
    current_recipient_y_pred <- y_pred_recipients[i]
    current_donors_x_df <- filtered_x_donors_list[[i]]
    current_donors_y <- filtered_y_donors_list[[i]]
    current_donors_y_pred <- filtered_y_pred_donors_list[[i]]

    if (nrow(current_donors_x_df) == 0) {
      warning("No donors available for a recipient. Using random sample from all original observed y.", call. = FALSE)
      imputed_val <- sample(y_donors, 1)
    } else {
      pred_diffs <- abs(current_donors_y_pred - current_recipient_y_pred)
      initial_pool_size <- min(max(donors * pmm_pool_factor, 10), nrow(current_donors_x_df))
      ordered_donors_idx_by_pred <- order(pred_diffs, decreasing = FALSE, na.last = TRUE)
      initial_donors_local_idx <- head(ordered_donors_idx_by_pred, initial_pool_size)
      x_initial_donors_pool <- current_donors_x_df[initial_donors_local_idx, , drop = FALSE]
      y_initial_donors_pool <- current_donors_y[initial_donors_local_idx]
      combined_data_for_gower <- rbind(current_recipient_x, x_initial_donors_pool)
      valid_cols_for_gower <- sapply(combined_data_for_gower, function(col) length(unique(stats::na.omit(col))) > 1)
      
      if (sum(valid_cols_for_gower) == 0) {
        warning("No valid columns for Gower's distance. Falling back to simple PMM.", call. = FALSE)
        selected_donor_local_idx <- sample(initial_donors_local_idx, 1)
        imputed_val <- current_donors_y[selected_donor_local_idx]
      } else {
        # *** FIX: Suppress benign warnings from daisy() ***
        gower_dist_matrix <- suppressWarnings(cluster::daisy(
            combined_data_for_gower[, valid_cols_for_gower, drop = FALSE], 
            metric = "gower", 
            stand = TRUE)
        )
        gower_distances <- as.matrix(gower_dist_matrix)[1, -1]
        gower_distances[is.na(gower_distances)] <- max(gower_distances, na.rm = TRUE) + 1
        num_final_donors <- min(donors, length(gower_distances))
        final_donors_local_idx_in_pool <- head(order(gower_distances, decreasing = FALSE, na.last = TRUE), num_final_donors)
        selected_donor_local_idx <- sample(final_donors_local_idx_in_pool, 1)
        imputed_val <- y_initial_donors_pool[selected_donor_local_idx]
      }
    }
    
    imputed_values[[i]] <- imputed_val
  }
  
  return(do.call(c, imputed_values))
}

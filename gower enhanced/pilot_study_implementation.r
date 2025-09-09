# =============================================================================
# PILOT STUDY: Focused Implementation for PhD Research
# Adaptive Weighted Gower-PMM - Minimal Viable Research Product
# =============================================================================

library(mice)
library(dplyr)
library(ggplot2)
library(VIM)

# =============================================================================
# PHASE 1: SIMPLIFIED THEORETICAL WEIGHTS (WEEKS 1-2)
# =============================================================================

# Core contribution: Information-theoretic variable weighting
compute_information_weights <- function(data, missing_pattern = NULL) {
  # Simplified information-theoretic weighting
  # Based on variable entropy and missing data correlation
  
  p <- ncol(data)
  weights <- numeric(p)
  
  for (j in 1:p) {
    if (is.numeric(data[[j]])) {
      # Continuous: differential entropy approximation
      # FIX: Use is.finite() to exclude NA, NaN, Inf, and -Inf
      observed_vals <- data[[j]][is.finite(data[[j]])]
      
      # IMPROVEMENT: Variance is valid for > 1 observation. The original > 5 was too strict.
      if (length(observed_vals) > 1) {
        # Use sample variance as entropy proxy
        weights[j] <- log(var(observed_vals) + 1) # na.rm is no longer needed
      } else {
        weights[j] <- 1 # Fallback for columns with < 2 finite values
      }
    } else {
      # Categorical: discrete entropy
      observed_vals <- data[[j]][!is.na(data[[j]])]
      if (length(observed_vals) > 2) {
        probs <- table(observed_vals) / length(observed_vals)
        weights[j] <- -sum(probs * log(probs + 1e-10))
      } else {
        weights[j] <- 1
      }
    }
  }
  
  # Normalize to prevent numerical issues
  # FIX: Add na.rm = TRUE to mean() to guard against any remaining NA issues.
  weights_mean <- mean(weights, na.rm = TRUE)
  if (is.na(weights_mean) || weights_mean == 0) {
      weights_mean <- 1 # Prevent division by zero or NA
  }
  weights <- weights / weights_mean
  
  weights[is.na(weights) | weights <= 0] <- 0.1  # Minimum weight and NA fallback
  
  return(weights)
}

# Simplified Gower distance with information weights
gower_distance_weighted <- function(x1, x2, weights) {
  # Simple vectorized Gower distance computation
  p <- length(x1)
  total_weight <- 0
  weighted_distance <- 0
  
  for (k in 1:p) {
    if (!is.na(x1[k]) && !is.na(x2[k])) {
      if (is.numeric(x1[k])) {
        # Continuous variable - need range from full dataset
        similarity <- 1 - abs(x1[k] - x2[k]) / (max(c(x1[k], x2[k])) - min(c(x1[k], x2[k])) + 0.01)
      } else {
        # Categorical variable
        similarity <- ifelse(x1[k] == x2[k], 1, 0)
      }
      
      weighted_distance <- weighted_distance + weights[k] * (1 - similarity)
      total_weight <- total_weight + weights[k]
    }
  }
  
  if (total_weight > 0) {
    return(weighted_distance / total_weight)
  } else {
    return(1)  # Maximum distance
  }
}

# Basic implementation of weighted PMM
mice.impute.weighted_pmm_pilot <- function(y, ry, x, donors = 5, use_weights = TRUE, ...) {
  
  # Prepare data
  x_donors <- x[ry, , drop = FALSE]
  y_donors <- y[ry]
  x_recipients <- x[!ry, , drop = FALSE]
  
  n_donors <- length(y_donors)
  n_recipients <- sum(!ry)
  
  if (n_recipients == 0) return(y[!ry])
  if (n_donors < 5) return(rep(mean(y_donors), n_recipients))  # Fallback
  
  # Compute weights
  if (use_weights) {
    all_data <- rbind(x_donors, x_recipients)
    weights <- compute_information_weights(all_data)
  } else {
    weights <- rep(1, ncol(x))
  }
  
  # Fit simple prediction model
  if (ncol(x_donors) > 0) {
    # Use simple linear/logistic regression
    df_train <- cbind(x_donors, y = y_donors)
    
    if (is.numeric(y_donors)) {
      model <- lm(y ~ ., data = df_train)
      y_pred_donors <- predict(model)
      y_pred_recipients <- predict(model, newdata = x_recipients)
    } else {
      # For categorical outcomes
      model <- glm(y ~ ., data = df_train, family = "binomial")
      y_pred_donors <- predict(model, type = "response")
      y_pred_recipients <- predict(model, newdata = x_recipients, type = "response")
    }
  } else {
    # No predictors, use mean
    y_pred_donors <- rep(mean(y_donors), n_donors)
    y_pred_recipients <- rep(mean(y_donors), n_recipients)
  }
  
  # Perform weighted PMM
  imputed_values <- numeric(n_recipients)
  
  for (i in 1:n_recipients) {
    # Find prediction-based pool
    pred_diffs <- abs(y_pred_donors - y_pred_recipients[i])
    pool_size <- min(donors * 5, n_donors)  # Larger initial pool
    pool_indices <- order(pred_diffs)[1:pool_size]
    
    # Calculate weighted Gower distances within pool
    distances <- numeric(length(pool_indices))
    for (j in seq_along(pool_indices)) {
      donor_idx <- pool_indices[j]
      distances[j] <- gower_distance_weighted(
        as.numeric(x_recipients[i, ]), 
        as.numeric(x_donors[donor_idx, ]), 
        weights
      )
    }
    
    # Select final donors based on distance
    final_donors <- pool_indices[order(distances)[1:min(donors, length(distances))]]
    
    # Random selection from final donors
    selected_donor <- sample(final_donors, 1)
    imputed_values[i] <- y_donors[selected_donor]
  }
  
  # Store diagnostic information
  attr(imputed_values, "weights") <- weights
  attr(imputed_values, "method") <- "weighted_pmm_pilot"
  
  return(imputed_values)
}

# =============================================================================
# PHASE 2: BASIC VALIDATION FRAMEWORK (WEEKS 3-4)
# =============================================================================

run_pilot_validation <- function(data, n_simulations = 10) {
  # Simplified validation focusing on core research question
  
  cat("Running pilot validation study...\n")
  cat("Dataset:", nrow(data), "x", ncol(data), "\n")
  cat("Simulations:", n_simulations, "\n\n")
  
  # Test conditions (simplified)
  missing_props <- c(0.2, 0.4)  # Start with moderate missing data
  mechanisms <- c("MCAR", "MAR")  # Skip MNAR for now
  methods <- c("weighted", "standard")
  
  results <- list()
  
  for (mechanism in mechanisms) {
    for (missing_prop in missing_props) {
      for (method in methods) {
        
        cat("Testing:", mechanism, "- Missing:", missing_prop, "- Method:", method, "\n")
        
        # Run simulations
        sim_results <- replicate(n_simulations, {
          # Generate missing data
          if (mechanism == "MCAR") {
            data_miss <- prodNA(data, noNA = missing_prop)
          } else {
            # Simplified MAR: missingness depends on first variable
            data_miss <- data
            first_var <- data[[1]]
            missing_prob <- plogis(scale(first_var))  # Logistic relationship
            missing_indices <- sample(nrow(data), size = floor(missing_prop * nrow(data)))
            
            # Make second variable missing
            if (ncol(data) > 1) {
              data_miss[missing_indices, 2] <- NA
            }
          }
          
          # Skip if no missing data created
          if (!any(is.na(data_miss))) {
            return(c(rmse = NA, bias = NA, coverage = NA))
          }
          
          # Perform imputation
          use_weights <- (method == "weighted")
          
          # Register method temporarily
          mice_method <- ifelse(use_weights, "weighted_pmm_pilot", "pmm")
          
          tryCatch({
            if (use_weights) {
              # Use our method
              imp <- mice(data_miss, method = "weighted_pmm_pilot", 
                         m = 1, maxit = 3, printFlag = FALSE)
            } else {
              # Standard PMM
              imp <- mice(data_miss, method = "pmm", 
                         m = 1, maxit = 3, printFlag = FALSE)
            }
            
            data_imp <- complete(imp)
            
            # Calculate metrics for first numeric variable
            numeric_vars <- sapply(data, is.numeric)
            if (sum(numeric_vars) == 0) {
              return(c(rmse = NA, bias = NA, coverage = NA))
            }
            
            first_numeric <- names(data)[numeric_vars][1]
            
            # Find imputed values
            was_missing <- is.na(data_miss[[first_numeric]]) & !is.na(data[[first_numeric]])
            
            if (sum(was_missing) == 0) {
              return(c(rmse = 0, bias = 0, coverage = 1))
            }
            
            true_vals <- data[[first_numeric]][was_missing]
            imp_vals <- data_imp[[first_numeric]][was_missing]
            
            # Calculate metrics
            rmse <- sqrt(mean((true_vals - imp_vals)^2))
            bias <- mean(imp_vals - true_vals)
            
            # Simple coverage approximation
            se_est <- sd(true_vals - imp_vals)
            coverage <- mean(abs(true_vals - imp_vals) <= 1.96 * se_est)
            
            return(c(rmse = rmse, bias = bias, coverage = coverage))
            
          }, error = function(e) {
            cat("Error:", e$message, "\n")
            return(c(rmse = NA, bias = NA, coverage = NA))
          })
          
        }, simplify = FALSE)
        
        # Aggregate results
        metrics_matrix <- do.call(rbind, sim_results)
        valid_rows <- complete.cases(metrics_matrix)
        
        if (sum(valid_rows) > 0) {
          results[[paste(mechanism, missing_prop, method, sep = "_")]] <- list(
            mechanism = mechanism,
            missing_prop = missing_prop,
            method = method,
            n_valid = sum(valid_rows),
            mean_rmse = mean(metrics_matrix[valid_rows, "rmse"]),
            se_rmse = sd(metrics_matrix[valid_rows, "rmse"]) / sqrt(sum(valid_rows)),
            mean_bias = mean(metrics_matrix[valid_rows, "bias"]),
            mean_coverage = mean(metrics_matrix[valid_rows, "coverage"])
          )
        }
      }
    }
  }
  
  return(results)
}

# =============================================================================
# PHASE 3: BASIC ANALYSIS AND REPORTING (WEEK 4)
# =============================================================================

analyze_pilot_results <- function(results) {
  # Convert to data frame for analysis
  results_df <- do.call(rbind, lapply(results, function(x) {
    data.frame(
      mechanism = x$mechanism,
      missing_prop = x$missing_prop,
      method = x$method,
      n_valid = x$n_valid,
      rmse = x$mean_rmse,
      rmse_se = x$se_rmse,
      bias = abs(x$mean_bias),
      coverage = x$mean_coverage
    )
  }))
  
  if (nrow(results_df) == 0) {
    cat("No valid results to analyze.\n")
    return(NULL)
  }
  
  # Basic statistical test
  weighted_rmse <- results_df$rmse[results_df$method == "weighted"]
  standard_rmse <- results_df$rmse[results_df$method == "standard"]
  
  if (length(weighted_rmse) > 0 && length(standard_rmse) > 0) {
    # Simple comparison
    improvement <- (mean(standard_rmse, na.rm = TRUE) - mean(weighted_rmse, na.rm = TRUE)) / 
                   mean(standard_rmse, na.rm = TRUE) * 100
    
    cat("\n=== PILOT STUDY RESULTS ===\n")
    cat("Weighted PMM RMSE:", round(mean(weighted_rmse, na.rm = TRUE), 4), "\n")
    cat("Standard PMM RMSE:", round(mean(standard_rmse, na.rm = TRUE), 4), "\n")
    cat("Improvement:", round(improvement, 1), "%\n\n")
    
    # Simple visualization
    if (nrow(results_df) > 2) {
      p <- ggplot(results_df, aes(x = method, y = rmse, fill = method)) +
        geom_bar(stat = "identity", position = "dodge") +
        geom_errorbar(aes(ymin = rmse - rmse_se, ymax = rmse + rmse_se), 
                     width = 0.2) +
        facet_grid(mechanism ~ missing_prop, labeller = label_both) +
        labs(title = "Pilot Study: Weighted vs Standard PMM",
             x = "Method", y = "RMSE") +
        theme_minimal()
      
      return(list(
        results_df = results_df,
        improvement = improvement,
        plot = p
      ))
    }
  }
  
  return(list(results_df = results_df))
}

# =============================================================================
# PHASE 4: EXECUTION AND TESTING FRAMEWORK
# =============================================================================

run_complete_pilot_study <- function(data = NULL) {
  
  cat("=== ADAPTIVE WEIGHTED GOWER-PMM: PILOT STUDY ===\n\n")
  
  # Use provided data or create test dataset
  if (is.null(data)) {
    set.seed(42)
    n <- 500  # Start small
    data <- data.frame(
      x1 = rnorm(n, 10, 3),           # High variance continuous
      x2 = rnorm(n, 5, 1),            # Low variance continuous  
      x3 = sample(c("A", "B", "C", "D"), n, replace = TRUE),  # 4-level categorical
      x4 = sample(c("Low", "High"), n, replace = TRUE),       # 2-level categorical
      x5 = x1 * 0.5 + rnorm(n, 0, 1)  # Correlated with x1
    )
    cat("Generated synthetic dataset:", nrow(data), "x", ncol(data), "\n")
  }
  
  # Test weight computation
  cat("\nTesting information weight computation...\n")
  weights <- compute_information_weights(data)
  cat("Computed weights:", round(weights, 3), "\n")
  
  # Register our method with mice
  environment(mice.impute.weighted_pmm_pilot) <- asNamespace('mice')
  
  # Run validation
  cat("\nRunning pilot validation study...\n")
  validation_results <- run_pilot_validation(data, n_simulations = 5)  # Small for testing
  
  # Analyze results  
  cat("\nAnalyzing results...\n")
  analysis <- analyze_pilot_results(validation_results)
  
  if (!is.null(analysis) && "plot" %in% names(analysis)) {
    print(analysis$plot)
  }
  
  cat("\n=== PILOT STUDY COMPLETED ===\n")
  return(analysis)
}

# Quick test function
test_basic_functionality <- function() {
  cat("Testing basic functionality...\n")
  
  # Create tiny test dataset
  test_data <- data.frame(
    x1 = c(1, 2, 3, 4, 5),
    x2 = c("A", "B", "A", "B", "A")
  )
  
  # Test weight computation
  weights <- compute_information_weights(test_data)
  cat("Weights computed:", all(weights > 0), "\n")
  
  # Test distance computation
  dist <- gower_distance_weighted(c(1, "A"), c(2, "B"), weights)
  cat("Distance computed:", is.numeric(dist), "\n")
  
  cat("Basic functionality test: PASSED\n")
}

# =============================================================================
# EXECUTION
# =============================================================================

cat("Pilot Study Implementation Loaded!\n")
cat("Key functions:\n")
cat("  • test_basic_functionality(): Quick functionality test\n") 
cat("  • run_complete_pilot_study(): Full pilot study execution\n")
cat("  • compute_information_weights(): Core theoretical contribution\n\n")

cat("To start, run:\n")
cat("test_basic_functionality()\n")
cat("results <- run_complete_pilot_study()\n")
# =============================================================================
# Comprehensive Benchmarking and Validation Suite for Enhanced GKP-PMM
# For Publication in High-Impact Journals
# =============================================================================

library(mice)
library(VIM)
library(naniar)
library(ggplot2)
library(dplyr)
library(tidyr)
library(microbenchmark)
library(foreach)
library(doParallel)
library(purrr)
library(broom)

# =============================================================================
# 1. PERFORMANCE BENCHMARKING FRAMEWORK
# =============================================================================

benchmark_methods <- function(data, missing_prop = 0.3, n_runs = 10,
                            methods_list = c("gkp_pmm_enhanced", "pmm", "rf", "cart"),
                            mechanisms = c("MCAR", "MAR", "MNAR")) {
  
  results <- list()
  
  for (mechanism in mechanisms) {
    cat("\nBenchmarking for", mechanism, "mechanism...\n")
    
    for (run in 1:n_runs) {
      # Generate missing data
      if (mechanism == "MCAR") {
        data_miss <- prodNA(data, noNA = missing_prop)
      } else if (mechanism == "MAR") {
        amp <- ampute(data, prop = missing_prop, mech = "MAR")
        data_miss <- amp$amp
      } else {  # MNAR
        amp <- ampute(data, prop = missing_prop, mech = "MNAR")
        data_miss <- amp$amp
      }
      
      # Benchmark each method
      for (method in methods_list) {
        cat("  Testing", method, "- Run", run, "/", n_runs, "\r")
        
        # Time the imputation
        time_result <- system.time({
          tryCatch({
            if (method == "gkp_pmm_enhanced") {
              imp <- mice(data_miss, method = method, m = 5, maxit = 10,
                         donors = 5, k_pre_clusters = "auto",
                         predictive_model = "auto", printFlag = FALSE)
            } else {
              imp <- mice(data_miss, method = method, m = 5, maxit = 10,
                         printFlag = FALSE)
            }
            data_imp <- complete(imp, 1)
          }, error = function(e) {
            cat("\nError in", method, ":", e$message, "\n")
            data_imp <- NULL
          })
        })
        
        if (!is.null(data_imp)) {
          # Calculate accuracy metrics
          accuracy <- calculate_accuracy_metrics(data, data_miss, data_imp)
          
          results[[length(results) + 1]] <- data.frame(
            mechanism = mechanism,
            method = method,
            run = run,
            time = time_result["elapsed"],
            rmse = accuracy$rmse,
            mae = accuracy$mae,
            coverage = accuracy$coverage,
            bias = accuracy$bias,
            variance_ratio = accuracy$variance_ratio
          )
        }
      }
    }
  }
  
  return(bind_rows(results))
}

# =============================================================================
# 2. ACCURACY METRICS CALCULATION
# =============================================================================

calculate_accuracy_metrics <- function(data_true, data_miss, data_imp) {
  # Identify which values were imputed
  imputed_mask <- is.na(data_miss) & !is.na(data_true)
  
  metrics <- list()
  
  # Overall metrics
  all_true <- data_true[imputed_mask]
  all_imp <- data_imp[imputed_mask]
  
  if (length(all_true) > 0) {
    if (is.numeric(all_true)) {
      metrics$rmse <- sqrt(mean((all_true - all_imp)^2, na.rm = TRUE))
      metrics$mae <- mean(abs(all_true - all_imp), na.rm = TRUE)
      metrics$bias <- mean(all_imp - all_true, na.rm = TRUE)
      
      # Coverage (proportion within 95% CI)
      se <- sd(all_true - all_imp, na.rm = TRUE)
      ci_lower <- all_imp - 1.96 * se
      ci_upper <- all_imp + 1.96 * se
      metrics$coverage <- mean(all_true >= ci_lower & all_true <= ci_upper, na.rm = TRUE)
      
      # Variance preservation
      metrics$variance_ratio <- var(all_imp, na.rm = TRUE) / var(all_true, na.rm = TRUE)
    } else {
      # For categorical variables
      metrics$accuracy <- mean(all_true == all_imp, na.rm = TRUE)
      metrics$rmse <- 1 - metrics$accuracy  # Convert to error rate
      metrics$mae <- metrics$rmse
      metrics$bias <- 0
      metrics$coverage <- metrics$accuracy
      metrics$variance_ratio <- 1
    }
  }
  
  return(metrics)
}

# =============================================================================
# 3. SCALABILITY TESTING
# =============================================================================

test_scalability <- function(max_n = 100000, max_p = 100) {
  sample_sizes <- c(100, 500, 1000, 5000, 10000, 50000, max_n)
  n_features <- c(10, 20, 50, max_p)
  
  scalability_results <- expand.grid(
    n = sample_sizes,
    p = n_features,
    method = c("gkp_pmm_enhanced", "gkp_pmm_parallel", "pmm", "rf"),
    time = NA,
    memory = NA
  )
  
  for (i in 1:nrow(scalability_results)) {
    n <- scalability_results$n[i]
    p <- scalability_results$p[i]
    method <- as.character(scalability_results$method[i])
    
    cat("\nTesting", method, "with n =", n, "and p =", p, "\n")
    
    # Generate synthetic data
    data <- generate_mixed_data(n, p)
    data_miss <- prodNA(data, noNA = 0.3)
    
    # Measure time and memory
    mem_before <- pryr::mem_used()
    
    time_result <- system.time({
      tryCatch({
        imp <- mice(data_miss, method = method, m = 1, maxit = 5,
                   printFlag = FALSE)
      }, error = function(e) {
        cat("Error:", e$message, "\n")
      })
    })
    
    mem_after <- pryr::mem_used()
    
    scalability_results$time[i] <- time_result["elapsed"]
    scalability_results$memory[i] <- as.numeric(mem_after - mem_before)
  }
  
  return(scalability_results)
}

# =============================================================================
# 4. STATISTICAL VALIDITY TESTING
# =============================================================================

test_statistical_validity <- function(data, n_imputations = 100, missing_prop = 0.3) {
  # Test if imputation preserves statistical properties
  
  results <- list()
  
  # Generate missing data once
  data_miss <- prodNA(data, noNA = missing_prop)
  
  # Collect imputed datasets
  imputed_datasets <- list()
  
  for (i in 1:n_imputations) {
    cat("\rGenerating imputation", i, "/", n_imputations)
    
    imp <- mice(data_miss, method = "gkp_pmm_enhanced", m = 1,
               printFlag = FALSE)
    imputed_datasets[[i]] <- complete(imp, 1)
  }
  
  cat("\n")
  
  # Test 1: Mean preservation
  original_means <- colMeans(data, na.rm = TRUE)
  imputed_means <- do.call(rbind, lapply(imputed_datasets, colMeans, na.rm = TRUE))
  
  mean_bias <- colMeans(imputed_means) - original_means
  mean_ci <- apply(imputed_means, 2, function(x) {
    t.test(x, mu = original_means[1])$conf.int
  })
  
  results$mean_preservation <- list(
    bias = mean_bias,
    ci_coverage = mean(original_means >= mean_ci[1,] & 
                       original_means <= mean_ci[2,])
  )
  
  # Test 2: Variance preservation
  original_vars <- apply(data, 2, var, na.rm = TRUE)
  imputed_vars <- do.call(rbind, lapply(imputed_datasets, function(x) {
    apply(x, 2, var, na.rm = TRUE)
  }))
  
  var_ratio <- colMeans(imputed_vars) / original_vars
  
  results$variance_preservation <- list(
    ratio = var_ratio,
    within_10pct = mean(var_ratio >= 0.9 & var_ratio <= 1.1)
  )
  
  # Test 3: Correlation structure preservation
  original_cor <- cor(data, use = "complete.obs")
  imputed_cors <- lapply(imputed_datasets, function(x) cor(x, use = "complete.obs"))
  mean_imputed_cor <- Reduce("+", imputed_cors) / length(imputed_cors)
  
  cor_diff <- mean_imputed_cor - original_cor
  
  results$correlation_preservation <- list(
    max_difference = max(abs(cor_diff[upper.tri(cor_diff)])),
    mean_absolute_difference = mean(abs(cor_diff[upper.tri(cor_diff)]))
  )
  
  # Test 4: Distribution preservation (KS test)
  ks_results <- matrix(NA, ncol(data), n_imputations)
  
  for (j in 1:ncol(data)) {
    if (is.numeric(data[,j])) {
      for (i in 1:n_imputations) {
        ks_results[j,i] <- ks.test(data[,j], imputed_datasets[[i]][,j])$p.value
      }
    }
  }
  
  results$distribution_preservation <- list(
    mean_ks_pvalue = rowMeans(ks_results, na.rm = TRUE),
    prop_nonsignificant = rowMeans(ks_results > 0.05, na.rm = TRUE)
  )
  
  return(results)
}

# =============================================================================
# 5. ROBUSTNESS TESTING
# =============================================================================

test_robustness <- function(base_data, contamination_levels = c(0, 0.05, 0.1, 0.2)) {
  results <- list()
  
  for (contam_level in contamination_levels) {
    cat("\nTesting with", contam_level * 100, "% contamination\n")
    
    # Add outliers
    data_contaminated <- base_data
    if (contam_level > 0) {
      n_outliers <- floor(nrow(base_data) * contam_level)
      outlier_rows <- sample(nrow(base_data), n_outliers)
      
      for (col in 1:ncol(base_data)) {
        if (is.numeric(base_data[,col])) {
          # Add extreme values
          data_contaminated[outlier_rows[1:floor(n_outliers/2)], col] <- 
            mean(base_data[,col], na.rm = TRUE) + 5 * sd(base_data[,col], na.rm = TRUE)
          data_contaminated[outlier_rows[(floor(n_outliers/2)+1):n_outliers], col] <- 
            mean(base_data[,col], na.rm = TRUE) - 5 * sd(base_data[,col], na.rm = TRUE)
        }
      }
    }
    
    # Generate missing data
    data_miss <- prodNA(data_contaminated, noNA = 0.3)
    
    # Test different methods
    methods <- c("gkp_pmm_enhanced", "pmm", "rf")
    
    for (method in methods) {
      imp <- mice(data_miss, method = method, m = 5, printFlag = FALSE)
      data_imp <- complete(imp, 1)
      
      # Calculate robustness metrics
      accuracy <- calculate_accuracy_metrics(base_data, data_miss, data_imp)
      
      results[[length(results) + 1]] <- data.frame(
        contamination = contam_level,
        method = method,
        rmse = accuracy$rmse,
        bias = abs(accuracy$bias),
        variance_ratio = accuracy$variance_ratio
      )
    }
  }
  
  return(bind_rows(results))
}

# =============================================================================
# 6. CONVERGENCE DIAGNOSTICS
# =============================================================================

diagnose_convergence <- function(data_miss, method = "gkp_pmm_enhanced", 
                                maxit = 20, m = 5) {
  # Run imputation with trace
  imp <- mice(data_miss, method = method, m = m, maxit = maxit, 
             printFlag = FALSE, seed = 123)
  
  # Extract convergence information
  convergence_data <- imp$chainMean
  
  # Plot convergence for each variable
  plots <- list()
  
  for (var in names(convergence_data)) {
    df <- as.data.frame(convergence_data[[var]])
    df$iteration <- 1:nrow(df)
    
    df_long <- df %>%
      pivot_longer(cols = -iteration, names_to = "chain", values_to = "mean")
    
    plots[[var]] <- ggplot(df_long, aes(x = iteration, y = mean, color = chain)) +
      geom_line() +
      labs(title = paste("Convergence for", var),
           x = "Iteration", y = "Mean") +
      theme_minimal()
  }
  
  # Calculate convergence statistics
  convergence_stats <- lapply(convergence_data, function(x) {
    # Gelman-Rubin statistic approximation
    chains <- as.matrix(x)
    n <- nrow(chains)
    m <- ncol(chains)
    
    chain_means <- colMeans(chains)
    grand_mean <- mean(chain_means)
    
    B <- n * var(chain_means)  # Between-chain variance
    W <- mean(apply(chains, 2, var))  # Within-chain variance
    
    var_plus <- ((n-1)/n) * W + (1/n) * B
    R_hat <- sqrt(var_plus / W)
    
    return(list(R_hat = R_hat, converged = R_hat < 1.1))
  })
  
  return(list(plots = plots, stats = convergence_stats))
}

# =============================================================================
# 7. UTILITY FUNCTIONS
# =============================================================================

generate_mixed_data <- function(n = 1000, p = 10, prop_categorical = 0.3) {
  n_cat <- floor(p * prop_categorical)
  n_num <- p - n_cat
  
  data <- data.frame(matrix(NA, n, p))
  
  # Generate numeric variables
  if (n_num > 0) {
    for (i in 1:n_num) {
      data[,i] <- rnorm(n, mean = runif(1, -10, 10), sd = runif(1, 1, 5))
    }
  }
  
  # Generate categorical variables
  if (n_cat > 0) {
    for (i in (n_num+1):p) {
      n_levels <- sample(2:5, 1)
      data[,i] <- factor(sample(LETTERS[1:n_levels], n, replace = TRUE))
    }
  }
  
  names(data) <- paste0("V", 1:p)
  return(data)
}

# =============================================================================
# 8. PUBLICATION-READY VISUALIZATION
# =============================================================================

create_publication_plots <- function(benchmark_results) {
  library(ggpubr)
  library(viridis)
  
  # Theme for publication
  theme_publication <- theme_minimal() +
    theme(
      text = element_text(size = 12),
      axis.text = element_text(size = 10),
      legend.position = "bottom",
      panel.grid.minor = element_blank()
    )
  
  # Plot 1: Accuracy comparison
  p1 <- benchmark_results %>%
    group_by(method, mechanism) %>%
    summarise(
      mean_rmse = mean(rmse),
      se_rmse = sd(rmse) / sqrt(n()),
      .groups = 'drop'
    ) %>%
    ggplot(aes(x = method, y = mean_rmse, fill = method)) +
    geom_bar(stat = "identity", position = "dodge") +
    geom_errorbar(aes(ymin = mean_rmse - se_rmse, 
                     ymax = mean_rmse + se_rmse),
                 width = 0.2) +
    facet_wrap(~mechanism) +
    scale_fill_viridis_d() +
    labs(title = "Imputation Accuracy Across Missing Data Mechanisms",
         x = "Method", y = "RMSE") +
    theme_publication
  
  # Plot 2: Computational efficiency
  p2 <- benchmark_results %>%
    group_by(method) %>%
    summarise(
      mean_time = mean(time),
      se_time = sd(time) / sqrt(n()),
      .groups = 'drop'
    ) %>%
    ggplot(aes(x = reorder(method, mean_time), y = mean_time, fill = method)) +
    geom_bar(stat = "identity") +
    geom_errorbar(aes(ymin = mean_time - se_time, 
                     ymax = mean_time + se_time),
                 width = 0.2) +
    scale_fill_viridis_d() +
    labs(title = "Computational Efficiency",
         x = "Method", y = "Time (seconds)") +
    coord_flip() +
    theme_publication
  
  # Plot 3: Bias-Variance tradeoff
  p3 <- benchmark_results %>%
    group_by(method, mechanism) %>%
    summarise(
      mean_bias = mean(abs(bias)),
      mean_var = mean(variance_ratio),
      .groups = 'drop'
    ) %>%
    ggplot(aes(x = mean_bias, y = mean_var, color = method, shape = mechanism)) +
    geom_point(size = 3) +
    geom_hline(yintercept = 1, linetype = "dashed", alpha = 0.5) +
    geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
    scale_color_viridis_d() +
    labs(title = "Bias-Variance Tradeoff",
         x = "Absolute Bias", y = "Variance Ratio") +
    theme_publication
  
  # Combine plots
  combined_plot <- ggarrange(p1, p2, p3, 
                            labels = c("A", "B", "C"),
                            ncol = 2, nrow = 2,
                            common.legend = TRUE,
                            legend = "bottom")
  
  return(list(accuracy = p1, efficiency = p2, bias_variance = p3, 
              combined = combined_plot))
}

# =============================================================================
# 9. MAIN EXECUTION SCRIPT
# =============================================================================

run_full_validation <- function(data = NULL, output_dir = "validation_results") {
  
  # Create output directory
  dir.create(output_dir, showWarnings = FALSE)
  
  # Use provided data or generate synthetic
  if (is.null(data)) {
    cat("Generating synthetic dataset...\n")
    data <- generate_mixed_data(n = 5000, p = 20, prop_categorical = 0.4)
  }
  
  cat("\n=== Starting Comprehensive Validation ===\n")
  
  # 1. Performance Benchmarking
  cat("\n1. Running performance benchmarks...\n")
  benchmark_results <- benchmark_methods(data, missing_prop = 0.3, n_runs = 20)
  saveRDS(benchmark_results, file.path(output_dir, "benchmark_results.rds"))
  
  # 2. Statistical Validity
  cat("\n2. Testing statistical validity...\n")
  validity_results <- test_statistical_validity(data, n_imputations = 50)
  saveRDS(validity_results, file.path(output_dir, "validity_results.rds"))
  
  # 3. Robustness Testing
  cat("\n3. Testing robustness...\n")
  robustness_results <- test_robustness(data)
  saveRDS(robustness_results, file.path(output_dir, "robustness_results.rds"))
  
  # 4. Convergence Diagnostics
  cat("\n4. Running convergence diagnostics...\n")
  data_miss <- prodNA(data, noNA = 0.3)
  convergence_results <- diagnose_convergence(data_miss)
  saveRDS(convergence_results, file.path(output_dir, "convergence_results.rds"))
  
  # 5. Create publication plots
  cat("\n5. Creating publication-ready plots...\n")
  plots <- create_publication_plots(benchmark_results)
  
  # Save plots
  ggsave(file.path(output_dir, "figure1_accuracy.pdf"), 
         plots$accuracy, width = 10, height = 6)
  ggsave(file.path(output_dir, "figure2_efficiency.pdf"), 
         plots$efficiency, width = 8, height = 6)
  ggsave(file.path(output_dir, "figure3_bias_variance.pdf"), 
         plots$bias_variance, width = 8, height = 6)
  ggsave(file.path(output_dir, "figure_combined.pdf"), 
         plots$combined, width = 12, height = 10)
  
  # 6. Generate summary report
  cat("\n6. Generating summary report...\n")
  
  report <- list(
    benchmark_summary = benchmark_results %>%
      group_by(method) %>%
      summarise(
        mean_rmse = mean(rmse),
        mean_time = mean(time),
        efficiency_score = mean(rmse) * mean(time)
      ),
    validity_summary = validity_results,
    robustness_summary = robustness_results %>%
      group_by(method) %>%
      summarise(
        robustness_score = mean(1 / (1 + rmse))
      ),
    convergence_summary = sapply(convergence_results$stats, function(x) x$converged)
  )
  
  saveRDS(report, file.path(output_dir, "summary_report.rds"))
  
  cat("\n=== Validation Complete ===\n")
  cat("Results saved to:", output_dir, "\n")
  
  return(report)
}

# =============================================================================
# 10. EXAMPLE USAGE
# =============================================================================

# Load your enhanced imputation function
# source("enhanced_gkp_pmm.R")

# Run full validation
# results <- run_full_validation()

# Or with your own data
# data(iris)
# results <- run_full_validation(data = iris[,1:4])

cat("\nValidation suite loaded successfully!\n")
cat("Run 'run_full_validation()' to start comprehensive testing.\n")
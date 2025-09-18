#' Comprehensive Benchmark Framework for Gower-PMM Evaluation
#'
#' This module provides a complete benchmarking system to evaluate the Gower-PMM
#' imputation method against other distance measures and imputation techniques.
#' Designed for thesis research and methodological comparison.
#'
#' @author Yaji Timothy T.
#' @references 
#' Gower, J. C. (1971). A general coefficient of similarity and some of its properties. Biometrics, 857-874.
#' Van Buuren, S., & Groothuis-Oudshoorn, K. (2011). mice: Multivariate imputation by chained equations in R. Journal of statistical software, 45(3), 1-67.

# --- Benchmark Configuration Class ---

#' S4 Class for Benchmark Configuration
#' 
#' @slot scenarios List of simulation scenarios
#' @slot methods List of methods to compare
#' @slot metrics List of evaluation metrics
#' @slot n_replications Number of Monte Carlo replications
#' @slot missing_patterns List of missing data patterns
#' @slot output_dir Directory for saving results
#' @export
setClass("BenchmarkConfig",
         representation(
           scenarios = "list",
           methods = "list", 
           metrics = "list",
           n_replications = "numeric",
           missing_patterns = "list",
           output_dir = "character"
         ),
         prototype(
           scenarios = list(),
           methods = list(),
           metrics = list(),
           n_replications = 100,
           missing_patterns = list(),
           output_dir = "benchmark_results"
         ))

#' S4 Class for Benchmark Results
#' 
#' @slot config The benchmark configuration used
#' @slot results Data frame containing all benchmark results
#' @slot summary_stats Summary statistics by method and scenario
#' @slot execution_time Total execution time
#' @slot timestamp When the benchmark was run
#' @export
setClass("BenchmarkResults",
         representation(
           config = "BenchmarkConfig",
           results = "data.frame",
           summary_stats = "list",
           execution_time = "numeric",
           timestamp = "POSIXct"
         ))

# --- Default Benchmark Configuration ---

#' Create Default Benchmark Configuration
#' 
#' Sets up a comprehensive benchmark comparing Gower-PMM against standard methods
#' 
#' @param n_replications Number of Monte Carlo replications (default: 100)
#' @param output_dir Directory for results (default: "benchmark_results")
#' @return BenchmarkConfig object
#' @export
create_default_benchmark_config <- function(n_replications = 100, output_dir = "benchmark_results") {
  
  # Define simulation scenarios
  scenarios <- list(
    # Small datasets
    small_balanced = list(n = 100, p_num = 3, p_cat = 2, p_ord = 1, 
                         correlation = "moderate", noise_level = "low"),
    small_mixed = list(n = 100, p_num = 2, p_cat = 3, p_ord = 2,
                      correlation = "high", noise_level = "medium"),
    
    # Medium datasets  
    medium_balanced = list(n = 500, p_num = 5, p_cat = 3, p_ord = 2,
                          correlation = "moderate", noise_level = "low"),
    medium_complex = list(n = 500, p_num = 4, p_cat = 4, p_ord = 3,
                         correlation = "mixed", noise_level = "high"),
    
    # Large datasets
    large_simple = list(n = 1000, p_num = 6, p_cat = 2, p_ord = 1,
                       correlation = "low", noise_level = "medium"),
    large_complex = list(n = 1000, p_num = 5, p_cat = 5, p_ord = 4,
                        correlation = "high", noise_level = "high"),
    
    # High-dimensional scenarios
    wide_data = list(n = 200, p_num = 10, p_cat = 5, p_ord = 3,
                    correlation = "moderate", noise_level = "medium"),
    
    # Extreme scenarios
    highly_correlated = list(n = 300, p_num = 4, p_cat = 2, p_ord = 1,
                           correlation = "very_high", noise_level = "low"),
    noisy_data = list(n = 300, p_num = 3, p_cat = 3, p_ord = 2,
                     correlation = "moderate", noise_level = "very_high")
  )
  
  # Define comparison methods
  methods <- list(
    # Gower-PMM variants
    gowerpmm_auto = list(
      name = "Gower-PMM (Auto Weights)",
      type = "gowerpmm",
      params = list(weights = "auto", k = 5)
    ),
    gowerpmm_equal = list(
      name = "Gower-PMM (Equal Weights)", 
      type = "gowerpmm",
      params = list(weights = "equal", k = 5)
    ),
    gowerpmm_range = list(
      name = "Gower-PMM (Range Scaling)",
      type = "gowerpmm", 
      params = list(weights = "auto", k = 5, scaling = "range")
    ),
    gowerpmm_iqr = list(
      name = "Gower-PMM (IQR Scaling)",
      type = "gowerpmm",
      params = list(weights = "auto", k = 5, scaling = "iqr")
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
      name = "MICE Random Forest",
      type = "mice_standard",
      params = list(method = "rf")
    ),
    
    # Distance-based comparisons
    fd_gowdis = list(
      name = "FD::gowdis + PMM",
      type = "fd_gowdis",
      params = list(k = 5)
    ),
    
    # VIM package methods
    vim_knn = list(
      name = "VIM k-NN",
      type = "vim",
      params = list(method = "kNN", k = 5)
    ),
    vim_irmi = list(
      name = "VIM IRMI", 
      type = "vim",
      params = list(method = "irmi")
    )
  )
  
  # Define evaluation metrics
  metrics <- list(
    # Distance measure evaluation
    distance_metrics = list(
      "rank_correlation_balance",
      "distance_preservation", 
      "computational_efficiency"
    ),
    
    # Imputation quality metrics
    imputation_metrics = list(
      "rmse",           # Root Mean Square Error
      "mae",            # Mean Absolute Error  
      "bias",           # Bias
      "coverage",       # Confidence interval coverage
      "pfc",            # Proportion of Falsely Classified
      "kolmogorov_smirnov", # Distribution preservation
      "correlation_preservation", # Correlation structure
      "variance_preservation"     # Variance preservation
    ),
    
    # Computational metrics
    computational_metrics = list(
      "execution_time",
      "memory_usage",
      "convergence_rate"
    )
  )
  
  # Define missing data patterns
  missing_patterns <- list(
    mcar_10 = list(type = "MCAR", rate = 0.10, description = "10% MCAR"),
    mcar_20 = list(type = "MCAR", rate = 0.20, description = "20% MCAR"),
    mcar_30 = list(type = "MCAR", rate = 0.30, description = "30% MCAR"),
    
    mar_10 = list(type = "MAR", rate = 0.10, description = "10% MAR"),
    mar_20 = list(type = "MAR", rate = 0.20, description = "20% MAR"),
    mar_30 = list(type = "MAR", rate = 0.30, description = "30% MAR"),
    
    mnar_10 = list(type = "MNAR", rate = 0.10, description = "10% MNAR"),
    mnar_20 = list(type = "MNAR", rate = 0.20, description = "20% MNAR"),
    
    # Complex patterns
    monotone = list(type = "MONOTONE", rate = 0.20, description = "20% Monotone"),
    block = list(type = "BLOCK", rate = 0.15, description = "15% Block pattern")
  )
  
  new("BenchmarkConfig",
      scenarios = scenarios,
      methods = methods,
      metrics = metrics, 
      n_replications = n_replications,
      missing_patterns = missing_patterns,
      output_dir = output_dir)
}

# --- Simulation Data Generation ---

#' Generate Simulation Dataset
#' 
#' Creates synthetic mixed-type datasets with specified characteristics
#' 
#' @param scenario List containing dataset parameters
#' @param seed Random seed for reproducibility
#' @return data.frame with mixed-type variables
#' @export
generate_simulation_data <- function(scenario, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  n <- scenario$n
  p_num <- scenario$p_num
  p_cat <- scenario$p_cat  
  p_ord <- scenario$p_ord
  correlation <- scenario$correlation
  noise_level <- scenario$noise_level
  
  # Generate correlation matrix based on correlation level
  p_total <- p_num + p_cat + p_ord
  if (correlation == "very_high") {
    base_corr <- 0.8
  } else if (correlation == "high") {
    base_corr <- 0.6
  } else if (correlation == "moderate") {
    base_corr <- 0.4
  } else if (correlation == "low") {
    base_corr <- 0.2
  } else if (correlation == "mixed") {
    # Mixed correlation structure
    base_corr <- NULL
  } else {
    base_corr <- 0.3
  }
  
  # Generate numeric variables
  if (p_num > 0) {
    if (is.null(base_corr)) {
      # Mixed correlation structure
      num_data <- generate_mixed_correlation_data(n, p_num)
    } else {
      # Uniform correlation structure
      sigma <- matrix(base_corr, p_num, p_num)
      diag(sigma) <- 1
      num_data <- MASS::mvrnorm(n, mu = rep(0, p_num), Sigma = sigma)
    }
    colnames(num_data) <- paste0("num_", 1:p_num)
  } else {
    num_data <- NULL
  }
  
  # Generate categorical variables
  if (p_cat > 0) {
    cat_data <- matrix(NA, n, p_cat)
    for (i in 1:p_cat) {
      n_levels <- sample(3:5, 1)  # 3-5 levels per categorical variable
      if (!is.null(num_data) && correlation != "low") {
        # Make categorical variables correlated with numeric ones
        probs <- plogis(num_data[,1] * 0.5)  # Use first numeric variable
        cat_data[,i] <- sample(1:n_levels, n, replace = TRUE, 
                              prob = c(probs[1], rep((1-probs[1])/(n_levels-1), n_levels-1)))
      } else {
        cat_data[,i] <- sample(1:n_levels, n, replace = TRUE)
      }
      cat_data[,i] <- factor(cat_data[,i], labels = paste0("Cat", i, "_", LETTERS[1:n_levels]))
    }
    colnames(cat_data) <- paste0("cat_", 1:p_cat)
  } else {
    cat_data <- NULL
  }
  
  # Generate ordinal variables  
  if (p_ord > 0) {
    ord_data <- matrix(NA, n, p_ord)
    for (i in 1:p_ord) {
      n_levels <- sample(4:7, 1)  # 4-7 levels per ordinal variable
      if (!is.null(num_data) && correlation != "low") {
        # Make ordinal variables correlated with numeric ones
        linear_pred <- num_data[,min(i, p_num)] * 0.8
        ord_data[,i] <- cut(linear_pred, breaks = n_levels, labels = FALSE, include.lowest = TRUE)
      } else {
        ord_data[,i] <- sample(1:n_levels, n, replace = TRUE)
      }
      ord_data[,i] <- ordered(ord_data[,i], levels = 1:n_levels, 
                             labels = paste0("Ord", i, "_", 1:n_levels))
    }
    colnames(ord_data) <- paste0("ord_", 1:p_ord)
  } else {
    ord_data <- NULL
  }
  
  # Combine all variables
  sim_data <- cbind(num_data, cat_data, ord_data)
  
  # Add noise based on noise level
  if (noise_level == "very_high") {
    noise_factor <- 2.0
  } else if (noise_level == "high") {
    noise_factor <- 1.5
  } else if (noise_level == "medium") {
    noise_factor <- 1.0
  } else if (noise_level == "low") {
    noise_factor <- 0.5
  } else {
    noise_factor <- 1.0
  }
  
  # Add noise to numeric variables
  if (p_num > 0) {
    for (i in 1:p_num) {
      sim_data[,i] <- sim_data[,i] + rnorm(n, 0, noise_factor * 0.1)
    }
  }
  
  as.data.frame(sim_data)
}

#' Generate Mixed Correlation Structure
#' @keywords internal
generate_mixed_correlation_data <- function(n, p) {
  # Create blocks of correlated variables
  if (p <= 2) {
    sigma <- diag(p)
  } else {
    sigma <- diag(p)
    # Create correlation blocks
    block_size <- max(2, floor(p/2))
    for (i in 1:(p-1)) {
      for (j in (i+1):p) {
        if (abs(i - j) <= block_size) {
          sigma[i,j] <- sigma[j,i] <- 0.6 * exp(-0.1 * abs(i-j))
        } else {
          sigma[i,j] <- sigma[j,i] <- 0.1
        }
      }
    }
  }
  MASS::mvrnorm(n, mu = rep(0, p), Sigma = sigma)
}

# --- Missing Data Pattern Generation ---

#' Introduce Missing Data Patterns
#' 
#' @param data Complete dataset
#' @param pattern Missing data pattern specification
#' @param seed Random seed
#' @return Dataset with missing values
#' @export
introduce_missing_data <- function(data, pattern, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  n <- nrow(data)
  p <- ncol(data)
  missing_rate <- pattern$rate
  
  if (pattern$type == "MCAR") {
    # Missing Completely At Random
    missing_indices <- sample(1:(n*p), size = floor(n*p*missing_rate))
    data_with_missing <- data
    data_with_missing[missing_indices] <- NA
    
  } else if (pattern$type == "MAR") {
    # Missing At Random - missingness depends on observed variables
    data_with_missing <- data
    # Use first variable to determine missingness in others
    if (is.numeric(data[,1])) {
      prob_missing <- plogis((data[,1] - mean(data[,1], na.rm = TRUE)) / sd(data[,1], na.rm = TRUE))
    } else {
      prob_missing <- runif(n, 0.1, 0.3)
    }
    
    for (j in 2:p) {
      n_missing <- floor(n * missing_rate / (p-1))
      missing_idx <- sample(1:n, size = n_missing, prob = prob_missing)
      data_with_missing[missing_idx, j] <- NA
    }
    
  } else if (pattern$type == "MNAR") {
    # Missing Not At Random - missingness depends on unobserved values
    data_with_missing <- data
    for (j in 1:p) {
      if (is.numeric(data[,j])) {
        # Higher values more likely to be missing
        prob_missing <- plogis((data[,j] - quantile(data[,j], 0.7)) / sd(data[,j]))
        n_missing <- floor(n * missing_rate / p)
        missing_idx <- sample(1:n, size = n_missing, prob = prob_missing)
        data_with_missing[missing_idx, j] <- NA
      }
    }
    
  } else if (pattern$type == "MONOTONE") {
    # Monotone missing pattern
    data_with_missing <- data
    for (j in 2:p) {
      n_missing <- floor(n * missing_rate * j / p)
      missing_idx <- sample(1:n, size = n_missing)
      data_with_missing[missing_idx, j:p] <- NA
    }
    
  } else if (pattern$type == "BLOCK") {
    # Block missing pattern
    data_with_missing <- data
    block_size <- max(1, floor(p/3))
    start_col <- sample(1:(p-block_size+1), 1)
    end_col <- start_col + block_size - 1
    n_missing <- floor(n * missing_rate)
    missing_idx <- sample(1:n, size = n_missing)
    data_with_missing[missing_idx, start_col:end_col] <- NA
    
  } else {
    stop("Unknown missing data pattern type: ", pattern$type)
  }
  
  data_with_missing
}
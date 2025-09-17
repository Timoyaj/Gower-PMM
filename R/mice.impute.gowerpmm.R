#' Optimized Gower-PMM Imputation for `mice`
#'
#' Imputes missing data using an optimized Gower's distance combined with
#' Predictive Mean Matching (PMM). This function is designed to be called
#' from the main `mice()` function.
#'
#' @param y A numeric vector with missing values to be imputed.
#' @param ry A logical vector indicating which values in `y` are observed (`TRUE`)
#'           and which are missing (`FALSE`).
#' @param x A `data.frame` or `matrix` of predictor variables.
#' @param k The number of nearest neighbors to use as donors for PMM.
#' @param weights A character string specifying the weighting method.
#'                `"auto"` (default) uses the optimized weights.
#'                `"equal"` uses equal weights.
#'                Can also be a user-provided numeric vector of weights.
#' @param control A list of control parameters from `gowerpmmControl()` for the
#'                optimization algorithm.
#' @param ... Additional arguments (not used).
#'
#' @return A numeric vector of imputed values.
#' @export
mice.impute.gowerpmm <- function(y, ry, x, k = 5, weights = "auto", control = gowerpmmControl(), ...) {
  
  # --- 1. Caching Mechanism for Optimal Weights ---
  # This is a critical performance optimization to avoid re-calculating weights
  # on every iteration and for every variable within a single mice run.
  
  # The `mice` algorithm operates in an environment where we can store cached results.
  # We navigate up two frames to find the main `mice` function's environment.
  mice_env <- parent.frame(n = 2)
  
  # Initialize the cache if it doesn't exist
  if (!exists("gower_cache", envir = mice_env)) {
    assign("gower_cache", new.env(parent = emptyenv()), envir = mice_env)
  }
  cache <- get("gower_cache", envir = mice_env)
  
  # --- 2. Determine Weights (with Caching) ---
  w <- NULL
  if (is.character(weights) && weights == "auto") {
    # Create a unique key for the current set of predictors and control settings
    cache_key <- digest::digest(list(colnames(x), control))
    
    if (exists(cache_key, envir = cache)) {
      # Cache Hit: Use the stored weights
      w <- get(cache_key, envir = cache)
    } else {
      # Cache Miss: Run the expensive optimization
      message("Optimizing Gower weights for predictor set (this happens once per set)...")
      opt_results <- calculate_optimal_gower_weights(as.data.frame(x), control = control)
      w <- opt_results@weights
      
      # Store the result in the cache for subsequent calls
      assign(cache_key, w, envir = cache)
    }
  } else if (is.numeric(weights)) {
    # Use user-provided weights (after validation)
    if (length(weights) != ncol(x)) stop("Length of weights must equal the number of predictors.")
    w <- weights / sum(weights) # Normalize
  } else {
    # Default to equal weights
    w <- rep(1 / ncol(x), ncol(x))
  }
  
  # --- 3. Perform PMM with Weighted Gower Distance ---
  x_obs <- as.data.frame(x[ry, , drop = FALSE])
  y_obs <- y[ry]
  x_mis <- as.data.frame(x[!ry, , drop = FALSE])
  
  if (nrow(x_mis) == 0) {
    return(numeric(0))
  }
  
  # Calculate the weighted Gower distance from each missing case to all observed cases
  gower_mat <- gower_dist_engine(data1 = x_mis, data2 = x_obs, weights = w, scaling = control$scaling %||% "range")
  
  # For each missing case, find the k nearest neighbors
  donors <- apply(gower_mat, 1, function(row_dists) {
    donor_indices <- order(row_dists)[1:k]
    # Randomly sample one of the k-nearest donors
    sample(donor_indices, 1)
  })
  
  # Return the observed y-values from the selected donors
  return(y_obs[donors])
}

# Helper for NULL default in lists
`%||%` <- function(a, b) {
  if (is.null(a)) b else a
}


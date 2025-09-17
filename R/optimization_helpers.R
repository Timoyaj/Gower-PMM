# This script contains the helper functions for the Gower weight optimization,
# including the main wrapper, the fitness function, control parameter setup,
# and the S4 class definition for results, following the best practices of the GA package.


#' @useDynLib gowerpmm, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL

# --- S4 Class Definition for Optimization Results ---

#' An S4 class to store the results of the Gower weight optimization.
#'
#' @slot weights The final optimized weight vector.
#' @slot objective_value The final value of the objective function.
#' @slot rank_correlations A vector of the final rank correlations.
#' @slot ga_summary A summary of the Genetic Algorithm run.
#' @export
setClass("optimal_gower_weights",
          representation(
            weights = "numeric",
            objective_value = "numeric",
            rank_correlations = "numeric",
            ga_summary = "ANY"
          ),
          prototype(
            weights = numeric(0),
            objective_value = numeric(0),
            rank_correlations = numeric(0),
            ga_summary = list()
          ))

# --- Control Function ---

#' Control parameters for the gowerpmm imputation method
#'
#' @param k Number of donors for PMM.
#' @param scaling Scaling method for numeric variables ("range" or "iqr").
#' @param ga_params A list of parameters to pass to the `GA::ga` function.
#' @param monitor A function to monitor the GA's progress.
#' @return A list of control parameters.
#' @export
gowerpmmControl <- function(k = 5, scaling = "range", ga_params = list(), monitor = NULL) {
  # Default GA parameters
  default_ga <- list(
    popSize = 50,
    maxiter = 100,
    run = 20,
    pmonitor = NULL,
    optim = TRUE # Use hybrid optimization
  )
  # Merge user-provided params with defaults
  ga_params <- utils::modifyList(default_ga, ga_params)

  list(k = k, scaling = scaling, ga_params = ga_params, monitor = monitor)
}


# --- Main Optimization Function ---

#' Calculate Optimal Gower Weights
#'
#' @param x The predictor `data.frame`.
#' @param control A list of control parameters from `gowerpmmControl`.
#' @return An object of class `optimal_gower_weights`.
#' @keywords internal
calculate_optimal_gower_weights <- function(x, control = gowerpmmControl()) {
  # --- 1. Pre-computation (in C++) ---
  var_types <- .get_var_types(x)
  dissim_matrices <- compute_dissimilarities_cpp(x, x, var_types, control$scaling %||% "range")

  # --- 2. Parallel Backend Setup ---
  n_cores <- parallel::detectCores()
  n_cores <- max(1, n_cores - 1)
  # Disable parallel for testing to avoid cluster issues
  n_cores <- 1
  if (n_cores > 1) {
    cl <- parallel::makeCluster(n_cores)
    doParallel::registerDoParallel(cl)
    on.exit(parallel::stopCluster(cl))
    # Load the package and export the C++ function to the cluster
    parallel::clusterEvalQ(cl, library(gowerpmm))
    parallel::clusterExport(cl, "objective_function_cpp", envir = environment())
  }
  
  # --- 3. Optimization using GA ---
  ga_results <- GA::ga(
    type = "real-valued",
    fitness = function(w) -objective_function_cpp(w, dissim_matrices),
    lower = rep(0, ncol(x)),
    upper = rep(1, ncol(x)),
    popSize = control$ga_params$popSize,
    maxiter = control$ga_params$maxiter,
    run = control$ga_params$run,
    monitor = if (is.null(control$monitor)) GA::gaMonitor else control$monitor,
    parallel = if(n_cores > 1) cl else FALSE,
    optim = control$ga_params$optim
  )

  # --- 4. Post-processing ---
  optimal_weights <- abs(ga_results@solution[1, ])
  optimal_weights <- optimal_weights / sum(optimal_weights)
  names(optimal_weights) <- colnames(x)

  # Recalculate final correlations and objective value for reporting
  final_obj <- objective_function_cpp(optimal_weights, dissim_matrices)
  final_corrs <- get_rank_correlations_cpp(optimal_weights, dissim_matrices)
  names(final_corrs) <- colnames(x)
  
  # Create S4 results object
  new("optimal_gower_weights",
      weights = optimal_weights,
      objective_value = -final_obj,
      rank_correlations = final_corrs,
      ga_summary = summary(ga_results))
}

# --- Utility Functions ---

#' Get variable types as integers for C++
#' 0: numeric, 1: factor, 2: ordered
#' @keywords internal
.get_var_types <- function(data) {
  sapply(data, function(col) {
    if (is.ordered(col)) return(2)
    if (is.factor(col)) return(1)
    if (is.numeric(col)) return(0)
    return(1) # Default to factor for other types like character
  })
}

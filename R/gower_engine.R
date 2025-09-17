#' @useDynLib gowerpmm, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL

#' Enhanced Gower Distance Engine
#'
#' Computes a Gower dissimilarity matrix for mixed-type data using a high-performance C++ backend.
#' This function serves as the core dissimilarity engine for the package.
#'
#' @param data1 A data.frame with mixed-type data (numeric, factor, ordered).
#' @param data2 An optional second data.frame with the same structure as data1.
#'        If provided, computes dissimilarities between rows of data1 and data2.
#'        If NULL, computes dissimilarities within data1.
#' @param weights A numeric vector of weights for each variable. If NULL (default), equal weights are used.
#' @param scaling A string specifying the scaling method for numeric variables.
#'        Can be "range" (default) or "iqr" for robust scaling.
#'
#' @return If data2 is NULL, a standard R `dist` object representing the Gower dissimilarity matrix.
#'         If data2 is provided, a matrix of dissimilarities between data1 rows and data2 rows.
#'
#' @export
gower_dist_engine <- function(data1, data2 = NULL, weights = NULL, scaling = "range") {
  # Constants for numerical stability
  TOLERANCE <- 1e-6
  MIN_SCALE <- 1e-8
  MEMORY_WARNING_MB <- 1000

  # --- Input Validation ---
  if (!is.data.frame(data1)) {
    stop("Input 'data1' must be a data.frame.", call. = FALSE)
  }
  if (nrow(data1) == 0 || ncol(data1) == 0) {
    stop("Input 'data1' must not be empty.", call. = FALSE)
  }
  if (!is.null(data2)) {
    if (!is.data.frame(data2) || ncol(data2) != ncol(data1)) {
      stop("Input 'data2' must be a data.frame with the same number of columns as 'data1'.", call. = FALSE)
    }
    if (nrow(data2) == 0) {
      stop("Input 'data2' must not be empty.", call. = FALSE)
    }
  }
  if (!scaling %in% c("range", "iqr")) {
    stop("'scaling' must be either 'range' or 'iqr'.", call. = FALSE)
  }

  # --- Memory Check for Large Datasets ---
  n1 <- nrow(data1)
  n2 <- if (is.null(data2)) n1 else nrow(data2)
  p <- ncol(data1)
  estimated_memory_mb <- (n1 * n2 * p * 8) / (1024^2)  # Rough estimate for double precision
  if (estimated_memory_mb > MEMORY_WARNING_MB) {
    warning(sprintf("Estimated memory usage: %.1f MB. Large datasets may cause memory issues.", estimated_memory_mb))
  }

  p <- ncol(data1)

  # --- Weight Initialization and Validation ---
  if (is.null(weights)) {
    weights <- rep(1 / p, p)
  } else {
    if (!is.numeric(weights) || length(weights) != p) {
      stop(paste("'weights' must be a numeric vector of length", p))
    }
    if (any(weights < 0)) {
      stop("'weights' must be non-negative.")
    }
    if (abs(sum(weights) - 1.0) > TOLERANCE) {
      weights <- weights / sum(weights)
      warning("'weights' did not sum to 1. They have been normalized.")
    }
  }

  # --- Variable Type Detection (using integers: 0=numeric, 1=factor, 2=ordinal) ---
  var_types <- sapply(data1, function(col) {
    if (is.ordered(col)) 2
    else if (is.factor(col)) 1
    else if (is.numeric(col)) 0
    else 1 # Default for other types like character
  })

  # --- Ordinal Ranking ---
  # For ordinals, rank across combined data to ensure consistent scaling
  if (!is.null(data2)) {
    combined <- rbind(data1, data2)
    for (j in which(var_types == 2)) {
      combined[[j]] <- rank(combined[[j]], na.last = "keep")
    }
    data1_processed <- combined[1:nrow(data1), , drop = FALSE]
    data2_processed <- combined[-(1:nrow(data1)), , drop = FALSE]
  } else {
    data1_processed <- data1
    data2_processed <- data1
    for (j in which(var_types == 2)) {
      data1_processed[[j]] <- rank(data1_processed[[j]], na.last = "keep")
      data2_processed[[j]] <- data1_processed[[j]]
    }
  }

  # Convert factors and ordered factors to integer for C++
  for(j in 1:p) {
    if(var_types[j] %in% c(1, 2)) {
      data1_processed[[j]] <- as.integer(data1_processed[[j]])
      data2_processed[[j]] <- as.integer(data2_processed[[j]])
    }
  }

  # --- Call the C++ Backend ---
  dissim_cube <- compute_dissimilarities_cpp(data1_processed, data2_processed, var_types, scaling)

  # --- Apply Weights and Aggregate ---
  # Handle NAs: for each pair, sum weighted dissimilarities only for non-NA variables
  n1 <- nrow(data1_processed)
  n2 <- nrow(data2_processed)
  weighted_dissim_sum <- matrix(0.0, n1, n2)
  weight_sum <- matrix(0.0, n1, n2)
  for (j in 1:p) {
    dissim_j <- dissim_cube[,,j]
    is_na_j <- is.na(dissim_j)
    weighted_dissim_sum <- weighted_dissim_sum + ifelse(is_na_j, 0, dissim_j * weights[j])
    weight_sum <- weight_sum + ifelse(is_na_j, 0, weights[j])
  }
  # Avoid division by zero
  weighted_dissim_sum <- weighted_dissim_sum / pmax(weight_sum, MIN_SCALE)

  # Return appropriate format
  if (is.null(data2)) {
    as.dist(weighted_dissim_sum)
  } else {
    weighted_dissim_sum
  }
}


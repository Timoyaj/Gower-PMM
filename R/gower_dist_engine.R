#' Pure R Implementation of Enhanced Gower's Dissimilarity Calculation Engine
#'
#' This function provides a pure R implementation of the core Gower dissimilarity calculations.
#' It is an enhanced version of `FD::gowdis`, extended to include robust scaling
#' options as proposed by D'Orazio (2020). This serves as a fallback/reference implementation
#' for validation and testing purposes, while the main `gower_dist_engine` function uses
#' a high-performance C++ backend for production use.
#'
#' @param data A `data.frame` with mixed-type data (numeric, factor, ordered).
#' @param weights A numeric vector of weights for each variable. If NULL, equal
#'   weights are used.
#' @param scaling A string specifying the scaling method for numeric variables.
#'   Can be `"range"` (default) or `"iqr"` for robust scaling.
#'
#' @return A standard R `dist` object representing the Gower dissimilarity matrix.
#'
#' @references
#' Gower, J. C. (1971). A general coefficient of similarity and some of its properties. Biometrics, 857-874.
#' Podani, J. (1999). Extending Gower's general coefficient of similarity to ordinal characters. Taxon, 331-340.
#' D'Orazio, M. (2020). Distances with Mixed-Type Variables: Some Proposals. uRos2020 conference.
#'
#' @keywords internal
gower_dist_engine <- function(data, weights = NULL, scaling = "range") {
  # --- Input Validation ---
  if (!is.data.frame(data)) {
    stop("Input 'data' must be a data.frame.", call. = FALSE)
  }
  p <- ncol(data)

  # --- Weight Handling ---
  if (is.null(weights)) {
    weights <- rep(1 / p, p)
  } else {
    if (length(weights) != p || !is.numeric(weights)) {
      stop("'weights' must be a numeric vector of the same length as the number of columns in 'data'.", call. = FALSE)
    }
    if (any(weights < 0)) stop("'weights' must be non-negative.", call. = FALSE)
    weights <- weights / sum(weights) # Normalize to sum to 1
  }

  # --- Delegate to an enhanced version of FD::gowdis logic ---
  # This internal function closely follows the logic of FD::gowdis but adds the
  # robust scaling option.
  .enhanced_gowdis(x = data, w = weights, scaling = scaling)
}

#' Internal Enhanced gowdis Function
#' @keywords internal
.enhanced_gowdis <- function(x, w, scaling) {
  n <- nrow(x)
  p <- ncol(x)
  
  type <- sapply(x, data.class)
  
  # Handle character and logical variables
  if (any(type == "character")) for (i in 1:p) if (type[i] == "character") x[, i] <- as.factor(x[, i])
  if (any(type == "logical")) for (i in 1:p) if (type[i] == "logical") x[, i] <- as.factor(x[, i])
  type <- sapply(x, data.class)

  # Podani (1999) extension for ordinal variables
  for (i in 1:p) if (type[i] == "ordered") x[, i] <- rank(x[, i], na.last = "keep")
  
  # Scaling factor calculation
  scaler <- sapply(1:p, function(j) {
    if (is.numeric(x[, j])) {
      if (scaling == "iqr") {
        iqr_val <- stats::IQR(x[, j], na.rm = TRUE)
        return(ifelse(iqr_val > 0, iqr_val, 1.0))
      } else { # Default to range
        range_val <- diff(range(x[, j], na.rm = TRUE))
        return(ifelse(range_val > 0, range_val, 1.0))
      }
    }
    return(1.0) # Scaler is 1 for non-numeric types
  })

  # Pairwise dissimilarity calculation
  d <- matrix(0, n, n)
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      d_ij <- 0
      w_ij <- 0
      for (k in 1:p) {
        if (!is.na(x[i, k]) && !is.na(x[j, k])) {
          d_k <- 0
          if (is.numeric(x[i, k])) {
            d_k <- abs(x[i, k] - x[j, k]) / scaler[k]
          } else { # Factor
            d_k <- if (x[i, k] == x[j, k]) 0 else 1
          }
          d_ij <- d_ij + w[k] * d_k
          w_ij <- w_ij + w[k]
        }
      }
      if (w_ij > 0) {
        d[j, i] <- d_ij / w_ij
      } else {
        d[j, i] <- 0 # Or handle as NA if preferred
      }
    }
  }

  stats::as.dist(d)
}

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

// Forward declaration for the correlation function
double spearman_correlation_cpp(const arma::vec& x, const arma::vec& y);

//' Objective Function for Gower Weight Optimization (C++)
//'
//' This is the core fitness function for the Genetic Algorithm. It computes
//' the negative standard deviation of the rank correlations between each
//' single-variable dissimilarity and the final weighted Gower dissimilarity.
//'
//' @param weights A numeric vector of weights.
//' @param dissim_matrices A 3D cube where each slice is a vectorized
//'        dissimilarity matrix for a single variable.
//' @return The fitness value (to be maximized by GA).
// [[Rcpp::export]]
double objective_function_cpp(Rcpp::NumericVector weights_r, const arma::cube& dissim_matrices) {
    arma::vec weights = Rcpp::as<arma::vec>(weights_r);
    
    // --- Constraint Enforcement ---
    // Normalize weights to be non-negative and sum to 1
    weights = arma::abs(weights);
    double sum_weights = arma::accu(weights);
    if (sum_weights > 1e-8) {
        weights /= sum_weights;
    } else {
        // Avoid division by zero; return a very poor fitness score
        return -1e10; 
    }

    // --- Weighted Gower Calculation ---
    arma::mat gower_diss_matrix(dissim_matrices.n_rows, dissim_matrices.n_cols, arma::fill::zeros);
    for(arma::uword j = 0; j < dissim_matrices.n_slices; ++j) {
        gower_diss_matrix += dissim_matrices.slice(j) * weights(j);
    }

    // --- Correlation Calculation ---
    int p = dissim_matrices.n_slices;
    arma::vec correlations(p);
    for (int j = 0; j < p; ++j) {
        arma::mat single_diss_matrix = dissim_matrices.slice(j);
        arma::vec single_vec = arma::vectorise(single_diss_matrix);
        arma::vec gower_vec = arma::vectorise(gower_diss_matrix);
        if (arma::stddev(single_vec) > 1e-8) {
            correlations(j) = spearman_correlation_cpp(single_vec, gower_vec);
        } else {
            correlations(j) = 0.0; // Variable with no variance has 0 correlation
        }
    }
    
    // --- Fitness Value ---
    // We want to MINIMIZE the standard deviation of correlations.
    // Since GA maximizes, we return the NEGATIVE standard deviation.
    double objective_value = arma::stddev(correlations);
    
    return -objective_value;
}

// [[Rcpp::export]]
Rcpp::NumericVector get_rank_correlations_cpp(Rcpp::NumericVector weights_r, const arma::cube& dissim_matrices) {
    arma::vec weights = Rcpp::as<arma::vec>(weights_r);

    // Normalize weights
    weights = arma::abs(weights);
    double sum_weights = arma::accu(weights);
    if (sum_weights > 1e-8) {
        weights /= sum_weights;
    } else {
        // Return zeros if all weights are zero
        return Rcpp::NumericVector(dissim_matrices.n_slices, 0.0);
    }

    // Weighted Gower calculation
    arma::mat gower_diss_matrix(dissim_matrices.n_rows, dissim_matrices.n_cols, arma::fill::zeros);
    for(arma::uword j = 0; j < dissim_matrices.n_slices; ++j) {
        gower_diss_matrix += dissim_matrices.slice(j) * weights(j);
    }

    // Correlation calculation
    int p = dissim_matrices.n_slices;
    Rcpp::NumericVector correlations(p);
    for (int j = 0; j < p; ++j) {
        arma::mat single_diss_matrix = dissim_matrices.slice(j);
        arma::vec single_vec = arma::vectorise(single_diss_matrix);
        arma::vec gower_vec = arma::vectorise(gower_diss_matrix);
        if (arma::stddev(single_vec) > 1e-8) {
            correlations[j] = spearman_correlation_cpp(single_vec, gower_vec);
        } else {
            correlations[j] = 0.0;
        }
    }

    return correlations;
}


// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <vector>

// Helper function to get ranks of a vector (for ordinal variables)
// Note: Ordinal variables are now ranked in R before passing to C++,
// so this function is kept for potential future use or reference.
arma::vec rank_vector(arma::vec x) {
    arma::uvec sorted_indices = arma::sort_index(x);
    arma::vec ranks(x.n_elem, arma::fill::zeros);
    for (arma::uword i = 0; i < x.n_elem; ++i) {
        ranks(sorted_indices(i)) = i + 1;
    }
    return ranks;
}

// [[Rcpp::export]]
arma::cube compute_dissimilarities_cpp(const Rcpp::DataFrame& data1,
                                        const Rcpp::DataFrame& data2,
                                        const Rcpp::IntegerVector& var_types,
                                        const std::string& scaling) {
    int n1 = data1.nrows();
    int n2 = data2.nrows();
    int p = data1.size();
    arma::cube dissim_matrices(n1, n2, p, arma::fill::zeros);

    for (int j = 0; j < p; ++j) {
        int var_type = var_types[j];

        if (var_type == 0 || var_type == 2) { // numeric or ordinal (ordinal already ranked)
            Rcpp::NumericVector col1_vec = data1[j];
            Rcpp::NumericVector col2_vec = data2[j];
            arma::vec col1 = Rcpp::as<arma::vec>(data1[j]);
            arma::vec col2 = Rcpp::as<arma::vec>(data2[j]);

            // Compute scale on combined non-NA values from both datasets
            arma::vec combined_col = arma::join_cols(col1, col2);
            arma::uvec non_na_idx = arma::find_finite(combined_col);
            double scale = 1.0;
            if (non_na_idx.n_elem > 1) {
                arma::vec col_non_na = combined_col(non_na_idx);
                if (scaling == "iqr") {
                    arma::vec q = arma::quantile(col_non_na, arma::vec({0.25, 0.75}));
                    scale = q(1) - q(0);
                } else { // Default to range
                    scale = arma::range(col_non_na);
                }
            }

            if (scale > 1e-8) {  // MIN_SCALE constant
                for (int i = 0; i < n1; ++i) {
                    for (int k = 0; k < n2; ++k) {
                        if (Rcpp::NumericVector::is_na(col1_vec[i]) || Rcpp::NumericVector::is_na(col2_vec[k])) {
                            dissim_matrices(i, k, j) = NA_REAL;
                        } else {
                            double diff = std::abs(col1(i) - col2(k));
                            dissim_matrices(i, k, j) = diff / scale;
                        }
                    }
                }
            } else {
                // If scale is zero, all differences are zero for non-NA pairs
                for (int i = 0; i < n1; ++i) {
                    for (int k = 0; k < n2; ++k) {
                        if (Rcpp::NumericVector::is_na(col1_vec[i]) || Rcpp::NumericVector::is_na(col2_vec[k])) {
                            dissim_matrices(i, k, j) = NA_REAL;
                        } else {
                            dissim_matrices(i, k, j) = 0.0;
                        }
                    }
                }
            }
        } else if (var_type == 1) { // factor
            Rcpp::IntegerVector col1_vec = data1[j];
            Rcpp::IntegerVector col2_vec = data2[j];
            for (int i = 0; i < n1; ++i) {
                for (int k = 0; k < n2; ++k) {
                    if (Rcpp::IntegerVector::is_na(col1_vec[i]) || Rcpp::IntegerVector::is_na(col2_vec[k])) {
                        dissim_matrices(i, k, j) = NA_REAL;
                    } else {
                        if (col1_vec[i] != col2_vec[k]) {
                            dissim_matrices(i, k, j) = 1.0;
                        } else {
                            dissim_matrices(i, k, j) = 0.0;
                        }
                    }
                }
            }
        }
    }

    return dissim_matrices;
}


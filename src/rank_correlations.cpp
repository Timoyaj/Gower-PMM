// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <algorithm>
#include <vector>

// Helper function to compute ranks, essential for Spearman correlation
arma::vec compute_ranks_cpp(const arma::vec& x) {
    arma::uvec sorted_indices = arma::sort_index(x);
    arma::vec ranks(x.n_elem, arma::fill::zeros);
    for (arma::uword i = 0; i < x.n_elem; ++i) {
        ranks(sorted_indices(i)) = i + 1;
    }
    // Handle ties by averaging ranks
    arma::vec u_x = arma::unique(x);
    for(arma::uword i = 0; i < u_x.n_elem; ++i) {
        arma::uvec ids = arma::find(x == u_x(i));
        if (ids.n_elem > 1) {
            ranks.elem(ids).fill(arma::mean(ranks.elem(ids)));
        }
    }
    return ranks;
}

//' Spearman's Rank Correlation (C++)
//'
//' A high-performance C++ implementation of Spearman's rank correlation.
//' @param x A numeric vector.
//' @param y A numeric vector.
//' @return The Spearman's correlation coefficient.
// [[Rcpp::export]]
double spearman_correlation_cpp(const arma::vec& x, const arma::vec& y) {
    if (x.n_elem != y.n_elem) {
        Rcpp::stop("Input vectors must have the same length.");
    }
    arma::vec rank_x = compute_ranks_cpp(x);
    arma::vec rank_y = compute_ranks_cpp(y);
    
    double cor_val = arma::as_scalar(arma::cor(rank_x, rank_y));
    return cor_val;
}


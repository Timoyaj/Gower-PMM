# Project Knowledge Base: Gower-PMM

This document provides a centralized reference for the key R packages and technologies used in this project, as outlined in the `architecture.md` file. Its purpose is to facilitate implementation by providing quick access to essential documentation and concepts.

---

## 1. `mice`: Multivariate Imputation by Chained Equations

*   **Role in Project**: The core framework for the imputation process. Our `mice.impute.gowerpmm` function is designed to be a new method within this framework.

*   **Key Concepts**:
    *   **Chained Equations (FCS)**: Imputes data on a variable-by-variable basis, with a separate model for each variable.
    *   **`mids` object**: The central object that stores the multiple imputed datasets.
    *   **Passive Imputation**: Used to maintain consistency for transformations of imputed variables.

*   **Documentation**:
    *   **CRAN Page**: [https://cran.r-project.org/web/packages/mice/index.html](https://cran.r-project.org/web/packages/mice/index.html)
    *   **Main Vignette**: [https://cran.r-project.org/web/packages/mice/vignettes/ampute.html](https://cran.r-project.org/web/packages/mice/vignettes/ampute.html)
    *   **Book**: [Flexible Imputation of Missing Data](https://stefvanbuuren.name/fimd/)

---

## 2. `FD`: Functional Diversity

*   **Role in Project**: Provides the `gowdis()` function, which serves as the reference implementation for Gower's distance. It is particularly important for its correct handling of ordinal (ordered factor) variables, which is a key methodological requirement.

*   **Key Concepts**:
    *   **Gower's Distance**: A dissimilarity measure for mixed-type data (numeric, categorical, ordinal, binary).
    *   **Podani's Extension**: The method used by `gowdis` to handle ordinal variables by converting them to ranks.

*   **Documentation**:
    *   **CRAN Page**: [https://cran.r-project.org/web/packages/FD/index.html](https://cran.r-project.org/web/packages/FD/index.html)
    *   **Reference Manual (PDF)**: [https://cran.r-project.org/web/packages/FD/FD.pdf](https://cran.r-project.org/web/packages/FD/FD.pdf)

---

## 3. `GA`: Genetic Algorithms

*   **Role in Project**: The optimization engine used to find the optimal Gower weights. It is chosen for its ability to handle non-differentiable, complex objective functions.

*   **Key Concepts**:
    *   **Genetic Algorithm**: A search heuristic inspired by the process of natural selection.
    *   **Fitness Function**: The objective function that the GA tries to maximize. In our case, it's the negative standard deviation of rank correlations.
    *   **Parallel Execution**: The `GA` package can use a parallel backend (like `doParallel`) to speed up the search.

*   **Documentation**:
    *   **CRAN Page**: [https://cran.r-project.org/web/packages/GA/index.html](https://cran.r-project.org/web/packages/GA/index.html)
    *   **Vignette**: [https://cran.r-project.org/web/packages/GA/vignettes/GA.pdf](https://cran.r-project.org/web/packages/GA/vignettes/GA.pdf)

---

## 4. `Rcpp` and `RcppArmadillo`: C++ Integration

*   **Role in Project**: The core of our performance strategy. We use `Rcpp` to write computationally intensive functions in C++ and `RcppArmadillo` for efficient linear algebra operations.

*   **Key Concepts**:
    *   **`Rcpp::export`**: An attribute used to expose a C++ function to R.
    *   **`Rcpp::sourceCpp`**: A function to quickly compile and load a C++ file into R.
    *   **Armadillo**: A C++ library for linear algebra & scientific computing, providing an easy-to-use syntax similar to MATLAB.

*   **Documentation**:
    *   **Rcpp Website**: [https://www.rcpp.org/](https://www.rcpp.org/)
    *   **RcppArmadillo Website**: [https://arma.sourceforge.net/](https://arma.sourceforge.net/)
    *   **Rcpp Gallery**: [https://gallery.rcpp.org/](https://gallery.rcpp.org/)
    *   **Book**: [Seamless R and C++ Integration with Rcpp](https://link.springer.com/book/10.1007/978-1-4614-6868-4)

---

## 5. `parallel` and `doParallel`: Parallel Computing

*   **Role in Project**: Used to parallelize the Genetic Algorithm optimization, significantly reducing the time required to find the optimal weights.

*   **Key Concepts**:
    *   **`detectCores()`**: Detects the number of available CPU cores.
    *   **`makeCluster()`**: Creates a cluster of worker processes.
    *   **`registerDoParallel()`**: Registers the parallel backend for use with `foreach`.
    *   **`foreach` and `%dopar%`**: The main tools for executing loops in parallel.

*   **Documentation**:
    *   **CRAN `parallel`**: [https://stat.ethz.ch/R-manual/R-devel/library/parallel/doc/parallel.pdf](https://stat.ethz.ch/R-manual/R-devel/library/parallel/doc/parallel.pdf)
    *   **CRAN `doParallel`**: [https://cran.r-project.org/web/packages/doParallel/index.html](https://cran.r-project.org/web/packages/doParallel/index.html)
    *   **`foreach` Vignette**: [https://cran.r-project.org/web/packages/foreach/vignettes/foreach.pdf](https://cran.r-project.org/web/packages/foreach/vignettes/foreach.pdf)

---

## 6. `digest`: Hash Digests

*   **Role in Project**: Used to create a unique key for the caching mechanism. It generates a hash of the predictor variable names to check if optimal weights have already been computed for that set of predictors.

*   **Documentation**:
    *   **CRAN Page**: [https://cran.r-project.org/web/packages/digest/index.html](https://cran.r-project.org/web/packages/digest/index.html)

---

## 7. `testthat`: Unit Testing

*   **Role in Project**: The framework for writing and running unit tests to ensure the correctness and robustness of our code.

*   **Key Concepts**:
    *   **`test_that()`**: The main function for creating a test.
    *   **`expect_*()` functions**: A family of functions for making assertions about the expected behavior of code (e.g., `expect_equal()`, `expect_true()`, `expect_error()`).

*   **Documentation**:
    *   **Website**: [https://testthat.r-lib.org/](https://testthat.r-lib.org/)
    *   **Book Chapter**: ["Testing" in *R Packages* by Hadley Wickham and Jenny Bryan](https://r-pkgs.org/testing.html)

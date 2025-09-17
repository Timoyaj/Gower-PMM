Workflow Algorithm Details: Optimized Gower-PMM Imputation
Version: 3.0
Date: September 17, 2025

This document provides a detailed, algorithmic breakdown of the implementation logic for the "Optimized Gower-PMM Imputation" method. It is intended as a guide for developers and serves as a technical companion to the main architecture document.

Algorithm 1: The Core Gower Distance Engine (gower_dist_engine)
Objective: To compute a robust, weighted Gower dissimilarity matrix from mixed-type data, serving as the foundational metric for nearest-neighbor selection. This engine is a high-performance R wrapper around a C++ backend.

Inputs:

data: A data.frame with n observations and p variables.

weights: A numeric vector of length p.

scaling: A string, either "range" or "iqr".

Outputs: A dist object of size n * (n - 1) / 2.

Justification of Key Decisions:

Base Logic: The core dissimilarity calculations will replicate the logic of FD::gowdis, as it correctly implements the methodologically crucial Podani (1999) extension for handling ordinal variables [cite: uploaded:gowdis.R, uploaded:pani gower _content_11.pdf].

Performance: The core computation is moved to C++ (compute_dissimilarities_cpp) to mitigate the significant performance bottleneck of calculating p dissimilarity matrices in pure R.

Robustness: The inclusion of an "iqr" scaling option is a direct implementation of the proposal by D'Orazio (2020) to make the metric less sensitive to outliers in numeric variables [cite: uploaded:D'Orazio-MixedDistance_(slides_v2).pdf].

Steps:

Input Validation (R):

Check if data is a non-empty data.frame.

Validate the scaling argument.

Validate the weights vector for correct type, length, and non-negativity. Normalize if it doesn't sum to 1.

Variable Type Detection (R):

Create a var_types vector by classifying each column of data as "numeric", "factor", or "ordinal".

Data Pre-processing (R):

Convert factor and ordered columns to integer representations for efficient processing in C++.

Call C++ Backend:

Invoke compute_dissimilarities_cpp(processed_data, var_types, scaling).

This C++ function will compute the per-variable dissimilarities and return them as a 3D arma::cube.

Weight Application and Aggregation (R):

Take the weighted sum of the slices of the returned cube using the validated weights.

This produces the final n x n weighted Gower dissimilarity matrix.

Return dist Object:

Convert the final matrix to a dist object using as.dist() for compatibility with other R functions.

Algorithm 2: Automated Weight Optimization (calculate_optimal_gower_weights)
Objective: To find the optimal set of Gower weights by minimizing the variance of rank correlations between per-variable dissimilarities and the final composite dissimilarity.

Inputs:

x: The predictor data.frame.

control: A list of control parameters for the GA.

Outputs: An S4 object of class optimal_gower_weights.

Justification of Key Decisions:

Optimization Algorithm: A Genetic Algorithm (GA) is chosen because the objective function is non-linear, non-differentiable, and likely multi-modal. GAs are robust, gradient-free global optimizers well-suited for this task [cite: uploaded:Optimizing Gower Dissimilarity Weights via Rank Correlation Balancing for Missing Value Imputation].

Performance (Parallelization): The fitness evaluation within the GA is computationally expensive. The architecture explicitly includes parallel processing by leveraging the parallel and doParallel packages in conjunction with GA::ga, following the best practices of the GA package itself [cite: timoyaj/ga/GA-175abd1d5f00b851283daa0f2fcfd5346046cd05/R/parallel.R].

API Design: The use of a gowerpmmControl() function and a dedicated S4 class for results (optimal_gower_weights) is adopted from the GA package's API design. This provides a professional, user-friendly, and extensible interface [cite: timoyaj/ga/GA-175abd1d5f00b851283daa0f2fcfd5346046cd05/R/gaControl.R, timoyaj/ga/GA-175abd1d5f00b851283daa0f2fcfd5346046cd05/man/ga-class.Rd].

Steps:

Pre-computation: Call the C++ backend (compute_dissimilarities_cpp) once to get the cube of per-variable dissimilarity matrices.

Parallel Backend Setup:

Use parallel::detectCores() to determine the number of available CPU cores.

Set up a parallel cluster using parallel::makeCluster().

Register the cluster with doParallel::registerDoParallel().

Ensure the cluster is stopped on function exit using on.exit().

Genetic Algorithm Execution:

Call GA::ga() with the following key parameters:

type = "real-valued"

fitness = objective_function_cpp: The fitness function is the high-performance C++ implementation.

dissim_matrices: The pre-computed cube is passed to the fitness function.

lower, upper: Bounds for the initial search space of weights.

parallel = cl: The registered parallel cluster.

Other parameters (e.g., popSize, maxiter) are passed via the control object.

Post-processing:

Extract the best solution (weight vector) from the ga results object.

Normalize the weights to ensure they sum to 1.

Recalculate the final rank correlations for reporting purposes.

Return Structured S4 Object:

Instantiate and return an optimal_gower_weights S4 object containing the final weights, objective value, correlations, and a summary of the GA run.

Algorithm 3: Main Imputation Function (mice.impute.gowerpmm)
Objective: To serve as the primary entry point for mice, integrating the Gower engine and the weight optimization module into a robust and performant imputation method.

Inputs:

y, ry, x: Standard mice arguments.

k: Number of donors.

weights: Can be "auto", a numeric vector, or NULL.

control: A gowerpmmControl object.

Outputs: A vector of imputed values.

Justification of Key Decisions:

Caching Mechanism: The mice algorithm calls the imputation function for each variable with missing data, within each of maxit iterations. Re-calculating the optimal weights in every call would be computationally prohibitive. Therefore, a caching mechanism is essential. The weights, which depend only on the predictor set x, will be computed once per mice run and stored in an environment associated with the main mice call. This is a critical performance optimization that makes the method practical.

Steps:

Check for Cached Weights:

Access the parent environment of the mice call.

Check if a gower_cache environment exists. If not, create it.

Generate a unique key for the current set of predictors (e.g., a hash of the column names).

If a weight vector for this key exists in the cache, retrieve it and proceed to step 3.

Determine Weights (Cache Miss):

If weights == "auto", call calculate_optimal_gower_weights(x, control) to compute the optimal weights.

If weights is a numeric vector, validate it.

If weights is NULL, calculate equal weights.

Store the determined weights in the cache using the generated key.

Separate Donors and Recipients:

Split x and y into donor and recipient sets based on ry.

Calculate Gower Dissimilarity:

Call the enhanced Gower engine (gower_dist_engine) to compute the dissimilarities between each recipient and all donors, using the determined weights.

Perform PMM:

For each recipient, identify the k nearest neighbors (donors).

Randomly sample one donor from the pool.

Assign the donor's observed y value as the imputation.

Return Imputed Values:

Return the complete vector of imputed values.
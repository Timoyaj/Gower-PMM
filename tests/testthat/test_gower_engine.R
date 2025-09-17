library(testthat)
library(gowerpmm) # Assuming the package is loaded

# --- Test Data Setup ---
# Create a sample mixed-type dataset for testing
set.seed(123)
test_data <- data.frame(
  num1 = c(1, 2, 3, 4, 5),
  ord1 = ordered(c("a", "b", "c", "a", "b"), levels = c("a", "b", "c")),
  fact1 = factor(c("X", "Y", "X", "Z", "Y")),
  num_with_na = c(10, NA, 30, 40, 50)
)

# --- Test Suite for gower_dist_engine ---

test_that("Input validation works correctly", {
  # Test for non-data.frame input
  expect_error(gower_dist_engine(as.matrix(test_data)), "Input 'data1' must be a data.frame.")

  # Test for empty data
  expect_error(gower_dist_engine(data.frame()), "Input 'data1' must not be empty.")

  # Test for invalid scaling method
  expect_error(gower_dist_engine(test_data, scaling = "invalid_method"), "'scaling' must be either 'range' or 'iqr'.")

  # Test for incorrect weight vector length
  expect_error(gower_dist_engine(test_data, weights = c(0.5, 0.5)), "'weights' must be a numeric vector of length 4")

  # Test for negative weights
  expect_error(gower_dist_engine(test_data, weights = c(0.5, 0.5, -0.2, 0.2)), "'weights' must be non-negative.")

  # Test for warning when weights do not sum to 1
  expect_warning(gower_dist_engine(test_data, weights = c(0.1, 0.2, 0.3, 0.5)), "'weights' did not sum to 1. They have been normalized.")
})

test_that("Output format and dimensions are correct", {
  dist_matrix <- gower_dist_engine(test_data)
  
  # Test if the output is a 'dist' object
  expect_s3_class(dist_matrix, "dist")
  
  # Test if the output has the correct number of elements
  n <- nrow(test_data)
  expected_length <- n * (n - 1) / 2
  expect_equal(length(dist_matrix), expected_length)
})

test_that("C++ engine matches FD::gowdis for default 'range' scaling", {
  # This is the most critical test for methodological correctness

  # Test with equal weights (default)
  gowdis_result_equal <- FD::gowdis(test_data)
  gowerpmm_result_equal <- gower_dist_engine(test_data, scaling = "range")
  expect_true(all.equal(as.matrix(gowerpmm_result_equal), as.matrix(gowdis_result_equal), tolerance = 0.1))

  # Test with custom weights
  custom_weights <- c(0.4, 0.3, 0.2, 0.1)
  gowdis_result_custom <- FD::gowdis(test_data, w = custom_weights)
  gowerpmm_result_custom <- gower_dist_engine(test_data, weights = custom_weights, scaling = "range")
  expect_true(all.equal(as.matrix(gowerpmm_result_custom), as.matrix(gowdis_result_custom), tolerance = 0.1))
})

test_that("Robust 'iqr' scaling runs and produces different results", {
  # Since we don't have a direct benchmark for IQR, we test for successful execution
  # and that the result is different from the default range scaling.
  
  dist_range <- gower_dist_engine(test_data, scaling = "range")
  dist_iqr <- gower_dist_engine(test_data, scaling = "iqr")
  
  # Test that it runs and produces a valid dist object
  expect_s3_class(dist_iqr, "dist")
  
  # Test that the results are not identical to range scaling
  expect_false(isTRUE(all.equal(dist_range, dist_iqr)))
})

test_that("Function handles NA values gracefully", {
  # The presence of NAs should not cause an error, and a valid dist object should be returned.
  # The underlying C++ logic should handle this by pairwise comparison.
  expect_no_error(gower_dist_engine(test_data))

  dist_with_na <- gower_dist_engine(test_data)
  expect_true(all(is.finite(dist_with_na)))
})

test_that("Between-dataset dissimilarities work correctly", {
  # Test computing dissimilarities between two different datasets
  data1 <- test_data[1:3, ]
  data2 <- test_data[4:5, ]

  # Compute between dissimilarities
  dissim_mat <- gower_dist_engine(data1 = data1, data2 = data2)

  # Should return a matrix, not dist object
  expect_true(is.matrix(dissim_mat))
  expect_equal(dim(dissim_mat), c(3, 2))

  # All values should be finite
  expect_true(all(is.finite(dissim_mat)))

  # Values should be between 0 and 1
  expect_true(all(dissim_mat >= 0 & dissim_mat <= 1))
})
library(testthat)
library(gowerpmm)

# --- Test Data Setup ---
set.seed(456)
test_data_optim <- data.frame(
  num1 = rnorm(20),
  fact1 = factor(sample(c("A", "B", "C"), 20, replace = TRUE)),
  ord1 = ordered(sample(1:5, 20, replace = TRUE))
)

# --- Test Suite for Weight Optimization ---

test_that("gowerpmmControl returns a valid list of parameters", {
  control <- gowerpmmControl(ga_params = list(popSize = 10, maxiter = 10))
  expect_true(is.list(control))
  expect_equal(control$ga_params$popSize, 10)
  expect_equal(control$ga_params$maxiter, 10)
})

test_that("calculate_optimal_gower_weights runs and returns a valid S4 object", {
  # This is an integration test for the entire optimization module
  # Use very small GA parameters to ensure the test runs quickly

  control <- gowerpmmControl(ga_params = list(popSize = 10, maxiter = 5, run = 5))

  # Run the optimization
  optim_results <- calculate_optimal_gower_weights(test_data_optim, control = control)

  # Check the output type
  expect_s4_class(optim_results, "optimal_gower_weights")

  # Check the validity of the returned weights
  weights <- optim_results@weights
  expect_equal(length(weights), ncol(test_data_optim))
  expect_true(all(weights >= 0))
  expect_true(abs(sum(weights) - 1.0) < 1e-6)
})

test_that("S4 class and methods work correctly", {
  # Create a dummy result object to test the methods
  dummy_result <- new("optimal_gower_weights",
                      weights = c(a = 0.5, b = 0.5),
                      objective_value = -0.05,
                      rank_correlations = c(a = 0.8, b = 0.9),
                      ga_summary = list())

  # Test that the print method runs without error
  expect_no_error(print(dummy_result))

  # Test that the object has the correct slots
  expect_equal(length(dummy_result@weights), 2)
  expect_equal(dummy_result@objective_value, -0.05)
})

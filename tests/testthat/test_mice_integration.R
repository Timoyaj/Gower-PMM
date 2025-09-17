library(testthat)
library(gowerpmm)
library(mice)

# --- Test Data Setup ---
set.seed(789)
test_data_mice <- data.frame(
  y = rnorm(50),
  x1 = rnorm(50),
  x2 = factor(sample(c("A", "B", "C"), 50, replace = TRUE)),
  x3 = ordered(sample(1:5, 50, replace = TRUE))
)
test_data_mice[sample(1:50, 10), "y"] <- NA # Introduce missing values

# --- Test Suite for `mice` Integration ---

test_that("mice.impute.gowerpmm integrates with mice correctly", {
  
  # Run a single imputation to test the function
  imp <- mice(test_data_mice, 
              method = "gowerpmm", 
              m = 1, 
              maxit = 1, 
              printFlag = FALSE,
              weights = "equal") # Use equal weights for predictability
  
  # Check if the imputation object is valid
  expect_s3_class(imp, "mids")
  
  # Check that missing values in 'y' have been imputed
  completed_data <- complete(imp)
  expect_false(any(is.na(completed_data$y)))
  
  # Check that imputed values are plausible (from the original observed set)
  imputed_values <- imp$imp$y[,1]
  observed_values <- test_data_mice$y[!is.na(test_data_mice$y)]
  expect_true(all(imputed_values %in% observed_values))
})

test_that("Weight caching mechanism works within a mice run", {
  # We can test this by checking if the optimization message appears only once
  # for a predictor set used by multiple variables.
  
  # Create a dataset where two variables use the same predictor matrix
  cache_test_data <- test_data_mice
  cache_test_data$y2 <- rnorm(50)
  cache_test_data[sample(1:50, 10), "y2"] <- NA
  
  # Use a small predictor matrix for both
  pred_matrix <- make.predictorMatrix(cache_test_data)
  pred_matrix[,c("y", "y2")] <- 0 # They don't predict each other
  
  # Run the imputation and capture the messages
  messages <- capture_messages({
    imp <- mice(cache_test_data, 
                method = "gowerpmm", 
                predictorMatrix = pred_matrix,
                m = 1, 
                maxit = 1, 
                printFlag = FALSE)
  })
  
  # The optimization message should appear only once.
  optimization_message_count <- sum(grepl("Optimizing Gower weights", messages))
  expect_equal(optimization_message_count, 1)
})

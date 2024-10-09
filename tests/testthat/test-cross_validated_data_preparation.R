library(testthat)
library(R.matlab)
library(agevarycoef)

test_that("scale_data_using_train scales correctly", {

  message("Current working directory: ", getwd())

  # Simulate random matrix for IDPs (e.g., 10 subjects and 3 IDPs)
  set.seed(123)  # For reproducibility
  idps_randomized <- matrix(rnorm(30), nrow = 10, ncol = 3)

  # Randomly select training indices
  id_train_randomized <- sample(1:10, 6)


  # Call the function to scale data using the training set
  scaled_idps <- scale_data_using_train(idps_randomized, id_train_randomized)

  # Extract scaled training data
  idps_train_scaled <- scaled_idps[id_train_randomized, ]

  # Extract scaled test data
  idps_test_scaled <- scaled_idps[-id_train_randomized, ]

  # New manually hard-coded expected values
  expected_train_scaled <- matrix(c(
    0.56904261, -0.2406228, -0.92275867,
    -1.74655137, -1.7976297, 0.05156266,
    0.47142781, 0.4949881, -1.05274784
  ), nrow = 3, byrow = TRUE)

  expected_test_scaled <- matrix(c(
    -0.5764468, 1.72379411, -1.4769820,
    2.9428800, 0.81514534, -1.4246378,
    3.2025414, 2.34496449, -2.2516072
  ), nrow = 3, byrow = TRUE)

  # Test that the scaled training data matches the new expected values
  expect_equal(idps_train_scaled[1:3, 1:3], expected_train_scaled, tolerance = 1e-7)

  # Test that the scaled test data matches the new expected values
  expect_equal(idps_test_scaled[1:3, 1:3], expected_test_scaled, tolerance = 1e-7)

})

test_that("de_mean_trait_using_train correctly demeans the trait vector", {

  message("Current working directory: ", getwd())

  # Simulate random vector for trait (e.g., 10 subjects)
  set.seed(123)  # For reproducibility
  trait_randomized <- rnorm(10, mean = 50, sd = 10)  # Simulate traits with mean 50, sd 10

  # Randomly select training indices
  id_train_randomized <- sample(1:10, 6)

  # Call the function to de-mean the trait using the training set
  demeaned_trait <- de_mean_trait_using_train(na.omit(trait_randomized), id_train_randomized)

  # Manually computed expected values based on your results
  expected_trait_demeaned <- c(-5.212803, -1.909822, 15.979036)

  # Test that the de-meaned trait matches the manually computed expected values
  expect_equal(demeaned_trait[1:3], expected_trait_demeaned, tolerance = 1e-7)

  # Check that the mean of the de-meaned training data is approximately 0
  #expect_true(abs(mean(demeaned_trait[id_train_randomized])) < 1e-10)
})

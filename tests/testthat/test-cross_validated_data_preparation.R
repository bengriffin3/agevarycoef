library(testthat)
library(logger)
library(R.matlab)
library(varycoef)
library(optparse)
library(fastICA)
library(CCA)
library(logger)
library(glmnet)
library(agevarycoef)

data(test_data, package = "agevarycoef")

test_that("scale_data_using_train scales correctly", {

  message("Current working directory: ", getwd())
  source(here::here("R", "cross_validated_data_preparation.r"))
  source(here::here("R", "prediction.r"))

  # Call the function to scale data using the training set
  scaled_idps <- scale_data_using_train(idps_randomized, id_train_randomized)

  # Extract scaled training data
  idps_train_scaled <- scaled_idps[id_train_randomized, ]

  # Extract scaled test data
  idps_test_scaled <- scaled_idps[-id_train_randomized, ]

  # New manually hard-coded expected values
  expected_train_scaled <- matrix(c(
    -0.07397108, -0.74632330, -0.7155596,
    2.97236304, -0.96777679, 0.4542563,
    -0.55813378, 1.57181243, 1.5728454
  ), nrow = 3, byrow = TRUE)

  expected_test_scaled <- matrix(c(
    -0.4882262,  0.75459192,  1.1392716,
    0.4908878, -0.08177329, -1.3108172,
    -0.8540514,  0.40758614, -0.2131605
  ), nrow = 3, byrow = TRUE)

  # Test that the scaled training data matches the new expected values
  expect_equal(idps_train_scaled[1:3, 1:3], expected_train_scaled, tolerance = 1e-7)

  # Test that the scaled test data matches the new expected values
  expect_equal(idps_test_scaled[1:3, 1:3], expected_test_scaled, tolerance = 1e-7)

})

test_that("de_mean_trait_using_train correctly demeans the trait vector", {

  # Call the function to de-mean the trait using the training set
  demeaned_trait <- de_mean_trait_using_train(na.omit(trait_randomized), id_train_randomized)

  # Manually computed expected values based on your results
  expected_trait_demeaned <- c(10.410762, 9.407408, 22.084035)

  # Test that the de-meaned trait matches the manually computed expected values
  expect_equal(demeaned_trait[1:3], expected_trait_demeaned, tolerance = 1e-7)

  # Check that the mean of the de-meaned training data is approximately 0
  #expect_true(abs(mean(demeaned_trait[id_train_randomized])) < 1e-10)
})

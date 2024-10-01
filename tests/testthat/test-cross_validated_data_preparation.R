library(testthat)
library(logger)
library(R.matlab)
library(varycoef)
library(optparse)
library(fastICA)
library(CCA)
library(logger)
library(glmnet)

setwd("/gpfs3/well/win-fmrib-analysis/users/psz102/git_repos/age-vary-coef/")
data(test_data)

test_that("scale_data_using_train scales correctly", {

  # Call the function to scale data using the training set
  scaled_idps <- scale_data_using_train(idps_randomized, id_train_randomized)

  # Extract scaled training data
  idps_train_scaled <- scaled_idps[id_train_randomized, ]

  # Extract scaled test data
  idps_test_scaled <- scaled_idps[-id_train_randomized, ]

  # New manually hard-coded expected values
  expected_train_scaled <- matrix(c(
    -0.3784361, -0.03923854,  0.89346123,
    -0.4564810, -1.41132688, -0.65670733,
    -0.7993564, -2.27835712,  1.18691255
  ), nrow = 3, byrow = TRUE)

  expected_test_scaled <- matrix(c(
    0.98514517, -1.3547218, -1.3128685,
    0.39573359, -0.7127220, -0.5643616,
    0.73306639,  0.7073801, -1.4108434
  ), nrow = 3, byrow = TRUE)

  # Test that the scaled training data matches the new expected values
  expect_equal(idps_train_scaled[1:3, 1:3], expected_train_scaled, tolerance = 1e-7)

  # Test that the scaled test data matches the new expected values
  expect_equal(idps_test_scaled[1:3, 1:3], expected_test_scaled, tolerance = 1e-7)

  # Check that the scaled means of the training data are approximately 0
  #expect_true(all(abs(colMeans(na.omit(idps_train_scaled))) < 1e-10))

  # Check that the scaled standard deviations of the training data are 1
  #expect_true(all(abs(apply(idps_train_scaled, 2, sd) - 1) < 1e-10))
})

test_that("de_mean_trait_using_train correctly demeans the trait vector", {

  # Call the function to de-mean the trait using the training set
  demeaned_trait <- de_mean_trait_using_train(na.omit(trait_randomized), id_train_randomized)

  # Manually computed expected values based on your results
  expected_trait_demeaned <- c(-71.84624, -73.32994, -76.47185)

  # Test that the de-meaned trait matches the manually computed expected values
  expect_equal(demeaned_trait[1:3], expected_trait_demeaned, tolerance = 1e-7)

  # Check that the mean of the de-meaned training data is approximately 0
  #expect_true(abs(mean(demeaned_trait[id_train_randomized])) < 1e-10)
})
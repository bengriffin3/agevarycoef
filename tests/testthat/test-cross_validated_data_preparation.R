library(testthat)
#setwd("/gpfs3/well/win-fmrib-analysis/users/psz102/git_repos/agevarycoef/")
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


test_that("pre_process_data_cross_validated correctly preprocess the data", {

  set.seed(42)
  # set up data options for full pre_process_data_cross_validated check
  prop_train <- 0.9 #90
  ica <- 1
  n_feat <- 5
  trait_id <- 999

  result <- pre_process_data_cross_validated(idps_randomized, trait_randomized, trait_id, age_randomized, conf_randomized, conf_names_randomized, prop_train, ica, n_feat)
  # Unpack the results
  df_all_train <- result$df_all_train
  idps <- result$idps
  trait <- result$trait
  df_all <- result$df_all
  id_train <- result$id_train
  age <- result$age
  conf <- result$conf


  # Update the expected values for idps to match the actual values
  expected_idps <- matrix(c(
    0.43600666, -0.9894784, -0.90391048,
    0.01770737,  0.4254333, -0.14037602,
    -0.09228352,  1.1677400, -0.28021882
  ), nrow = 3, byrow = TRUE)

  # Test if the idps matches the updated expected values
  expect_equal(idps[1:3, 1:3], expected_idps, tolerance = 1e-7)

  # Manually hard-coded expected values for trait
  expected_trait <- c(0.2395, 0.2690, 1.2543)

  # Check if trait matches expected values
  expect_equal(trait[1:3], expected_trait, tolerance = 1e-3)

  # Manually hard-coded expected values for id_train
  expected_id_train <- c(1, 2, 3)

  # Check if id_train matches expected values
  expect_equal(id_train[1:3], expected_id_train)

  # Manually hard-coded expected values for age
  expected_age <- c(65.93970, 71.53632, 75.87305)

  # Check if age matches expected values
  expect_equal(age[1:3], expected_age, tolerance = 1e-7)

})

#library(testthat)
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
    0.2668068,  1.63889145, 0.03133779,
    -0.7986769,  0.08736424, 2.52321454,
    0.6055321, -1.76666590, 1.32474920
  ), nrow = 3, byrow = TRUE)

  expected_test_scaled <- matrix(c(
    -0.04830547,  1.2637210, -0.331411,
    -1.37831883, -0.5182675, -1.422242,
    1.22937218, -1.6635561, -1.802373
  ), nrow = 3, byrow = TRUE)

    # Test that the scaled training data matches the new expected values
  expect_equal(idps_train_scaled[1:3, 1:3], expected_train_scaled, tolerance = 1e-7)

  # Test that the scaled test data matches the new expected values
  expect_equal(idps_test_scaled[1:3, 1:3], expected_test_scaled, tolerance = 1e-7)

  # Check that the scaled means of the training data are approximately 0
  expect_true(all(abs(colMeans(idps_train_scaled)) < 1e-10))

  # Check that the scaled standard deviations of the training data are 1
  expect_true(all(abs(apply(idps_train_scaled, 2, sd) - 1) < 1e-10))
})




test_that("de_mean_trait_using_train correctly demeans the trait vector", {

  # Call the function to de-mean the trait using the training set
  demeaned_trait <- de_mean_trait_using_train(trait_randomized, id_train_randomized)

  # Manually computed expected values based on your results
  expected_trait_demeaned <- c(-61.5349406, -0.2990617, 22.4904990)

  # Test that the de-meaned trait matches the manually computed expected values
  expect_equal(demeaned_trait[1:3], expected_trait_demeaned, tolerance = 1e-7)

  # Check that the mean of the de-meaned training data is approximately 0
  expect_true(abs(mean(demeaned_trait[id_train_randomized])) < 1e-10)
})

# set up data for full pre_process_data_cross_validated check
prop_train <- 0.9 #90
ica <- 1
n_feat <- 5
trait_id <- 999
#list[df_all_train, idps, trait, df_all, id_train, age, conf] <- 
#pre_process_data_cross_validated(idps_randomized, trait_randomized, trait_id, age_randomized, conf_randomized, conf_names_randomized, prop_train, id_train_randomized, ica, n_feat)
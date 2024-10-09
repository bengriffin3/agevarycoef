library(testthat)
library(logger)
library(R.matlab)
# library(varycoef)
library(glmnet)
library(agevarycoef)

test_that("elastic net predictions run correctly", {

  message("Current working directory: ", getwd())

  set.seed(42)

  # Simulate data
  n_samples <- 10000  # Number of samples
  n_idps <- 10        # Number of IDPs

  # Simulated IDPs
  set.seed(42)  # For reproducibility
  idps_randomized <- matrix(rnorm(n_samples * n_idps), nrow = n_samples, ncol = n_idps)

  # Simulated trait (response variable) with correlation to the first 5 IDPs
  # Create a trait that is a linear combination of the first 5 IDPs
  trait_randomized <- rowSums(idps_randomized[, 1:5]) * 2 + rnorm(n_samples, mean = 0, sd = 5)  # Linear combination with noise

  # Simulated ages
  age_randomized <- sample(20:70, n_samples, replace = TRUE)  # Random ages between 20 and 70

  # Randomly assign training indices
  id_train_randomized <- sample(1:n_samples, size = round(0.7 * n_samples))  # 70% training data

  # Prepare the training data frame with the response variable
  df_train_randomized <- data.frame(
    y = trait_randomized  # response variable
  )

  # Add IDPs to the training data frame
  for (i in 1:n_idps) {
    df_train_randomized[[paste0("x.", i)]] <- idps_randomized[, i]
  }

  model_age <- 1

  # determine best features (linear)
  n_feat <- 5
  best_features_linear <- 1:n_feat

  # prepare linear model data
  linear_data <- prepare_linear_data(df_train_randomized, idps_randomized, best_features_linear)
  df_train_linear <- linear_data$df_train_linear
  idps_linear <- linear_data$idps_linear

  ################ RUN ELASTIC NET ################

  elastic_net_results <- run_elastic_net_model(idps_linear, trait_randomized, id_train_randomized, age_randomized, model_age)
  yhat_train <- elastic_net_results$yhat_train
  yhat_test <- elastic_net_results$yhat_test

  model_results <- get_model_stats(yhat_train, yhat_test, trait_randomized[id_train_randomized], trait_randomized[-id_train_randomized])
  corr_train <- model_results$corr_train
  corr_test <- model_results$corr_test

  expected_corr_train <- 0.66772
  expected_corr_test <- 0.64483

  # Actual values (assume these are calculated in your function)
  actual_corr_train <- corr_train
  actual_corr_test <- corr_test

  # Check that the actual in-sample correlation matches the expected value
  expect_equal(round(actual_corr_train, 5), round(expected_corr_train, 5), tolerance = 1e-7)

  # Check that the actual out-of-sample correlation matches the expected value
  expect_equal(round(actual_corr_test, 5), round(expected_corr_test, 5), tolerance = 1e-7)

})

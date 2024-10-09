library(testthat)
library(logger)
library(R.matlab)
# library(varycoef)
library(glmnet)
library(agevarycoef)
library(here)

test_that("elastic net predictions run correctly", {

  message("Current working directory: ", getwd())
  source(here::here("R", "data_preparation.r"))
  source(here::here("R", "prediction.r"))
  # load(here::here("data", "test_data.rda"))
  # load(here::here("data", "test_data_pred.rda"))

  set.seed(42)

  # Simulate data
  n_samples <- 10000  # Number of samples
  n_idps <- 10      # Number of IDPs

  # Simulated IDPs (e.g., 100 samples with 10 IDPs)
  idps_randomized <- matrix(rnorm(n_samples * n_idps), nrow = n_samples, ncol = n_idps)

  # Simulated trait (response variable)
  trait_randomized <- rnorm(n_samples)

  # Simulated ages
  age_randomized <- sample(20:70, n_samples, replace = TRUE)  # Random ages between 20 and 70

  # Randomly assign training indices
  id_train_randomized <- sample(1:n_samples, size = round(0.7 * n_samples))  # 70% training data

  df_train_randomized <- data.frame(
    y = trait_randomized + rnorm(n_samples, mean = 0, sd = 10)  # response variable with noise
  )
  
  # Add 1435 random predictors
  for (i in 1:n_idps) {
    df_train_randomized[[paste0("x.", i)]] <- rnorm(n_samples)
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

  expected_corr_train <- 0.13833
  expected_corr_test <- 0.11169

  # Actual values (assume these are calculated in your function)
  actual_corr_train <- corr_train
  actual_corr_test <- corr_test

  # Check that the actual in-sample correlation matches the expected value
  expect_equal(round(actual_corr_train, 5), round(expected_corr_train, 5), tolerance = 1e-7)

  # Check that the actual out-of-sample correlation matches the expected value
  expect_equal(round(actual_corr_test, 5), round(expected_corr_test, 5), tolerance = 1e-7)

})

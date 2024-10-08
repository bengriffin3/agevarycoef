library(testthat)
library(logger)
library(R.matlab)
library(varycoef)
library(glmnet)
library(agevarycoef)

library(here)

source(here("R", "prediction.r"))
source(here("R", "data_preparation.r"))

# load('/gpfs3/well/win-fmrib-analysis/users/psz102/git_repos/agevarycoef/data/test_data.rda')
# load('/gpfs3/well/win-fmrib-analysis/users/psz102/git_repos/agevarycoef/data/test_data_pred.rda')

load('data/test_data.rda')
load('test_data_pred.rda')

n_feat <- 5
model_age <- 1

test_that("elastic net predictions run correctly", {
  set.seed(42)
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

  expected_corr_train <- 0.1383271
  expected_corr_test <- 0.1116852

  # Actual values (assume these are calculated in your function)
  actual_corr_train <- corr_train
  actual_corr_test <- corr_test

  # Check that the actual in-sample correlation matches the expected value
  expect_equal(round(actual_corr_train, 5), round(expected_corr_train, 5), tolerance = 1e-7)

  # Check that the actual out-of-sample correlation matches the expected value
  expect_equal(round(actual_corr_test, 5), round(expected_corr_test, 5), tolerance = 1e-7)

})

library(testthat)
library(logger)
library(varycoef)
library(glmnet)

#setwd("/gpfs3/well/win-fmrib-analysis/users/psz102/git_repos/agevarycoef/")
data(test_data_pred)

n_feat <- 5
model_age <- 1

test_that("linear predictions run correctly", {
  set.seed(42)
  # determine best features (linear)
  best_features_linear <- 1:n_feat

  # prepare linear model data
  #list[df_train_linear, idps_linear] <- prepare_linear_data(df_train_randomized, idps_randomized, best_features_linear)
  linear_data <- prepare_linear_data(df_train_randomized, idps_randomized, best_features_linear)
  df_train_linear <- linear_data$df_train_linear
  idps_linear <- linear_data$idps_linear



  # prepare age data if we are investigating different age options
  #list[df_train_linear, df] <- prepare_age_data(df_train_linear, df_randomized, age_train_randomized, model_age, age_randomized)
  age_data <- prepare_age_data(df_train_linear, df_randomized, age_train_randomized, model_age, age_randomized)
  df_train_linear <- age_data$df_train_linear
  df <- age_data$df


  ################ RUN LINEAR MODEL ################
  #list[fit_lm, se_lm, corr_lm_in, corr_lm_out, lm_yhat] <- run_linear_model(df_train_linear, id_train, df)
  linear_results <- run_linear_model(df_train_linear, id_train_randomized, df)
  fit_lm <- linear_results$fit_lm
  se_lm <- linear_results$se_lm
  corr_lm_in <- linear_results$corr_lm_in
  corr_lm_out <- linear_results$corr_lm_out
  lm_yhat <- linear_results$lm_yhat



  # Expected values (these are the values you've provided)
  expected_se_lm <- c(1086.25718, 304.61554, 29.18595)
  expected_corr_lm_in <- 0.1492521
  expected_corr_lm_out <- 0.1751417
  expected_lm_yhat <- c(-2.8320490, 0.8914371, 2.0784881)

  # Actual values (assume these are calculated in your function)
  actual_se_lm <- se_lm[1:3]
  actual_corr_lm_in <- corr_lm_in
  actual_corr_lm_out <- corr_lm_out
  actual_lm_yhat <- lm_yhat[1:3]

  # Check that the actual values match the expected values
  expect_equal(actual_se_lm, expected_se_lm, tolerance = 1e-7)
  expect_equal(round(actual_corr_lm_in, 5), round(expected_corr_lm_in, 5), tolerance = 1e-7)
  expect_equal(round(actual_corr_lm_out, 5), round(expected_corr_lm_out, 5), tolerance = 1e-7)
  expect_equal(actual_lm_yhat, expected_lm_yhat, tolerance = 1e-7)

})

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
  #list[enet_coefficients, cvfit_glmnet, enet_yhat, se_enet, corr_enet_in, corr_enet_out] <- run_elastic_net_model(idps_randomized, idps_linear, trait_randomized, id_train_randomized, age_randomized, model_age)
  #list[enet_coefficients, cvfit_glmnet, enet_yhat, se_enet, corr_enet_in, corr_enet_out] <- run_elastic_net_model(idps_randomized, idps_linear, trait_randomized, id_train_randomized, age_randomized, model_age)
  elastic_net_results <- run_elastic_net_model(idps_randomized, idps_linear, trait_randomized, id_train_randomized, age_randomized, model_age)
  enet_coefficients <- elastic_net_results$enet_coefficients
  cvfit_glmnet <- elastic_net_results$cvfit_glmnet
  enet_yhat <- elastic_net_results$enet_yhat
  se_enet <- elastic_net_results$se_enet
  corr_enet_in <- elastic_net_results$corr_enet_in
  corr_enet_out <- elastic_net_results$corr_enet_out


  expected_enet_yhat <- c(-3.2944718, -0.2900734, -0.1323354)
  expected_se_enet <- c(1116.95246, 347.25387, 57.96121)
  #expected_se_enet <- c(6.41061408, 0.02410067, 0.30670160)
  expected_corr_enet_in <- 0.1383271
  expected_corr_enet_out <- 0.1116852

  # Actual values (assume these are calculated in your function)
  actual_enet_yhat <- enet_yhat[1:3]
  actual_se_enet <- se_enet[1:3]
  actual_corr_enet_in <- corr_enet_in
  actual_corr_enet_out <- corr_enet_out


  # Check that the actual predicted values (yhat) match the expected values
  expect_equal(actual_enet_yhat, expected_enet_yhat, tolerance = 1e-7)

  # Check that the actual standard errors match the expected values
  expect_equal(actual_se_enet, expected_se_enet, tolerance = 1e-7)

  # Check that the actual in-sample correlation matches the expected value
  expect_equal(round(actual_corr_enet_in, 5), round(expected_corr_enet_in, 5), tolerance = 1e-7)

  # Check that the actual out-of-sample correlation matches the expected value
  expect_equal(round(actual_corr_enet_out, 5), round(expected_corr_enet_out, 5), tolerance = 1e-7)

})


test_that("SVC predictions run correctly", {

  set.seed(42)
  prof <- FALSE
  #perc_train <- 90
  tap <- 0
  taper <- 0
  cov <- "exp"
  model_age <- 1
  #ica <- 0
  #n_feat <- 5
  #trait_id <- 999

  # determine best svc features
  #best_features_svc <- determine_best_svc_features(trait_id, n_feat, n_sub, perc_train, prof, tap, cov, ica)
  best_features_svc <- c(1, 5, 6, 9, 10)

  # prepare svc data
  #list[df_train_svc, idps_svc, idps_svc_train] <- prepare_svc_data(df_train_randomized, best_features_svc, id_train_randomized, idps_randomized, age_randomized)
  svc_data <- prepare_svc_data(df_train_randomized, best_features_svc, id_train_randomized, idps_randomized, age_randomized, df_all_randomized)
  df_train_svc <- svc_data$df_train_svc
  idps_svc <- svc_data$idps_svc
  idps_svc_train <- svc_data$idps_svc_train
  df_svc <- svc_data$df_svc


  # fit svc model and use for prediction
  #list[df_svc_pred, se_svc, corr_svc_in, corr_svc_out] <- fit_svc_model(best_features_svc, df_all_train_randomized, df_all_randomized, df_train_svc, cov, prof, taper, id_train_randomized, age_randomized, age_train_randomized, model_age)
  svc_results <- fit_svc_model(best_features_svc, df_all_train_randomized, df_all_randomized, df_train_svc, cov, prof, taper, id_train_randomized, age_randomized, age_train_randomized, model_age, df_svc)
  # unpack results
  se_svc <- svc_results$se_svc
  corr_svc_in <- svc_results$corr_svc_in
  corr_svc_out <- svc_results$corr_svc_out

  # Actual values (as calculated by your function)
  actual_se_svc <- se_svc[1:3]
  actual_corr_svc_in <- corr_svc_in
  actual_corr_svc_out <- corr_svc_out

  # Expected values
  #expected_se_svc <- c(196.787282, 64.133961, 2.294036)
  expected_se_svc <- c(209.3331375, 319.5240798, 0.1692568)
  expected_corr_svc_in <- -0.07991901 # -0.01855146
  expected_corr_svc_out <- 0.1345724 #-0.1258423

  # Check that the actual values match the expected values
  expect_equal(actual_se_svc, expected_se_svc, tolerance = 1e-7)
  expect_equal(actual_corr_svc_in, expected_corr_svc_in, tolerance = 1e-7)
  expect_equal(round(actual_corr_svc_out, 5), round(expected_corr_svc_out, 5), tolerance = 1e-7)

})
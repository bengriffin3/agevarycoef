library(testthat)
library(logger)
library(varycoef)
library(glmnet)

#setwd("/gpfs3/well/win-fmrib-analysis/users/psz102/git_repos/agevarycoef/")
data(test_data_pred)

n_feat <- 5
model_age <- 1

# test_that("linear predictions run correctly", {
#   set.seed(42)
#   # determine best features (linear)
#   best_features_linear <- 1:n_feat

#   # prepare linear model data
#   #list[df_train_linear, idps_linear] <- prepare_linear_data(df_train_randomized, idps_randomized, best_features_linear)
#   linear_data <- prepare_linear_data(df_train_randomized, idps_randomized, best_features_linear)
#   df_train_linear <- linear_data$df_train_linear
#   idps_linear <- linear_data$idps_linear



#   # prepare age data if we are investigating different age options
#   #list[df_train_linear, df] <- prepare_age_data(df_train_linear, df_randomized, age_train_randomized, model_age, age_randomized)
#   age_data <- prepare_age_data(df_train_linear, df_randomized, age_train_randomized, model_age, age_randomized)
#   df_train_linear <- age_data$df_train_linear
#   df <- age_data$df


#   ################ RUN LINEAR MODEL ################
#   linear_model <- run_linear_model(df_train_linear, id_train_randomized, df)
#   lm_yhat <- linear_model$lm_yhat
#   yhat_train <- lm_yhat[id_train_randomized]
#   yhat_test <- lm_yhat[-id_train_randomized]

#   model_results <- get_model_stats(yhat_train, yhat_test, trait_randomized[id_train_randomized], trait_randomized[-id_train_randomized])
#   corr_train <- model_results$corr_train
#   corr_test <- model_results$corr_test

#   # Expected values (these are the values you've provided)
#   expected_corr_train <- 0.1492521
#   expected_corr_test <- 0.1751417
#   expected_lm_yhat <- c(-2.8320490, 0.8914371, 2.0784881)

#   # Actual values (assume these are calculated in your function)
#   actual_corr_train <- corr_train
#   actual_corr_test <- corr_test
#   actual_lm_yhat <- lm_yhat[1:3]

#   # Check that the actual values match the expected values
#   expect_equal(round(actual_corr_train, 5), round(expected_corr_train, 5), tolerance = 1e-7)
#   expect_equal(round(actual_corr_test, 5), round(expected_corr_test, 5), tolerance = 1e-7)
#   expect_equal(actual_lm_yhat, expected_lm_yhat, tolerance = 1e-7)

# })

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
  elastic_net_results <- run_elastic_net_model(idps_randomized, idps_linear, trait_randomized, id_train_randomized, age_randomized, model_age)  
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


# test_that("SVC predictions run correctly", {

#   set.seed(42)
#   prof <- FALSE
#   #perc_train <- 90
#   tap <- 0
#   taper <- 0
#   cov <- "exp"
#   model_age <- 1
#   #ica <- 0
#   #n_feat <- 5
#   #trait_id <- 999

#   # determine best svc features
#   #best_features_svc <- determine_best_svc_features(trait_id, n_feat, n_sub, perc_train, prof, tap, cov, ica)
#   best_features_svc <- c(1, 5, 6, 9, 10)

#   # prepare svc data
#   #list[df_train_svc, idps_svc, idps_svc_train] <- prepare_svc_data(df_train_randomized, best_features_svc, id_train_randomized, idps_randomized, age_randomized)
#   svc_data <- prepare_svc_data(df_train_randomized, best_features_svc, id_train_randomized, idps_randomized, age_randomized, df_all_randomized)
#   df_train_svc <- svc_data$df_train_svc
#   idps_svc <- svc_data$idps_svc
#   idps_svc_train <- svc_data$idps_svc_train
#   df_svc <- svc_data$df_svc


#   # fit svc model and use for prediction
#   #list[df_svc_pred, se_svc, corr_svc_in, corr_svc_out] <- fit_svc_model(best_features_svc, df_all_train_randomized, df_all_randomized, df_train_svc, cov, prof, taper, id_train_randomized, age_randomized, age_train_randomized, model_age)
#   svc_results <- fit_svc_model(best_features_svc, df_all_train_randomized, df_all_randomized, df_train_svc, cov, prof, taper, id_train_randomized, age_randomized, age_train_randomized, model_age, df_svc)
#   # unpack results
#   yhat_train <- svc_results$yhat_train
#   yhat_test <- svc_results$yhat_test

#   model_results <- get_model_stats(yhat_train, yhat_test, trait_randomized[id_train_randomized], trait_randomized[-id_train_randomized])
#   corr_train <- model_results$corr_train
#   corr_test <- model_results$corr_test

#   # Actual values (as calculated by your function)
#   actual_corr_train <- corr_train
#   actual_corr_test <- corr_test

#   # Expected values
#   expected_corr_train <- 0.97172365 # -0.07991901 # -0.01855146
#   expected_corr_test <- -0.18884077 # 0.1345724 #-0.1258423

#   # Check that the actual values match the expected values
#   expect_equal(actual_corr_train, expected_corr_train, tolerance = 1e-7)
#   expect_equal(round(actual_corr_test, 5), round(expected_corr_test, 5), tolerance = 1e-7)

# })
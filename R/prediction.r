library(R.matlab)
library(varycoef)
library(optparse)
library(fastICA)
library(CCA)
library(logger)
library(glmnet)
library(lightgbm)
library(xgboost)
library(data.table)
library(caret)


run_linear_model <- function(df, id_train) {

  df_train <- df[id_train, ]
  df_test <- df[-id_train, ]

  print("Fitting linear model to test data...")
  fit_lm <- lm(y ~ ., data = df_train)

  # ridge regression (it doesn't really make sense to constrain the sum of the single coefficient to equal something)
  # so we don't do this for univariate predictions (could be worth for multivariate but just use glmnet)
  #fit_ridge <-  lm.ridge(y ~ ., data = df_train_linear, lambda = seq(0, .4, 1e-3))
  lm_yhat <- numeric(length(df$y))
  lm_yhat[id_train] <- predict(fit_lm, newdata = df_train)
  lm_yhat[-id_train] <- predict(fit_lm, newdata = df_test)

  return(list(fit_lm = fit_lm, lm_yhat = lm_yhat))
}


run_linear_model_cv <- function(df) {
  # We manually do the cross-validation so we can look at the predictions and compare
  # training / test set accuracies
  print("Fitting linear model with 10-fold cross-validation...")

  # Initialize vectors to store correlations for each fold
  corr_train <- numeric(10)
  corr_test <- numeric(10)
  mse_train <- numeric(10)
  mse_test <- numeric(10)
  r_squared_train <- numeric(10)
  r_squared_test <- numeric(10)


  # Initialize a vector to store predictions
  lm_yhat <- numeric(length(df$y))

  # Perform 10-fold cross-validation manually to calculate correlations
  folds <- createFolds(df$y, k = 10)

  #models <- list()
  models <- vector("list", 10)  # 10 folds, so 10 models

  for (i in seq_along(folds)) {
    print(paste("Fold", i))
    test_idx <- folds[[i]]
    train_idx <- setdiff(seq_along(df$y), test_idx)

    # Train the model on the training fold
    fit_lm <- lm(y ~ ., data = df[train_idx, ])
    #models <- c(models, fit_lm)
    models[[i]] <- fit_lm

    # Predictions for the training and test set
    lm_yhat_train <- predict(fit_lm, newdata = df[train_idx, ])
    lm_yhat_test <- predict(fit_lm, newdata = df[test_idx, ])

    # Note accuracies
    corr_train[i] <- cor(df$y[train_idx], lm_yhat_train)
    corr_test[i] <- cor(df$y[test_idx], lm_yhat_test)
    mse_train[i] <- mean((df$y[train_idx] - lm_yhat_train) ^ 2)
    mse_test[i] <- mean((df$y[test_idx] - lm_yhat_test) ^ 2)
    r_squared_train[i] <- 1 - (sum((df$y[train_idx] - lm_yhat_train) ^ 2) / sum((df$y[train_idx] - mean(df$y[train_idx])) ^ 2))
    r_squared_test[i] <- 1 - (sum((df$y[test_idx] - lm_yhat_test) ^ 2) / sum((df$y[test_idx] - mean(df$y[test_idx])) ^ 2))

    # Store predictions in the original index positions
    #lm_yhat[train_idx] <- lm_yhat_train
    lm_yhat[test_idx] <- lm_yhat_test
  }

  # Return the fitted model, predictions, and accuracies
  return(list(
    lm_yhat = lm_yhat,
    corr_train = corr_train,
    corr_test = corr_test,
    mse_train = mse_train,
    mse_test = mse_test,
    r_squared_train = r_squared_train,
    r_squared_test = r_squared_test,
    models = models
  ))
}

# Function to calculate performance metrics
calculate_metrics <- function(actual_train, predicted_train, actual_test, predicted_test) {
  # Initialize results list
  metrics <- list()
  
  # Calculate correlations
  metrics$corr_train <- cor(actual_train, predicted_train)
  metrics$corr_test <- cor(actual_test, predicted_test)
  
  # Calculate Mean Squared Error
  metrics$mse_train <- mean((actual_train - predicted_train) ^ 2)
  metrics$mse_test <- mean((actual_test - predicted_test) ^ 2)
  
  # Calculate R-squared
  metrics$r_squared_train <- 1 - (sum((actual_train - predicted_train) ^ 2) / sum((actual_train - mean(actual_train)) ^ 2))
  metrics$r_squared_test <- 1 - (sum((actual_test - predicted_test) ^ 2) / sum((actual_test - mean(actual_test)) ^ 2))
  
  return(metrics)
}


run_elastic_net_model <- function(idps, idps_linear, trait, id_train, age, model_age) {

  print("Fitting elastic net model...")
  trait_train <- trait[id_train]
  idps_linear <- idps_linear
  idps_linear_train <- idps_linear[id_train, ]
  age_train <- matrix(age[id_train])

  if (model_age == 1 ||  model_age == 3) {
    cvfit_glmnet <- cv.glmnet(cbind(idps_linear_train, age_train), trait_train)
  } else if (model_age == 0) {
    cvfit_glmnet <- cv.glmnet(idps_linear_train, trait_train)
  } else if (model_age == 2) {
    cvfit_glmnet <- cv.glmnet(cbind(rep(1, length(age_train)), age_train), trait_train)
  }

  # enet_coefficients <- coef(cvfit_glmnet, s = "lambda.min") # note betas
  # print(enet_coefficients)
  # print(enet_coefficients != 0)
  #idx_non_zero_coeff <- determine_non_zero_coeff(enet_coefficients) # get best features
  # idx_non_zero_coeff <- which(enet_coefficients != 0)
  # print(paste0("Number of non-zero features: ", length(idx_non_zero_coeff)))

  enet_yhat <- numeric(length(trait))
  if (model_age == 1 ||  model_age == 3) {
    enet_yhat[id_train] <- predict(cvfit_glmnet, newx = cbind(idps_linear_train, age_train), s = "lambda.min")
    enet_yhat[-id_train] <- predict(cvfit_glmnet, newx = cbind(idps_linear[-id_train, ], age[-id_train]), s = "lambda.min")
  } else if (model_age == 0) {
    enet_yhat[id_train] <- predict(cvfit_glmnet, newx = idps_linear_train, s = "lambda.min")
    enet_yhat[-id_train] <- predict(cvfit_glmnet, newx = idps_linear[-id_train, ], s = "lambda.min")
  } else if (model_age == 2) {
    enet_yhat[id_train] <- predict(cvfit_glmnet, newx = cbind(rep(1, length(age_train)), age_train), s = "lambda.min")
    enet_yhat[-id_train] <- predict(cvfit_glmnet, newx = cbind(rep(1, length(age[-id_train])), age[-id_train]) , s = "lambda.min")
  }

  yhat_train <- enet_yhat[id_train]
  yhat_test <- enet_yhat[-id_train]

  return(list(cvfit_glmnet = cvfit_glmnet, yhat_train = yhat_train, yhat_test = yhat_test))

}


run_svc_model <- function(best_features, df, cov, prof, taper, age, model_age, id_train) {

  svc_config <- configure_svc_model(cov, prof, taper)

  # create data frame for SVC
  df_svc <- df[, c(1, best_features + 1)]
  df_svc$age <- age
  fit_svc <- fit_svc_model(best_features, df, df_svc, id_train, age, model_age, svc_config) 
  
  # Predict using the trained model
  df_svc_pred <- predict(fit_svc, newdata = df_svc, newlocs = age, control = svc_config)
  yhat_train <- df_svc_pred$y.pred[id_train]
  yhat_test <- df_svc_pred$y.pred[-id_train]

  return(list(
    fit_svc = fit_svc,
    yhat_train = yhat_train,
    yhat_test = yhat_test,
    df_svc_pred = df_svc_pred
  ))
}

configure_svc_model <- function(cov, prof, taper) {
  # Configure SVC
  svc_config <- SVC_mle_control(
    cov.name = c(cov),
    profileLik = prof,
    tapering = taper,
    hessian = TRUE
  )

  return(svc_config)
}

fit_svc_model <- function(best_features, df, df_svc, id_train, age, model_age, svc_config) {

  print("Fitting SVC")

  # Generate feature names
  xnam <- paste("x.", best_features, sep = "")
  fmla <- as.formula(paste("y ~ ", paste(xnam, collapse = "+")))

  # get training data
  df_train_svc <- df_svc[id_train, ]

  # Fit the SVC model based on model_age
  # Prepare formula for SVC model
  if (model_age == 0) {
    # Model without intercept
    fmla_mi <- as.formula(paste("y ~ ", paste(xnam, collapse = "+"), " - 1"))
    fit_svc <- SVC_mle(
      formula = fmla, 
      data = df_train_svc, 
      locs = df_train_svc$age, 
      control = svc_config, 
      RE_formula = fmla_mi
    )
  } else if (model_age == 1) {
    # Model with age as a covariate
    fit_svc <- SVC_mle(
      formula = fmla, 
      data = df_train_svc, 
      locs = df_train_svc$age, 
      control = svc_config
    )
  } else if (model_age == 2) {
    # Model with only age
    fit_svc <- SVC_mle(
      formula = as.formula("y ~ age"), 
      data = df_train_svc, 
      locs = df_train_svc$age, 
      control = svc_config
    )
  } else if (model_age == 3) {
    # Model with additional features
    load(sprintf("/gpfs3/well/win-fmrib-analysis/users/psz102/age_varying_coefficients/results/univariate/ranking_idps_trait_id_%i_run_svc_0_nsub_46471_ra_0.RData", trait_id))
    best_features_linear_add <- idp_ranking_trait[1:length(best_features)]

    xnam <- paste("x.", c(best_features, best_features_linear_add), sep = "")
    fmla_fix <- as.formula(paste("y ~ ", paste(xnam, collapse = "+")))
    
    # Add additional features to training and test data
    extracted_columns <- df[id_train, paste0("x.", best_features_linear_add)]
    df_train_svc <- cbind(df_train_svc, extracted_columns)
    
    extracted_columns_df <- df[, paste0("x.", best_features_linear_add)]
    df_svc <- cbind(df, extracted_columns_df)
    
    fit_svc <- SVC_mle(
      formula = fmla_fix, 
      data = df_train_svc, 
      locs = df_train_svc$age, 
      control = svc_config, 
      RE_formula = fmla
    )
  }

    return(fit_svc)
}



run_svc_model_cv <- function(best_features, df, cov, prof, taper, age, model_age, n_folds = 10) {


  # Create cross-validation folds
  folds <- createFolds(df$y, k = n_folds)
  
  # Initialize vectors to store results
  yhat_test_all <- numeric(nrow(df))
  corr_train <- numeric(n_folds)
  corr_test <- numeric(n_folds)
  mse_train <- numeric(n_folds)
  mse_test <- numeric(n_folds)
  r_squared_train <- numeric(n_folds)
  r_squared_test <- numeric(n_folds)


   svc_config <- configure_svc_model(cov, prof, taper)

  # create data frame for SVC
  df_svc <- df[, c(1, best_features + 1)]
  df_svc$age <- age

  models <- vector("list", 10)  # 10 folds, so 10 models

  # Perform cross-validation
  for (i in seq_along(folds)) {
    print(paste("Processing fold", i))
    
    # Define train and test indices
    test_idx <- folds[[i]]
    train_idx <- setdiff(seq_len(nrow(df)), test_idx)
    
    # Split data for this fold
    df_train <- df[train_idx, ]
    df_test <- df[test_idx, ]

    fit_svc <- fit_svc_model(best_features, df, df_svc, train_idx, age, model_age, svc_config)
    models[[i]] <- fit_svc
  
    # Predict using the trained model
    df_svc_pred <- predict(fit_svc, newdata = df_svc, newlocs = age, control = svc_config)
    yhat_train <- df_svc_pred$y.pred[train_idx]
    yhat_test <- df_svc_pred$y.pred[test_idx]
    
    # Store predictions
    yhat_test_all[test_idx] <- yhat_test
    
    # Calculate and store metrics
    corr_train[i] <- cor(df$y[train_idx], yhat_train)
    corr_test[i] <- cor(df$y[test_idx], yhat_test)
    mse_train[i] <- mean((df$y[train_idx] - yhat_train) ^ 2)
    mse_test[i] <- mean((df$y[test_idx] - yhat_test) ^ 2)
    r_squared_train[i] <- 1 - (sum((df$y[train_idx] - yhat_train) ^ 2) / sum((df$y[train_idx] - mean(df$y[train_idx])) ^ 2))
    r_squared_test[i] <- 1 - (sum((df$y[test_idx] - yhat_test) ^ 2) / sum((df$y[test_idx] - mean(df$y[test_idx])) ^ 2))
  }
  
  # Return results
  return(list(
    yhat_test_all = yhat_test_all,
    corr_train = corr_train,
    corr_test = corr_test,
    mse_train = mse_train,
    mse_test = mse_test,
    r_squared_train = r_squared_train,
    r_squared_test = r_squared_test,
    models=models
  ))
}

determine_non_zero_coeff <- function(enet_coefficients) {

  idx_best_features <- which(enet_coefficients != 0)
  # to save index: save(idx_best_features, file = "idx_best_features.RData")
  return(idx_best_features)
}

penalized_svc <- function(idps, trait, age, configl) {
  # Penalized MLE for SVC (i.e., L1 feature selection)
  # first alter config: SVC_control_config$extract_fun <- TRUE
  # run SVC_mle(y = trait, X = idps, W = idps, locs = age, control = config)
}

extract_best_features <- function(corr, n) {

  best_features <- sort(corr, index.return = TRUE, decreasing = TRUE)
  best_features_idx <- best_features$ix
  best_n_features <- sort(best_features_idx[1:n])

  return(best_n_features)
}


run_lgboost_model <- function(df_all_train_x, df_all_train_y, df_all_test_x, df_all_test_y, params = 0) {

  train_lgb <- lgb.Dataset(data = df_all_train_x, label = df_all_train_y)
  test_lgb <- lgb.Dataset(data = df_all_test_x, label = df_all_test_y)
  if (any(params)) {
    params <- params
  } else {
  params <- list(
    objective = "regression",
    metric = "l2",
    learning_rate = 0.01,
    num_leaves = 15,
    min_data_in_leaf = 10,
    max_depth = 10,
    lambda_l1 = 0.5,           # L1 regularization
    lambda_l2 = 0.5,            # L2 regularization
    bagging_fraction = 0.8,
    bagging_freq = 1
  )
}
  print("Running LGBoost")
  # Train the model
  lgb_model <- lgb.train(
                         params = params,
                         data = train_lgb,
                         nrounds = 1000,
                         early_stopping_rounds = 10,
                         valids = list(train = train_lgb, test = test_lgb),
                         verbose = 1)


  # Predict using the trained model
  predictions_test <- predict(lgb_model, df_all_test_x)
  predictions_train <- predict(lgb_model, df_all_train_x)


 return(list(predictions_train, predictions_test, lgb_model))
 }


 run_lgboost_model_cv <- function(df_boost, params = NULL, n_folds = 10) {
  print("Running LGBoost with cross-validation...")
  
  # Prepare vectors to store metrics for each fold
  corr_train <- numeric(n_folds)
  corr_test <- numeric(n_folds)
  mse_train <- numeric(n_folds)
  mse_test <- numeric(n_folds)
  r_squared_train <- numeric(n_folds)
  r_squared_test <- numeric(n_folds)
  
  # Initialize a matrix to store predictions
  predictions_test_all <- numeric(nrow(df_boost))
  
  # Generate fold indices
  folds <- createFolds(df_boost$y, k = n_folds)
  
  # Default LightGBM parameters if none are provided
  if (is.null(params)) {
    params <- list(
      objective = "regression",
      metric = "l2",
      learning_rate = 0.01,
      num_leaves = 15,
      min_data_in_leaf = 10,
      max_depth = 10,
      lambda_l1 = 0.5,
      lambda_l2 = 0.5,
      bagging_fraction = 0.8,
      bagging_freq = 1
    )
  }

  for (i in seq_along(folds)) {
    print(paste("Fold", i))

    test_idx <- folds[[i]]
    train_idx <- setdiff(seq_len(nrow(df_boost)), test_idx)

    df_train_x <- as.matrix(df_boost[train_idx, ][, !names(df_boost[train_idx, ]) %in% c("y")])
    df_train_y <- as.matrix(df_boost[train_idx, ]$y)

    df_test_x <- as.matrix(df_boost[test_idx, ][, !names(df_boost[test_idx, ]) %in% c("y")])
    df_test_y <- as.matrix(df_boost[test_idx, ]$y)

    train_lgb <- lgb.Dataset(data = df_train_x, label = df_train_y)
    test_lgb <- lgb.Dataset(data = df_test_x, label = df_test_y)

    # Train the LightGBM model
    lgb_model <- lgb.train(
      params = params,
      data = train_lgb,
      nrounds = 1000,
      early_stopping_rounds = 100,
      valids = list(train = train_lgb, test = test_lgb),
      verbose = 1
    )
    
    # Predict using the trained model
    predictions_train <- predict(lgb_model, df_train_x)
    predictions_test <- predict(lgb_model, df_test_x)
    
    # Store predictions
    predictions_test_all[test_idx] <- predictions_test
    
    # Note accuracies
    corr_train[i] <- cor(df_boost$y[train_idx], predictions_train)
    corr_test[i] <- cor(df_boost$y[test_idx], predictions_test)
    mse_train[i] <- mean((df_boost$y[train_idx] - predictions_train) ^ 2)
    mse_test[i] <- mean((df_boost$y[test_idx] - predictions_test) ^ 2)
    r_squared_train[i] <- 1 - (sum((df_boost$y[train_idx] - predictions_train) ^ 2) / sum((df_boost$y[train_idx] - mean(df_boost$y[train_idx])) ^ 2))
    r_squared_test[i] <- 1 - (sum((df_boost$y[test_idx] - predictions_test) ^ 2) / sum((df_boost$y[test_idx] - mean(df_boost$y[test_idx])) ^ 2))
  }

  # Return results
  return(list(
    predictions_test = predictions_test_all,
    corr_train = corr_train,
    corr_test = corr_test,
    mse_train = mse_train,
    mse_test = mse_test,
    r_squared_train = r_squared_train,
    r_squared_test = r_squared_test
  ))
}

get_model_stats <- function(predictions_train, predictions_test, df_all_train_y, df_all_test_y) {

  # Calculate Mean Squared Error
  mse_test <- mean((predictions_test - df_all_test_y) ^ 2)
  mse_train <- mean((predictions_train - df_all_train_y) ^ 2)

  # Calculate R-squared
  r_squared_test <- 1 - (sum((df_all_test_y - predictions_test) ^ 2) / sum((df_all_test_y - mean(df_all_test_y)) ^ 2))
  r_squared_train <- 1 - (sum((df_all_train_y - predictions_train) ^ 2) / sum((df_all_train_y - mean(df_all_train_y)) ^ 2))

  # Calculate correlation
  corr_test <- cor(predictions_test, df_all_test_y)
  corr_train <- cor(predictions_train, df_all_train_y)

  print(paste("Mean Squared Error (TRAIN):", mse_train))
  print(paste("Mean Squared Error (TEST):", mse_test))
  print(paste("R-squared (TRAIN):", r_squared_train))
  print(paste("R-squared (TEST):", r_squared_test))
  print(paste("Correlation (TRAIN):", corr_train))
  print(paste("Correlation (TEST):", corr_test))

return(list(mse_test = mse_test, mse_train = mse_train, r_squared_test = r_squared_test, r_squared_train = r_squared_train, corr_test = corr_test, corr_train = corr_train))}

run_xgboost_model <- function(df_all_train_x, df_all_train_y, df_all_test_x, df_all_test_y) {

  dtrain <- xgb.DMatrix(data = df_all_train_x, label = df_all_train_y)
  dtest <- xgb.DMatrix(data = df_all_test_x, label = df_all_test_y)

  # Set XGBoost parameters
  params <- list(
    objective = "reg:squarederror",  # Regression objective
    eval_metric = "rmse",            # Evaluation metric
    #booster = "gblinear",
    max_depth = 4,                   # Maximum depth of a tree
    eta = 0.1,                       # Learning rate
    subsample = 0.8,                 # Subsample ratio of the training instances
    colsample_bytree = 0.8           # Subsample ratio of columns when constructing each tree
  )

  xgb_model <- xgb.train(
    params = params,
    data = dtrain,
    nrounds = 40,                    # Number of boosting rounds # 100
    #watchlist = list(train = dtrain), # Watchlist to see performance on training data
    watchlist = list(train = dtrain, test = dtest),
    verbose = 1                       # Print out training progress
  )

  predictions_test <- predict(xgb_model, newdata = df_all_test_x)
  predictions_train <- predict(xgb_model, newdata = df_all_train_x)

  return(list(predictions_train, predictions_test, xgb_model))
}

get_importance_matrix <- function(model, run_svc) {

  if (run_svc == 2) {
    importance_matrix <- lgb.importance(model = model, percentage = TRUE)
  } else if (run_svc == 3) {
    importance_matrix <- xgb.importance(model = model)
  }

  return(importance_matrix)
}

get_mean_accuracies_across_cv <- function(corr_train_cv, corr_test_cv, mse_train_cv, mse_test_cv, r_squared_train_cv, r_squared_test_cv) {
    corr_train <- mean(corr_train_cv)
    corr_test <- mean(corr_test_cv)
    mse_train <- mean(mse_train_cv)
    mse_test <- mean(mse_test_cv)
    r_squared_train <- mean(r_squared_train_cv)
    r_squared_test <- mean(r_squared_test_cv)

    print(paste("Mean Squared Error (TRAIN):", mse_train))
    print(paste("Mean Squared Error (TEST):", mse_test))
    print(paste("R-squared (TRAIN):", r_squared_train))
    print(paste("R-squared (TEST):", r_squared_test))
    print(paste("Correlation (TRAIN):", corr_train))
    print(paste("Correlation (TEST):", corr_test))

return(list(corr_train=corr_train, corr_test=corr_test, mse_train=mse_train, mse_test=mse_test, r_squared_train=r_squared_train, r_squared_test=r_squared_test))

}

# Define cross-validation function
# cv_model_splines <- function(data, k = 10) {
#     # Create the cross-validation folds
#     folds <- createFolds(data$cognition, k = k, list = TRUE, returnTrain = TRUE)
    
#     # Initialize vectors to store predictions and actual values
#     all_predictions_linear <- numeric(nrow(data))
#     all_predictions_varycoef <- numeric(nrow(data))
#     all_predictions_elastic_net <- numeric(nrow(data))
#     all_actuals <- data$cognition
    
#     cor_train_linear <- numeric(k)
#     cor_test_linear <- numeric(k)
#     cor_train_varycoef <- numeric(k)
#     cor_test_varycoef <- numeric(k)
#     cor_train_elastic_net <- numeric(k)
#     cor_test_elastic_net <- numeric(k)
    
#     # Perform cross-validation
#     for (i in seq_along(folds)) {
#         # Split data into training and test sets
#         train_data <- data[folds[[i]], ]
#         test_data <- data[-folds[[i]], ]
        
#         # Fit models
#         model_linear <- lm(cognition ~ Age + structural_IDP, data = train_data)
#         model_varycoef <- lm(cognition ~ ns(Age, df = 4) * structural_IDP, data = train_data)
        
#         # Prepare data for Elastic Net
#         x_train <- model.matrix(cognition ~ Age + structural_IDP, data = train_data)[, -1]
#         y_train <- train_data$cognition
#         x_test <- model.matrix(cognition ~ Age + structural_IDP, data = test_data)[, -1]
#         y_test <- test_data$cognition
        
#         # Fit Elastic Net model with cross-validation for lambda
#         model_elastic_net <- cv.glmnet(x_train, y_train, alpha = 0.5)
        
#         # Make predictions
#         train_data$predict_linear <- predict(model_linear, newdata = train_data)
#         train_data$predict_varycoef <- predict(model_varycoef, newdata = train_data)
#         test_data$predict_linear <- predict(model_linear, newdata = test_data)
#         test_data$predict_varycoef <- predict(model_varycoef, newdata = test_data)
        
#         # Predict with Elastic Net model
#         test_data$predict_elastic_net <- predict(model_elastic_net, newx = x_test, s = "lambda.min")
        
#         # Compute correlations and store them
#         cor_train_linear[i] <- cor(train_data$cognition, train_data$predict_linear)
#         cor_test_linear[i] <- cor(test_data$cognition, test_data$predict_linear)
#         cor_train_varycoef[i] <- cor(train_data$cognition, train_data$predict_varycoef)
#         cor_test_varycoef[i] <- cor(test_data$cognition, test_data$predict_varycoef)
#         cor_train_elastic_net[i] <- cor(train_data$cognition, predict(model_elastic_net, newx = x_train, s = "lambda.min"))
#         cor_test_elastic_net[i] <- cor(test_data$cognition, test_data$predict_elastic_net)
        
#         # Store predictions
#         all_predictions_linear[which(!1:nrow(data) %in% folds[[i]])] <- test_data$predict_linear
#         all_predictions_varycoef[which(!1:nrow(data) %in% folds[[i]])] <- test_data$predict_varycoef
#         all_predictions_elastic_net[which(!1:nrow(data) %in% folds[[i]])] <- test_data$predict_elastic_net
#     }
    
#     # Calculate correlation
#     cor_linear <- cor(all_actuals, all_predictions_linear)
#     cor_varycoef <- cor(all_actuals, all_predictions_varycoef)
#     cor_elastic_net <- cor(all_actuals, all_predictions_elastic_net)
    
#     # Calculate R-squared
#     ss_total <- sum((all_actuals - mean(all_actuals))^2)
#     ss_residual_linear <- sum((all_actuals - all_predictions_linear)^2)
#     r_squared_linear <- 1 - (ss_residual_linear / ss_total)
    
#     ss_residual_varycoef <- sum((all_actuals - all_predictions_varycoef)^2)
#     r_squared_varycoef <- 1 - (ss_residual_varycoef / ss_total)
    
#     ss_residual_elastic_net <- sum((all_actuals - all_predictions_elastic_net)^2)
#     r_squared_elastic_net <- 1 - (ss_residual_elastic_net / ss_total)
    
#     list(
#         cor_linear = cor_linear,
#         cor_varycoef = cor_varycoef,
#         cor_elastic_net = cor_elastic_net,
#         r_squared_linear = r_squared_linear,
#         r_squared_varycoef = r_squared_varycoef,
#         r_squared_elastic_net = r_squared_elastic_net,
#         cor_train_linear = cor_train_linear,
#         cor_test_linear = cor_test_linear,
#         cor_train_varycoef = cor_train_varycoef,
#         cor_test_varycoef = cor_test_varycoef,
#         cor_train_elastic_net = cor_train_elastic_net,
#         cor_test_elastic_net = cor_test_elastic_net,
#         all_predictions_linear = all_predictions_linear,
#         all_predictions_varycoef = all_predictions_varycoef,
#         all_predictions_elastic_net = all_predictions_elastic_net
#     )
# }
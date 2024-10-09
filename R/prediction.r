library(R.matlab)
# library(varycoef)
library(fastICA)
library(CCA)
library(glmnet)
library(lightgbm)
library(xgboost)
library(data.table)
library(caret)
library(tvem)
library(gsubfn)

#' Run Linear Model on Train and Test Data
#'
#' This function fits a linear model on the training data and uses it to predict 
#' outcomes for both the training and test data. The predictions are returned alongside 
#' the fitted linear model object.
#'
#' @param df_linear A data frame containing the dependent variable `y` and the independent variables (features).
#' @param id_train A numeric vector indicating the indices of the training data within the data frame.
#'
#' @return A list containing:
#' \item{fit_lm}{The fitted linear model object.}
#' \item{lm_yhat}{A numeric vector of predictions for both training and test data.}
#' @importFrom stats lm
#' @importFrom stats predict
#' @export
run_linear_model <- function(df_linear, id_train) {

  df_train <- df_linear[id_train, ]
  df_test <- df_linear[-id_train, ]

  print("Fitting linear model to test data...")
  fit_lm <- lm(y ~ ., data = df_train)

  # ridge regression (it doesn't really make sense to constrain the sum of the single coefficient to equal something)
  # so we don't do this for univariate predictions (could be worth for multivariate but just use glmnet)
  #fit_ridge <-  lm.ridge(y ~ ., data = df_train_linear, lambda = seq(0, .4, 1e-3))
  lm_yhat <- numeric(length(df_linear$y))
  lm_yhat[id_train] <- predict(fit_lm, newdata = df_train)
  lm_yhat[-id_train] <- predict(fit_lm, newdata = df_test)

  return(list(fit_lm = fit_lm, lm_yhat = lm_yhat))
}


#' Run Linear Model with 10-Fold Cross-Validation
#'
#' This function fits a linear model using 10-fold cross-validation. It returns 
#' predictions, correlations for both training and test sets, and the fitted models 
#' for each fold.
#'
#' @param idps A matrix or data frame of imaging-derived phenotypes (IDPs).
#' @param trait A numeric vector representing the trait (target) variable.
#' @param age A numeric vector representing age.
#' @param conf A matrix or data frame of confounding variables.
#' @param conf_names A character vector of confounder names.
#' @param trait_id A numeric identifier for the trait type; 999 indicates CCA, 0 indicates PCA.
#' @param remove_age A binary flag (0 or 1) indicating whether to remove age from confounding variables.
#' @param model_age A binary flag indicating whether to include age in the model.
#'
#' @return A list containing:
#' \item{lm_yhat}{A numeric vector of predictions for each fold (test set).}
#' \item{corr_train}{A numeric vector of correlations between true and predicted values on the training set for each fold.}
#' \item{corr_test}{A numeric vector of correlations between true and predicted values on the test set for each fold.}
#' \item{models}{A list of fitted linear models for each fold.}
#' \item{trait_transformed}{A numeric vector of the transformed trait used for prediction.}
#' @importFrom stats lm
#' @importFrom stats predict
#' @export
run_linear_model_cv <- function(idps, trait, age, conf, conf_names, trait_id, remove_age, model_age) {
  # We manually do the cross-validation so we can look at the predictions and compare
  # training / test set accuracies
  print("Fitting linear model with 10-fold cross-validation...")

  # Initialize vectors to store correlations for each fold
  corr_train <- numeric(10)
  corr_test <- numeric(10)

  # Perform 10-fold cross-validation manually to calculate correlations
  lm_yhat <- numeric(length(age)) # to store predictions
  trait_transformed <- numeric(length(age)) # to store predictions
  folds <- caret::createFolds(age, k = 10)
  models <- vector("list", 10)  # 10 folds, so 10 models

  for (i in seq_along(folds)) {
    print(paste("Fold", i))
    test_idx <- folds[[i]]
    train_idx <- setdiff(seq_along(age), test_idx)

    preprocessed_data <- pre_process_data_cross_validated(idps, trait, age, conf, conf_names, trait_id, train_idx, remove_age, ica=0, n_feat=0)
    idps_fold <- preprocessed_data$idps
    trait_fold <- preprocessed_data$trait
    age_fold <- preprocessed_data$age

    df <- data.frame(y = trait_fold, idps_fold)

    df <- prepare_age_data(df, age_fold, model_age)

    # Train the model on the training fold (and save model info)
    fit_lm <- lm(y ~ ., data = df[train_idx, ])
    models[[i]] <- fit_lm

    # Predictions for the training and test set
    lm_yhat_train <- predict(fit_lm, newdata = df[train_idx, ])
    lm_yhat_test <- predict(fit_lm, newdata = df[test_idx, ])

    # Note train/test accuracies
    corr_train[i] <- stats::cor(df$y[train_idx], lm_yhat_train)
    corr_test[i] <- stats::cor(df$y[test_idx], lm_yhat_test)

    # Store predictions in the original index positions
    lm_yhat[test_idx] <- lm_yhat_test
    trait_transformed[test_idx] <- trait_fold[test_idx]
  }

  # Return the fitted model, predictions, and accuracies
  return(list(
    lm_yhat = lm_yhat,
    corr_train = corr_train,
    corr_test = corr_test,
    models = models,
    trait_transformed = trait_transformed
  ))
}

# Function to calculate performance metrics
calculate_metrics <- function(actual_train, predicted_train, actual_test, predicted_test) {
  # Initialize results list
  metrics <- list()

  # Calculate correlations
  metrics$corr_train <- stats::cor(actual_train, predicted_train)
  metrics$corr_test <- stats::cor(actual_test, predicted_test)

  # Calculate Mean Squared Error
  metrics$mse_train <- mean((actual_train - predicted_train) ^ 2)
  metrics$mse_test <- mean((actual_test - predicted_test) ^ 2)

  # Calculate R-squared
  metrics$r_squared_train <- 1 - (sum((actual_train - predicted_train) ^ 2) / sum((actual_train - mean(actual_train)) ^ 2))
  metrics$r_squared_test <- 1 - (sum((actual_test - predicted_test) ^ 2) / sum((actual_test - mean(actual_test)) ^ 2))

  return(metrics)
}

#' Run Elastic Net Model
#'
#' Fits an Elastic Net model with or without age as a covariate and returns the fitted model
#' and predicted values for both training and test sets.
#'
#' @param idps_linear A matrix of IDP (input data) features.
#' @param trait A vector of target variable (traits) to predict.
#' @param id_train Indices of training data.
#' @param age A vector of ages corresponding to the observations.
#' @param model_age Integer to specify how age is treated in the model:
#' 0 = No age, 1 = Age as a covariate, 2 = Age only, 3 = IDPs + Age.
#'
#' @return A list containing:
#' \item{cvfit_glmnet}{The fitted Elastic Net model.}
#' \item{yhat_train}{Predictions for the training data.}
#' \item{yhat_test}{Predictions for the test data.}
#' @import glmnet
#' @export
run_elastic_net_model <- function(idps_linear, trait, id_train, age, model_age) {

  print("Fitting elastic net model...")
  trait_train <- trait[id_train]
  idps_linear <- idps_linear
  idps_linear_train <- idps_linear[id_train, ]
  age_train <- matrix(age[id_train])

  if (model_age == 1 ||  model_age == 3) {
    cvfit_glmnet <- glmnet::cv.glmnet(cbind(idps_linear_train, age_train), trait_train)
  } else if (model_age == 0) {
    cvfit_glmnet <- glmnet::cv.glmnet(idps_linear_train, trait_train)
  } else if (model_age == 2) {
    cvfit_glmnet <- glmnet::cv.glmnet(cbind(rep(1, length(age_train)), age_train), trait_train)
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

#' Run Elastic Net Model with Cross-Validation
#'
#' Performs k-fold cross-validation using Elastic Net, predicting the trait from IDP data
#' with optional inclusion of age as a covariate.
#'
#' @param idps A matrix of IDP (input data) features.
#' @param trait A vector of target variable (traits) to predict.
#' @param age A vector of ages corresponding to the observations.
#' @param conf A matrix of confounding variables to be removed (optional).
#' @param conf_names Column names for the confound matrix.
#' @param trait_id Trait identifier used in pre-processing.
#' @param remove_age Logical, whether to remove age effect during pre-processing.
#' @param model_age Integer to specify how age is treated in the model:
#' 0 = No age, 1 = Age as a covariate, 2 = Age only, 3 = IDPs + Age.
#' @param alpha Elastic net mixing parameter (0 = ridge, 1 = lasso).
#' @param n_folds Number of cross-validation folds.
#'
#' @return A list containing:
#' \item{yhat}{Predicted values for all samples across the test folds.}
#' \item{corr_train}{Training set correlations for each fold.}
#' \item{corr_test}{Test set correlations for each fold.}
#' \item{models}{List of fitted Elastic Net models for each fold.}
#' \item{trait_transformed}{Transformed trait used during cross-validation.}
#' @import glmnet
#' @export
run_elastic_net_model_cv <- function(idps, trait, age, conf, conf_names, trait_id, remove_age, model_age, alpha = 1, n_folds = 10) {
  print("Fitting elastic net model with cross-validation...")

  # Initialize vectors to store predictions
  yhat <- numeric(length(age))
  trait_transformed <- numeric(length(age))

  # Create folds for cross-validation
  folds <- caret::createFolds(age, k = n_folds)

  # Initialize lists to store models and metrics
  models <- vector("list", n_folds)
  corr_train <- numeric(n_folds)
  corr_test <- numeric(n_folds)

  for (i in seq_along(folds)) {
    print(paste("Fold", i))

    # Get train and test indices for the current fold
    test_idx <- folds[[i]]
    train_idx <- setdiff(seq_along(age), test_idx)

    preprocessed_data <- pre_process_data_cross_validated(idps, trait, age, conf, conf_names, trait_id, train_idx, remove_age, ica=0, n_feat=0)
    idps_fold <- preprocessed_data$idps
    trait_fold <- preprocessed_data$trait
    age_fold <- preprocessed_data$age

    # Prepare training and testing data
    trait_train <- trait_fold[train_idx]
    idps_train <- idps_fold[train_idx, ]
    age_train <- matrix(age_fold[train_idx])

    trait_test <- trait_fold[test_idx]
    idps_test <- idps_fold[test_idx, ]
    age_test <- matrix(age_fold[test_idx])

    # Fit the elastic net model
    if (model_age == 1 || model_age == 3) {
      cvfit_glmnet <- glmnet::cv.glmnet(cbind(idps_train, age_train), trait_train, alpha = alpha)
    } else if (model_age == 0) {
      cvfit_glmnet <- glmnet::cv.glmnet(idps_train, trait_train, alpha = alpha)
      # cvfit_glmnet <- cva.glmnet(idps_linear_train, trait_train, alpha = seq(0, 1, 0.05))
    } else if (model_age == 2) {
      cvfit_glmnet <- glmnet::cv.glmnet(cbind(rep(1, length(age_train)), age_train), trait_train, alpha = alpha)
    }

    # Store the model
    models[[i]] <- cvfit_glmnet

    # Make predictions
    if (model_age == 1 || model_age == 3) {
      yhat_train <- predict(cvfit_glmnet, newx = cbind(idps_train, age_train), s = "lambda.min")
      yhat_test <- predict(cvfit_glmnet, newx = cbind(idps_test, age_test), s = "lambda.min")
    } else if (model_age == 0) {
      yhat_train <- predict(cvfit_glmnet, newx = idps_train, s = "lambda.min")
      yhat_test <- predict(cvfit_glmnet, newx = idps_test, s = "lambda.min")
    } else if (model_age == 2) {
      yhat_train <- predict(cvfit_glmnet, newx = cbind(rep(1, length(age_train)), age_train), s = "lambda.min")
      yhat_test <- predict(cvfit_glmnet, newx = cbind(rep(1, length(age_test)), age_test), s = "lambda.min")
    }

    # Store predictions in the original index positions
    yhat[test_idx] <- yhat_test
    trait_transformed[test_idx] <- trait_fold[test_idx]

    # Calculate and store performance metrics
    corr_train[i] <- stats::cor(trait_train, yhat_train)
    corr_test[i] <- stats::cor(trait_test, yhat_test)
  }

  # Return results
  return(list(
    yhat = yhat,
    corr_train = corr_train,
    corr_test = corr_test,
    models = models,
    trait_transformed = trait_transformed
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


#' Run LightGBM Model
#'
#' Fits a LightGBM model on training data and evaluates it on test data.
#'
#' @param df_all_train_x Training data matrix (features).
#' @param df_all_train_y Training data vector (target).
#' @param df_all_test_x Test data matrix (features).
#' @param df_all_test_y Test data vector (target).
#' @param params List of LightGBM parameters (optional).
#'
#' @return A list containing:
#' \item{predictions_train}{Predicted values for the training data.}
#' \item{predictions_test}{Predicted values for the test data.}
#' \item{lgb_model}{The trained LightGBM model.}
#' @import lightgbm
#' @export
run_lgboost_model <- function(df_all_train_x, df_all_train_y, df_all_test_x, df_all_test_y, params = 0) {

  train_lgb <- lightgbm::lgb.Dataset(data = df_all_train_x, label = df_all_train_y)
  test_lgb <- lightgbm::lgb.Dataset(data = df_all_test_x, label = df_all_test_y)
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
  lgb_model <- lightgbm::lgb.train(
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


#' Run LightGBM with Cross-Validation
#'
#' Performs k-fold cross-validation using LightGBM, fitting the model on the training data
#' and evaluating it on the test data.
#'
#' @param df_boost A data frame containing features and target variable ("y").
#' @param params A list of LightGBM parameters (optional).
#' @param n_folds Number of cross-validation folds.
#'
#' @return A list containing:
#' \item{predictions_test}{Predicted values for all samples across the test folds.}
#' \item{corr_train}{Training set correlations for each fold.}
#' \item{corr_test}{Test set correlations for each fold.}
#' \item{mse_train}{Mean squared error on the training set for each fold.}
#' \item{mse_test}{Mean squared error on the test set for each fold.}
#' \item{r_squared_train}{R-squared values for the training set for each fold.}
#' \item{r_squared_test}{R-squared values for the test set for each fold.}
#' @import lightgbm
#' @export
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
  folds <- caret::createFolds(df_boost$y, k = n_folds)

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

    train_lgb <- lightgbm::lgb.Dataset(data = df_train_x, label = df_train_y)
    test_lgb <- lightgbm::lgb.Dataset(data = df_test_x, label = df_test_y)

    # Train the LightGBM model
    lgb_model <- lightgbm::lgb.train(
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
    corr_train[i] <- stats::cor(df_boost$y[train_idx], predictions_train)
    corr_test[i] <- stats::cor(df_boost$y[test_idx], predictions_test)
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
  corr_test <- stats::cor(predictions_test, df_all_test_y)
  corr_train <- stats::cor(predictions_train, df_all_train_y)

  print(paste("Mean Squared Error (TRAIN):", mse_train))
  print(paste("Mean Squared Error (TEST):", mse_test))
  print(paste("R-squared (TRAIN):", r_squared_train))
  print(paste("R-squared (TEST):", r_squared_test))
  print(paste("Correlation (TRAIN):", corr_train))
  print(paste("Correlation (TEST):", corr_test))

return(list(mse_test = mse_test, mse_train = mse_train, r_squared_test = r_squared_test, r_squared_train = r_squared_train, corr_test = corr_test, corr_train = corr_train))
}


get_cv_model_stats <- function(yhat, trait) {
  # Calculate metrics
  corr <- stats::cor(yhat, trait)
  mse <- mean((yhat - trait) ^ 2)
  r_squared <- 1 - (sum((trait - yhat) ^ 2) / sum((trait - mean(trait)) ^ 2))

  print(paste("Correlation:", corr))
  print(paste("Mean Squared Error:", mse))
  print(paste("R-squared:", r_squared))

  return(list(corr = corr, mse = mse, r_squared = r_squared))
}

#' Run XGBoost Model
#'
#' Trains an XGBoost model using the provided training data and tests the model on the provided test data.
#'
#' @param df_all_train_x A matrix or data frame containing the training data features.
#' @param df_all_train_y A vector containing the training data labels.
#' @param df_all_test_x A matrix or data frame containing the test data features.
#' @param df_all_test_y A vector containing the test data labels.
#'
#' @return A list containing:
#'   - predictions_train: Predicted values for the training data.
#'   - predictions_test: Predicted values for the test data.
#'   - xgb_model: The fitted XGBoost model object.
#'
#' @import xgboost
#' @export
run_xgboost_model <- function(df_all_train_x, df_all_train_y, df_all_test_x, df_all_test_y) {

  dtrain <- xgboost::xgb.DMatrix(data = df_all_train_x, label = df_all_train_y)
  dtest <- xgboost::xgb.DMatrix(data = df_all_test_x, label = df_all_test_y)

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

  xgb_model <- xgboost::xgb.train(
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
    importance_matrix <- lightgbm::lgb.importance(model = model, percentage = TRUE)
  } else if (run_svc == 3) {
    importance_matrix <- xgboost::xgb.importance(model = model)
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

#' Run Spline Model with 10-fold Cross-Validation
#'
#' Fits a spline model using 10-fold cross-validation and calculates correlations for both training and test sets.
#'
#' @param idps A matrix or data frame of input features (IDPs).
#' @param trait A vector of the trait to predict.
#' @param age A vector of subject ages.
#' @param conf A data frame of additional confounding variables.
#' @param conf_names A character vector of confound variable names.
#' @param trait_id An integer representing the trait ID.
#' @param remove_age A boolean indicating whether to remove the age variable from the model.
#' @param model_age A boolean indicating whether to include age as a linear feature in the model.
#'
#' @return A list containing:
#'   - spline_yhat: Predicted trait values for the test data.
#'   - corr_train: Correlations for the training folds.
#'   - corr_test: Correlations for the test folds.
#'   - models: List of fitted models for each fold.
#' @importFrom stats lm
#' @importFrom stats predict
#' @import splines
#' @export
run_spline_model_cv <- function(idps, trait, age, conf, conf_names, trait_id, remove_age, model_age) {
  # We manually do the cross-validation so we can look at the predictions and compare
  # training / test set accuracies
  print("Fitting spline model with 10-fold cross-validation...")

  # Initialize vectors to store correlations for each fold
  corr_train <- numeric(10)
  corr_test <- numeric(10)


  # Initialize a vector to store predictions
  spline_yhat <- numeric(length(age))

  # Perform 10-fold cross-validation manually to calculate correlations
  folds <- caret::createFolds(age, k = 10)
  models <- vector("list", 10)  # 10 folds, so 10 models

  for (i in seq_along(folds)) {
    print(paste("Fold", i))
    test_idx <- folds[[i]]
    train_idx <- setdiff(seq_along(df_spline$y), test_idx)

    preprocessed_data <- pre_process_data_cross_validated(idps, trait, age, conf, conf_names, trait_id, train_idx, remove_age, ica=0, n_feat=0)
    idps_fold <- preprocessed_data$idps
    trait_fold <- preprocessed_data$trait
    age_fold <- preprocessed_data$age

    df_spline <- data.frame(y = trait_fold, x = idps_fold, age = age_fold)

    # rename x. column to idp
    idp_column <- grep("^x\\.", names(df_spline), value = TRUE)
    names(df_spline)[names(df_spline) == idp_column] <- "structural_IDP"

    # determine best df for spline
    spline_df <- determine_optimal_df_spline(df_spline, train_idx)

    print(paste("Spline df: ", spline_df))

    # Train the model on the training data
    fit_spline <- lm(y ~ splines::ns(age, df = spline_df, intercept = TRUE) * structural_IDP, data = df_spline[train_idx, ])
    # option 2: fit_spline <- lm(y ~ splines::ns(age[train_idx], df = spline_df, intercept = TRUE) * structural_IDP, data = df_spline[train_idx, ])
    models[[i]] <- fit_spline

    # Predictions for the training and test set
    spline_yhat_train <- predict(fit_spline, newdata = df_spline[train_idx, ])
    spline_yhat_test <- predict(fit_spline, newdata = df_spline[test_idx, ])

    # Note accuracies
    corr_train[i] <- stats::cor(df_spline$y[train_idx], spline_yhat_train)
    corr_test[i] <- stats::cor(df_spline$y[test_idx], spline_yhat_test)

    # Store predictions in the original index positions
    spline_yhat[test_idx] <- spline_yhat_test
  }

  # Return the fitted model, predictions, and accuracies
  return(list(
    spline_yhat = spline_yhat,
    corr_train = corr_train,
    corr_test = corr_test,
    models = models
  ))
}


#' Run Spline Model with Elastic Net Regularization (10-fold Cross-Validation)
#'
#' Fits a spline model with elastic net regularization using 10-fold cross-validation and calculates correlations for both training and test sets.
#'
#' @param idps A matrix or data frame of input features (IDPs).
#' @param trait A vector of the trait to predict.
#' @param age A vector of subject ages.
#' @param conf A data frame of additional confounding variables.
#' @param conf_names A character vector of confound variable names.
#' @param trait_id An integer representing the trait ID.
#' @param remove_age A boolean indicating whether to remove the age variable from the model.
#' @param model_age A boolean indicating whether to include age as a linear feature in the model.
#'
#' @return A list containing:
#'   - spline_yhat: Predicted trait values for the test data.
#'   - corr_train_cv: Correlations for the training folds.
#'   - corr_test_cv: Correlations for the test folds.
#'   - models: List of fitted models for each fold.
#'   - trait_transformed: Transformed trait values for the test data.
#' @import splines
#' @export
run_spline_model_cv_enet <- function(idps, trait, age, conf, conf_names, trait_id, remove_age, model_age) {
  print("Fitting spline model with 10-fold cross-validation...")

  # Initialize vectors to store correlations and MSE for each fold
  corr_train_cv <- numeric(10)
  corr_test_cv <- numeric(10)
  spline_yhat <- numeric(length(age))
  trait_transformed <- numeric(length(age))

  # Perform 10-fold cross-validation manually
  folds <- caret::createFolds(age, k = 10)
  models <- vector("list", 10)

  for (i in seq_along(folds)) {
    print(paste("Fold", i))
    test_idx <- folds[[i]]
    train_idx <- setdiff(seq_along(age), test_idx)

    # preprocess data
    preprocessed_data <- pre_process_data_cross_validated(idps, trait, age, conf, conf_names, trait_id, train_idx, remove_age, ica=0, n_feat=0)
    idps_fold <- preprocessed_data$idps
    trait_fold <- preprocessed_data$trait
    age_fold <- preprocessed_data$age

    # create data frame for spline model
    if (is.vector(idps_fold)) {
      df_spline <- data.frame(y = trait_fold, x.1 = idps_fold, age = age_fold)
      idp_columns <- "x.1"
    } else {
      df_spline <- data.frame(y = trait_fold, x = idps_fold, age = age_fold)
      idp_columns <- grep("^x\\.", names(df_spline), value = TRUE)
    }

    # Determine the best df for the spline
    ########### need to sort this out with regularization I think
    spline_deg_free <- determine_optimal_df_spline(df_spline, train_idx)
    print(paste("Spline df: ", spline_deg_free))

    # Create the spline basis matrix for the age variable for training data
    spline_basis_train <- splines::ns(df_spline$age[train_idx], df = spline_deg_free, intercept = TRUE)

    interaction_terms <- paste0("spline_basis_train * ", idp_columns, collapse = " + ")

    # add age as a linear features
    if (model_age == 1) {interaction_terms <- paste0("age + ", interaction_terms)}

    # Create design matrix for training data
    X_train <- stats::model.matrix(stats::as.formula(paste("~", interaction_terms, "- 1")), data = df_spline[train_idx, ])


    # Create the spline basis matrix for test data (using training data)
    spline_basis_test <- predict(spline_basis_train, newx = df_spline$age[test_idx])
    interaction_terms <- paste0("spline_basis_test * ", idp_columns, collapse = " + ")
    if (model_age == 1) {interaction_terms <- paste0("age + ", interaction_terms)}
    X_test <- stats::model.matrix(stats::as.formula(paste("~", interaction_terms, "- 1")), data = df_spline[test_idx, ])


    # Fit the model using Ridge regularization (L2)
    # WE USE RIDGE HERE BECAUSE WE ARE DOING UNIVARIATE PREDICTIONS, SO L1 REGULARIZATION IS FUTILE I ASSUME
    fit_spline <- glmnet::cv.glmnet(X_train, df_spline$y[train_idx], alpha = 0, nfolds = 10)
    models[[i]] <- fit_spline

    # Predictions for the training and test set
    spline_yhat_train <- predict(fit_spline, newx = X_train, s = "lambda.min")
    spline_yhat_test <- predict(fit_spline, newx = X_test, s = "lambda.min")

    # Note accuracies
    corr_train_cv[i] <- stats::cor(df_spline$y[train_idx], spline_yhat_train)
    corr_test_cv[i] <- stats::cor(df_spline$y[test_idx], spline_yhat_test)

    # Store predictions in the original index positions
    spline_yhat[test_idx] <- spline_yhat_test
    trait_transformed[test_idx] <- trait_fold[test_idx]

  }

  # Return the fitted model, predictions, and accuracies
  return(list(
    spline_yhat = spline_yhat,
    corr_train_cv = corr_train_cv,
    corr_test_cv = corr_test_cv,
    models = models,
    trait_transformed = trait_transformed
  ))
}

determine_optimal_df_spline <- function(df_spline, train_idx) {

  train_control <- caret::trainControl(method = "cv", number = 10)

  idp_columns <- grep("^x\\.", names(df_spline), value = TRUE)

  # Initialize a list to store results
  results <- list()

  # Define the range of degrees of freedom you want to test
  df_values <- 2:5 # 2 is linear i think # no it's not it's how many knots to use (I think 2 is just boundary knots so not using a spline basically (other than at the boundaries))

  # Loop over each df value
  for (df in df_values) {
    # Create the formula dynamically
    # Generate the spline terms
    spline_term <- paste0("splines::ns(age, df = ", df, ", intercept = TRUE)")

    # Create interaction terms with idp_columns
    interaction_terms <- paste0(spline_term, " * ", paste(idp_columns, collapse = " + "))
    
  # Construct the formula
   formula <- stats::as.formula(paste("y ~", interaction_terms))

    # Train the model using the dynamically created formula
    model <- caret::train(formula,
                   data = df_spline[train_idx, ],
                   method = "lm",
                   trControl = train_control)

  # Store the model and associated df in the results list
    results[[paste("df", df, sep = "_")]] <- model
  }

  # Compare models by their performance
  model_performance <- sapply(results, function(model) {
    min(model$results$RMSE)  # or another performance metric
  })

  # Find the best df value
  best_df <- df_values[which.min(model_performance)]
  #best_model <- results[[paste("df", best_df, sep = "_")]]

  return(best_df)

}

interpolate_predictions <- function(model, ages, idp_values) {
  # Extract the time grid and fitted coefficients from the model
  time_grid <- model$time_grid
  beta_idp <- model$grid_fitted_coefficients$idp1$estimate
  beta_intercept <- model$grid_fitted_coefficients$`(Intercept)`$estimate

  # Interpolate to find the coefficients for the given ages
  interpolated_idp <- stats::approx(time_grid, beta_idp, xout = ages, rule = 2)$y
  interpolated_intercept <- stats::approx(time_grid, beta_intercept, xout = ages, rule = 2)$y

  # Calculate the predictions
  predictions <- interpolated_intercept + interpolated_idp * idp_values

  return(predictions)
}

#' Run 10-Fold Cross-Validation for Time-Varying Effects Model (TVEM)
#'
#' This function performs a 10-fold cross-validation to fit a time-varying effects model
#' (TVEM) and predict a given trait using IDP (imaging-derived phenotypes) and age. 
#' It supports optional confounder adjustment, alignment of trait directions across folds,
#' and provides both training and test set correlations.
#'
#' @param idp A matrix or dataframe of imaging-derived phenotypes (IDPs).
#' @param trait A numeric vector representing the trait to be predicted.
#' @param age A numeric vector representing the age of the subjects.
#' @param conf A matrix or dataframe of confounding variables to adjust for (optional).
#' @param conf_names A vector of strings with the names of confounding variables (optional).
#' @param trait_id An identifier for the trait being modeled.
#' @param remove_age Logical; if TRUE, removes age from the set of predictors.
#' @param model_age Logical; if TRUE, includes age as a linear effect in the model.
#'
#' @return A list containing the following components:
#'   \item{tvem_yhat}{A numeric vector of predictions for the trait based on the TVEM model.}
#'   \item{corr_train}{A numeric vector of correlations between predicted and actual trait values in the training set for each fold.}
#'   \item{corr_test}{A numeric vector of correlations between predicted and actual trait values in the test set for each fold.}
#'   \item{models}{A list of fitted TVEM models for each fold.}
#'   \item{trait_transformed}{A numeric vector of the transformed trait, aligned across folds.}
#'
#' @export
run_tvem_model_cv <- function(idp, trait, age, conf, conf_names, trait_id, remove_age, model_age) {
  # We manually do the cross-validation so we can look at the predictions and compare
  # training / test set accuracies
  print("Fitting time-varying effects model with 10-fold cross-validation...")

  # Initialize vectors to store correlations for each fold
  corr_train <- numeric(10)
  corr_test <- numeric(10)

  # Perform 10-fold cross-validation manually to calculate correlations
  tvem_yhat <- numeric(length(age)) # to store predictions
  trait_transformed <- numeric(length(age)) # to store predictions
  folds <- caret::createFolds(age, k = 10)
  models <- vector("list", 10)  # 10 folds, so 10 models

  # Use the first fold as a reference for aligning the trait direction
  reference_trait <- NULL

  for (i in seq_along(folds)) {
    print(paste("Fold", i))
    test_idx <- folds[[i]]
    train_idx <- setdiff(seq_along(age), test_idx)

    preprocessed_data <- pre_process_data_cross_validated(idp, trait, age, conf, conf_names, trait_id, train_idx, remove_age, ica=0, n_feat=0)
    idps_fold <- preprocessed_data$idps
    trait_fold <- preprocessed_data$trait
    age_fold <- preprocessed_data$age

    # Use the first fold's trait direction as the reference
    if (is.null(reference_trait)) {
      reference_trait <- trait_fold
    } else {
      # Align current fold's trait direction with the reference trait direction
      if (stats::cor(trait_fold, reference_trait) < 0) {
        trait_fold <- -trait_fold
      }
    }

    # Generate subject IDs based on the number of subjects in the fold
    subject_id <- seq_along(trait_fold)

    # Create dataframe with subject_id included
    df <- data.frame(subject_id = subject_id, y = trait_fold, idps_fold)
    df <- prepare_age_data(df, age_fold, 1)

    # # Train the model on the training fold (and save model info)
    if (model_age == 0) {
      fit_tvem <- tvem::select_tvem(data = df[train_idx, ],
                     formula = y ~ idps_fold,
                     id = subject_id, # Assuming subject_id is available
                     time = age,
                     grid = length(trait_fold),
                     max_knots = 5)
    } else if (model_age == 1) {
      fit_tvem <- tvem::select_tvem(data = df[train_idx, ],
                     formula = y ~ idps_fold,
                     #invar_effect = ~age,
                     id = subject_id, # Assuming subject_id is available
                     time = age,
                     #num_knots = 1,
                     spline_order = 3,
                     grid = length(trait_fold),
                     max_knots = 1)
    }

    models[[i]] <- fit_tvem

    # Predictions for the training and test set
    tvem_yhat_train <- interpolate_predictions(fit_tvem, df$age[train_idx], df$idps_fold[train_idx])
    tvem_yhat_test <- interpolate_predictions(fit_tvem, df$age[test_idx], df$idps_fold[test_idx])

    # Note train/test accuracies
    corr_train[i] <- stats::cor(df$y[train_idx], tvem_yhat_train)
    corr_test[i] <- stats::cor(df$y[test_idx], tvem_yhat_test)

    print(paste("corr fold (train):", corr_train[i]))
    print(paste("corr fold (test):", corr_test[i]))

    # Store predictions in the original index positions
    tvem_yhat[test_idx] <- tvem_yhat_test
    trait_transformed[test_idx] <- trait_fold[test_idx]
  }

  # Return the fitted model, predictions, and accuracies
  return(list(
    tvem_yhat = tvem_yhat,
    corr_train = corr_train,
    corr_test = corr_test,
    models = models,
    trait_transformed = trait_transformed
  ))
}

# run_tvem_model_cv_2 <- function(idp, trait, age, conf, conf_names, trait_id, remove_age, model_age) {
#   print("Fitting time-varying effects model with 10-fold cross-validation...")

#   # Initialize vectors to store correlations for each fold
#   corr_train <- numeric(10)
#   corr_test <- numeric(10)

#   # Perform 10-fold cross-validation manually to calculate correlations
#   tvem_yhat <- numeric(length(age)) # to store predictions
#   trait_transformed <- numeric(length(age)) # to store predictions
#   folds <- caret::createFolds(age, k = 10)
#   models <- vector("list", 10)  # 10 folds, so 10 models

#   # Use the first fold as a reference for aligning the trait direction
#   reference_trait <- NULL

#   for (i in seq_along(folds)) {
#     print(paste("Fold", i))
#     test_idx <- folds[[i]]
#     train_idx <- setdiff(seq_along(age), test_idx)

#     preprocessed_data <- pre_process_data_cross_validated(idp, trait, age, conf, conf_names, trait_id, train_idx, remove_age, ica=0, n_feat=0)
#     idps_fold <- preprocessed_data$idps
#     trait_fold <- preprocessed_data$trait
#     age_fold <- preprocessed_data$age

#     # Use the first fold's trait direction as the reference
#     if (is.null(reference_trait)) {
#       reference_trait <- trait_fold
#     } else {
#       if (stats::cor(trait_fold, reference_trait) < 0) {
#         trait_fold <- -trait_fold
#       }
#     }

#     # Generate subject IDs based on the number of subjects in the fold
#     subject_id <- seq_along(trait_fold)

#     # Create dataframe with subject_id included
#     df <- data.frame(subject_id = subject_id, y = trait_fold, idps_fold)
#     df <- prepare_age_data(df, age_fold, 1)

#     df_train <- df[train_idx, ]

#     print(head(df_train))

#     # Train the model on the training fold (and save model info)
#     if (model_age == 0) {
#       print("test")
#       fit_tvem <- run_tvem(df_train)
#       # fit_tvem <- tvem::select_tvem(data = df_train,
#       #                formula = y ~ idps_fold,
#       #                id = subject_id,
#       #                time = age)
#       print("test2")
#                      #grid = length(trait_fold),
#                      #max_knots = 5)
#       # fit_tvem <- tvem::tvem(data = df_train,
#       #                formula = y ~ idps_fold,
#       #                id = subject_id,
#       #                time = age,
#       #                grid = length(trait_fold),
#       #                num_knots = 1)
#     } else if (model_age == 1) {
#       fit_tvem <- tvem::select_tvem(data = df_train,
#                      formula = y ~ idps_fold,
#                      id = subject_id,
#                      time = age,
#                      spline_order = 3,
#                      grid = length(trait_fold),
#                      max_knots = 1)
#     }

#     models[[i]] <- fit_tvem

#     # Predictions for the training and test set
#     tvem_yhat_train <- interpolate_predictions(fit_tvem, df$age[train_idx], df$idps_fold[train_idx])
#     tvem_yhat_test <- interpolate_predictions(fit_tvem, df$age[test_idx], df$idps_fold[test_idx])

#     # Note train/test accuracies
#     corr_train[i] <- stats::cor(df$y[train_idx], tvem_yhat_train)
#     corr_test[i] <- stats::cor(df$y[test_idx], tvem_yhat_test)

#     print(paste("corr fold (train):", corr_train[i]))
#     print(paste("corr fold (test):", corr_test[i]))

#     # Store predictions in the original index positions
#     tvem_yhat[test_idx] <- tvem_yhat_test
#     trait_transformed[test_idx] <- trait_fold[test_idx]
#   }

#   # Return the fitted model, predictions, and accuracies
#   return(list(
#     tvem_yhat = tvem_yhat,
#     corr_train = corr_train,
#     corr_test = corr_test,
#     models = models,
#     trait_transformed = trait_transformed
#   ))
# }
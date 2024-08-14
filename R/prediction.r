library(R.matlab)
library(varycoef)
library(optparse)
library(fastICA)
library(CCA)
library(logger)
library(glmnet)


get_best_feat_svc <- function(trait_id, n_feat, n_sub, perc_train, prof, tap, cov, n, ica) {
  dir <- '/gpfs3/well/win-fmrib-analysis/users/psz102/age_varying_coefficients/results/univariate'
  perc_train <- 90 #75 # might want to update this to 90
  n_idps <- 1436 # 10
  n_sub <- 20000 # 10000
  n_feat <- 1436 # 10
  model_intercept <- 0
  corr <- rep(NA, n_idps)


  for (idp in 1:n_idps) {
    print(idp)
    pred <- try(load(sprintf("%s/yhat_t_%i_n_f_%i_f_1436_n_%i_p_%i_l_%i_ta_%i_c_%s_idp_%i_ica_%i_mi_%i.RData",
                  dir, trait_id, n_feat, n_sub, perc_train, prof, tap, cov, idp, ica, model_intercept)))
    if (inherits(pred, "try-error")) {
      #error handling code, maybe just skip this iteration using
      next
    }
    # pred <- mget(load(sprintf("%s/yhat_t_%i_f_1436_n_%i_p_%i_l_%i_ta_%i_c_%s_idp_%i_ica_%i.RData",
    #              dir, trait_id, n_sub, perc_train, prof, tap, cov, idp, ica), envir=(NE. <- new.env())), envir=NE.)

    # y <- extractorRData(sprintf("%s/yhat_t_%i_n_f_%i_f_1436_n_%i_p_%i_l_%i_ta_%i_c_%s_idp_%i_ica_%i_mi_%i.RData",
    #               dir, trait_id, n_feat, n_sub, perc_train, prof, tap, cov, idp, ica, model_intercept), 'y') 

    # yhat <- extractorRData(sprintf("%s/yhat_t_%i_n_f_%i_f_1436_n_%i_p_%i_l_%i_ta_%i_c_%s_idp_%i_ica_%i_mi_%i.RData",
    #               dir, trait_id, n_feat, n_sub, perc_train, prof, tap, cov, idp, ica, model_intercept), 'yhat') 

    corr_svc <- extractorRData(sprintf("%s/yhat_t_%i_n_f_%i_f_1436_n_%i_p_%i_l_%i_ta_%i_c_%s_idp_%i_ica_%i_mi_%i.RData",
                   dir, trait_id, n_feat, n_sub, perc_train, prof, tap, cov, idp, ica, model_intercept), 'corr_svc')
                                      
    corr[idp] <- corr_svc #cor(y, yhat)
    #corr[idp] <- cor(pred$y, pred$yhat)
  }
  # get idx of best 'n' features
  idx_best_features <- extract_best_features(corr, n)

  return(idx_best_features)
}

run_univ_pred_svc <- function(df_train, id_train_inner, idp, age_train, svc_config, model_vary_intercept) {

  # get training subjects dataframe (including age to act as modulator)
  age_train_inner <- age_train[id_train_inner]
  age_test_inner <- age_train[-id_train_inner]

  # create dataframes
  df_train_inner <- data.frame(df_train[id_train_inner, c(1, idp + 1)], age_train_inner)
  df_test_inner <- data.frame(df_train[-id_train_inner, c(1, idp + 1)], age_test_inner)

  colnames(df_train_inner)[2] <- "idp"
  colnames(df_test_inner)[2] <- "idp"

  # fit model (depends on if we want to model intercept)
  if (model_vary_intercept == 1) {
    fit_svc <- SVC_mle(formula = y ~ idp, data = df_train_inner, locs = df_train_inner$age_train_inner)
  } else if (model_vary_intercept == 0) {
    fit_svc <- SVC_mle(formula = y ~ idp, data = df_train_inner, locs = df_train_inner$age_train_inner, RE_formula = y ~ idp - 1)
  }

  # make prediction
  df_svc_pred <- predict(fit_svc, newdata = df_test_inner, newlocs = df_test_inner$age_test_inner, control = svc_config)
  df_svc_pred_insample <- predict(fit_svc, newdata = df_train_inner, newlocs = df_train_inner$age_train_inner, control = svc_config)

  # note down predictions
  yhat <- df_svc_pred$y.pred
  y <- df_train$y[-id_train_inner]
  yhat_insample <- df_svc_pred_insample$y.pred
  y_insample <- df_train$y[id_train_inner]

  # make prediction
  #df_svc_pred <- predict(fit_svc, newlocs = df_test_inner$age_test_inner)


  # create list of prediction and target
  prediction_list <- list("yhat" = yhat, "y" = y, "fit_svc" = fit_svc, "yhat_insample" = yhat_insample, "y_insample" = y_insample)

  return(prediction_list)
}


run_univ_pred_lm <- function(df_train_all, id_train_inner, idp) {

  # get training subjects dataframe
  df_train_inner <- data.frame(df_train_all[id_train_inner, c(1, idp + 1)])
  df_test_inner <- data.frame(df_train_all[-id_train_inner, c(1, idp + 1)])

  # fit model
  fit_lm <- lm(y ~ ., df_train_inner)

  # make predictions
  yhat_insample <- predict(fit_lm, df_train_inner)
  y_insample <- df_train_inner$y
  yhat <- predict(fit_lm, df_test_inner)
  y <- df_test_inner$y

  prediction_list <- list("yhat" = yhat, "y" = y, "fit_lm" = fit_lm, "yhat_insample" = yhat_insample, "y_insample" = y_insample)

  return(prediction_list)
}

get_best_feat_lin <- function(df_train, prop_train_inner, extract_n_feat, age_train, model) {

  n_feat <- dim(df_train)[2] - 1
  corr <- numeric(n_feat)

  id_train_inner <- get_idp_train_inner(df_train, prop_train_inner)

  for (idp in 1:n_feat) {
    print(paste("IDP number: ", idp))

    if (model == "linear") {

      # run univariate predictions
      prediction_list <- run_univ_pred_lm(df_train, id_train_inner, idp)

      yhat <- prediction_list$yhat

      # test for significance (of betas or accuracy?)
      corr[idp] <- cor(yhat, df_train$y[-id_train_inner])

    } else if (model == "svc") {
      # run univariate predictions
      prediction_list <- run_univ_pred_svc(df_train, id_train_inner, idp, age_train, svc_config)

      # test for significance (of betas or accuracy?)
      corr[idp] <- cor(prediction_list$yhat, df_train$y[-id_train_inner])
    }
  }

  # select best features
  best_n_features <- extract_best_features(corr, extract_n_feat)

  return(best_n_features)
}


run_linear_model <- function(df_train_linear, id_train, df) {
  print("Fitting linear model...") 
  fit_lm <- lm(y ~ ., data = df_train_linear)
  lm_yhat <- numeric(length(df$y))
  lm_yhat[id_train] <- predict(fit_lm, newdata = df_train_linear)
  lm_yhat[-id_train] <- predict(fit_lm, newdata = df[-id_train, ])
  se_lm <- (lm_yhat - df$y)^2
  corr_lm_in <- cor(lm_yhat[id_train], df_train_linear$y)
  corr_lm_out <- cor(lm_yhat[-id_train], df[-id_train, ]$y)

  # linear_model <- list()
  # linear_model$fit_lm <- fit_lm
  # linear_model$se_lm <- se_lm
  # linear_model$corr_lm_in <- corr_lm_in
  # linear_model$corr_lm_out <- corr_lm_out
  # linear_model$lm_yhat <- lm_yhat

  # return(linear_model)

  return(list(fit_lm = fit_lm, se_lm = se_lm, corr_lm_in = corr_lm_in, corr_lm_out = corr_lm_out, lm_yhat = lm_yhat))
}

run_elastic_net_model <- function(idps, idps_linear, trait, id_train, age, model_age) {

  print("Fitting elastic net model...")
  trait_train <- trait[id_train]
  idps_linear_train <- idps_linear[id_train, ]
  age_train <- matrix(age[id_train])

  if (model_age == 1 ||  model_age == 3) {
    cvfit_glmnet <- cv.glmnet(cbind(idps_linear_train, age_train), trait_train)
  } else if (model_age == 0) {
    cvfit_glmnet <- cv.glmnet(idps_linear_train, trait_train)
  } else if (model_age == 2) {
    cvfit_glmnet <- cv.glmnet(cbind(rep(1, length(age_train)), age_train), trait_train)
  }

  print("this is the right function")
  enet_coefficients <- coef(cvfit_glmnet, s = "lambda.min") # note betas
  print(enet_coefficients)
  print(enet_coefficients != 0)
  #idx_non_zero_coeff <- determine_non_zero_coeff(enet_coefficients) # get best features
  # idx_non_zero_coeff <- which(enet_coefficients != 0)
  # print(paste0("Number of non-zero features: ", length(idx_non_zero_coeff)))

  enet_yhat <- numeric(length(df$y))
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

  se_enet <- (enet_yhat - df$y)^2
  corr_enet_in <- cor(enet_yhat[id_train], trait_train)
  corr_enet_out <- cor(enet_yhat[-id_train], df[-id_train, ]$y)


  # elastic_net_model <- list()
  # elastic_net_model$enet_coefficients <- enet_coefficients
  # elastic_net_model$cvfit_glmnet <- cvfit_glmnet
  # elastic_net_model$enet_yhat <- enet_yhat
  # elastic_net_model$se_enet <- se_enet
  # elastic_net_model$corr_enet_in <- corr_enet_in
  # elastic_net_model$corr_enet_out <- corr_enet_out


  # return(elastic_net_model)

  return(list(enet_coefficients=enet_coefficients, cvfit_glmnet=cvfit_glmnet, enet_yhat=enet_yhat, se_enet=se_enet, corr_enet_in=corr_enet_in, corr_enet_out=corr_enet_out))

}

fit_svc_model <- function (best_features_svc, df_all_train, df_all, df_train_svc, cov, prof, taper, id_train, age, age_train, model_age, df_svc, prop_train_inner=0.9) {

  # configure SVC
  svc_config <- SVC_mle_control(cov.name = c(cov),
                                profileLik = prof,
                                tapering = taper,
                                hessian = TRUE)

  # fit svc model
  print("Fitting SVC")

  # to run penalized MLE (not currently working), call:
  # e.g., 'penalized_svc(idps_train, trait_train, age_train, svc_config)'
  #xnam <- paste("x.", 1:n_feat, sep = "")
  xnam <- paste("x.", best_features_svc, sep = "")
  fmla <- as.formula(paste("y ~ ", paste(xnam, collapse = "+")))
  if (model_age == 0) {
    fmla_mi <-  as.formula(paste("y ~ ", paste(xnam, collapse = "+"), " - 1"))
    fit_svc <- SVC_mle(formula = fmla, data = df_train_svc, locs = age_train, control = svc_config, RE_formula = fmla_mi)
  } else if (model_age == 1) {
    fit_svc <- SVC_mle(formula = fmla, data = df_train_svc, locs = age_train, control = svc_config)
  } else if (model_age == 2) {
    fmla <- as.formula("y ~ age")
    fit_svc <- SVC_mle(formula = fmla, data = df_train_svc, locs = age_train, control = svc_config)
  } else if (model_age == 3) {
    best_features_linear_add <- get_best_feat_lin(df_all_train, prop_train_inner, n_feat, age_train, "linear")
    fmla_vary <- fmla
    xnam <- paste("x.", c(best_features_svc, best_features_linear_add), sep = "")
    fmla_fix <-  as.formula(paste("y ~ ", paste(xnam, collapse = "+")))
    # add the columns to our training dataframe
    column_names <- paste0("x.", best_features_linear_add)
    extracted_columns <- df_all_train[, column_names]
    df_train_svc <- cbind(df_train_svc, extracted_columns)
    extracted_columns_df <- df_all[, column_names]
    #df <- cbind(df_all, extracted_columns_df)
    df_svc <- cbind(df_all, extracted_columns_df)
    fit_svc <- SVC_mle(formula = fmla_fix, data = df_train_svc, locs = age_train, control = svc_config, RE_formula = fmla_vary)
  }

  # run predictions (for train and test at same time)
  df_svc_pred <- predict(fit_svc, newdata = df_svc, newlocs = age, control = svc_config)

  # calculate squared errors
  se_svc <- (df_svc_pred$y.pred - df_svc$y)^2
  corr_svc_in <- cor(df_svc_pred$y.pred[id_train], df_svc$y[id_train])
  corr_svc_out <- cor(df_svc_pred$y.pred[-id_train], df_svc$y[-id_train])


  return(list(df_svc_pred = df_svc_pred, se_svc = se_svc, corr_svc_in = corr_svc_in, corr_svc_out = corr_svc_out))

}

determine_best_svc_features <- function(trait_id, n_feat, n_sub, perc_train, prof, tap, cov, ica) {

  # determine best features (svc)
  print("Getting best svc features")
  if (ica == 3 || ica == 5) {
    best_features_svc <- 1:n_feat
  } else {
    best_features_svc <- get_best_feat_svc(trait_id, n_feat, n_sub, perc_train, prof, tap, cov, n_feat, ica)
  }

  return(best_features_svc)

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

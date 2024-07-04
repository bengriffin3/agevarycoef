get_best_feat_lin <- function(df_train, prop_train_inner, extract_n_feat, age_train, model, svc_config) {

  n_feat_lin <- dim(df_train)[2] - 1
  corr <- numeric(n_feat_lin)

  id_train_inner <- get_idp_train_inner(df_train, prop_train_inner)

  for (idp in 1:n_feat_lin) {
    print(paste("IDP number: ", idp))

    if (model == "linear") {

      # run univariate predictions
      prediction_list <- run_univ_pred_lm(df_train, id_train_inner, idp)

      y_hat <- prediction_list$yhat

      # test for significance (of betas or accuracy?)
      corr[idp] <- stats::cor(y_hat, df_train$y[-id_train_inner])

    } else if (model == "svc") {
      # run univariate predictions
      prediction_list <- run_univ_pred_svc(df_train, id_train_inner, idp, age_train, svc_config)

      # test for significance (of betas or accuracy?)
      corr[idp] <- stats::cor(prediction_list$yhat, df_train$y[-id_train_inner])
    }
  }

  # select best features
  best_n_features <- extract_best_features(corr, extract_n_feat)

  return(best_n_features)
}


extract_best_features <- function(corr, n) {

  best_features <- sort(corr, index.return = TRUE, decreasing = TRUE)
  best_features_idx <- best_features$ix
  best_n_features <- sort(best_features_idx[1:n])

  return(best_n_features)
}

get_best_feat_svc <- function(dir, trait_id, n_feat, n_sub, perc_train, prof, tap, cov, n, ica) {
  #dir <- '/gpfs3/well/win-fmrib-analysis/users/psz102/age_varying_coefficients/results/univariate/all_idps'
  dir <- '/Users/bengriffin/OneDrive - Aarhus Universitet/Dokumenter/git_repos/vary_coef_bg/results/univariate/all_idps'
  perc_train <- 90 # might want to update this to 90
  n_idps <- 1436
  n_sub <- 20000
  n_feat <- 1436
  corr <- rep(NA, n_idps)


  #for (idp in 1:n_idps) {
  for (idp in 1:5) {
    print(idp)
    pred <- try(load(sprintf("%s/yhat_t_%i_n_f_%i_f_1436_n_%i_p_%i_l_%i_ta_%i_c_%s_idp_%i_ica_%i_mi_0.RData",
                             dir, trait_id, n_feat, n_sub, perc_train, prof, tap, cov, idp, ica)))
    if (inherits(pred, "try-error")) {
      #error handling code, maybe just skip this iteration using
      next
    }
    # pred <- mget(load(sprintf("%s/yhat_t_%i_f_1436_n_%i_p_%i_l_%i_ta_%i_c_%s_idp_%i_ica_%i.RData",
    #              dir, trait_id, n_sub, perc_train, prof, tap, cov, idp, ica), envir=(NE. <- new.env())), envir=NE.)

    # y <- extractorRData(sprintf("%s/yhat_t_%i_n_f_%i_f_1436_n_%i_p_%i_l_%i_ta_%i_c_%s_idp_%i_ica_%i_mi_0.RData",
    #                             dir, trait_id, n_feat, n_sub, perc_train, prof, tap, cov, idp, ica), 'y_svc')
    #
    # yhat <- extractorRData(sprintf("%s/yhat_t_%i_n_f_%i_f_1436_n_%i_p_%i_l_%i_ta_%i_c_%s_idp_%i_ica_%i_mi_0.RData",
    #                                dir, trait_id, n_feat, n_sub, perc_train, prof, tap, cov, idp, ica), 'yhat_svc')

    corr_svc <- 0
    yhat <- extractorRData(sprintf("%s/yhat_t_%i_n_f_%i_f_1436_n_%i_p_%i_l_%i_ta_%i_c_%s_idp_%i_ica_%i_mi_0.RData",
                                   dir, trait_id, n_feat, n_sub, perc_train, prof, tap, cov, idp, ica), 'corr_svc')

    #corr[idp] <- stats::cor(y, yhat)
    corr[idp] <- corr_svc
    #corr[idp] <- stats::cor(pred$y, pred$yhat)
  }
  # get idx of best 'n' features
  idx_best_features <- extract_best_features(corr, n)

  return(idx_best_features)
}

run_univ_pred_lm <- function(df_train_all, id_train_inner, idp) {

  # get training subjects dataframe
  df_train_inner <- data.frame(df_train_all[id_train_inner, c(1, idp + 1)])
  df_test_inner <- data.frame(df_train_all[-id_train_inner, c(1, idp + 1)])

  # fit model
  fit_lm <- stats::lm(y ~ ., df_train_inner)

  # make predictions
  yhat_insample <- stats::predict(fit_lm, df_train_inner)
  y_insample <- df_train_inner$y
  yhat <- stats::predict(fit_lm, df_test_inner)
  y <- df_test_inner$y

  prediction_list <- list("yhat" = yhat, "y" = y, "fit_lm" = fit_lm, "yhat_insample" = yhat_insample, "y_insample" = y_insample)

  return(prediction_list)
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
    fit_svc <- varycoef::SVC_mle(formula = y ~ idp, data = df_train_inner, locs = df_train_inner$age_train_inner)
  } else if (model_vary_intercept == 0) {
    #fit_svc <- varycoef::SVC_mle(formula = y ~ idp - 1, data = df_train_inner, locs = df_train_inner$age_train_inner)
    fit_svc <- varycoef::SVC_mle(formula = y ~ idp, data = df_train_inner, locs = df_train_inner$age_train_inner, RE_formula = y ~ idp - 1)
  }

  # make prediction
  #df_svc_pred <- varycoef::predict(fit_svc, newdata = df_test_inner, newlocs = df_test_inner$age_test_inner, control = svc_config)
  #df_svc_pred_insample <- varycoef::predict(fit_svc, newdata = df_train_inner, newlocs = df_train_inner$age_train_inner, control = svc_config)
  df_svc_pred <- predict.SVC_mle(fit_svc, newdata = df_test_inner, newlocs = df_test_inner$age_test_inner, control = svc_config)
  df_svc_pred_insample <- predict.SVC_mle(fit_svc, newdata = df_train_inner, newlocs = df_train_inner$age_train_inner, control = svc_config)

  # note down predictions
  yhat <- df_svc_pred$y.pred
  y <- df_train$y[-id_train_inner]
  yhat_insample <- df_svc_pred_insample$y.pred
  y_insample <- df_train$y[id_train_inner]

  # make prediction
  #df_svc_pred <- varycoef::SVC_mle.predict(fit_svc, newlocs = df_test_inner$age_test_inner)


  # create list of prediction and target
  prediction_list <- list("yhat" = yhat, "y" = y, "fit_svc" = fit_svc, "yhat_insample" = yhat_insample, "y_insample" = y_insample)

  return(prediction_list)
}


calculate_stats <- function(yhat, df, trait_train, id_train) {
  stats <- list()
  stats$se <- (yhat - df$y)^2
  stats$corr_in <- stats::cor(yhat[id_train], trait_train)
  stats$corr_out <- stats::cor(yhat[-id_train], df[-id_train, ]$y)

  return(stats)
}


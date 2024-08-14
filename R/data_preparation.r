library(R.matlab)
library(varycoef)
library(optparse)
library(fastICA)
library(CCA)
library(logger)
library(glmnet)



load_prediction_variables <- function(n_subjects, proj_dir) {
  # Load the mat file
  if (n_subjects <= 5000) {
    print(n_subjects)
    data <- readMat(paste0(proj_dir, "/data/prediction_vars_5000.mat"))
  } else if (n_subjects > 5000 && n_subjects <= 20000) {
    data <- readMat(paste0(proj_dir, "/data/prediction_vars_20000.mat"))
  } else if (n_subjects > 20000) {
    data <- readMat(paste0(proj_dir, "/data/prediction_vars.mat"))
  }
  # note index of best 30 traits
  idx30 <- readMat(paste0(proj_dir, "/data/idx_30.mat"))
  data$idx30 <- idx30$idx

  return(data)
}

load_cog_variables <- function(n_subjects, proj_dir) {
  # load cognitive vars
  if (n_subjects <= 1000) {
    data <- readMat(paste0(proj_dir, "/data/cognitive_vars_1000.mat"))
  } else if (n_subjects <= 20000) {
    data <- readMat(paste0(proj_dir, "/data/cognitive_vars_20000.mat"))
  } else {
    data <- readMat(paste0(proj_dir, "/data/cognitive_vars.mat"))
  }
  return(data)
}

get_top_30_cog_vars <- function(data, n_sub, proj_dir) {

  cog_data <- load_cog_variables(n_sub, proj_dir)
  y <- cog_data$vars.cog
  y_top_30 <- y[1:n_sub, data$idx30 + 1]
  return(y_top_30)
}

remove_nan_sub <- function(idps, trait, age, conf) {

  idx_keep <- complete.cases(idps, trait, age, conf)

  idps <- idps[idx_keep, ]
  if (is.array(trait)) {
    trait <- trait[idx_keep, ]
  } else {
    trait <- trait[idx_keep]
  }
  age <- age[idx_keep]
  conf <- conf[idx_keep, ]

  list_vars <- mget(c("idps", "trait", "age", "conf"))

  log_warn(paste0("Number of subjects removed: ", sum(!idx_keep)))
  log_info(paste0("Number of subjects remaining: ", dim(idps)[1]))
  log_info(paste0("Number of features/IDPs: ", dim(idps)[2]))


  return(list_vars)
}


get_training_sample_id <- function(trait, prop_train) {
  # sample training data
  if (is.vector(trait)) {
    id_train <- sort(sample(length(trait), length(trait) * prop_train))
  } else if (is.matrix(trait)) {
    id_train <- sort(sample(dim(trait)[1], dim(trait)[1] * prop_train))
  }
  return(id_train)
}

get_traits <- function(data, trait_id, proj_dir, n_sub) {

  # load traits
  if (trait_id == 1 || trait_id == 2 || trait_id == 3) {
    trait <- data$vars.cog.best.3[1:n_sub, trait_id]
  } else {
    # if performing PC or CCA, let's use the top 30 cog traits
    trait <- get_top_30_cog_vars(data, n_sub, proj_dir)
  }

  return(trait)
}


perform_pca <- function(idps, n_feat) {

  pca_result <- prcomp(idps, center = TRUE, scale. = TRUE)
  idps <- pca_result$x[, 1:n_feat]
  return(idps)
}



deconfoundPhen <- function(data, conf, beta = NULL, my = NULL) {
  if (is.null(beta)) {
    # De-mean the response
    my <- mean(data)
    data <- data - my
    # Compute the coefficients (betaY) for the confounds
    beta <- solve(t(conf) %*% conf + 0.00001 * diag(ncol(conf))) %*% t(conf) %*% data
  }
  # Compute the residuals
  residuals <- data - conf %*% beta

  # organise into list
  my_list <- list("beta" = beta, "my" = my, "residuals" = residuals)

  return(my_list)
}

prepare_options <- function(tap, prof, perc_train) {
  taper <- if (tap == 0) NULL else tap
  prof <- ifelse(prof == 1, TRUE, FALSE)
  prop_train <- perc_train / 100


  return(list(taper, prof, prop_train))
}


organise_data <- function(trait, idps, id_train, age) {
  df <- data.frame(y = trait, x = idps)
  df_train <- df[id_train, ]
  df_test <- df[-id_train, ]

  # determine best features using univariate predictions
  age_train <- matrix(age[id_train])

  return(list(df, df_train, df_test, age_train))
}

prepare_age_data <- function(df_train_linear, df, age_train, model_age, age) {
  if (model_age == 1 ||  model_age == 3) {
    df_train_linear$age <- as.vector(age_train)
    df$age <- as.vector(age)
  } else if (model_age == 2) {
    # keep only y (remove all IDPs)
    df_train_linear <- df_train_linear[c("y")]
    df <- df[c("y")]

    # add age which will be the only predictor
    df_train_linear$age <- as.vector(age_train)
    df$age <- as.vector(age)
  }

  return(list(df_train_linear = df_train_linear, df = df))
}

prepare_linear_data <- function(df_train, idps, best_features_linear) {
  df_train_linear <- df_train[, c(1, best_features_linear + 1)]
  idps_linear <- idps[, best_features_linear]

  return(list(df_train_linear = df_train_linear, idps_linear = idps_linear))
}

prepare_svc_data <- function(df_train, best_features_svc, id_train, idps, age, df) {


  df_train_svc <- df_train[, c(1, best_features_svc + 1)]
  df_train_svc$age <- age[id_train]
  idps_svc <- idps[, best_features_svc]
  idps_svc_train <- idps_svc[id_train, ]

  df_svc <- df[, c(1, best_features_svc + 1)]
  df_svc$age <- age

  # svc_data <- list()
  # svc_data$df_train_svc <- df_train_svc
  # svc_data$idps_svc <- idps_svc
  # svc_data$idps_svc_train <- idps_svc_train

  # return(svc_data)

  return(list(df_train_svc = df_train_svc, idps_svc = idps_svc, idps_svc_train = idps_svc_train, df_svc = df_svc))

}

unpack_data <- function(proj_dir, data, trait_id, n_sub) {

  # Extract the variables from the mat file and subsample
  idps <- data$structural.idps[1:n_sub, ]
  age <- data$age1[1:n_sub]
  trait <- get_traits(data, trait_id, proj_dir, n_sub)
  conf <- data$conf1[1:n_sub, ]
  conf_names <- data$conf.names


  return(list(idps, age, trait, conf, conf_names))
}

# bootstrap_pca <- function(cognitive_traits, idx_trait_keep, idps) {

#   n_keep <- length(idx_trait_keep)

#   # randomly selected 30 traits (with replacement)
#   y <- cognitive_traits[, idx_trait_keep + 1] # + 1 for python indexing
#   idx_bootstrap <- sort(sample(seq_len(n_keep), n_keep))
#   y <- y[, idx_bootstrap]

#   # scale y
#   y <- scale(y)

#   # remove subjects with NaNs
#   subject_keep <- !apply(is.na(y), 1, any)
#   y <- y[subject_keep, ] # remove rows where any NaNs

#   # perform PCA
#   pca <- prcomp(y, scale. = TRUE, center = TRUE)
#   principal_components <- pca$x
#   y <- principal_components[, 1] # take the first PC

#   x <- idps[subject_keep, ]

#   my_list <- list(x, y, subject_keep)

#   return(my_list)
# }


get_x_pcs <- function(x, n_pc_comp) {
  # perform PCA including de-meaning and setting unit variance
  pca_result_x <- prcomp(x, center = TRUE, scale. = TRUE)

  # note down the PCs
  x_pca <- pca_result_x$x[, 1:n_pc_comp]
  return(x_pca)
}

# de_mean_trait <- function(trait) {
#   if (length(dim(trait)) == 1) {
#     my <- mean(trait)
#     trait <- trait - my
#   } else if (length(dim(trait)) == 2) {
#     my <- colMeans(trait)
#     trait <- trait - rep(my, rep.int(nrow(trait), ncol(trait)))
#   }
#   return(trait)
# }
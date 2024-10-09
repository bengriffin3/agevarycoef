library(R.matlab)
library(gsubfn)
library(logger)

load_prediction_variables <- function(n_subjects, proj_dir) {
  # Load the mat file
  if (n_subjects <= 5000) {
    print(n_subjects)
    data <- R.matlab::readMat(paste0(proj_dir, "/data/prediction_vars_5000_NEW.mat"))
  } else if (n_subjects > 5000 && n_subjects <= 20000) {
    data <- R.matlab::readMat(paste0(proj_dir, "/data/prediction_vars_20000_NEW.mat"))
  } else if (n_subjects > 20000) {
    data <- R.matlab::readMat(paste0(proj_dir, "/data/prediction_vars_NEW.mat"))
  }
  # note index of best 30 traits
  idx30 <- R.matlab::readMat(paste0(proj_dir, "/data/idx_30.mat"))
  data$idx30 <- idx30$idx

  return(data)
}

load_cog_variables <- function(n_subjects, proj_dir) {
  # load cognitive vars
  if (n_subjects <= 1000) {
    data <- R.matlab::readMat(paste0(proj_dir, "/data/cognitive_vars_1000.mat"))
  } else if (n_subjects <= 20000) {
    data <- R.matlab::readMat(paste0(proj_dir, "/data/cognitive_vars_20000.mat"))
  } else {
    data <- R.matlab::readMat(paste0(proj_dir, "/data/cognitive_vars.mat"))
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

  idx_keep <- stats::complete.cases(age, conf) #idx_keep <- stats::complete.cases(idps, trait, age, conf)

  if (is.array(idps)) {
    idps <- idps[idx_keep, ]
  } else {
    idps <- idps[idx_keep]
  }

  if (is.array(trait)) {
    trait <- trait[idx_keep, ]
  } else {
    trait <- trait[idx_keep]
  }
  age <- age[idx_keep]
  conf <- conf[idx_keep, ]

  list_vars <- mget(c("idps", "trait", "age", "conf"))

  logger::log_warn(paste0("Number of subjects removed due to missing age: ", sum(!idx_keep)))
  logger::log_info(paste0("Number of subjects remaining: ", dim(idps)[1]))
  logger::log_info(paste0("Number of features/IDPs: ", dim(idps)[2]))


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

#' Get traits based on the specified trait_id.
#'
#' @name get_traits
#' @param data The input data containing traits.
#' @param trait_id The identifier for the trait to load.
#' @param proj_dir The project directory path.
#' @param n_sub The number of subjects.
#' @return A vector of traits.
#' @importFrom utils globalVariables
utils::globalVariables("simulated_trait")
get_traits <- function(data, trait_id, proj_dir, n_sub) {

  # load traits
  if (trait_id == 1 || trait_id == 2 || trait_id == 3) {
    idx_trait <- data$idx30[trait_id]
    trait <- data$vars.cog[1:n_sub, idx_trait]
  } else if (trait_id == 0 || trait_id == 999) {
    # if performing PC or CCA, let's use the top 30 cog traits
    trait <- get_top_30_cog_vars(data, n_sub, proj_dir)
  } else if (trait_id == 1000) {
    # load simulated trait
    load("/gpfs3/well/win-fmrib-analysis/users/psz102/age_varying_coefficients/data/simulated_trait.Rdata")
    trait <- simulated_trait[1:n_sub]
  }

  return(trait)
}


perform_pca <- function(idps, n_feat) {

  pca_result <- stats::prcomp(idps, center = TRUE, scale. = TRUE)
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


prepare_age_data <- function(df_linear, age, model_age) {
  if (model_age == 1) {
    # df_train_linear$age <- as.vector(age_train)
    df_linear$age <- as.vector(age)
  } else if (model_age == 2) {
    # keep only y (remove all IDPs)
    # df_train_linear <- df_train_linear[c("y")]
    df_linear <- df_linear[c("y")]

    # add age which will be the only predictor
    # df_train_linear$age <- as.vector(age_train)
    df_linear$age <- as.vector(age)
  } else if (model_age == 3) {
    stop("Error: model_age == 3 is only appropriate with SVC.")
  }

  return(df_linear)
}

prepare_linear_data <- function(df, idps, best_features) {
  df_linear <- df[, c(1, best_features + 1)]
  idps_linear <- idps[, best_features]

  return(list(df_linear = df_linear, idps_linear = idps_linear))
}


prepare_svc_data <- function(df_train, best_features_svc, id_train, idps, age, df) {


  df_train_svc <- df_train[, c(1, best_features_svc + 1)]
  df_train_svc$age <- age[id_train]
  idps_all <- idps[, best_features_svc]

  df_svc <- df[, c(1, best_features_svc + 1)]
  df_svc$age <- age

  # svc_data <- list()
  # svc_data$df_train_svc <- df_train_svc
  # svc_data$idps_svc <- idps_svc
  # svc_data$idps_svc_train <- idps_svc_train

  # return(svc_data)

  return(list(df_train_svc = df_train_svc, idps_all = idps_all, df_svc = df_svc))

}

unpack_data <- function(proj_dir, data, trait_id, n_sub) {

  # Extract the variables from the mat file and subsample
  #idps <- data$structural.idps[1:n_sub, ]
  idps <- data$IDPs1KEEP[1:n_sub, ]
  age <- data$age1[1:n_sub]
  trait <- get_traits(data, trait_id, proj_dir, n_sub)
  conf <- data$conf1[1:n_sub, ]
  conf_names <- data$conf.names


  return(list(idps, age, trait, conf, conf_names))
}


get_x_pcs <- function(x, n_pc_comp) {
  # perform PCA including de-meaning and setting unit variance
  pca_result_x <- stats::prcomp(x, center = TRUE, scale. = TRUE)

  # note down the PCs
  x_pca <- pca_result_x$x[, 1:n_pc_comp]
  return(x_pca)
}


add_idp_age_interaction_term <- function(df_all, age) {

  # Create new features by multiplying 'age' with each feature (excluding 'y')
  # List all feature columns except 'y'
  feature_columns <- names(df_all)[!names(df_all) %in% c("y", "age", "age_squared")]

  # Multiply 'age' with each feature column and add to the dataframe
  for (feature in feature_columns) {
    new_feature_name <- paste("age_", feature, sep = "")
    df_all[[new_feature_name]] <- df_all$age * df_all[[feature]]
  }


  return(df_all)
}

perform_pca_after_feature_engineering <- function(df_all) {

  features_for_pca <- names(df_all)[!names(df_all) %in% c("y", "age", "age_squared")]

  pca_result_x <- stats::prcomp(df_all[, features_for_pca], center = TRUE, scale. = TRUE)  # Scaling is recommended
  pca_components <- as.data.frame(pca_result_x$x[, 1:600])
  df_all <- df_all[, !(names(df_all) %in% features_for_pca)]
  df_all <- cbind(df_all, pca_components)

  return(df_all)
}

prepare_boost_data <- function(df, id_train, trait_all, age_all, model_age, n_feat, best_features=0) {

  # things I've previously tried
  #idps <- perform_pca_using_train(idps, id_train, n_PCs = 1000)
  #df_all <- data.frame(y = trait, x = idps[[1]]) # perform PCA before feature engineering
  #df_all <- data.frame(y = trait, x = age, x1 = age^2)
  #df_all$age_squared <- age^2 # add age^2 as a feature
  #df_all <- add_idp_age_interaction_term(df_all, age) # add interaction terms  
  #df_all <- perform_pca_after_feature_engineering(df_all) # perform PCA after feature engineering


  # df_all$age_cubed <- age^3 # add age as a feature
  # df_all$age_x1 <- age * df_all$x.1 # add age as a feature
  # df_all$age_squared_x1 <- age^2 * df_all$x.1 # add age as a feature
  # df_all$age_cubed_x1 <- age^3 * df_all$x.1 # add age as a feature

  # # Separate the target column (y) and the features
  # y_column <- df_all$y
  # features <- df_all[, !(names(df_all) %in% "y")]
  # scaled_features <- scale(features) # Scale the features (mean = 0, standard deviation = 1)  
  # scaled_features <- as.data.frame(scaled_features) # Convert the scaled features back to a dataframe
  # df_scaled <- cbind(y = y_column, scaled_features) # Combine the y column back with the scaled features
  # df_scaled <- as.data.frame(df_scaled) # Ensure df_scaled is still a dataframe
  # df_all <- df_scale



  # if model_age == 0, we just leave it as it is, then 1 = model age, 2 = only model age, 3 = fixed with varying?
  if (model_age == 1) {
    df$age <- (age_all - mean(age_all)) / stats::sd(age_all) # add age as a feature
  } else if (model_age == 2) {
    df <- data.frame(y = trait_all, x = age_all)
  } else if (model_age == 3) {
    stop("Error: model_age == 3 is only appopriate with SVC model.")
  }

  #   # prepare data for boosting model
  # df_all_save <- df_all
  # best_features <- c(1)
  # df_all <- df_all[, c(1, best_features + 1)]

  df_train_x <- as.matrix(df[id_train, ][, !names(df[id_train, ]) %in% c("y")])
  df_train_y <- as.matrix(df[id_train, ]$y)

  df_test_x <- as.matrix(df[-id_train, ][, !names(df[-id_train, ]) %in% c("y")])
  df_test_y <- as.matrix(df[-id_train, ]$y)

  return(list(df_train_x, df_train_y, df_test_x, df_test_y))
}

determine_best_boost_features <- function(df_all, id_train, n_feat, model_age) {

  best_features <- 0
  if (n_feat < 1436 && model_age != 2) {

    # for now we load the best features that have been determined elsewhere - in the future this will have to be done within fold
    best_features <- c(1) # load('')

    #best_features <- get_best_boost_features(df_al, id_train, n_feat)
  }

  return(best_features)
}

simulate_linear_trait <- function(age, features) {
  # Set seed for reproducibility
  set.seed(5)

  # Number of samples
  n <- length(age)
  n_feat <- dim(features)[2]

  features_scaled <- scale(features)
  age_scaled <- scale(age)


  # Impute missing values with the mean of each column
  # Calculate the mean of each column, excluding NAs
  feature_means <- colMeans(features_scaled, na.rm = TRUE)
  age_mean <- mean(age, na.rm = TRUE)

  # Replace NA values with column means
  features_imputed <- features_scaled
  for (i in 1:n_feat) {
    features_imputed[is.na(features_scaled[, i]), i] <- feature_means[i]
  }

  # Replace NA values with column means
  age_imputed <- age_scaled
  age_imputed[is.na(age)] <- age_mean


  # Assign coefficients
  coeff_age <- 0.6
  #coeff_X <- stats::runif(n_feat, min = -0.01, max = 0.0001)  # Random coefficients between 0.05 and 0.3
  coeff_X <- stats::runif(n_feat, min = 0, max = 0.0)  # Random coefficients between 0.05 and 0.3

  # Simulate noise
  noise <- stats::rnorm(n, mean = 0, sd = 1)

  trait <- coeff_age * age_imputed + features_imputed %*% coeff_X + noise

  return(trait)
}

impute_and_filter_data <- function(idps, trait, age, conf, max_nans = 1000) {

  # If univariate prediction then just impute the variable
  if (is.vector(idps)) {
    keep_mask <- rep(TRUE, length(idps))
    idps_filtered <- idps[keep_mask]
    # Calculate the means for each IDP, ignoring NA values
    idps_means <- mean(idps_filtered, na.rm = TRUE)
  } else {
    # Compute the number of NaNs for each subject
    nan_counts <- rowSums(is.na(idps))
    # Create a mask for subjects with <= max_nans NaNs
    keep_mask <- nan_counts <= max_nans
    # filter idps and calc mean of filtered
    idps_filtered <- idps[keep_mask, ]
    # Calculate the means for each IDP, ignoring NA values
    idps_means <- colMeans(idps_filtered, na.rm = TRUE)
  }

  # Filter out subjects in trait and age based on the mask
  if (is.vector(trait)) {
    trait_filtered <- trait[keep_mask]
  } else if (is.matrix(trait)) {
    trait_filtered <- trait[keep_mask, ]
  } else {
    stop("Trait must be either a 1D vector or a 2D matrix.")
  }

  age_filtered <- age[keep_mask]
  conf_filtered <- conf[keep_mask, ]

  # Impute the missing values in the filtered dataset
  idps_imputed <- impute_with_mean(idps_filtered, idps_means)

  #logger::log_warn(paste0("Number of subjects removed due to missing idp/trait: ", sum(!keep_mask)))

  return(list(idps = idps_imputed, trait = trait_filtered, age = age_filtered, conf = conf_filtered))
}


# Impute missing values with the means
impute_with_mean <- function(data, means) {

  if (length(means) == 1) {
    data[is.na(data)] <- means

  } else {
    for (i in seq_along(means)) {
      data[is.na(data[, i]), i] <- means[i]
    }
  }

  return(data)
}
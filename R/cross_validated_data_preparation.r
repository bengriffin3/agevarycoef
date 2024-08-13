library(R.matlab)
library(varycoef)
library(optparse)
library(fastICA)
library(CCA)
library(logger)
library(glmnet)


pre_process_data_cross_validated <- function(idps, trait, trait_id, age, conf, conf_names, prop_train, id_train, ica, n_feat) {

  # Remove subjects with NaNs in idps and age
  idps_with_nans <- idps
  list_of_vars <- remove_nan_sub(idps, trait, age, conf)
  idps <- list_of_vars$idps
  trait <- list_of_vars$trait
  age <- list_of_vars$age
  conf <- list_of_vars$conf

  # note training subjects
  id_train <- get_training_sample_id(trait, prop_train)

  # Select confounds
  conf <- confound_selection_using_train(conf, conf_names, id_train)

  # impute idps (remove NaNs instead as this is better)
  # idps <- impute_data_using_train(idps, id_train)

  # scale idps
  idps <- scale_data_using_train(idps, id_train)

  # perform CCA or PC
  if (trait_id == 999) {
    cca_list <- perform_cca_using_train(idps, id_train, trait)
    trait <- cca_list[[1]]
    #cca_object_idps_trait <- cca_list[[2]]
    #pca_object_idps_for_cca <- cca_list[[3]]
    #pca_object_trait <- cca_list[[4]] # we have already performed PCA so ignore this
  } else if (trait_id == 0) {
    # perform pca
    trait_pca_list <- perform_pca_using_train(trait, id_train, 1)
    trait_pca <- trait_pca_list[[1]]
    #pca_object_trait <- trait_pca_list[[2]]
    trait <- trait_pca
  }

  # perform ICA on idps
  df_all <- data.frame(y = trait, x = idps)
  df_all_train <- df_all[id_train, ]
  if (ica == 1) {
    # choose whether to upweight idps (this seems to do very very little)
    upweight <- 0
    if (upweight == 1) {
      idps_upweight <- upweight_idps_using_train(idps, id_train, trait)
      idps <- idps_upweight
    }
    idps_ica_list <- perform_ica_using_train(idps, id_train, n_feat)
    idps_ica <- idps_ica_list[[1]]
    idps <- idps_ica
    #ica_object_idps <- idps_ica_list[[2]]
  } else if (ica == 4) {
    idps_pca_list <- perform_pca_using_train(idps, id_train)
    idps <- idps_pca_list[[1]]
  } else if (ica == 5) {
    load("/gpfs3/well/win-fmrib-analysis/users/psz102/age_varying_coefficients/data/idx_1436_keep_16_ordered.Rdata")
    # select best idps
    idps_best_16 <- idps[, idx_1436_keep_16_ordered]
    idps <- idps_best_16[, 1:n_feat]

  }

  # de-mean target
  trait <- de_mean_trait_using_train(trait, id_train)

  # deconfound data
  idps <- deconfound_data_using_train(idps, conf, id_train)
  trait <- deconfound_data_using_train(trait, conf, id_train)



#   pre_processed_data <- list()
# #  pre_processed_data$idps_with_nans <- idps_with_nans
#   pre_processed_data$df_all_train <- df_all_train
#   pre_processed_data$idps <- idps
#   pre_processed_data$trait <- trait
#   pre_processed_data$df_all <- df_all
#   pre_processed_data$id_train <- id_train
#   pre_processed_data$n_feat <- n_feat
#   pre_processed_data$age <- age

#   return(pre_processed_data)

  return(list(df_all_train, idps, trait, df_all, id_train, age, conf))
}


scale_data_using_train <- function(idps, id_train) {

  # note train and test idps
  idps_train <- idps[id_train, ]
  idps_test <- idps[-id_train, ]

  # Scale the training IDPs
  idps_train_scaled <- scale(idps_train)

  # Save the scaling parameters
  train_mean <- attr(idps_train_scaled, "scaled:center")
  train_sd <- attr(idps_train_scaled, "scaled:scale")

  # Scale the test IDPs using the training scaling parameters
  idps_test_scaled <- sweep(idps_test, 2, train_mean, "-")
  idps_test_scaled <- sweep(idps_test_scaled, 2, train_sd, "/")

  # remake idps variable
  idps[id_train, ] <- idps_train_scaled
  idps[-id_train, ] <- idps_test_scaled

  return(idps)
}

de_mean_trait_using_train <- function(trait, id_train) {

  if (is.vector(trait)) {
    # note train and test trait
    trait_train <- trait[id_train]
    trait_test <- trait[-id_train]

    # De-mean the trait
    trait_mean <- mean(trait_train)
    trait_train_demeaned <- trait_train - trait_mean
    trait_test_demeaned <- trait_test - trait_mean

    # Recreate the trait variable with de-meaned data
    trait[id_train] <- trait_train_demeaned
    trait[-id_train] <- trait_test_demeaned

    return(trait)
  }
  # If trait is a 2D array
  else if (is.matrix(trait)) {
    # De-mean each column separately
    trait_mean <- colMeans(trait)
    trait_demeaned <- sweep(trait, 2, trait_mean)

    return(trait_demeaned)
  }
  else {
    stop("Input trait must be either a vector or a matrix.")
  }
}

#### we no longer impute because it's mostly just subjects without
#### any recorded IDPs so better to remove than have 'mean' subjects
# impute_data_using_train <- function(idps, id_train) {

#   # note train and test idps
#   idps_train <- idps[id_train, ]
#   idps_test <- idps[-id_train, ]

#   # Compute the mean of each column in the training set, ignoring NA values
#   train_means <- colMeans(idps_train, na.rm = TRUE)

#   # Function to impute missing values with the provided means
#   impute_with_mean <- function(data, means) {
#     for (i in seq_along(means)) {
#       data[is.na(data[, i]), i] <- means[i]
#     }
#     return(data)
#   }

#   # Impute missing values in the training and test sets using the training means
#   idps_train_imputed <- impute_with_mean(idps_train, train_means)
#   idps_test_imputed <- impute_with_mean(idps_test, train_means)

#   # remake idps variable
#   idps[id_train, ] <- idps_train_imputed
#   idps[-id_train, ] <- idps_test_imputed

#   return(idps)
# }


confound_selection_using_train <- function(conf, conf_names, id_train) {

  ### first we do stuff not based on training/test split
  log_info(paste0("Number of total possible confounds: ", dim(conf)[2]))

  # first remove essential confounds
  conf_list <- remove_essential_confounds(conf, conf_names)
  conf_ess <- conf_list$conf_ess
  conf_remove_ess <- conf_list$conf_remove_ess
  conf_names_remove_ess <- conf_list$conf_names_remove_ess

  # now remove 'age' confounds
  conf_remove_age <- remove_confounds_by_name(conf_remove_ess, conf_names_remove_ess, "Age")

  # Dimensionality reduction
  conf_svd <- svd_reduce_conf_using_train(conf_remove_age, id_train)

  # Combine essential and reduced confounds
  conf <- cbind(conf_ess, conf_svd)

  log_info(paste0("Number of confounds remaining: ", dim(conf)[2]))

  return(conf)
}

# Function to perform PCA using SVD to dimensionality reduce confounds
svd_reduce_conf_using_train <- function(conf, id_train, variance_threshold = 0.80) {

  # de-mean confounds (important for SVD)
  conf <- de_mean_trait_using_train(conf, id_train)

  # note training and test confound data
  conf_train <- conf[id_train, ]
  conf_test <- conf[-id_train, ]

  # perform SVD
  svd_result <- svd(conf_train)

  # note cumulative explained variance of principal components
  explained_variance <- cumsum(svd_result$d^2) / sum(svd_result$d^2)

  # note how many components required to reach EV threshold
  num_components <- which(explained_variance > variance_threshold)[1]

  # note the new confounds based on PCs (left singular vectors * singular values)
  # or could do original data * right singular vectors
  conf_train_pc <- svd_result$u[, 1:num_components] %*% diag(svd_result$d[1:num_components])

  # perform same transformation to test data
  conf_test_pc <- conf_test %*% svd_result$v[, 1:num_components]

  # create reduced confound matrix
  conf_reduced <- matrix(NA, nrow = nrow(conf), ncol = num_components)
  conf_reduced[id_train, ] <- conf_train_pc
  conf_reduced[-id_train, ] <- conf_test_pc

  return(conf_reduced)
}


deconfound_data_using_train <- function(data, conf, id_train) {
  # Note train and test confounds
  conf_train <- conf[id_train, ]
  conf_test <- conf[-id_train, ]

  # Separate train and test data
  if (is.vector(data)) {
    data_train <- data[id_train]
    data_test <- data[-id_train]
  } else {
    data_train <- data[id_train, ]
    data_test <- data[-id_train, ]
  }

  # Deconfound data (training)
  my_list <- deconfoundPhen(data_train, conf_train)
  beta <- my_list$beta
  intercept_train <- my_list$my
  residuals_train <- my_list$residuals

  # Apply to test data
  if (is.vector(data)) {
    data_test <- data_test - intercept_train
    residuals_test <- data_test - conf_test %*% beta
  } else {
    data_test <- data_test - intercept_train
    residuals_test <- data_test - conf_test %*% beta
  }

  # Recombine the deconfounded data
  if (is.vector(data)) {
    data[id_train] <- residuals_train
    data[-id_train] <- residuals_test
  } else {
    data[id_train, ] <- residuals_train
    data[-id_train, ] <- residuals_test
  }

  return(data)
}

perform_ica_using_train <- function(idps, id_train, n_feat, upweight = 0) {

  n_ICs <- n_feat

  # # perform pca
  # idps_pca_list <- perform_pca_using_train(idps, id_train)
  # idps_pca <- idps_pca_list[[1]]
  # pca_object_idps <- idps_pca_list[[2]]


  # # split idps into train and test
  # idps_train <- idps_pca[id_train, ]
  # idps_test <- idps_pca[-id_train, ]

  idps_train <- idps[id_train, ]
  idps_test <- idps[-id_train, ]

  # now apply ICA to the PCA'd idps
  ica_object_idps <- fastICA(idps_train, n_ICs)
  idps_train_ica <- ica_object_idps$S

  # idps
  idps_test_ica <- idps_test %*% ica_object_idps$K %*% ica_object_idps$W

  # now create new idps based on the independent component data
  idps_ic <- matrix(numeric(), nrow = nrow(idps), ncol = n_feat)
  idps_ic[id_train, ] <- as.matrix(idps_train_ica)
  idps_ic[-id_train, ] <-  as.matrix(idps_test_ica)

  # final IDP data based on independent components
  idps <- idps_ic

  #ica_list <- list(idps, ica_object_idps, pca_object_idps)
  ica_list <- list(idps, ica_object_idps)

  return(ica_list)
}

perform_pca_using_train <- function(idps, id_train, n_PCs = 100) {

  # split idps into train and test
  idps_train <- idps[id_train, ]
  idps_test <- idps[-id_train, ]

  # perform pca on training data
  pca_result_x <- prcomp(idps_train, center = TRUE, scale. = TRUE)
  idps_train_pca <- pca_result_x$x[, 1:n_PCs]

  # now apply that same transformation to the test subjects
  idps_test_scale <- scale(idps_test, center = pca_result_x$center, scale = pca_result_x$scale)

  # apply PCA coefficients to data to project to PCs (first 100 PCs)
  idps_test_pca <- idps_test_scale %*% pca_result_x$rotation[, 1:n_PCs]

  # now create new idps based on the independent component data
  idps_pca <- matrix(numeric(), nrow = nrow(idps), ncol = n_PCs)
  idps_pca[id_train, ] <- as.matrix(idps_train_pca)
  idps_pca[-id_train, ] <-  as.matrix(idps_test_pca)

  idps <- idps_pca

  pca_list <- list(idps, pca_result_x)

  return(pca_list)
}

upweight_idps_using_train <- function(idps, id_train, trait) {

  # split idps and trait into train and test
  idps_train <- idps[id_train, ]
  idps_test <- idps[-id_train, ]
  trait_train <- trait[id_train]
  trait_test <- trait[-id_train]


  idp_weightings <- cor(idps_train, trait_train)^2

  # weight the features (UDPs)
  idps_weighted <- sweep(idps, 2, idp_weightings, `*`)

  idps <- idps_weighted

  return(idps)
}

perform_cca_using_train <- function(idps, id_train, trait) {

  # perform PCA on X first
  idps_pca_list <- perform_pca_using_train(idps, id_train)
  idps_pca <- idps_pca_list[[1]]
  pca_object_idps <- idps_pca_list[[2]]

  # scale target trait
  trait <- scale_data_using_train(trait, id_train)
  # trait_pca_list <- perform_pca_using_train(trait, id_train)
  # trait_pca <- trait_pca_list[[1]]
  # pca_object_trait <- trait_pca_list[[2]]
  # we have already performed PCA so don't need to perform it again here


  # split idps and trait into train and test
  idps_pca_train <- idps_pca[id_train, ]
  # idps_pca_test <- idps_pca[-id_train, ]
  trait_train <- trait[id_train, ]
  trait_test <- trait[-id_train, ]

  # perform CCA
  cc1_train <- cc(idps_pca_train, as.matrix(trait_train))

  # note coefficients
  y_coef <- cc1_train$ycoef[, 1]

  # Apply the CCA transformation to both the training and test data
  cca_train_scores <- as.matrix(trait_train) %*% y_coef
  cca_test_scores <- as.matrix(trait_test) %*% y_coef

  # note new trait
  new_trait <- numeric(dim(trait)[1])
  new_trait[id_train] <- cca_train_scores
  new_trait[-id_train] <- cca_test_scores

  #cca_list <- list(new_trait, cc1_train, pca_object_idps, pca_object_trait)
  cca_list <- list(new_trait, cc1_train, pca_object_idps)

  return(cca_list)

  # note first mode
  # cca_mode_1_train <- cc1_train$scores$yscores[, 1]  
  # trait <- cca_mode_1_train
}

get_idp_train_inner <- function(df_train, prop_train_inner) {

  n_train <- dim(df_train)[1]

  # get training subjects
  id_train_inner <- sort(sample(n_train, n_train * prop_train_inner))

  return(id_train_inner)

}

remove_essential_confounds <- function(conf, conf_names) {
  # Main function to select confounds
  
  # note essential confounds
  ess_confounds <- c("Sex_1_Site_1", "Sex_1_Site_2", "Sex_1_Site_3",
                     "Site_1_vs_2", "Site_1_vs_3", "HeadSize_Site_1",
                     "HeadSize_Site_2", "HeadSize_Site_3",
                     "HeadMotion_mean_rfMRI_rel_Site_1",
                     "HeadMotion_mean_rfMRI_rel_Site_2",
                     "HeadMotion_mean_rfMRI_rel_Site_3")

  # get values of essential confounds
  ess_confounds_idx <- sapply(ess_confounds, function(name) which(unlist(conf_names) == name))
  conf_ess <- conf[, ess_confounds_idx]

  # Remove essential confounds from the list
  idx_remove_ess <- setdiff(seq_len(ncol(conf)), ess_confounds_idx)
  conf_remove_ess <- conf[, idx_remove_ess]
  conf_names_remove_ess <- conf_names[idx_remove_ess]

  # create list of confound data
  conf_list <- list("conf_ess" = conf_ess, "conf_remove_ess" = conf_remove_ess, "conf_names_remove_ess" = conf_names_remove_ess)

  return(conf_list)
}



# Function to remove confounds containing a specific substring in their names
remove_confounds_by_name <- function(conf, conf_names, name_substring) {
  idx_to_remove <- which(grepl(name_substring, conf_names))
  idx_to_keep <- setdiff(seq_len(ncol(conf)), idx_to_remove)
  conf <- conf[, idx_to_keep]
  return(conf)
}

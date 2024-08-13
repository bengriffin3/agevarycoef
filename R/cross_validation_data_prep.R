get_idp_train_inner <- function(df_train, prop_train_inner) {

  n_train <- dim(df_train)[1]

  # get training subjects
  id_train_inner <- sort(sample(n_train, n_train * prop_train_inner))

  return(id_train_inner)

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
  ica_object_idps <- fastICA::fastICA(idps_train, n_ICs)
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
  pca_result_x <- stats::prcomp(idps_train, center = TRUE, scale. = TRUE)
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

confound_selection_using_train <- function(conf, conf_names, id_train) {

    ### first we do stuff not based on training/test split
    logger::log_info(paste0("Number of total possible confounds: ", dim(conf)[2]))

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

    logger::log_info(paste0("Number of confounds remaining: ", dim(conf)[2]))

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



upweight_idps_using_train <- function(idps, id_train, trait) {

  # split idps and trait into train and test
  idps_train <- idps[id_train, ]
  idps_test <- idps[-id_train, ]
  trait_train <- trait[id_train]
  trait_test <- trait[-id_train]


  idp_weightings <- stats::cor(idps_train, trait_train)^2

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
  cc1_train <- CCA::cc(idps_pca_train, as.matrix(trait_train))

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


preprocess_traits <- function(trait_code, idps, id_train, trait, n_feat) {

  #Options
  # 0 = PCA
  # 999 = perform CCA

  preprocessing_list = list()

  if (trait_code == 0) {
    trait_pca_list <- perform_pca_using_train(trait, id_train, 1)
    trait <- trait_pca_list[[1]] # pca_object_trait <- trait_pca_list[[2]]
    preprocessing_list$trait <- trait
  } else if (trait_code == 999) {
    cca_list <- perform_cca_using_train(idps, id_train, trait)
    trait <- cca_list[[1]]
    cca_object_idps_trait <- cca_list[[2]]
    preprocessing_list$trait <- trait
    preprocessing_list$cca_object_idps_trait <- cca_object_idps_trait
  } else  {
    trait = trait
  }

  return(preprocessing_list)
}

preprocess_idps <- function(idp_preproc, idps, id_train, trait, n_feat) {

  #Options
  # 0 = leave it as it is
  # 1 = perform ICA
  # 2 = perform weighted-ICA
  # 3 = perform CCA
  # 4 = perform PCA
  # 5 = select best features


  preprocessing_list = list()

  if (idp_preproc == 0) {
    idps <- idps
    preprocessing_list$idps = idps
  } else if (idp_preproc %in% c(1, 2)) {

    # this seems to do very little
    if (idp_preproc == 2) {
      idps <- upweight_idps_using_train(idps, id_train, trait)
    }

    idps_ica_list <- perform_ica_using_train(idps, id_train, n_feat)
    idps <- idps_ica_list[[1]]
    ica_object_idps <- idps_ica_list[[2]]
    preprocessing_list$idps = idps
    preprocessing_list$ica_object_idps = ica_object_idps

  } else if (idp_preproc == 3) {
    # this is CCA (done elsewhere?)
    preprocessing_list$idps = idps
  } else if (idp_preproc == 4) {
    idps_pca_list <- perform_pca_using_train(idps, id_train)
    idps <- idps_pca_list[[1]]
    # pca_object_trait <- trait_pca_list[[2]]
    preprocessing_list$idps = idps
  } else if (idp_preproc == 5) {
    load(paste0(proj_dir, "/data/idx_1436_keep_16_ordered.Rdata"))
    # select best idp
    idps_best_16 <- idps[, idx_1436_keep_16_ordered]
    idps <- idps_best_16[, 1:n_feat]
    preprocessing_list$idps = idps
  }

    return(preprocessing_list)
}







preprocessing_pipeline <- function(data, n_sub, trait_code, prop_train, idp_preproc) {

  # Extract the variables from the mat file and subsample
  idps <- data$structural.idps[1:n_sub, ]
  age <- data$age1[1:n_sub]
  trait <- get_traits(data, trait_code, n_sub)
  conf <- data$conf1[1:n_sub, ]
  conf_names <- data$conf.names


  # Remove subjects with NaNs (can impute but removing is better)
  clean_vars <- remove_nan_sub(idps, trait, age, conf)
  age <- clean_vars$age

  # note training subjects
  id_train <- get_training_sample_id(clean_vars$trait, prop_train)

  # Select confounds
  conf <- confound_selection_using_train(clean_vars$conf, conf_names, id_train)

  # scale idps
  idps <- scale_data_using_train(clean_vars$idps, id_train)

  # perform preprocessing on traits (e.g., PCA/CCA)
  trait_preproc_list <- preprocess_traits(trait_code, idps, id_train, clean_vars$trait, n_feat)
  trait <- trait_preproc_list$trait

  # perform preprocessing on IDPs (e.g., ICA/PCA/CCA)
  idp_preproc_list <- preprocess_idps(idp_preproc, idps, id_train, trait, n_feat)
  idps <- idp_preproc_list$idps

  # de-mean target
  trait <- de_mean_trait_using_train(trait, id_train)

  # deconfound data
  idps <- deconfound_data_using_train(idps, conf, id_train)
  trait <- deconfound_data_using_train(trait, conf, id_train)

  pre_processed_data <- list()
  pre_processed_data$idps <- idps
  pre_processed_data$trait <- trait
  pre_processed_data$age <- age
  pre_processed_data$conf <- conf
  pre_processed_data$id_train <- id_train

  return(pre_processed_data)

}

remove_nan_sub <- function(idps, trait, age, conf) {

  idx_keep <- stats::complete.cases(idps, trait, age, conf)

  idps <- idps[idx_keep, ]
  if (is.array(trait)) {
    trait <- trait[idx_keep, ]
  } else {
    trait <- trait[idx_keep]
  }
  age <- age[idx_keep]
  conf <- conf[idx_keep, ]

  list_vars <- mget(c("idps", "trait", "age", "conf"))

  logger::log_warn(paste0("Number of subjects removed: ", sum(!idx_keep)))
  logger::log_info(paste0("Number of subjects remaining: ", dim(idps)[1]))
  logger::log_info(paste0("Number of features/IDPs: ", dim(idps)[2]))


  return(list_vars)
}


perform_pca <- function(idps, n_feat) {

  pca_result <- stats::prcomp(idps, center = TRUE, scale. = TRUE)
  idps <- pca_result$x[, 1:n_feat]
  return(idps)
}

get_x_pcs <- function(x, n_pc_comp) {
  # perform PCA including de-meaning and setting unit variance
  pca_result_x <- stats::prcomp(x, center = TRUE, scale. = TRUE)

  # note down the PCs
  x_pca <- pca_result_x$x[, 1:n_pc_comp]
  return(x_pca)
}


# Main function to select confounds
remove_essential_confounds <- function(conf, conf_names) {

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

get_training_sample_id <- function(trait, prop_train) {
  # sample training data
  if (is.vector(trait)) {
    id_train <- sort(sample(length(trait), length(trait) * prop_train))
  } else if (is.matrix(trait)) {
    id_train <- sort(sample(dim(trait)[1], dim(trait)[1] * prop_train))
  }
  return(id_train)
}

get_traits <- function(data, trait_id, n_sub) {

  # load traits
  if (trait_id == 1 || trait_id == 2 || trait_id == 3) {
    trait <- data$vars.cog.best.3[1:n_sub, trait_id]
  } else {
    # if performing PC or CCA, let's use the top 30 cog traits
    trait <- get_top_30_cog_vars(data, n_sub)
  }

  return(trait)
}

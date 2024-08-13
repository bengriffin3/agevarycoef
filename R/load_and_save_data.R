library(R.matlab)

load_prediction_variables <- function(n_subjects) {
  # Load the mat file
  if (n_subjects <= 5000) {
    print(n_subjects)
    data <- R.matlab::readMat(paste0("../vary_coef_bg/data/prediction_vars_5000.mat"))
  } else if (n_subjects > 5000 && n_subjects <= 20000) {
    data <- R.matlab::readMat(paste0("../vary_coef_bg/data/prediction_vars_20000.mat"))
  } else if (n_subjects > 20000) {
    data <- R.matlab::readMat(paste0("../vary_coef_bg/data/prediction_vars.mat"))
  }
  # note index of best 30 traits
  idx30 <- R.matlab::readMat(paste0("../vary_coef_bg/data/idx_30.mat"))
  data$idx30 <- idx30$idx

  return(data)
}

load_cog_variables <- function(n_subjects) {
  # load cognitive vars
  if (n_subjects <= 1000) {
    data <- R.matlab::readMat(paste0("../vary_coef_bg/data/cognitive_vars_1000.mat"))
  } else if (n_subjects <= 20000) {
    data <- R.matlab::readMat(paste0("../vary_coef_bg/data/cognitive_vars_20000.mat"))
  } else {
    data <- R.matlab::readMat(paste0("../vary_coef_bg/data/cognitive_vars.mat"))
  }
  return(data)
}

# Save functions used in the script for better organization
save_training_data <- function(dir, df_train, age, id_train, idp_preproc) {
  print("Saving training data then exiting...")
  age_train <- age[id_train]
  save(list = intersect(ls(), c("df_train", "age_train", "id_train", "ica_object_idps")),
       file = sprintf("%s/train_data_t_%i_n_f_%i_f_1436_n_%i_p_%i_idp_preproc_%i.RData",
                      dir, trait_code, n_feat, n_sub, perc_train, idp_preproc))
  quit(save="ask")
  #Sys.sleep(5)
}


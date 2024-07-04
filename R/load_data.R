library(R.matlab)

load_prediction_variables <- function(n_subjects, proj_dir) {
  # Load the mat file
  if (n_subjects <= 5000) {
    print(n_subjects)
    data <- R.matlab::readMat(paste0(proj_dir, "/data/prediction_vars_5000.mat"))
  } else if (n_subjects > 5000 && n_subjects <= 20000) {
    data <- R.matlab::readMat(paste0(proj_dir, "/data/prediction_vars_20000.mat"))
  } else if (n_subjects > 20000) {
    data <- R.matlab::readMat(paste0(proj_dir, "/data/prediction_vars.mat"))
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

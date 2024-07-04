get_top_30_cog_vars <- function(data, n_sub, proj_dir) {

  cog_data <- load_cog_variables(n_sub, proj_dir)
  y <- cog_data$vars.cog
  y_top_30 <- y[1:n_sub, data$idx30 + 1]
  return(y_top_30)
}

#' Determining non-zero elastic net coefficients
#'
#' @param enet_coefficients coefficients from a fitted elastic net model
#'
#' @return A vector of the best features
#' @export
#'
#' @examples
#' enet_coefficients <- c(0,0,0.1,0.2,0,0)
#' determine_non_zero_coeff(enet_coefficients)
determine_non_zero_coeff <- function(enet_coefficients) {
  idx_best_features <- which(enet_coefficients != 0)
  return(idx_best_features)
}

determine_non_zero_coeff2 <- function(enet_coefficients) {

  idx_best_features <- determine_non_zero_coeff(enet_coefficients)
  # to save index: save(idx_best_features, file = "idx_best_features.RData")
  return(idx_best_features)
}

extractorRData <- function(file, object) {

  E <- new.env()
  load(file, envir = E)
  # Function for extracting an object from a .RData file created by R's save() command
  # Inputs: RData file, object name
  return(get(object, envir = E, inherits = F))
}

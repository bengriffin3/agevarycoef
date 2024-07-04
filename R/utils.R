get_top_30_cog_vars <- function(data, n_sub) {

  cog_data <- load_cog_variables(n_sub)
  y <- cog_data$vars.cog
  y_top_30 <- y[1:n_sub, data$idx30 + 1]
  return(y_top_30)
}


extractorRData <- function(file, object) {

  E <- new.env()
  load(file, envir = E)
  # Function for extracting an object from a .RData file created by R's save() command
  # Inputs: RData file, object name
  return(get(object, envir = E, inherits = F))
}

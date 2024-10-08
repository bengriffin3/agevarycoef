% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/prediction.r
\name{run_linear_model_cv}
\alias{run_linear_model_cv}
\title{Run Linear Model with 10-Fold Cross-Validation}
\usage{
run_linear_model_cv(
  idps,
  trait,
  age,
  conf,
  conf_names,
  trait_id,
  remove_age,
  model_age
)
}
\arguments{
\item{idps}{A matrix or data frame of imaging-derived phenotypes (IDPs).}

\item{trait}{A numeric vector representing the trait (target) variable.}

\item{age}{A numeric vector representing age.}

\item{conf}{A matrix or data frame of confounding variables.}

\item{conf_names}{A character vector of confounder names.}

\item{trait_id}{A numeric identifier for the trait type; 999 indicates CCA, 0 indicates PCA.}

\item{remove_age}{A binary flag (0 or 1) indicating whether to remove age from confounding variables.}

\item{model_age}{A binary flag indicating whether to include age in the model.}
}
\value{
A list containing:
\item{lm_yhat}{A numeric vector of predictions for each fold (test set).}
\item{corr_train}{A numeric vector of correlations between true and predicted values on the training set for each fold.}
\item{corr_test}{A numeric vector of correlations between true and predicted values on the test set for each fold.}
\item{models}{A list of fitted linear models for each fold.}
\item{trait_transformed}{A numeric vector of the transformed trait used for prediction.}
}
\description{
This function fits a linear model using 10-fold cross-validation. It returns
predictions, correlations for both training and test sets, and the fitted models
for each fold.
}

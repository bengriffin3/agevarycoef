% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/prediction.r
\name{run_lgboost_model}
\alias{run_lgboost_model}
\title{Run LightGBM Model}
\usage{
run_lgboost_model(
  df_all_train_x,
  df_all_train_y,
  df_all_test_x,
  df_all_test_y,
  params = 0
)
}
\arguments{
\item{df_all_train_x}{Training data matrix (features).}

\item{df_all_train_y}{Training data vector (target).}

\item{df_all_test_x}{Test data matrix (features).}

\item{df_all_test_y}{Test data vector (target).}

\item{params}{List of LightGBM parameters (optional).}
}
\value{
A list containing:
\item{predictions_train}{Predicted values for the training data.}
\item{predictions_test}{Predicted values for the test data.}
\item{lgb_model}{The trained LightGBM model.}
}
\description{
Fits a LightGBM model on training data and evaluates it on test data.
}

#' Prediction information printing
#'
#' @param trait Target trait to predict
#' @param n_feat Number of features to use to predict
#' @param n_sub Number of subjects
#' @param perc_tn Percentage of training subjects
#' @param prof Use profile likelihood
#' @param tap Use covariance tapering
#' @param cov Which covariance function to use
#' @param run_svc Whether to run varying-coefficient model
#'
#' @return Print information
#' @export
#'
#' @examples
#' trait <- 1; n_feat <- 1436; n_sub <- 5000; perc_tn <- 90; prof <- 1;
#' tap <- 5; cov = "exp"; run_svc <- 1
#' print_info(trait, n_feat, n_sub, perc_tn, prof, tap, cov, run_svc)
print_info <- function(trait, n_feat, n_sub, perc_tn, prof, tap, cov, run_svc) {
  print(paste0("Trait: ", trait))
  print(paste0("Number of features: ", n_feat))
  print(paste0("Number of subjects: ", n_sub))
  print(paste0("Proportion of training subjects: ", perc_tn))
  print(paste0("Profile likelihood: ", prof))
  print(paste0("Covariance tapering: ", tap))
  print(paste0("Cov func: ", cov))
  print(paste0("Run SVC: ", run_svc))
}


print_linear_accuracy_info <- function(se_lm, se_enet, corr_lm_in, corr_enet_in, corr_lm_out, corr_enet_out, id_train) {

  print(knitr::kable(data.frame(
    model = c("linear", "enet"),
    `in-sample RMSE` =
      round(sqrt(c(mean(se_lm[id_train]), mean(se_enet[id_train]))), 3),
    `out-of-sample RMSE` =
      round(sqrt(c(mean(se_lm[-id_train]), mean(se_enet[-id_train]))), 3)
  )))

  knitr::kable(data.frame(
    model = c("linear", "enet"),
    `in-sample corr` =
      c(corr_lm_in, corr_enet_in),
    `out-of-sample corr` =
      c(corr_lm_out, corr_enet_out)
  ))

}

print_accuracy_info <- function(se_lm, se_enet, se_svc, corr_lm_in, corr_enet_in, corr_svc_in, corr_lm_out, corr_enet_out, corr_svc_out, id_train) {

  print("Show dataframe of accuracy")
  knitr::kable(data.frame(
    model = c("linear", "SVC", "enet"),
    `in-sample RMSE` =
      round(sqrt(c(mean(se_lm[id_train]),
                   mean(se_svc[id_train]), mean(se_enet[id_train]))), 3),
    `out-of-sample RMSE` =
      round(sqrt(c(mean(se_lm[-id_train]),
                   mean(se_svc[-id_train]), mean(se_enet[-id_train]))), 3),
    `in-sample corr` =
      c(corr_lm_in, corr_svc_in, corr_enet_in),
    `out-of-sample corr` =
      c(corr_lm_out, corr_svc_out, corr_enet_out)
  ))
}

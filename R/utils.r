library(R.matlab)

save_training_data <- function (proj_dir, trait_id, n_feat, n_sub, perc_train, ica, df_train, age_train, id_train) {

  #, "ica_object_idps")),
  dir <- paste0(proj_dir, "/results/training_data")
  save(list = intersect(ls(), c("df_train", "age_train", "id_train")),
       file = sprintf("%s/train_data_t_%i_n_f_%i_f_1436_n_%i_p_%i_ica_%i.RData",
                      dir, trait_id, n_feat, n_sub, perc_train, ica))
  #Sys.sleep(5)
  quit()


}

check_for_existing_model <- function(proj_dir, trait_id, n_feat, n_sub, perc_train, prof, tap, cov, ica, run_svc, model_age, remove_age) {

  save_dir <- paste0(proj_dir, "/results/results_with_boost")
  if (remove_age == 0) {
    if (run_svc == 1) {
      file1 <- sprintf("%s/results_t_%i_f_%i_n_%i_p_%i_ica_%i_mi_%i_svc_%i_l_%i_ta_%i_c_%s.RData",
                        save_dir, trait_id, n_feat, n_sub, perc_train, ica, model_age, run_svc, prof, tap, cov)
    } else {
        file1 <- sprintf("%s/results_t_%i_f_%i_n_%i_p_%i_ica_%i_mi_%i_svc_%i.RData",
                          save_dir, trait_id, n_feat, n_sub, perc_train, ica, model_age, run_svc)
    }
  } else if (remove_age == 1) {
    if (run_svc == 1) {
        file1 <- sprintf("%s/results_t_%i_f_%i_n_%i_p_%i_ica_%i_mi_%i_svc_%i_l_%i_ta_%i_c_%s_ra_%i.RData",
                          save_dir, trait_id, n_feat, n_sub, perc_train, ica, model_age, run_svc, prof, tap, cov, remove_age)
    } else {
        file1 <- sprintf("%s/results_t_%i_f_%i_n_%i_p_%i_ica_%i_mi_%i_svc_%i_ra_%i.RData",
                          save_dir, trait_id, n_feat, n_sub, perc_train, ica, model_age, run_svc, remove_age)
    }
  }


  if (file.exists(file1)) {
    print("Model already run so exiting...")
    stop("Model already run")
    #Sys.sleep(5)
    #quit()
  }

}

#' Extract an object from a .RData file
#'
#' This function extracts an object from an RData file created by R's save() command.
#'
#' @param file A character string specifying the path to the RData file.
#' @param object A character string specifying the name of the object to extract.
#' @return The object extracted from the RData file.
#' @export
extractorRData <- function(file, object) {
  E <- new.env()
  load(file = file, envir = E)
  return(get(object, envir = E, inherits = F))
}
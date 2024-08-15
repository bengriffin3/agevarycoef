library(R.matlab)
library(varycoef)
library(optparse)
library(fastICA)
library(CCA)
library(logger)
library(glmnet)

save_training_data <- function (proj_dir, trait_id, n_feat, n_sub, perc_train, ica, df_train, age_train, id_train) {

  #, "ica_object_idps")),
  dir <- paste0(proj_dir, "/results/training_data")
  save(list = intersect(ls(), c("df_train", "age_train", "id_train")),
       file = sprintf("%s/train_data_t_%i_n_f_%i_f_1436_n_%i_p_%i_ica_%i.RData",
                      dir, trait_id, n_feat, n_sub, perc_train, ica))
  #Sys.sleep(5)
  quit()


}

check_for_existing_model <- function(proj_dir, trait_id, n_feat, n_sub, perc_train, prof, tap, cov, ica, run_svc, model_age) {

  dir <- paste0(proj_dir, "/results/results_GP")
  file1 <- sprintf("%s/rel_t_%i_f_%i_n_%i_p_%i_l_%i_ta_%i_c_%s_ica_%i_svc_%i_mi_%i.RData",
                      dir, trait_id, n_feat, n_sub, perc_train, prof, tap, cov, ica, run_svc, model_age)
  if (file.exists(file1)) {
    print("Model already run so exiting...")
    #Sys.sleep(5)
    #quit()
  }

}

extractorRData <- function(file, object) {
      #' Function for extracting an object from a .RData file created by R's save() command
      #' Inputs: RData file, object name
      E <- new.env()
      load(fil = file, envir = E)
      return(get(object, envir = E, inherits = F))
    }
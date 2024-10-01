
print_info <- function(trait_id, n_sub, run_svc, remove_age, model_age, n_feat=1) {
  print(paste0("Trait: ", trait_id))
  print(paste0("Number of features: ", n_feat))
  print(paste0("Number of subjects: ", n_sub))
  print(paste0("Run SVC: ", run_svc))
  print(paste0("Deconfound age: ", remove_age))
  print(paste0("Model age: ", model_age))
}
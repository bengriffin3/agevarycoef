library(testthat)
library(agevarycoef)

# Test 1: Check that NA values are correctly removed
test_that("Subjects with missing age or confounders are removed", {

  # Example data for testing
  idps <- matrix(1:12, nrow = 4, ncol = 3)  # 4 subjects, 3 IDPs
  trait <- c(1, 2, NA, 4)               # Trait vector with one NA value
  age <- c(25, 30, NA, 45)              # Age vector with one NA value
  conf <- matrix(c(1, 2, NA, 4, 1, 2, 3, 4), ncol=2) # Confounders with one NA value


  result <- remove_nan_sub(idps, trait, age, conf)

  # Expect 2 subjects to be kept (since 2 rows have NA)
  expect_equal(dim(result$idps)[1], 3)
  expect_equal(length(result$trait), 3)
  expect_equal(length(result$age), 3)
  expect_equal(dim(result$conf)[1], 3)
})

# Test 2: Check that IDPs and trait are correctly subsetted when no NA values
test_that("Correct behavior when no NA values", {
  idps_clean <- matrix(1:12, nrow = 4, ncol = 3)  # Clean IDP data
  trait_clean <- c(1, 2, 3, 4)
  age_clean <- c(25, 30, 35, 45)
  conf_clean <- matrix(1:8, nrow = 4, ncol = 2)  # Clean confounders

  result <- remove_nan_sub(idps_clean, trait_clean, age_clean, conf_clean)

  # All subjects should be kept
  expect_equal(dim(result$idps)[1], 4)
  expect_equal(length(result$trait), 4)
  expect_equal(length(result$age), 4)
  expect_equal(dim(result$conf)[1], 4)
})

# Test 3: Check handling of vector input for IDPs and trait
test_that("Handles non-array input for IDPs and trait", {
  idps_vec <- c(1, 2, NA, 4)  # Vector IDP with one NA
  trait_vec <- c(5, 6, 7, 8)
  age_vec <- c(25, 30, NA, 45)
  conf_vec <- matrix(c(1, 2, NA, 4, 5, 6, 7, 8), ncol = 2)  # Confounders

  result <- remove_nan_sub(idps_vec, trait_vec, age_vec, conf_vec)

  # Expect 2 subjects to be kept
  expect_equal(length(result$idps), 3)
  expect_equal(length(result$trait), 3)
  expect_equal(length(result$age), 3)
  expect_equal(dim(result$conf)[1], 3)
})
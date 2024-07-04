test_that("determine_non_zero_coeff() ignore 0s", {
  expect_equal(determine_non_zero_coeff(c(5,2,1,0,0,2,3)), c(1, 2, 3, 6, 7))
})

test_that("determine_non_zero_coeff() includes negative", {
  expect_equal(determine_non_zero_coeff(c(-5,-2,1,0,0,2,-3)), c(1, 2, 3, 6, 7))
})

test_that("determine_non_zero_coeff2() includes negative", {
  expect_equal(determine_non_zero_coeff2(c(-5,-2,1,0,0,2,-3)), c(1, 2, 3, 6, 7))
})

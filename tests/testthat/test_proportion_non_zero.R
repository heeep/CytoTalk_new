library(testthat) 
library(Matrix)


# Start writing tests
test_that("Basic Functionality", {
  mat <- matrix(c(1, 0, 2, 3, 0, 4), nrow = 2, byrow = TRUE)
  expect_equal(proportion_non_zero(mat), c(2/3, 2/3))
})


test_that("Single Row and Single Column Matrices", {
  mat_row <- matrix(c(1, 0, 2), nrow = 1)
  expect_equal(proportion_non_zero(mat_row), 2/3)
  
  mat_col <- matrix(c(0, 0, 2), ncol = 1)
  expect_equal(proportion_non_zero(mat_col), c(0,0,1))
})



test_that("mi_mat_parallel returns correct MI values", {
  mat_test <- matrix(c(1, 2, 3, 2, 3, 4), ncol=2)

  mi_values <- mi_mat_parallel(mat_test, method="emp", normalize=FALSE)
  expect_true(is.matrix(mi_values))
  expect_equal(dim(mi_values), c(2, 2))
  expect_true(mi_values[1, 2] > 0)
  mi_values_norm <- mi_mat_parallel(mat_test, method="emp", normalize=TRUE)
  expect_true(is.matrix(mi_values_norm))
  expect_equal(dim(mi_values_norm), c(2, 2))
  expect_true(all(mi_values_norm >= 0 & mi_values_norm <= 1))
})


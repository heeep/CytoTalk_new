test_that("Empty matrix returns an empty matrix", {
  expect_equal(subset_non_zero(matrix(numeric(0), nrow = 0, ncol = 0), 0.5), matrix(numeric(0), nrow = 0, ncol = 0))
})

test_that("All-zero matrix returns an empty matrix", {
  mat <- matrix(0, nrow = 3, ncol = 3)
  expect_equal(subset_non_zero(mat, 0.5), matrix(numeric(0), nrow = 0, ncol = 3))
})

test_that("Proportion of non-zero elements equals cutoff", {
  mat <- matrix(c(1, 0, 0, 1, 1, 0, 0, 0, 0), nrow = 3)
  expect_equal(subset_non_zero(mat, 1/3), mat[c(1,2), ])
})

test_matrix <- matrix(c(1, 2, 3, 4), nrow = 2)
expected_normalized <- log1p(t(t(test_matrix) / colSums(test_matrix) * 10000))


test_that("Matrix is normalized correctly", {
  result <- normalize_sparse(test_matrix)
  expect_equal(result, expected_normalized)
})
test_that("Empty matrix is handled", {
  empty_matrix <- matrix(numeric(0), nrow = 0, ncol = 0)
  result <- normalize_sparse(empty_matrix)
  expect_equal(dim(result), c(0, 0))
})

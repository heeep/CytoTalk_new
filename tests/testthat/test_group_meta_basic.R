test_that("Basic functionality with multiple groups", {
  mat <- matrix(1:12, nrow = 3, ncol = 4)
  labels <- c("A", "B", "A", "B")
  result <- group_meta_basic(mat, labels)
  expect_equal(length(result), 2)
  expect_equal(names(result), sort(unique(labels)))
  expect_true(all(c("A", "B") %in% names(result)))
  expect_equal(result$A, mat[, c(1, 3), drop = FALSE])
  expect_equal(result$B, mat[, c(2, 4), drop = FALSE])
})


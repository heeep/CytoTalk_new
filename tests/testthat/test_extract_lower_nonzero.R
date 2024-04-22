
test_that("Handles all-zero matrix correctly", {
  mat <- matrix(0, nrow = 3, ncol = 3)
  rownames(mat) <- colnames(mat) <- c("Row1", "Row2", "Row3")
  
  result <- extract_lower_nonzero(mat)
  
  expect_equal(nrow(result), 0)
})



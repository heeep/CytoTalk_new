library(testthat)

# Assuming the other functions (is_integer, normalize_sparse) are provided and correct.

# Test cases for ungroup_meta
test_that("ungroup_meta correctly combines matrices and assigns cell types", {
  # Create a list of matrices with names
  lst <- list(
    "type1" = matrix(1:4, ncol = 2),
    "type2" = matrix(5:8, ncol = 2)
  )
  
  # Expected result
  expected_mat <- matrix(c(1:4, 5:8), ncol = 4)
  expected_cell_types <- c("type1", "type1", "type2", "type2")
  
  # Apply the function
  result <- ungroup_meta(lst)
  
  # Check if the result is a list with correct elements
  expect_true(is.list(result))
  expect_true("mat" %in% names(result))
  expect_true("cell_types" %in% names(result))
  
  # Check if the combined matrix is correct
  expect_equal(result$mat, expected_mat)
  
  # Check if the cell types are correctly assigned
  expect_equal(result$cell_types, expected_cell_types)
})



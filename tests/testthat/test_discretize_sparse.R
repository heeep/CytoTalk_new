test_that("Test with normal input data", {

  input_matrix <- matrix(1:9, nrow = 3, ncol = 3)
  
  result <- discretize_sparse(input_matrix)
  
  expect_is(result, "Matrix")
  
  expect_equal(dim(result), c(3, 3))
  
})




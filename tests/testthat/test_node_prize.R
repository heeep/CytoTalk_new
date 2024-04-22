test_that("Function works correctly with standard input", {
 
  mat_pem <- matrix(c(1, 2, -1, 4), nrow = 2, ncol = 2)
  colnames(mat_pem) <- c("Type1", "Type2")
  cell_type <- "Type1"
  vec_relev <- c(Gene1 = 0.5, Gene2 = 0.6)
  

  result <- node_prize(mat_pem, cell_type, vec_relev)
  
 expect_equal(length(result), 2)
   

})
test_that("Handles unmatched column names gracefully", {
  mat_pem <- matrix(c(1, 2, 3, 4), nrow = 2, ncol = 2)
  colnames(mat_pem) <- c("TypeA", "TypeB")
  cell_type <- "TypeC"     #non-existent type
  vec_relev <- c(Gene1 = 0.5, Gene2 = 0.6)
  
  
  result <- node_prize(mat_pem, cell_type, vec_relev)
  

 
  expect_true(all(is.na(result))) 
})
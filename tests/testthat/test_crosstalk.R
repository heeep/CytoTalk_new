library(testthat)

test_that("crosstalk function processes inputs correctly with mat_type matrix", {
  mat_pem <- matrix(c(1, 2, -1, NA, 3, 4), nrow = 3, byrow = TRUE)
  rownames(mat_pem) <- c("L1", "R1", "L2")
  colnames(mat_pem) <- c("CellA", "CellB")
  
  vec_nst_a <- c(0.5, 0.8)  
  vec_nst_b <- c(0.6, 0.9)  
  
  mat_type <- matrix(c(1, 0, 0, 1, 1, 0), nrow = 3, byrow = TRUE)
  rownames(mat_type) <- c("L1", "R1", "L2")
  colnames(mat_type) <- c("CellA", "CellB")
  
  lrp <- data.frame(ligand = c("L1", "L2"), receptor = c("R1", "R1"))

  result <- crosstalk(mat_pem, "CellA", "CellB", vec_nst_a, vec_nst_b, mat_type, lrp)
  
  expect_true(all(c("ligand", "receptor", "ligand_type", "receptor_type", "nst", "expr") %in% names(result)),
              "The result should have all expected columns.")

  expect_true(nrow(result) >= nrow(lrp) && nrow(result) <= nrow(lrp) * 2,
              "The result should have row count greater than or equal to the number of lrp rows, and less than or equal to twice the number of lrp rows")
})




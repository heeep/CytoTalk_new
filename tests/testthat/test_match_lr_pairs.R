library(testthat)
mat <- matrix(1:4, nrow=2, dimnames = list(c("A", "B"), NULL))
lrp <- data.frame(ligand = c("A", "C"), receptor = c("B", "D"))

test_that("Existing pairs are matched correctly", {
  result <- match_lr_pairs(mat, lrp)
  expect_equal(nrow(result), 1) 
  expect_equal(result$ligand, "A")
  expect_equal(result$receptor, "B")
})


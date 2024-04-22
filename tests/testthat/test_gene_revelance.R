test_that("Empty matrix is handled properly", {
  mat_intra <- matrix(numeric(0), nrow = 0, ncol = 0)
  lrp <- character(0) 
  expect_error(gene_relevance(mat_intra, lrp))
})
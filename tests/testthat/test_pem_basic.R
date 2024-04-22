library(testthat)
library(Matrix)

mat_scrna <- matrix(c(1, 2, 3, 2, 3, 4, 3, 4, 5), nrow = 3, byrow = FALSE,
                    dimnames = list(c("gene1", "gene2", "gene3"), c("cell1", "cell2", "cell3")))

labels <- c("type1", "type1", "type2")

test_that("pem_basic calculates expected mathematical results", {
  
  result <- pem_basic(mat_scrna, labels)
  
 
  expected_result <- matrix(c(-0.04533511, -0.00550334, 0.00733284, 
                              0.01026995, 0.00131816, -0.0017889), 
                            ncol = 2,
                            dimnames = list(c("gene1", "gene2", "gene3"), c("type1", "type2")))
  
      8
  expect_equal(result, expected_result, tolerance = 1e-5)
})

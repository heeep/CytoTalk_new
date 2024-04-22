library(testthat)

# Define the extract_group_basic function
extract_group_basic <- function(group, mat, labels) {
  mat[, labels == group, drop = FALSE]
}

# Start testing
test_that("Basic Functionality Test", {
  mat <- matrix(1:9, nrow = 3)
  labels <- c("A", "B", "A")
  expect_equal(extract_group_basic("A", mat, labels), mat[, c(TRUE, FALSE, TRUE)])
})

test_that("Empty Matrix Test", {
  mat <- matrix(numeric(0), nrow = 0)
  labels <- character(0)
  expect_equal(extract_group_basic("A", mat, labels), matrix(numeric(0), nrow = 0))
})

test_that("Non-matching Labels Test", {
  mat <- matrix(1:4, nrow = 2)
  labels <- c("A", "B")
  expect_equal(dim(extract_group_basic("C", mat, labels)), c(2, 0))
})


                  

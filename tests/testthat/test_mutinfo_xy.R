library(testthat)
library(infotheo)

test_that("mutinfo_xy calculates mutual information correctly", {
  x <- c(1, 2, 3, 4, 5)
  y <- c(5, 4, 3, 2, 1)
  expect_true(is.numeric(mutinfo_xy(x, y)))
  expect_true(mutinfo_xy(x, y) > 0)
})

test_that("mutinfo_xy handles normalization correctly", {
  x <- c(1, 2, 3, 4, 5)
  y <- c(5, 4, 3, 2, 1)
  normalized_mutinfo <- mutinfo_xy(x, y, normalize = TRUE)
  expect_true(normalized_mutinfo >= 0 && normalized_mutinfo <= 1)
})

test_that("mutinfo_xy returns 0 for identical variables when normalized", {
  x <- c(1, 1, 1, 1, 1)
  y <- x
  expect_equal(mutinfo_xy(x, y, normalize = TRUE), 0)
})

test_that("mutinfo_xy handles completely independent variables correctly", {
  x <- c(1, 2, 3, 4, 5)
  y <- c(1, 1, 1, 1, 1) # y is constant, implying independence
  expect_equal(mutinfo_xy(x, y), 0)
})

# 可以添加更多测试用例以覆盖不同的输入和边缘情况

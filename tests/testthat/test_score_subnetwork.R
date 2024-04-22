test_that("Basic functionality with positive values", {
  prizes <- c(10, 20, 30)
  costs <- c(5, 5, 5)
  beta <- 2
  result <- score_subnetwork(prizes, costs, beta)
  expect_equal(result, c(mean(c(10, 20, 30)), mean(c(5, 5, 5)), 2 * sum(c(10, 20, 30)) - sum(c(5, 5, 5))))
})

test_that("Empty vectors return NaN for means and -sum(costs) for potential", {
  prizes <- numeric(0)
  costs <- numeric(0)
  beta <- 1
  result <- score_subnetwork(prizes, costs, beta)
  expect_true(all(is.nan(result[1:2])))
  expect_equal(result[3], -sum(costs))
})

test_that("Handling of negative values", {
  prizes <- c(-10, -20)
  costs <- c(-5, -5)
  beta <- 1
  result <- score_subnetwork(prizes, costs, beta)
  expect_equal(result, c(mean(prizes), mean(costs), beta * sum(prizes) - sum(costs)))
})

test_that("Behavior with zero values", {
  prizes <- c(0, 0, 0)
  costs <- c(0, 0, 0)
  beta <- 1
  result <- score_subnetwork(prizes, costs, beta)
  expect_equal(result, c(0, 0, 0))
})


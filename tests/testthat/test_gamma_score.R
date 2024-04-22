test_that("gamma_score returns values in the correct range", {
  set.seed(123) # For reproducibility of random results
  v <- rnorm(100)
  
  # Test without inversion
  score <- gamma_score(v)
  expect_is(score, "numeric")
  expect_true(score >= 0)
  expect_true(score <= 1)
  
  # Test with inversion
  score_inv <- gamma_score(v, inverse = TRUE)
  expect_is(score_inv, "numeric")
  expect_true(score_inv >= 0)
  expect_true(score_inv <= 1)
})

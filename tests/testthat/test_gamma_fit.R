test_that("gamma_fit works correctly", {
  # 创建测试数据
  x <- rnorm(100, 10, 2) 
  y <- 5
  
  result <- gamma_fit(x, y)
  expect_length(result, 2)
  expect_type(result$params, "double")
  expect_length(result$params, 2)
  expect_true(all(result$params > 0))
  expect_is(result$kstest, "htest")
  expect_true(result$kstest$p.value > 0.05) 
})
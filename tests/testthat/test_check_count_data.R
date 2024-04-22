
count_matrix <- matrix(c(1L, 2L, 3L, 4L), nrow = 2)
non_count_matrix <- matrix(c(1.5, 2.5, 3.5, 4.5), nrow = 2)




test_that("Non-integer matrix is not transformed", {
  result <- check_count_data(non_count_matrix, auto_transform = TRUE)
  expect_equal(result, non_count_matrix)
})

test_that("No transformation occurs when auto_transform is FALSE", {
  result <- check_count_data(count_matrix, auto_transform = FALSE)
  expect_equal(result, count_matrix)
})



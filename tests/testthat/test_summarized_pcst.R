library(testthat)

lst_pcst <- list(
  edges = data.frame(
    beta = c(1, 1, 2, 2, 2),
    omega = c(0.1, 0.1, 0.2, 0.2, 0.3),
    node1 = c("A", "B", "C", "D", "E"),
    node2 = c("B", "A", "D", "C", "F"),
    stringsAsFactors = FALSE
  )
)

filter_arti <- function(df) {
  return(df)
}

test_that("summarize_pcst correctly summarizes pcst results", {
  summarized <- summarize_pcst(lst_pcst)

  expect_true(all(summarized$n_edges == c(2, 2, 1)))
  
})



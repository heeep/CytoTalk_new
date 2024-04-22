node_names_a <- c("Gene1_CellA", "Gene2_CellA")
node_names_b <- c("Gene3_CellB", "Gene4_CellB")
vec_np_a <- c(0.8, 0.5)
vec_np_b <- c(0.7, 0.6)
vec_pem_match_a <- c(10, 8)
vec_pem_match_b <- c(9, 11)
vec_gr_a <- c(0.8, 0.6)
vec_gr_b <- c(0.7, 0.9)
edge_full <- rbind(c("Gene1_CellA", "Gene3_CellB"), c("Gene2_CellA", "Gene4_CellB"))
cost_full <- c(1, 2)

df_net_node <- data.frame(
  node = c(node_names_a, node_names_b),
  prize = c(vec_np_a, vec_np_b),
  pem = c(vec_pem_match_a, vec_pem_match_b),
  gene_relevance = c(vec_gr_a, vec_gr_b)
)

df_net_edge <- data.frame(
  node1 = edge_full[, 1],
  node2 = edge_full[, 2],
  cost = cost_full
)

lst_net <- list(nodes = df_net_node, edges = df_net_edge)
beta_max <- 5
omega_min <- 0.1
omega_max <- 1.0

# 编写测试用例
test_that("run_pcst returns a list with nodes and edges data frames", {
  result <- run_pcst(lst_net, beta_max, omega_min, omega_max)
  expect_true(is.list(result))
 expect_true(is.data.frame(result$nodes))
  expect_true(is.data.frame(result$edges))
  expect_false(is.null(nrow(result$nodes)))
  expect_false(is.null(nrow(result$edges)))
  expect_gt(nrow(result$nodes), 0)
  expect_gt(nrow(result$edges), 0)
})


test_that("integrate_network returns correct structure and types", {
  vec_nst_a <- c(Gene1 = 0.8, Gene2 = 0.5)
  vec_nst_b <- c(Gene3 = 0.7, Gene4 = 0.6)
  
  # 内部交互矩阵
  mat_intra_a <- matrix(c(0, 0.4, 0.4, 0), nrow = 2, byrow = TRUE,
                        dimnames = list(c("Gene1", "Gene2"), c("Gene1", "Gene2")))
  mat_intra_b <- matrix(c(0, 0.5, 0.5, 0), nrow = 2, byrow = TRUE,
                        dimnames = list(c("Gene3", "Gene4"), c("Gene3", "Gene4")))
  
  # 特异性表达矩阵
  mat_pem <- matrix(c(10, 8, 1, 2, 2, 3, 9, 11), nrow = 4, byrow = FALSE,
                    dimnames = list(c("Gene1", "Gene2", "Gene3", "Gene4"), c("CellA", "CellB")))
  
  # 配体-受体对列表（lrp）
  lrp <- data.frame(
    ligand = c("Gene1", "Gene2"),
    receptor = c("Gene3", "Gene4")
  )
  mat_type <- matrix(c(1, 0, 0, 1, 1, 0, 0, 1), nrow = 4, byrow = FALSE,
                     dimnames = list(c("Gene1", "Gene2", "Gene3", "Gene4"), c("CellA", "CellB")))
  cell_type_a = "CellA"
  cell_type_b = "CellB"
  result <- integrate_network(
    vec_nst_a, vec_nst_b, mat_intra_a, mat_intra_b,
    cell_type_a, cell_type_b, mat_pem, mat_type, lrp
  )
  test_that("integrate_network returns correct structure", {
    result <- integrate_network(vec_nst_a, vec_nst_b, mat_intra_a, mat_intra_b,
                                cell_type_a, cell_type_b, mat_pem, mat_type, lrp)
    
    
    expect_true(all(c("nodes", "edges") %in% names(result)))
    
    
    
  })
  
})

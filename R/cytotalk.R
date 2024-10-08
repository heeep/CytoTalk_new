#' Main CytoTalk Pipeline
#'
#' @param lst_scrna List containing scRNA-seq data; for example, lists returned
#'   from `read_matrix_folder` or `read_matrix_with_meta`
#'
#' @param cell_type_a Name of cell type A that matches scRNA-seq file; for
#'   example, `"Fibroblasts"`
#'
#' @param cell_type_b Name of cell type B that matches scRNA-seq file; for
#'   example, `"LuminalEpithelialCells"`
#'
#' @param dir_out Folder used for output; if not specified, a "CytoTalk-output"
#'   folder will be generated
#'
#' @param cutoff_a Proportional threshold for lowly expressed genes in cell
#'   type A (range of \[0-1\]); for example, 0.2 means genes with some
#'   expression in at least 20% of cells are retained
#'
#' @param cutoff_b Proportional expression threshold for cell type B (range of
#'   \[0-1\])
#'
#' @param pcg A character vector, contains the names of protein coding genes;
#'   by default, uses the `pcg_human` data. This package also includes
#'   `pcg_mouse`, but you can also use your own data
#'
#' @param lrp A dataframe or matrix object with two columns, ligands names and
#'   the names of their receptors; by default, uses the `lrp_human` data. This
#'   package also includes `lrp_mouse`, but you can also use your own data
#'
#' @param beta_max Upper limit of the test values of the PCSF objective
#'   function parameter $I^2$, which is inversely proportional to the total
#'   number of genes in a given cell-type pair; suggested to be 100 (default)
#'   if the total number of genes in a given cell-type pair is above 10,000; if
#'   the total number of genes is below 5,000, increase to 500
#'
#' @param omega_min Start point of omega range; omega represents the edge cost
#'   of the artificial network, but has been found to be less significant than
#'   beta. Recommended minimum of `0.5`
#'
#' @param omega_max End point of range between `omega_min` and `omega_max`,
#'   step size of `0.1`. Recommended maximum of `1.5`
#'
#' @param depth Starting at each ligand-receptor pair in the resultant network,
#'   how many steps out from that pair should be taken to generate each
#'   neighborhood?
#'
#' @param ntrial How many random network subsets shall be created to get an
#'   empirical p-value for node prize and edge cost?
#'
#' @param cores How many cores to use for parallel processing?
#'
#' @param echo Should update messages be printed?
#' 
#' @param silent Should warning messages be exported? 
#'
#' @examples 
#' \donttest{
#' cell_type_a <- "Fibroblasts"
#' cell_type_b <- "LuminalEpithelialCells"
#' cutoff_a <- 0.6
#' cutoff_b <- 0.6
#' dir_in <- system.file("extdata", package = "CytoTalknew")
#' lst_scrna <- CytoTalknew::read_matrix_folder(dir_in)
#' # result <- CytoTalknew::run_cytotalk(lst_scrna,
#' #                                  cell_type_a, cell_type_b,
#' #                                  cutoff_a, cutoff_b,
#' #                                  cores = 2)
#' }
#'
#' @return A list containing model parameters, prefential expression measure,
#' the integrated co-expression network, the results of the PCST, and resulting
#' pathways from the final extracted network
#'
#' @export
run_cytotalk <- function(
    lst_scrna, cell_type_a, cell_type_b,
    cutoff_a=0.2, cutoff_b=0.2,
    pcg=pcg_human, lrp=lrp_human,
    beta_max=100, omega_min=0.5, omega_max=0.5,
    depth=3, ntrial=1000,
    cores=NULL, echo=TRUE, dir_out=NULL,silent=TRUE) {

    # save numeric parameters
    params <- list(
        cell_type_a = cell_type_a, cell_type_b = cell_type_b,
        cutoff_a = cutoff_a, cutoff_b = cutoff_b,
        beta_max = beta_max, omega_min = omega_min, omega_max = omega_max,
        depth = depth, ntrial = ntrial
    )

    # register parallel backend
    if (is.null(cores) || 1 < cores) {
        unregister_parallel()
        register_parallel(cores)
    }

    # create directory
    if (!is.null(dir_out) && !dir.exists(dir_out)) {
        dir.create(dir_out, recursive = TRUE)
    }

    if (echo) {
        tick(1, "Preprocessing...")
    }

    mat_pem <- pem(lst_scrna)
    mat_a <- extract_group(cell_type_a, lst_scrna)
    mat_b <- extract_group(cell_type_b, lst_scrna)
    vec_nst_a <- nonselftalk(mat_a, lrp)
    vec_nst_b <- nonselftalk(mat_b, lrp)
    mat_filt_a <- subset_rownames(subset_non_zero_old(mat_a, cutoff_a), pcg)
    mat_filt_b <- subset_rownames(subset_non_zero_old(mat_b, cutoff_b), pcg)

    # write out PEM matrix
    if (!is.null(dir_out)) {
        fpath <- file.path(dir_out, "PEM.txt")
        vroom_write_silent(mat_pem, fpath, rownames = TRUE,silent=silent)
    }

    if (echo) {
        tick(2, "Mutual information matrix...")
    }

    mat_disc_a <- discretize_sparse(Matrix::t(mat_filt_a))
    mat_disc_b <- discretize_sparse(Matrix::t(mat_filt_b))
    mat_mi_a <- mi_mat_parallel(mat_disc_a, method = "mm")
    mat_mi_b <- mi_mat_parallel(mat_disc_b, method = "mm")
    dimnames(mat_mi_a) <- list(colnames(mat_disc_a), colnames(mat_disc_a))
    dimnames(mat_mi_b) <- list(colnames(mat_disc_b), colnames(mat_disc_b))

    if (echo) {
        tick(3, "Indirect edge-filtered network...")
    }

    mat_intra_a <- Matrix::Matrix(parmigene::aracne.m(zero_diag(mat_mi_a)))
    mat_intra_b <- Matrix::Matrix(parmigene::aracne.m(zero_diag(mat_mi_b)))

    if (echo) {
        tick(4, "Integrate network...")
    }

    lst_net <- integrate_network(
        vec_nst_a, vec_nst_b, mat_intra_a, mat_intra_b,
        cell_type_a, cell_type_b, mat_pem, mat_a, lrp
    )

    # write out integrated nodes and edges
    if (!is.null(dir_out)) {
        fpath <- file.path(dir_out, "IntegratedNodes.txt")
        vroom_write_silent(lst_net$nodes, fpath)
        fpath <- file.path(dir_out, "IntegratedEdges.txt")
        vroom_write_silent(lst_net$edges, fpath, silent=silent)
    }

    if (echo) {
        tick(5, "PCSF...")
    }

    lst_pcst <- run_pcst(lst_net, beta_max, omega_min, omega_max)

    # write out PCST nodes and edges
    if (!is.null(dir_out)) {
        fpath <- file.path(dir_out, "PCSTNodeOccurance.txt")
        vroom_write_silent(lst_pcst$nodes, fpath,silent)
        fpath <- file.path(dir_out, "PCSTEdgeOccurance.txt")
        vroom_write_silent(lst_pcst$edges, fpath,silent=silent)
    }

    if (echo) {
        tick(6, "Determine best signaling network...")
    }

    df_test <- ks_test_pcst(lst_pcst, silent=silent)
    index <- order(as.numeric(df_test[, "pval"]))[1]
    beta <- df_test[index, "beta"]
    omega <- df_test[index, "omega"]

    # write out PCST scores
    if (!is.null(dir_out)) {
        vroom_write_silent(df_test, file.path(dir_out, "PCSTScores.txt"), silent=silent)
    }

    if (echo) {
        tick(7, "Generate network output...")
    }

    df_net <- extract_network(lst_net, lst_pcst, mat_pem, beta, omega)
    lst_path <- extract_pathways(df_net, cell_type_a, depth)
    lst_graph <- lapply(lst_path, graph_pathway)

    # no pathways found
    if (is.null(lst_path)) {
        result <- list(
            params = params,
            pem = mat_pem,
            integrated_net = lst_net,
            pcst = list(
                occurances = lst_pcst,
                ks_test_pval = df_test,
                final_network = df_net
            ),
            pathways = NULL
        )

        # unregister parallel backend
        if (is.null(cores) || 1 < cores) {
            unregister_parallel()
        }

        tick(8, "NOTE: No pathways found, analysis skipped!")
        return(result)
    }

    # write out pathways
    if (!is.null(dir_out)) {
        dir_path <- file.path(dir_out, "pathways")
        if (!dir.exists(dir_path)) {
            dir.create(dir_path, recursive = TRUE)
        }
        fnames <- names(lst_path)
        for (fn in fnames) {
            fpath <- file.path(dir_path, sprintf("%s.txt", fn))
            vroom_write_silent(lst_path[[fn]], fpath, silent=silent)
        }

        dir_gv <- file.path(dir_out, "graphviz")
        if (!dir.exists(dir_gv)) {
            dir.create(dir_gv, recursive = TRUE)
        }
        fnames <- names(lst_graph)
        for (fn in fnames) {
            fpath <- file.path(dir_gv, sprintf("%s.svg", fn))
            content <- DiagrammeRsvg::export_svg(lst_graph[[fn]])
            write(content, fpath)
        }
    }

    # write out final network
    if (!is.null(dir_out)) {
        write_network_sif(df_net, cell_type_a, dir_out)
        vroom_write_silent(df_net, file.path(dir_out, "FinalNetwork.txt"), silent=silent)
    }

    if (echo) {
        tick(8, "Analyze pathways...")
    }

    lst_pval <- lapply(
        lst_path, analyze_pathway, lst_net,
        cell_type_a, cell_type_b, beta, ntrial
    )

    # format the ligand and receptor cell types
    nodes <- do.call(rbind, strsplit(names(lst_path), "--"))
    df_pval <- do.call(rbind, apply(nodes, 1, function(x) {
        t <- ifelse(endsWith(x, "_"), cell_type_b, cell_type_a)
        x <- gsub("_$", "", x)
        data.frame(
            ligand = x[1], receptor = x[2],
            ligand_type = t[1], receptor_type = t[2]
        )
    }))

    # combine with scores and sort
    df_pval <- cbind(df_pval, do.call(rbind, lst_pval))
    df_pval <- df_pval[order(as.numeric(df_pval$pval_potential)), ]

    # write out analysis
    if (!is.null(dir_out)) {
        fpath <- file.path(dir_out, "PathwayAnalysis.txt")
        vroom_write_silent(df_pval, fpath, silent=silent)
    }

    # return out
    result <- list(
        params = params,
        pem = mat_pem,
        integrated_net = lst_net,
        pcst = list(
            occurances = lst_pcst,
            ks_test_pval = df_test,
            final_network = df_net
        ),
        pathways = list(
            raw = lst_path,
            graphs = lst_graph,
            df_pval = df_pval
        )
    )

    # write out session
    if (!is.null(dir_out)) {
        fpath <- file.path(dir_out, "CytoTalkSession.rda")
        save(result, file = fpath, version = 2)
    }

    # unregister parallel backend
    if (is.null(cores) || 1 < cores) {
        unregister_parallel()
    }

    result
}

#' @noRd
new_named_list <- function(mat, cell_types) {
  list(mat = mat, cell_types = cell_types)
}

#' @noRd
vroom_silent <- function(..., silent = TRUE) {
  if (silent) {
    suppressMessages(vroom::vroom(..., progress = FALSE))
  } else {
    vroom::vroom(..., progress = FALSE)
  }
}

#' @noRd
vroom_write_silent <- function(x, file, rownames = FALSE, silent = TRUE) {
  x <- as.data.frame(x)
  if (rownames) { x <- tibble::rownames_to_column(x) }
  if (silent) {
    suppressMessages(vroom::vroom_write(x, file, progress = FALSE))
  } else {
    vroom::vroom_write(x, file, progress = FALSE)
  }
}


#' @noRd
vroom_with_rownames <- function(..., row_names=1, silent=TRUE) {
    dat <- vroom_silent(..., silent=silent)
    tibble::column_to_rownames(dat, names(dat)[row_names])
}

#' @noRd
vroom_sparse_with_rownames <- function(..., row_names=1,silent=TRUE) {
    dat <- vroom_with_rownames(..., row_names = row_names, silent=silent)
    Matrix::Matrix(Matrix::as.matrix(dat), sparse = TRUE)
}

from_single_cell_experiment <- function(sce) {
    count <- SingleCellExperiment::logcounts(sce)
    names <- colnames(sce)
    new_named_list(count, names)
}

#' Read Folder with scRNAseq Data
#'
#' @param dpath The path of a directory, which contains scRNAseq matrices
#'
#' @param pattern A regular expression, matches scRNAseq filenames
#'
#' @param auto_transform Should count data be transformed if detected?
#' 
#' @param silent Should warning messages be exported? 
#'
#' @examples {
#' dir_in <- system.file("extdata", package = "CytoTalknew")
#' lst_scrna <- CytoTalknew::read_matrix_folder(dir_in)
#' table(lst_scrna$cell_types)
#' }
#'
#' @return A named list containing a sparse data matrix and cell type metadata
#'
#' @export
read_matrix_folder <- function(
    dpath, pattern=".*scRNAseq_(.+?)\\..+", auto_transform=TRUE, silent=TRUE) {

    # initial parameters
    mat <- NULL
    msg <- "not all rownames identical between input files"

    # determine filepaths
    fpaths <- dir_full(dpath, pattern = pattern)

    # read in all files
    for (fpath in fpaths) {
        cell_type <- gsub(pattern, "\\1", fpath)
        if (is.null(mat)) {
            # read in
            mat <- vroom_sparse_with_rownames(fpath,silent=silent)
            # start cell types vector
            cell_types <- rep(cell_type, ncol(mat))
        } else {
            # read in
            new <- vroom_sparse_with_rownames(fpath)
            # check for identical rownames
            errorifnot(identical(rownames(mat), rownames(new)), msg)
            # accumulate cell type names
            cell_types <- c(cell_types, rep(cell_type, ncol(new)))
            # combine all data
            mat <- cbind(mat, new)
        }
    }

    # check for count data
    mat <- check_count_data(mat, auto_transform)

    # return sparse matrix
    new_named_list(mat, cell_types)
}

#' Read scRNAseq Data Matrix and Metadata
#'
#' @param fpath_mat The path of a file containing a scRNAseq data matrix
#'
#' @param fpath_meta The path of a file contianing column metadata (cell types)
#'
#' @param auto_transform Should count data be transformed if detected?
#'
#' @examples {
#' base_path <- CytoTalknew::get_example_data()
#' fpath_mat <- file.path(base_path, "sample_counts.txt")
#' fpath_meta <- file.path(base_path, "sample_meta.txt")
#' lst_scrna <- CytoTalknew::read_matrix_with_meta(fpath_mat, fpath_meta)
#' table(lst_scrna$cell_types)
#' }
#'
#' @return A named list containing a sparse data matrix and cell type metadata
#'
#' @export
read_matrix_with_meta <- function(fpath_mat, fpath_meta, auto_transform=TRUE
                                  ,silent=TRUE) {
  # read in
  if(silent){
    suppressMessages(mat<- vroom::vroom(fpath_mat, 
                skip = 1,col_names = FALSE))
  }
  else{
    mat <- vroom::vroom(fpath_mat, skip = 1,col_names = FALSE)
    }
    mat <- tibble::column_to_rownames(mat, var = colnames(mat)[1])
    col_names <- vroom::vroom(fpath_mat,col_names = FALSE, n_max = 1)
    colnames(mat) <- unlist(col_names)
    meta <- vroom::vroom(fpath_meta,col_names = FALSE)
    # match meta to matrix
    index <- match(meta[[2]], colnames(mat))
    # ensure index matches exactly
    errorifnot(!any(is.na(index)), "meta file does not match matrix colnames")
    # otherwise, reorder
    cell_types <- as.character(meta$X2)
    # check for count data
    mat <- check_count_data(mat, auto_transform)

    # return sparse matrix
    new_named_list(mat, cell_types)
}

#' "Download the cellphoneDB data in your current directory from this package."
#' 
#' @return A file path which contains the example cellphone DB files. 
#' 
#' @examples {
#' cat(CytoTalknew::get_example_data())
#' }
#' @export
get_example_data <- function() {
  # Get the current working directory
  current_dir <- getwd()
  
  # Create a path for a new folder
  new_folder_path <- file.path(current_dir, "scRNAseq-data-cpdb")
  if (dir.exists(new_folder_path)) {
    # If the directory already exists, return the path and stop the function
    return(new_folder_path)
  }
  
  # If the directory does not exist, create it
  dir.create(new_folder_path)
  
  # Load the data from the specified package
  data("scrna_cpdb", package = "mycytotalk")
  
  # Ensure the mat object retains row and column names
  mat <- Matrix::as.matrix(scrna_cpdb$mat, sparse = FALSE)
  rownames(mat) <- rownames(scrna_cpdb$mat)
  colnames(mat) <- colnames(scrna_cpdb$mat)
  
  # Write the matrix to a file, preserving row and column names
  write.table(mat, file.path(new_folder_path, "sample_counts.txt"), sep = "\t", 
              quote = FALSE, row.names = TRUE, col.names = TRUE)
  
  # Save the cell_types vector as sample_meta.txt
  write.table(scrna_cpdb$cell_types, file.path(new_folder_path, 
              "sample_meta.txt"), sep = "\t", quote = FALSE, col.names = FALSE)
  
  # Return the path to the new folder
  return(new_folder_path)
}

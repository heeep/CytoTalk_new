#' @noRd
is_integer <- function(x) {
    identical(x, round(x))
}

#' @noRd
zero_diag <- function(mat) {
    diag(mat) <- 0
    mat
}

#' @noRd
img <- function(mat, ..., scl=0) {
    mat <- as.matrix(mat)
    diag(mat) <- diag(mat) * scl
    graphics::image(mat, ...)
}

#' @noRd
proportion_non_zero <- function(mat) {
    row_props <- Matrix::rowSums(mat != 0) / ncol(mat)
    return(row_props)
}

subset_non_zero <- function(mat, cutoff) {
    index <- (cutoff <= proportion_non_zero(mat))
    mat[index, ]
}

#' @noRd
subset_non_zero_old <- function(mat, cutoff) {
    thresh <- floor(ncol(mat) * cutoff)
    index <- thresh <= Matrix::rowSums(mat != 0)
    mat[index, ]
}

subset_rownames <- function(mat, labels) {
    index <- which(!is.na(cmatch(rownames(mat), labels)))
    mat[index, ]
}

#' @noRd
add_noise <- function(mat) {
    n <- nrow(mat)
    m <- ncol(mat)
    dummy <- c(1e-20, rep(0, m - 1))
    mat <- mat + matrix(replicate(n, sample(dummy)), n, byrow = TRUE)
    Matrix::Matrix(mat)
}

#' @noRd
extract_group_basic <- function(group, mat, labels) {
    mat[, labels == group, drop = FALSE]
}

extract_group <- function(group, lst) {
    extract_group_basic(group, lst[[1]], lst[[2]])
}

#' @noRd
group_meta_basic <- function(mat, labels) {
    groups <- sort(unique(labels))
    lst <- lapply(groups, extract_group_basic, mat, labels)
    names(lst) <- groups
    lst
}

group_meta <- function(lst) {
    group_meta_basic(lst[[1]], lst[[2]])
}

ungroup_meta <- function(lst) {
    ncols <- vapply(lst, ncol, numeric(1))
    cell_types <- unlist(lapply(seq_len(length(ncols)), function(i) {
        rep(names(ncols)[i], ncols[i])
    }))
    new_named_list(do.call(cbind, lst), cell_types)
}

match_lr_pairs <- function(mat, lrp) {
    # save a copy of the rownames
    hold <- rownames(mat)

    # match lr_pairs to mat rownames
    index <- data.frame(apply(lrp, 2, function(x) {
        match(toupper(x), toupper(hold))
    }))
    index <- index[rowSums(is.na(index)) == 0, ]

    # return out
    data.frame(
        ligand = hold[index[, 1]],
        receptor = hold[index[, 2]]
    )
}

normalize_sparse <- function(mat, scale.factor=10000) {
    log1p(Matrix::t(Matrix::t(mat) / Matrix::colSums(mat) * scale.factor))
}

check_count_data <- function(mat, auto_transform=TRUE) {
    check <- is_integer(mat)
    if (check) {
        if (auto_transform) {
            mat <- normalize_sparse(mat)
            # warn that normalization was performed
            msg <- paste(
                "count data detected;",
                "auto-transformed (see `?check_count_data`)"
            )
            warnifnot(!check, msg)
        } else {
            # warn if counts detected
            msg <- "count data detected, make sure to transform it"
            warnifnot(!check, msg)
        }
    }
    # return matrix
    mat
}

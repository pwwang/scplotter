#' @keywords internal
#' @importFrom utils getFromNamespace
check_columns <- getFromNamespace("check_columns", "plotthis")

#' @keywords internal
#' @importFrom utils getFromNamespace
combine_plots <- getFromNamespace("combine_plots", "plotthis")

#' Try to read a H5Group as a data.frame
#' @keywords internal
#' @param x H5Group object
#' @return A data.frame if successful
h5group_to_dataframe <- function(x) {
    if (!inherits(x, "H5Group")) {
        stop("Input must be an H5Group object.")
    }

    keys <- names(x)
    metadf <- stats::setNames(
        lapply(keys, function(k) {
            # correct categorical data
            if (inherits(x[[k]], "H5Group") && setequal(names(x[[k]]), c("categories", "codes"))) {
                # Make it more robust?
                # https://github.com/mojaveazure/seurat-disk/blob/877d4e18ab38c686f5db54f8cd290274ccdbe295/R/Convert.R#L308
                codes <- x[[k]][['codes']]$read()
                categories <- x[[k]][['categories']]$read()
                factor(categories[codes + 1], levels = categories)
            } else if (inherits(x[[k]], "H5Group")) {
                stop("Complex data structure detected in metadata (obs) column:", k)
            } else {
                x[[k]]$read()
            }
        }),
        keys
    )

    as.data.frame(metadf)
}

#' Try to read a H5Group as a matrix
#' @keywords internal
#' @param x H5Group object
#' @return A sparseMatrix if successful
h5group_to_matrix <- function(x) {
    if (!inherits(x, "H5Group")) {
        stop("Input must be an H5Group object.")
    }

    if (!'data' %in% names(x) || !'indices' %in% names(x) || !'indptr' %in% names(x)) {
        stop("H5Group does not contain the required structure (data, indices, indptr) for a matrix.")
    }
    if ('shape' %in% names(x)) {
        shape <- x[['shape']]$read()
    } else if ('shape' %in% hdf5r::h5attr_names(x)) {
        shape <- hdf5r::h5attr(x, 'shape')
    } else {
        stop("H5Group does not contain the 'shape' information.")
    }

    if (length(shape) != 2) {
        stop("Shape must be a vector of length 2 for a matrix.")
    }

    data <- x[['data']]$read()
    indices <- x[['indices']]$read()
    indptr <- x[['indptr']]$read()

    # Construct sparse matrix (dgCMatrix is default)
    Matrix::sparseMatrix(
        i = indices + 1,  # Convert to 1-based indexing
        p = indptr,  # Convert to 1-based indexing
        x = data,
        dims = rev(shape)
    )
}

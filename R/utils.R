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

#' Subset function with automatic layer repair for Seurat v5 objects
#'
#' This function handles corrupted layers in Seurat v5 assays that can cause
#' "incorrect number of dimensions" errors during subsetting.
#'
#' Strategy:
#' 1. Evaluate subset normally
#' 2. Check for corrupted layers proactively
#' 3. Repair any corrupted layers before subsetting
#' 4. Perform the subset operation
#' 5. If subset still fails, re-throw the error
#'
#' Example: subset_seurat(obj, EFS %in% c("EFS_L", "EFS_S"), nCount_RNA > 1000)
#' @param object Seurat object to subset
#' @param ... Arguments passed to `subset()`
#' @return Subsetted Seurat object with repaired layers if needed
#' @keywords internal
subset_seurat <- function(object, ...) {
  # Try subsetting first
  result <- tryCatch({
    subset(object, ...)
  }, error = function(e) {
    if (grepl("incorrect number of dimensions", e$message)) {
      warning("Subsetting failed due to corrupted layers. Attempting to repair...", call. = FALSE)

      # Scan and repair assays with multiple layers
      for (assay_name in names(object@assays)) {
        assay <- object@assays[[assay_name]]

        # Check if it's a v5 assay with layers
        if (inherits(assay, "Assay5") && length(assay@layers) > 1) {
          cat(sprintf("Scanning assay '%s' with %d layers...\n", assay_name, length(assay@layers)))

          # Identify intact and corrupted layers
          intact_layers <- list()
          corrupted_layers <- character()

          for (layer_name in names(assay@layers)) {
            layer <- assay@layers[[layer_name]]
            if (is.matrix(layer) || inherits(layer, "dgCMatrix")) {
              intact_layers[[layer_name]] <- layer
            } else {
              corrupted_layers <- c(corrupted_layers, layer_name)
            }
          }

          # If all layers are corrupted, remove the assay
          if (length(intact_layers) == 0) {
            warning(sprintf("All layers in assay '%s' are corrupted. Removing assay.", assay_name), call. = FALSE)
            object@assays[[assay_name]] <- NULL
          } else {
            # Restore corrupted layers based on intact ones
            if (length(corrupted_layers) > 0) {
              # Get dimensions from first intact layer
              template_layer <- intact_layers[[1]]
              n_features <- nrow(template_layer)
              n_cells <- ncol(template_layer)

              for (layer_name in corrupted_layers) {
                corrupted_data <- assay@layers[[layer_name]]

                # Check if it's a 1D vector with correct number of features
                if (is.vector(corrupted_data) && length(corrupted_data) == n_features) {
                  # Restore as single-column sparse matrix
                  warning(sprintf("Restoring corrupted layer '%s' in assay '%s' (1D vector -> sparse matrix)",
                                layer_name, assay_name), call. = FALSE)
                  object@assays[[assay_name]]@layers[[layer_name]] <-
                    Matrix::Matrix(corrupted_data, nrow = n_features, ncol = 1, sparse = TRUE)
                } else {
                  # Cannot restore, remove it
                  warning(sprintf("Cannot restore layer '%s' in assay '%s' (incompatible dimensions). Removing.",
                                layer_name, assay_name), call. = FALSE)
                  object@assays[[assay_name]]@layers[[layer_name]] <- NULL
                }
              }
            }
          }
        }
      }

      # Try subsetting again after repairs
      subset(object, ...)
    } else {
      # Re-throw other errors
      stop(e)
    }
  })

  return(result)
}

#' Call a function with a list of arguments
#'
#' Different from `do.call`, this is a faster version, especially when there are
#' big objects in the list.
#'
#' @param fn A function to call
#' @param args A list of arguments to pass to the function
#' @param quote Whether to quote the arguments
#' @param envir The environment to evaluate the function in
#' @return The result of the function call
#' @keywords internal
do_call <- function(fn, args, quote = FALSE, envir = parent.frame()) {
    # source: Gmisc
    # author: Max Gordon <max@gforge.se>

    if (quote) {
        args <- lapply(args, enquote)  # nocov
    }

    if (is.null(names(args)) ||
        is.data.frame(args)) {
        argn <- args
        args <- list()
    } else {
        # Add all the named arguments
        argn <- lapply(names(args)[names(args) != ""], as.name)
        names(argn) <- names(args)[names(args) != ""]
        # Add the unnamed arguments
        argn <- c(argn, args[names(args) == ""])
        args <- args[names(args) != ""]
    }

    if (inherits(fn, "character")) {
        if (is.character(fn)) {
            fn <- strsplit(fn, "[:]{2,3}")[[1]]
            fn <- if (length(fn) == 1) {
                get(fn[[1]], envir = envir, mode = "function")
            } else {
                get(fn[[2]], envir = asNamespace(fn[[1]]), mode = "function")
            }
        }
        call <- as.call(c(list(fn), argn))
    } else if (inherits(fn, "function")) {
        f_name <- deparse(substitute(fn))
        call <- as.call(c(list(as.name(f_name)), argn))
        args[[f_name]] <- fn
    } else if (inherits(fn, "name")) {  # nocov
        call <- as.call(c(list(fn, argn)))  # nocov
    }

    eval(call,
        envir = args,
        enclos = envir
    )
}

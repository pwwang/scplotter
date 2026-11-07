#' Look up an object in the scRepertoire namespace
#'
#' @description
#' Returns `NULL` instead of erroring when `name` is absent. scRepertoire
#' renamed several internal helpers between v2.2.1 and v2.9.1, and users may
#' have either version installed, so callers need a non-throwing probe.
#'
#' @param name Name of the object to look up.
#' @return The object, or `NULL` if scRepertoire does not define it.
#' @keywords internal
.screp_ns <- function(name) {
    tryCatch(
        utils::getFromNamespace(name, "scRepertoire"),
        error = function(e) NULL
    )
}

#' Get the first available scRepertoire internal
#'
#' @description
#' Tries candidate names in order of preference and returns the first one that
#' exists. Resolution is by name rather than by installed version because the
#' changes were renames, which no version floor can express: requiring >= 2.9.1
#' would drop working older installs, while requesting the old name would still
#' fail on newer ones.
#'
#' @param ... Candidate names, in order of preference.
#' @return The first object found.
#' @keywords internal
.screp_get <- function(...) {
    candidates <- c(...)
    for (candidate in candidates) {
        obj <- .screp_ns(candidate)
        if (!is.null(obj)) {
            return(obj)
        }
    }
    stop(
        "Could not find any of ",
        paste0("'", candidates, "'", collapse = ", "),
        " in the installed scRepertoire ",
        utils::packageVersion("scRepertoire"),
        ".",
        call. = FALSE
    )
}

#' Test whether an object is a Seurat or SummarizedExperiment object
#'
#' @description
#' v2.2.1 exposed `is_seurat_object()` and `is_se_object()` separately; v2.9.1
#' merged them into `.is.seurat.or.se.object()`. The unified predicate is
#' preferred when present, otherwise the two old predicates are OR'ed.
#'
#' @param x An object.
#' @return `TRUE` if `x` is a Seurat or SummarizedExperiment object.
#' @keywords internal
.screp_is_object <- function(x) {
    is_obj <- .screp_ns(".is.seurat.or.se.object")
    if (!is.null(is_obj)) {
        return(is_obj(x))
    }
    is_seurat_object <- .screp_ns("is_seurat_object")
    is_se_object <- .screp_ns("is_se_object")
    (!is.null(is_seurat_object) && is_seurat_object(x)) ||
        (!is.null(is_se_object) && is_se_object(x))
}

#' Minimum number of unique clonotypes across samples
#'
#' @description
#' v2.2.1 had `.short.check()`; v2.9.1 removed it. The fallback reproduces its
#' semantics exactly: the minimum, over samples, of the number of unique non-NA
#' clonotypes. The empty-input guard keeps `min()` from warning and returning
#' `Inf`.
#'
#' @param df A named list of per-sample data frames.
#' @param clone_call Column holding the clonotype call.
#' @return The minimum number of unique clonotypes across samples.
#' @keywords internal
.screp_short_check <- function(df, clone_call) {
    if (length(df) == 0) {
        return(0L)
    }
    short_check <- .screp_ns(".short.check")
    if (!is.null(short_check)) {
        return(short_check(df, clone_call))
    }
    min(vapply(
        df,
        function(d) length(which(!is.na(unique(d[[clone_call]])))),
        integer(1)
    ))
}

#' Resolve a scplotter diversity method to its scRepertoire implementation
#'
#' @description
#' v2.9.1 replaced the per-method dotted helpers (`.shannon`, `.invsimpson`,
#' ...) with the `.div.registry` named list, whose keys differ from scplotter's
#' `method` values (`ACE` -> `ace`, `gini.coeff` -> `gini`). Resolution is by
#' name rather than by version for the same reason as [.screp_get()].
#' scplotter's own `gini.coeff`, `d50` and `dXX` methods are handled before
#' this function is called.
#'
#' @param method A scplotter diversity method name.
#' @return The corresponding diversity function.
#' @keywords internal
.screp_div_fn <- function(method) {
    registry <- .screp_ns(".div.registry")
    if (!is.null(registry)) {
        keys <- c(
            shannon = "shannon", inv.simpson = "inv.simpson",
            norm.entropy = "norm.entropy", gini.simpson = "gini.simpson",
            chao1 = "chao1", ACE = "ace", gini.coeff = "gini"
        )
        key <- if (method %in% names(keys)) keys[[method]] else method
        div_fn <- registry[[key]]
        if (is.null(div_fn)) {
            stop(
                "scRepertoire ", utils::packageVersion("scRepertoire"),
                " has no diversity method '", key, "' in '.div.registry'.",
                call. = FALSE
            )
        }
        return(div_fn)
    }
    # v2.2.1 and earlier: one dotted helper per method, e.g. ".invsimpson"
    .screp_get(paste0(".", gsub(".", "", method, fixed = TRUE)))
}

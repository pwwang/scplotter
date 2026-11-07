#!/usr/bin/env Rscript
# Runs a fixed matrix of scplotter calls over Giotto mini objects and records
# per-call status plus an md5 of the built ggplot data.
# usage: Rscript tools/giotto-matrix.R <outdir>

args <- commandArgs(trailingOnly = TRUE)
outdir <- if (length(args)) args[1] else "giotto-matrix-out"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

suppressMessages({
    library(scplotter)
    library(ggplot2)
})

cat("== versions ==\n")
for (p in c("Giotto", "GiottoClass", "GiottoData", "GiottoVisuals", "GiottoUtils", "plotthis", "scplotter")) {
    cat(sprintf("%-14s %s\n", p, tryCatch(as.character(packageVersion(p)), error = function(e) "NA")))
}
cat("Giotto URL:", tryCatch(packageDescription("Giotto")$URL, error = function(e) NA), "\n\n")

first_categorical <- function(g) {
    md <- GiottoClass::getCellMetadata(g, output = "data.table")
    cand <- names(md)[vapply(md, function(x) is.character(x) || is.factor(x), logical(1))]
    setdiff(cand, c("cell_ID", "fov"))[1]
}

run <- function(id, expr) {
    res <- tryCatch({
        p <- force(expr)
        d <- ggplot_build(p)$data
        csv <- file.path(outdir, paste0(id, ".csv"))
        utils::write.csv(
            do.call(rbind, lapply(seq_along(d), function(i) {
                x <- as.data.frame(d[[i]])
                x$layer <- i
                x
            })),
            csv,
            row.names = FALSE
        )
        ggsave(file.path(outdir, paste0(id, ".png")), p, width = 7, height = 6, dpi = 100)
        list(id = id, ok = TRUE, rows = sum(vapply(d, nrow, integer(1))),
             md5 = unname(tools::md5sum(csv)), error = NA_character_)
    }, error = function(e) {
        list(id = id, ok = FALSE, rows = NA_integer_, md5 = NA_character_,
             error = conditionMessage(e))
    })
    cat(sprintf("%-42s %-5s %s\n", res$id, res$ok,
                ifelse(res$ok, res$md5, substr(res$error, 1, 140))))
    res
}

results <- list()
add <- function(id, expr) {
    results[[length(results) + 1L]] <<- run(id, expr)
    invisible(NULL)
}

minis <- c("visium", "cosmx", "vizgen", "seqfish", "starmap", "spatialgenomics")
for (m in minis) {
    g <- tryCatch(GiottoData::loadGiottoMini(m), error = function(e) NULL)
    if (is.null(g)) {
        add(paste0(m, ":load"), stop("load failed"))
        next
    }
    grp <- first_categorical(g)
    feats <- head(GiottoClass::featIDs(g, feat_type = "rna"), 2)
    if (!is.na(grp)) {
        add(paste0(m, ":SpatDimPlot"), SpatDimPlot(g, group_by = grp))
        add(paste0(m, ":SpatDimPlot+image"), SpatDimPlot(g, group_by = grp, image = TRUE))
        add(paste0(m, ":SpatDimPlot+shapes"), SpatDimPlot(g, group_by = grp, shapes = TRUE))
        add(paste0(m, ":SpatDimPlot+graph"), SpatDimPlot(g, group_by = grp, graph = TRUE))
        add(paste0(m, ":CellDimPlot"), CellDimPlot(g, group_by = grp))
        add(paste0(m, ":CellStatPlot"), CellStatPlot(g, ident = grp, plot_type = "bar"))
    }
    if (length(feats)) {
        add(paste0(m, ":SpatFeaturePlot"), SpatFeaturePlot(g, features = feats))
        add(paste0(m, ":FeatureStatPlot"), FeatureStatPlot(g, features = feats, plot_type = "violin", ident = grp))
        add(paste0(m, ":FeatureStatPlot:dim"), FeatureStatPlot(g, features = feats, plot_type = "dim", reduction = "umap"))
    }
}

df <- do.call(rbind, lapply(results, as.data.frame))
utils::write.csv(df, file.path(outdir, "summary.csv"), row.names = FALSE)
cat("\nwrote", file.path(outdir, "summary.csv"), "-", sum(df$ok), "of", nrow(df), "calls ok\n")

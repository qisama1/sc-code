function (object, ident.1 = NULL, ident.2 = NULL, group.by = NULL,
    subset.ident = NULL, assay = NULL, slot = "data", reduction = NULL,
    features = NULL, logfc.threshold = 0.25, test.use = "wilcox",
    min.pct = 0.1, min.diff.pct = -Inf, verbose = TRUE, only.pos = FALSE,
    max.cells.per.ident = Inf, random.seed = 1, latent.vars = NULL,
    min.cells.feature = 3, min.cells.group = 3, mean.fxn = NULL,
    fc.name = NULL, base = 2, densify = FALSE, ...)
{
    if (!is.null(x = group.by)) {
        if (!is.null(x = subset.ident)) {
            object <- subset(x = object, idents = subset.ident)
        }
        Idents(object = object) <- group.by
    }
    if (!is.null(x = assay) && !is.null(x = reduction)) {
        stop("Please only specify either assay or reduction.")
    }
    if (length(x = ident.1) == 0) {
        stop("At least 1 ident must be specified in `ident.1`")
    }
    if (is.null(x = reduction)) {
        assay <- assay %||% DefaultAssay(object = object)
        data.use <- object[[assay]]
        cellnames.use <- colnames(x = data.use)
    }
    else {
        data.use <- object[[reduction]]
        cellnames.use <- rownames(x = data.use)
    }
    cells <- IdentsToCells(object = object, ident.1 = ident.1,
        ident.2 = ident.2, cellnames.use = cellnames.use)
    if (!is.null(x = latent.vars)) {
        latent.vars <- FetchData(object = object, vars = latent.vars,
            cells = c(cells$cells.1, cells$cells.2))
    }
    norm.command <- paste0("NormalizeData.", assay)
    norm.method <- if (norm.command %in% Command(object = object) &&
        is.null(x = reduction)) {
        Command(object = object, command = norm.command, value = "normalization.method")
    }
    else {
        NULL
    }
    de.results <- FindMarkers(object = data.use, slot = slot,
        cells.1 = cells$cells.1, cells.2 = cells$cells.2, features = features,
        logfc.threshold = logfc.threshold, test.use = test.use,
        min.pct = min.pct, min.diff.pct = min.diff.pct, verbose = verbose,
        only.pos = only.pos, max.cells.per.ident = max.cells.per.ident,
        random.seed = random.seed, latent.vars = latent.vars,
        min.cells.feature = min.cells.feature, min.cells.group = min.cells.group,
        mean.fxn = mean.fxn, base = base, fc.name = fc.name,
        densify = densify, norm.method = norm.method, ...)
    return(de.results)
}

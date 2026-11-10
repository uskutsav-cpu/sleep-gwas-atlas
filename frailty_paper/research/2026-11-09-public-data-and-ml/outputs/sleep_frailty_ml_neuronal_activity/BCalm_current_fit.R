fit_elements <- function(object, normalize = TRUE, block = NULL, endomorphic = FALSE, normalizeSize = 1e9, plot=FALSE, ...) {
    if (!("endomorphic" %in% names(formals(mpralm))) && endomorphic) {
        warning(
            "The 'endomorphic' argument is not available in the version ",
            "of mpra you have installed. Please update mpra to use this ",
            "argument. Proceeding with endomorphic = FALSE."
        )
        endomorphic <- FALSE
    }
    design <- data.frame(rep(1, ncol(object)))
    if ("endomorphic" %in% names(formals(mpralm))) {
        mpralm_fit <- mpralm(
            object = object, design = design, aggregate = "none",
            normalize = normalize, model_type = "corr_groups",
            block = block, endomorphic = endomorphic, normalizeSize = normalizeSize, ...
        )
    } else {
        mpralm_fit <- mpralm(
            object = object, design = design, aggregate = "none",
            normalize = normalize, model_type = "corr_groups",
            block = block, ...
        )
    }
    if (!endomorphic) {
        mpralm_fit$logFC <- mpralm_fit$coefficients
        mpralm_fit$label <- getLabel(object)
    }
    return(mpralm_fit)
}

compute_logratio <- function(object, aggregate = c("mean", "sum", "none")) {
    .is_mpra_or_stop(object)

    aggregate <- match.arg(aggregate)

    if (aggregate == "sum") {
        dna <- getDNA(object, aggregate = TRUE)
        rna <- getRNA(object, aggregate = TRUE)
        logr <- log2(rna + 1) - log2(dna + 1)
    } else if (aggregate == "none") {
        dna <- getDNA(object, aggregate = FALSE)
        rna <- getRNA(object, aggregate = FALSE)
        logr <- log2(rna + 1) - log2(dna + 1)
    } else if (aggregate == "mean") {
        dna <- getDNA(object, aggregate = FALSE)
        rna <- getRNA(object, aggregate = FALSE)
        eid <- getEid(object)
        logr <- log2(rna + 1) - log2(dna + 1)

        by_out <- by(logr, eid, colMeans, na.rm = TRUE)
        logr <- do.call("rbind", by_out)
        rownames(logr) <- names(by_out)
    }
    return(logr)
}

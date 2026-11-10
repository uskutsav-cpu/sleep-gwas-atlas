plot_groups <- function(object, percentile = NULL, neg_label = NULL, test_label = NULL) {
    if (is(object, "MPRASet")) {
        object <- rowData(object)
    } else {
        object <- as.data.frame(object)
    }
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("The 'ggplot2' package is required but not installed. Please install it.")
    }
    if (utils::packageVersion("ggplot2") < "3.3.0") {
        stop("The 'ggplot2' package version must be at least 3.3.0. Please update it.")
    }
    if (is.null(object$logFC)) {
        stop("Your MPRASet does not contain logFC values. Please run mpralm or fit_elements first.")
    }
    if (is.null(object$label) && (!is.null(neg_label) || !is.null(test_label))) {
        stop("Your MPRASet should contain labels.")
    }
    if (!is.null(neg_label) && !neg_label %in% unique(object$label)) {
        stop("The negative label you provided is not in the label column of the mpra fit object.")
    }
    if (!is.null(test_label) && !test_label %in% unique(object$label)) {
        stop("The test label you provided is not in the label column of the mpra fit object.")
    }
    if (is.null(percentile) && !is.null(neg_label)) {
        percentile <- 0.95
        print("No percentile provided, using 0.95.")
    }

    if (!is.null(test_label) && !is.null(neg_label)) {
        object <- object[object$label %in% c(test_label, neg_label), ]
    }

    plot <- ggplot(object, aes(x = logFC, fill = label, y = after_stat(density))) +
        geom_histogram(alpha = 0.5, position = "identity", binwidth = 0.1) +
        geom_density(alpha = 0.2, adjust = 1) +
        theme_minimal() +
        labs(title = "Normalized Histogram of logratio Values", x = "Logratio", y = "Density", color = NULL) +
        xlim(c(min(object$logFC), max(object$logFC)))

    if (!is.null(neg_label)) {
        percentile_up <- quantile(object$logFC[object$label == neg_label], percentile)
        up_label <- paste(percentile, "th percentile of negative controls", sep = "")

        percentile_down <- quantile(object$logFC[object$label == neg_label], 1 - percentile)
        down_label <- paste(1 - percentile, "th percentile of negative controls", sep = "")

        # ggplot2 >= 3.4.0 uses `linewidth` for lines; older versions only know `size`
        line_width <- if (utils::packageVersion("ggplot2") >= "3.4.0") list(linewidth = 1) else list(size = 1)

        plot <- plot +
            do.call(geom_vline, c(list(aes(xintercept = percentile_up, color = up_label), linetype = "dashed"), line_width)) +
            do.call(geom_vline, c(list(aes(xintercept = percentile_down, color = down_label), linetype = "dashed"), line_width)) +
            scale_color_manual(
                values = setNames(c("green", "orange"), c(up_label, down_label)),
                guide = guide_legend(override.aes = list(linetype = "dashed"))
            )
    }

    return(plot)
}

# 	Moderated t-statistics relative to a logFC threshold.
# 	Davis McCarthy, Gordon Smyth, adapted by Pia Keukeleire from original function in limma package.
# 	This version shifts the mean of the coefficients to the mean of
#   the negative controls in order to work with negative upper thresholds.
# 	percentile: The percentile of the negative controls to use as the threshold.
# 	25 November 2024.
#
# 	Null hypothesis: the true (centred) effect lies in [lfc_left, lfc_right], where lfc_right and
# 	lfc_left are the `percentile` and `1 - percentile` quantiles of the centred negative controls.
# 	For skewed controls |lfc_left| != lfc_right, so the two sides of the null are tested separately.
# 	With se = stdev.unscaled * sqrt(s2.post) and C the centred coefficient, the p-value is
# 	  C >= mid: P(T >= (C - lfc_right)/se) + P(T >= (C - lfc_left)/se)
# 	  C <  mid: P(T <= (C - lfc_left)/se)  + P(T <= (C - lfc_right)/se)
# 	with mid = (lfc_left + lfc_right)/2 (0 for symmetric thresholds) and T ~ t(df.total).
# 	The first term tests the side the coefficient lies on; the p-value is always <= 1.
# 	With symmetric thresholds (lfc_right = -lfc_left) this equals limma::treat.
# 	t: (|C| - threshold of its own side)/se, positive for up and negative for down, 0 inside the null.
# 	Fixed 3 November 2026: negative coefficients were previously tested against the wrong thresholds.
mpra_treat <- function(fit, percentile = 0.975, neg_label, trend = FALSE, robust = FALSE, winsor.tail.p = c(0.05, 0.1)) {
    # 	Check fit
    if (is(fit, "MPRASet")) {
        mpra <- attr(fit, "MArrayLM")
        mpra$logFC <- rowData(fit)$logFC
        mpra$label <- getLabel(fit)
        fit <- mpra
    }
    if (!is(fit, "MArrayLM")) stop("fit must be an MArrayLM object")
    if (is.null(fit$coefficients)) stop("coefficients not found in fit object")
    if (is.null(fit$stdev.unscaled)) stop("stdev.unscaled not found in fit object")
    if (is.null(fit$label)) stop("Your mpra fit object should contain a label column.")

    fit$lods <- NULL

    neg_mean <- mean(fit$coefficients[fit$label == neg_label])

    coefficients <- as.matrix(fit$coefficients - neg_mean)
    coefficients_neg <- as.matrix(fit$coefficients[fit$label == neg_label] - neg_mean)

    stdev.unscaled <- as.matrix(fit$stdev.unscaled)
    sigma <- fit$sigma
    df.residual <- fit$df.residual
    if (is.null(coefficients) || is.null(stdev.unscaled) || is.null(sigma) || is.null(df.residual)) {
        stop("No data, or argument is not a valid lmFit object")
    }
    if (max(df.residual) == 0) {
        stop("No residual degrees of freedom in linear model fits")
    }
    if (!any(is.finite(sigma))) {
        stop("No finite residual standard deviations")
    }
    if (trend) {
        covariate <- fit$Amean
        if (is.null(covariate)) stop("Need Amean component in fit to estimate trend")
    } else {
        covariate <- NULL
    }
    sv <- squeezeVar(sigma^2, df.residual, covariate = covariate, robust = robust, winsor.tail.p = winsor.tail.p)
    fit$df.prior <- sv$df.prior
    fit$s2.prior <- sv$var.prior
    fit$s2.post <- sv$var.post
    df.total <- df.residual + sv$df.prior
    df.pooled <- sum(df.residual, na.rm = TRUE)
    df.total <- pmin(df.total, df.pooled)
    fit$df.total <- df.total

    se <- stdev.unscaled * sqrt(fit$s2.post)
    lfc_right <- quantile(coefficients_neg, percentile)
    lfc_left <- quantile(coefficients_neg, 1 - percentile)
    # The null is lfc_left <= coef <= lfc_right, which is asymmetric for skewed controls.
    # As in limma::treat, the main term tests against the threshold on the coefficient's own
    # side and the second term adds the (tiny) probability mass of the opposite side.
    # "Own side" is decided relative to the midpoint of the null interval (0 for limma::treat);
    # this keeps p <= 1 and is also correct when lfc_right < 0 for extremely skewed controls.
    is_up <- coefficients >= (lfc_right + lfc_left) / 2
    thr_same <- ifelse(is_up, lfc_right, -lfc_left)
    thr_opp <- ifelse(is_up, -lfc_left, lfc_right)
    # coefficient expressed as distance in the direction of its own side
    side_coef <- ifelse(is_up, coefficients, -coefficients)
    tstat.same <- (side_coef - thr_same) / se
    tstat.opp <- (side_coef + thr_opp) / se
    fit$t <- array(0, dim(coefficients), dimnames = dimnames(coefficients))
    fit$p.value <- pt(tstat.same, df = df.total, lower.tail = FALSE) + pt(tstat.opp, df = df.total, lower.tail = FALSE)
    tstat.same <- pmax(tstat.same, 0)
    fc.up <- (coefficients > lfc_right)
    fc.down <- (coefficients < lfc_left)
    fit$t[fc.up] <- tstat.same[fc.up]
    fit$t[fc.down] <- -tstat.same[fc.down]
    fit$treat.lfc_right <- lfc_right
    fit$treat.lfc_left <- lfc_left
    fit
}

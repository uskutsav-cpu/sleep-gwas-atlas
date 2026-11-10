plot_groups <- function(object, percentile=NULL, neg_label=NULL, test_label=NULL) {
	if (is(object, "MPRASet")) {
		object <- rowData(object)
	} else {
		object <- as.data.frame(object)
	}
	if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("The 'ggplot2' package is required but not installed. Please install it.")
    }
	if (is.null(object$logFC)) {
		stop("Your MPRASet does not contain logFC values. Please run mpralm or fit_elements first.")
	}
	if (is.null(object$label) && (! is.null(neg_label) || ! is.null(test_label))) {
		stop("Your MPRASet should contain labels.")
	}
	if (! is.null(neg_label) && ! neg_label %in% unique(object$label)) {
		stop("The negative label you provided is not in the label column of the mpra fit object.")
	}
	if (! is.null(test_label) && ! test_label %in% unique(object$label)) {
		stop("The test label you provided is not in the label column of the mpra fit object.")
	}
	if (is.null(percentile) && ! is.null(neg_label)) {
		percentile <- 0.95
		print("No percentile provided, using 0.95.")
	}

	if (! is.null(test_label) && ! is.null(neg_label)) {
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
		up_label <- paste(percentile, "th percentile of negative controls", sep="")

		percentile_down <- quantile(object$logFC[object$label == neg_label], 1 - percentile)
		down_label <- paste(1 - percentile, "th percentile of negative controls", sep="")

		plot <- plot + 
			geom_vline(aes(xintercept = percentile_up, color = up_label), linetype = "dashed", size = 1) +
			geom_vline(aes(xintercept = percentile_down, color = down_label), linetype = "dashed", size = 1) +
			scale_color_manual(values = setNames(c("green", "orange"), c(up_label, down_label)), 
				guide = guide_legend(override.aes = list(linetype = "dashed"))) 
	}

	return(plot)
}


mpra_treat <- function(fit, percentile=0.975, neg_label, trend=FALSE, robust=FALSE, winsor.tail.p=c(0.05,0.1))
#	Moderated t-statistics relative to a logFC threshold.
#	Davis McCarthy, Gordon Smyth, adapted by Pia Keukeleire from original function in limma package.
#	This version shifts the mean of the coefficients to the mean of the negative controls in order to work with negative upper thresholds.
#	percentile: The percentile of the negative controls to use as the threshold.
#	25 November 2024.
{
#	Check fit
	if (is(fit, "MPRASet")) {
		mpra <- attr(fit, "MArrayLM")
		mpra$logFC <- rowData(fit)$logFC
		mpra$label <- getLabel(fit)
		fit <- mpra
	}
	if(!is(fit,"MArrayLM")) stop("fit must be an MArrayLM object")
	if(is.null(fit$coefficients)) stop("coefficients not found in fit object")
	if(is.null(fit$stdev.unscaled)) stop("stdev.unscaled not found in fit object")
	if (is.null(fit$label)) stop("Your mpra fit object should contain a label column.")
	
	fit$lods <- NULL

	neg_mean <- mean(fit$coefficients[fit$label == neg_label])

	coefficients <- as.matrix(fit$coefficients - neg_mean)
	coefficients_neg <- as.matrix(fit$coefficients[fit$label == neg_label] - neg_mean)

	stdev.unscaled <- as.matrix(fit$stdev.unscaled)
	sigma <- fit$sigma
	df.residual <- fit$df.residual
	if (is.null(coefficients) || is.null(stdev.unscaled) || is.null(sigma) || 
		is.null(df.residual)) 
		stop("No data, or argument is not a valid lmFit object")
	if (max(df.residual) == 0) 
		stop("No residual degrees of freedom in linear model fits")
	if (!any(is.finite(sigma))) 
		stop("No finite residual standard deviations")
	if(trend) {
		covariate <- fit$Amean
		if(is.null(covariate)) stop("Need Amean component in fit to estimate trend")
	} else {
		covariate <- NULL
	}
	sv <- squeezeVar(sigma^2, df.residual, covariate=covariate, robust=robust, winsor.tail.p=winsor.tail.p)
	fit$df.prior <- sv$df.prior
	fit$s2.prior <- sv$var.prior
	fit$s2.post <- sv$var.post
	df.total <- df.residual + sv$df.prior
	df.pooled <- sum(df.residual,na.rm=TRUE)
	df.total <- pmin(df.total,df.pooled)
	fit$df.total <- df.total

	acoef <- abs(coefficients)
	se <- stdev.unscaled*sqrt(fit$s2.post)
	lfc_right <- quantile(coefficients_neg, percentile)
	lfc_left <- quantile(coefficients_neg, 1-percentile)
	tstat.right <- (acoef-lfc_right)/se
	tstat.left <- (acoef-lfc_left)/se
	fit$t <- array(0,dim(coefficients),dimnames=dimnames(coefficients))
	fit$p.value <- pt(tstat.right, df=df.total,lower.tail=FALSE) + pt(tstat.left,df=df.total,lower.tail=FALSE)
	tstat.right <- pmax(tstat.right,0)
	tstat.left <- pmax(tstat.left,0)
	fc.up <- (coefficients > lfc_right)
	fc.down <- (coefficients < lfc_left)
	fit$t[fc.up] <- tstat.right[fc.up]
	fit$t[fc.down] <- tstat.left[fc.down]
	fit$treat.lfc_right <- lfc_right
	fit$treat.lfc_left <- lfc_left
	fit
}

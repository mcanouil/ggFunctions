#' qqplot create a qqplot
#'
#' @param pvalue [vector(numeric) or matrix(numeric)]: A matrix, data.frame or vector which contains pvalues.
#' @param lambdaNames [vector(character)]: Names of lambda GC, if several colums is given in 'pvalue'.
#' @param pt.size [numeric]: Point size.
#' @param bw [logical]: Turn the Grey theme into a a Black&White theme.
#' @param noGrid [logical]: A grid in the background of the plot should be drawn.
#' @param base_size [numeric]: Font size.
#' @return A qqplot in ggplot2 format.
#' @export
# @examples
# qqplot()
qqplot <- function (pvalue, lambdaNames = NULL, pt.size = 1, bw = TRUE, noGrid = FALSE, base_size = 12) {
    if (!(is.matrix(pvalue) | is.data.frame(pvalue))) {
        pvalue <- matrix(pvalue)
    } else {}

    .ggplotColours <- function(n = 6, h = c(0, 360) + 15) {
        if ((diff(h)%%360) < 1) {
            h[2] <- h[2] - 360/(n-1)
        } else {}
        return(c("black", hcl(h = (seq(h[1], h[2], length = (n-1))), c = 100, l = 65)))
    }
    if (is.null(lambdaNames) & any(colnames(pvalue)=="")) {
        lambdaNames <- paste0("lambda", seq(ncol(pvalue)))
    } else {}

    # labs <- NULL
    # res <- NULL
    # for (i in seq(ncol(pvalue))) {
        # pv <- pvalue[, i]
        # X2 <- qnorm(pv, lower.tail = FALSE)^2
        # gc <- median(X2, na.rm = TRUE)/qchisq(0.5, df = 1)
        # obspval <- sort(pv)
        # logobspval <- -(log10(obspval))
        # exppval <- c(1:length(obspval))
        # logexppval <- -(log10((exppval-0.5)/length(exppval)))
        # labs <- c(labs, bquote(lambda[gc]^.(lambdaNames[i]) == .(round(gc, 4))))
        # res <- rbind(res, data.frame(logexppval, logobspval, i))
    # }

    # res <- do.call("rbind", mclapply(seq(ncol(pvalue)), mc.cores = detectCores(), function (i) {
    res <- do.call("rbind", lapply(seq(ncol(pvalue)), function (i) {
        pv <- pvalue[, i]
        X2 <- qnorm(pv, lower.tail = FALSE)^2
        gc <- median(X2, na.rm = TRUE)/qchisq(0.5, df = 1)
        obspval <- sort(pv)
        logobspval <- -(log10(obspval))
        exppval <- c(1:length(obspval))
        logexppval <- -(log10((exppval-0.5)/length(exppval)))
        # labs <- bquote(lambda[gc]^.(lambdaNames[i]) == .(round(gc, 4)))
        labs <- paste0(lambdaNames[i], "=", round(gc, 4))
        tmp <- data.frame(logexppval, logobspval, i, labs)
        whichInfinite <- apply(res, 1, function (iRow) {
            return(any(sapply(as.numeric(iRow[-3]), is.infinite)))
        })
        return(tmp[!whichInfinite, ])
    }))
    p <- ggplot(data = res)
    if (bw) {
        blackwhite <- function (base_size = 12, base_family = "", noGrid = FALSE) {
            if (noGrid) {
                noGridColour <- c("transparent", "transparent") # "white"
            } else {
                noGridColour <- c("grey90", "grey98")
            }
            theme_grey(base_size = base_size, base_family = base_family) %+replace%
                theme(
                    axis.text = element_text(size = rel(0.8)),
                    axis.ticks = element_line(colour = "black"),
                    legend.background = element_rect(fill = "white", colour = "black"),
                    legend.key = element_rect(fill = "white", colour = "black"),
                    legend.position = "right",
                    legend.justification = "center",
                    legend.box = NULL,
                    panel.background = element_rect(fill = "white", colour = NA),
                    panel.border = element_rect(fill = NA, colour = "black"),
                    panel.grid.major = element_line(colour = noGridColour[1], size = 0.2),
                    panel.grid.minor = element_line(colour = noGridColour[2], size = 0.5),
                    strip.background = element_rect(fill = "grey80", colour = "black", size = 0.2)
                )
        }
        p <- p + blackwhite(base_size = base_size, noGrid = noGrid)
    } else {}
    p <- p + geom_abline(intercept = 0, slope = 1)
    p <- p + geom_point(aes_string(x = "logexppval", y = "logobspval", colour = "i"), size = pt.size, shape = 1)
    p <- p + labs(x = bquote(Expected -log[10](P[value])), y = bquote(Observed -log[10](P[value])))
    if (ncol(pvalue) > length(c("dodgerblue", "firebrick2", "springgreen3", "maroon2", "goldenrod2", "deepskyblue"))) {
        p <- p + scale_colour_manual(name = element_blank(), breaks = seq(ncol(pvalue)), labels = labs, values = .ggplotColours(ncol(pvalue)))
    } else {
        p <- p + scale_colour_manual(name = element_blank(), breaks = seq(ncol(pvalue)), labels = labs, values =  c("dodgerblue", "firebrick2", "springgreen3", "maroon2", "goldenrod2", "deepskyblue"))
    }
    p <- p + labs(title = "Q-Q plot") + theme(plot.title = element_text(lineheight = 0.8, face = "bold"))
    axisLim <- range(pretty_breaks()(range(unlist(lapply(seq(ncol(pvalue)), function (i) {
        range(res[res[, "i"]%in%i, c(1, 2)])
    })))))
    # axisLim <- range(pretty_breaks()(range(unlist(mclapply(seq(ncol(pvalue)), mc.cores = detectCores(), function (i) {
        # range(res[res[, "i"]%in%i, c(1, 2)])
    # })))))
    p <- p + xlim(axisLim) + ylim(axisLim)
    return(p)
}
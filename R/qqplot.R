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
    if (is.null(lambdaNames)) {
        lambdaNames <- colnames(pvalue)
    } else {}

    obsmax <- NA
    expmax <- NA
    labs <- NULL
    res <- NULL
    for (i in seq(ncol(pvalue))) {
        pv  <- pvalue[, i]
        X2  <- qnorm(pv, lower.tail = FALSE)^2
        gc  <- median(X2, na.rm = TRUE)/qchisq(0.5, df = 1)
        obspval <- sort(pv)
        logobspval <- -(log10(obspval))
        exppval <- c(1:length(obspval))
        logexppval <- -(log10((exppval-0.5)/length(exppval)))

        obsmax <- max(obsmax, trunc(max(logobspval, na.rm = TRUE))+1, na.rm = TRUE)
        expmax <- max(expmax, trunc(max(logexppval, na.rm = TRUE))+1, na.rm = TRUE)
        labs <- c(labs, bquote(lambda[gc]^.(lambdaNames[i]) == .(round(gc, 4))))
        res <- rbind(res, data.frame(logexppval, logobspval, i))
    }

    res[, "i"] <- factor(res[, "i"], levels = seq(ncol(pvalue)))
    whichInfinite <- apply(res, 1, function (iRow) {
        return(any(sapply(as.numeric(iRow[-3]), is.infinite)))
    })
    res <- res[!whichInfinite, ]
    p <- ggplot(data = res)
    if (bw) {
        bw <- function (base_size = 12, base_family = "", noGrid = FALSE) {
            if (noGrid) {
                noGridColour <- "transparent"
            } else {
                noGridColour <- c("gray90", "grey95")
            }
            theme(
                line = element_line(colour = "black", size = 0.5, linetype = 1, lineend = "butt"),
                rect = element_rect(fill = "white", colour = "black", size = 0.5, linetype = 1),
                text = element_text(family = base_family, face = "plain", colour = "black", size = base_size, hjust = 0.5, vjust = 0.5, angle = 0, lineheight = 0.9),
                axis.text = element_text(size = rel(0.8), colour = "black"),
                strip.text = element_text(size = rel(0.8)),
                axis.line = element_blank(),
                axis.text.x = element_text(vjust = 1),
                axis.text.y = element_text(hjust = 1),
                axis.ticks = element_line(colour = "black"),
                axis.title.x = element_text(),
                axis.title.y = element_text(angle = 90),
                axis.ticks.length = unit(0.15, "cm"),
                axis.ticks.margin = unit(0.1, "cm"),
                legend.background = element_rect(fill = "white", colour = "black"),
                legend.margin = unit(0.2, "cm"),
                legend.key = element_rect(fill = "white", colour = "black"),
                legend.key.size = unit(1.2, "lines"),
                legend.key.height = NULL,
                legend.key.width = NULL,
                legend.text = element_text(size = rel(0.8)),
                legend.text.align = NULL,
                legend.title = element_text(size = rel(0.8), face = "bold", hjust = 0),
                legend.title.align = NULL,
                legend.position = "right",
                legend.direction = NULL,
                legend.justification = "center",
                legend.box = NULL,
                panel.background = element_rect(fill = "white", colour = "black"),
                panel.border = element_blank(),
                panel.grid.major = element_line(colour = noGridColour[1]),
                panel.grid.minor = element_line(colour = noGridColour[length(noGridColour)], size = 0.25),
                panel.margin = unit(0.25, "lines"),
                strip.background = element_rect(fill = "black", colour = "black"),
                strip.text.x = element_text(colour = "white"),
                strip.text.y = element_text(angle = -90, colour = "white"),
                plot.background = element_rect(colour = "white"),
                plot.title = element_text(size = rel(1.2)),
                plot.margin = unit(c(1, 1, 0.5, 0.5), "lines"),
                complete = TRUE
            )
        }
        p <- p + bw(base_size = base_size, noGrid = noGrid)
    } else {}
    p <- p + geom_point(aes_string(x = "logexppval", y = "logobspval", colour = "i"), size = pt.size)
    p <- p + geom_abline(intercept = 0, slope = 1)
    p <- p + labs(x = bquote(Expected -log[10](P[value])), y = bquote(Observed -log[10](P[value])))
    if (ncol(pvalue) > length(c("dodgerblue", "firebrick2", "springgreen3", "maroon2", "goldenrod2", "deepskyblue"))) {
        p <- p + scale_colour_manual(name = element_blank(), breaks = seq(ncol(pvalue)), labels = labs, values = .ggplotColours(ncol(pvalue)))
    } else {
        p <- p + scale_colour_manual(name = element_blank(), breaks = seq(ncol(pvalue)), labels = labs, values =  c("dodgerblue", "firebrick2", "springgreen3", "maroon2", "goldenrod2", "deepskyblue"))
    }
    p <- p + labs(title = "Q-Q plot") + theme(plot.title = element_text(lineheight = 0.8, face = "bold"))
    axisLim <- c(0, max(p$data[, as.character(p$layers[[1]]$mapping[c(1, 2)])]))
    p <- p + xlim(axisLim) + ylim(axisLim)
    if (which.max(apply(p$data[, c("logexppval", "logobspval")], 2, max))==2) {
         p <- p + theme(legend.justification = c(1, 0), legend.position = c(1, 0))
    } else {
        p <- p + theme(legend.justification = c(0, 1), legend.position = c(0, 1))
    }
    return(p)
}
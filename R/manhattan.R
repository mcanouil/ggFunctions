#' manhattan create a Manhanttan plot
#'
#' @param data [data.frame]: A 'data.frame' which include some columns: chromosome, position and pvalues.
#' @param chr [character]: Name of the 'chr' (chromosome) column in 'data'.
#' @param position [character]: Name of the 'position' column in 'data'.
#' @param y [character]: Name of the 'p-value' column in 'data' (or anything you want to plot with a manhattan plot).
#' @param title [character]: The text for the title.
#' @param xlab [character]: The text for the x axis.
#' @param ylab [character]: The text for the y axis.
#' @param sep [numeric]: Width of the space between bars.
#' @return A Manhattan plot in ggplot2 format.
#' @export
# @examples
# manhattan()
manhattan <- function (data, chr, position, y, title = "Manhattan plot", xlab = "Chromosomes", ylab = "P-Value", sep = 0.02) {
    data[, chr] <- gsub("chr", "", tolower(data[, chr]))
    data[, chr] <- factor(toupper(data[, chr]), levels = c(seq(22), "X", "Y"))
    data <- data[order(data[, chr], data[, position]), ]
    notNA <- apply(data[, c(chr, position)], 1, function (iRow) {any(is.na(iRow))})
    data <- data[which(!notNA), ]
    chrSize <- table(data[, chr])
    if (length(chrSize)!=24 | any(chrSize==0)) {
        CHR <- c(seq(22), "X", "Y")
        equalZero <- names(which(chrSize==0))
        notIn <- CHR[which(!CHR %in% names(chrSize))]
        chr2Add <- unique(c(equalZero, notIn))
        newLines <- data.frame(do.call("rbind", lapply(chr2Add, function (iChr) {
            newLine <- matrix(as.numeric(rep(NA, ncol(data))), ncol = ncol(data), dimnames = list(NULL, colnames(data)))
            newLine <- data.frame(newLine)
            newLine[1, chr] <- iChr
            return(newLine)
        })))
        colnames(newLines) <- colnames(data)
        data <- rbind(data, newLines)
        data <- data[order(data[, chr], data[, position]), ]
    } else {}
    chrSize <- chrSize[chrSize!=0]
    data <- data[!is.na(data[, y]), ]
    chrSizeNew <- table(data[, chr])
    chrSizeNew <- chrSizeNew[chrSizeNew!=0]
    chrStep <- floor(sum(chrSizeNew) * sep)
    data[, "xPos"] <- unlist(sapply(seq_along(chrSizeNew), function (iSize) {
        if (chrSizeNew[iSize]!=0) {
            xPos <- seq_len(chrSizeNew[iSize])
            range(xPos)
            if (iSize>1) {
                xPos <- xPos + sum(chrSizeNew[seq(iSize-1)]) + chrStep*(iSize-1)
                range(xPos)
            } else {}
            return(xPos)
        } else {}
    }), use.names = FALSE)
    avoidZero <- rep(0, length(chrSize))
    avoidZero[which(chrSize==0)] <- chrStep
    whichIsCenter <- ceiling(c(cumsum(chrSizeNew) - diff(c(0, cumsum(chrSizeNew)))/2))
    xBreaks <- data[whichIsCenter, "xPos"]
    p <- ggplot(data = data, aes_string(x = "xPos", y = y, colour = chr)) +
        geom_hline(yintercept = 0) +
        geom_point(size = 1.5, shape = 1, na.rm = TRUE) +
        scale_x_continuous(
            breaks = xBreaks,
            labels = names(chrSize),
            limits = c(min(data[, "xPos"]), max(data[, "xPos"])+sum(avoidZero)),
            expand = c(0.01, 0.01)
        ) +
        labs(title = title, y = ylab, x = xlab)
    return(p)
}
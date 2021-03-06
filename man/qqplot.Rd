% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/qqplot.R
\name{qqplot}
\alias{qqplot}
\title{qqplot create a qqplot}
\usage{
qqplot(pvalue, lambdaNames = NULL, pt.size = 1)
}
\arguments{
\item{pvalue}{[vector(numeric) or matrix(numeric)]: A matrix, data.frame or vector which contains pvalues.}

\item{lambdaNames}{[vector(character)]: Names of lambda GC, if several colums is given in 'pvalue'.}

\item{pt.size}{[numeric]: Point size.}
}
\value{
A qqplot in ggplot2 format.
}
\description{
qqplot create a qqplot
}


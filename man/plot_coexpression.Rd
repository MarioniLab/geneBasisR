% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/visualization_functions.R
\name{plot_coexpression}
\alias{plot_coexpression}
\title{plot_coexpression}
\usage{
plot_coexpression(sce, genes, title = NULL, ...)
}
\arguments{
\item{sce}{SingleCellExperiment object containing gene counts matrix (stored in 'logcounts' assay).}

\item{genes}{Character vector containing gene names to evaluate for co-expression.}

\item{title}{String to be passed as a title. Default title=NULL.}

\item{...}{Additional arguments to pass.}
}
\value{
Heatmap for co-expression.
}
\description{
plot_coexpression
}

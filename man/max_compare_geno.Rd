% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/compare_geno.R
\name{max_compare_geno}
\alias{max_compare_geno}
\alias{max.compare_geno}
\title{Find pair with most similar genotypes}
\usage{
max_compare_geno(object, ...)

\method{max}{compare_geno}(object, ...)
}
\arguments{
\item{object}{A square matrix with genotype comparisons for pairs
of individuals, as output by \code{\link[=compare_geno]{compare_geno()}}.}

\item{...}{Ignored}
}
\value{
Data frame with individual pair, proportion matches, number
of mismatches, number of matches, and total markers genotyped.
}
\description{
From results of \code{\link[=compare_geno]{compare_geno()}}, show the pair with most similar genotypes.
}
\examples{
grav2 <- read_cross2(system.file("extdata", "grav2.zip", package="qtl2geno"))
cg <- compare_geno(grav2)
max(cg)
}
\keyword{utilities}

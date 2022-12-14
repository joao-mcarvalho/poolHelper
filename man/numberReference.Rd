% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions_pool.R
\name{numberReference}
\alias{numberReference}
\title{Compute the number of reference reads at multiple loci}
\usage{
numberReference(genotypes, indContribution, error)
}
\arguments{
\item{genotypes}{is a list, where each entry corresponds to a different
locus. Each entry should be a matrix containing the genotypes (coded as 0,
1 or 2). Each column of that matrix should be a different site and each row
a different individual.}

\item{indContribution}{either a list or a matrix (that the function will
convert to a list). Each list entry should be a matrix of individual
contributions.  Each row of that matrix is a different individual and each
column is a different site. Thus, each entry of the matrix should contain
the number of reads contributed by that individual at that particular site.}

\item{error}{a numeric value with error rate associated with the sequencing
and mapping process. This error rate is assumed to be symmetric:
error(reference -> alternative) = error(alternative -> reference). This
number should be between 0 and 1.}
}
\value{
a list with one entry per locus. Each of those entries is a matrix
with the number of reference allele reads contributed by each individual.
Each matrix row represents a different individual and each column is a
different site.
}
\description{
This function computes the number of reference reads over multiple loci and
for a single population.
}
\details{
Note that this function will also work on a single locus, provided that the
input is in the list format.
}
\keyword{internal}

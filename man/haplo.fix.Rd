% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions_pool.R
\name{haplo.fix}
\alias{haplo.fix}
\title{Create invariable sites}
\usage{
haplo.fix(haplotypes, nHap)
}
\arguments{
\item{haplotypes}{a list of haplotypes obtained from the simulations done
with \code{\link[scrm:scrm]{scrm::scrm()}}. Each entry of the list is a matrix that corresponds to
a given locus. At each matrix, each column is a different site and each row
is a different haplotype.}

\item{nHap}{an integer representing the total number of haplotypes simulated.}
}
\value{
a list of haplotypes identical to \code{haplotypes}, but without empty
loci.
}
\description{
This function applies a correction for the situations where \code{\link[scrm:scrm]{scrm::scrm()}}
does not produce a single polymorphic site for a given locus. In this
situation, two artificial sites are created at that locus. All individuals
are assumed to be homozygous for the reference allele at those sites.
}
\keyword{internal}

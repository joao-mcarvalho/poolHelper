% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions_pool.R
\name{Ifreqs}
\alias{Ifreqs}
\title{Compute allele frequencies from genotypes}
\usage{
Ifreqs(nDip, genotypes)
}
\arguments{
\item{nDip}{an integer representing the total number of diploid individuals
to simulate. Note that \code{\link[scrm:scrm]{scrm::scrm()}} actually simulates haplotypes, so the
number of simulated haplotypes is double of this.}

\item{genotypes}{a list of simulated genotypes, where each entry is a matrix
corresponding to a different locus. At each matrix, each column is a
different SNP and each row is a different individual.}
}
\value{
a list of allele frequencies. Each entry of the list corresponds to a
different locus.
}
\description{
Computes alternative allele frequencies from genotypes by dividing the total
number of alternative alleles by the total number of gene copies.
}
\examples{
genotypes <- run_scrm(nDip = 10, nloci = 10)
Ifreqs(nDip = 10, genotypes)

}

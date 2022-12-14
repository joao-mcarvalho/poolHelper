% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions_pool.R
\name{numberReferencePop}
\alias{numberReferencePop}
\title{Compute the number of reference reads for multiple populations}
\usage{
numberReferencePop(genotypes, indContribution, size, error)
}
\arguments{
\item{genotypes}{either a list with a single entry (one locus) or a matrix
(that the function will convert to a list) containing the genotypes (coded
as 0, 1 or 2). Each column of that matrix should be a different site and
each row a different individual.}

\item{indContribution}{a list where each entry contains the information for a
single population. Each entry should be a matrix, with as many rows as the
number of individuals of that population. Each row contains the number of
contributed reads for a given individual and across all sites.}

\item{size}{a list with one entry per population. Each entry should be a
vector containing the size (in number of diploid individuals) of each pool.
Thus, if a population was sequenced using a single pool, the vector should
contain only one entry. If a population was sequenced using two pools, each
with 10 individuals, this vector should contain two entries and both will
be 10.}

\item{error}{a numeric value with error rate associated with the sequencing
and mapping process. This error rate is assumed to be symmetric:
error(reference -> alternative) = error(alternative -> reference). This
number should be between 0 and 1.}
}
\value{
a list with one entry per population. Each entry contains the number
of reference allele reads for the individuals of that population and for
that locus. Different individuals are in different rows and each columns
represents a different site.
}
\description{
This function computes the number of reference reads over a single locus for
multiple populations.
}
\details{
Note that this function will not work as intended if the input consists of
multiple loci.
}
\examples{
# simulate coverage at 5 SNPs for two populations, assuming 20x mean coverage
reads <- simulateCoverage(mean = c(20, 20), variance = c(100, 100), nSNPs = 5, nLoci = 1)

# simulate the number of reads contributed by each individual
# for each population there are two pools, each with 5 individuals
indContribution <- popsReads(list_np = rep(list(rep(5, 2)), 2), coverage = reads, pError = 5)

# set seed and create a random matrix of genotypes for the 20 individuals - 10 per population
set.seed(10)
genotypes <- matrix(rpois(100, 0.5), nrow = 20)

# simulate the number of reference reads for the two populations
numberReferencePop(genotypes = genotypes, indContribution = indContribution,
size = rep(list(rep(5, 2)), 2), error = 0.01)

}

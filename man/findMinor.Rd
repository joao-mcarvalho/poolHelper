% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions_pool.R
\name{findMinor}
\alias{findMinor}
\title{Define major and minor alleles}
\usage{
findMinor(reference, alternative, coverage)
}
\arguments{
\item{reference}{is a matrix of reference allele reads. Each row of the
matrix should be a different population and each column a different site.
Thus, each entry of the matrix contains the number of observed reads with
the reference allele for that population at a given site.}

\item{alternative}{is a matrix of alternative allele reads. Each row of the
matrix should be a different population and each column a different site.
Thus, each entry of the matrix contains the number of observed reads with
the alternative allele for that population at a given site.}

\item{coverage}{is a matrix of total coverage. Each row of the matrix should
be a different population and each column a different site. Thus, each
entry of the matrix contains the total number of observed reads for that
population at a given site.}
}
\value{
a list with three names entries

\item{major}{a list with one entry per locus. Each entry is a matrix with
the number of major allele reads for each population. Each column
represents a different site and each row a different population.}

\item{minor}{a list with one entry per locus. Each entry is a matrix with
the number of minor allele reads for each population. Each column
represents a different site and each row a different population.}

\item{total}{a list with one entry per locus. Each entry is a matrix with
the coverage of each population. Each column represents a different site
and each row a different population.}
}
\description{
This function checks which of the two simulated alleles (reference or
alternative) corresponds to the minor allele. This function can also be used
to remove sites according to a minor-allele reads threshold.
}
\details{
More precisely, this function counts the number of reads with the reference
or alternative allele at each site and then sets the minor allele as the
least frequent of the two. This is done across all populations and so the
major and minor alleles are defined at a global level. Then if the
\code{min.minor} input is not NA, sites where the number of minor allele reads,
across all populations, is below the user-defined threshold are removed.
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
readsReference <- numberReferencePop(genotypes = genotypes, indContribution = indContribution,
size = rep(list(rep(5, 2)), 2), error = 0.01)

# create Pooled DNA sequencing data for these two populations and for a single locus
pools <- poolPops(nPops = 2, nLoci = 1, indContribution = indContribution,
readsReference = readsReference)

# define the major and minor alleles for this Pool-seq data
# we have to select the first entry of the pools list because this function works for matrices
findMinor(reference = pools$reference[[1]], alternative = pools$alternative[[1]],
coverage = pools$total[[1]])

}

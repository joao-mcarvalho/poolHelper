% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions_pool.R
\name{simReads}
\alias{simReads}
\title{Simulate coverage at a single locus}
\usage{
simReads(mean, variance, nSNPs = NA, genotypes = NA)
}
\arguments{
\item{mean}{an integer that defines the mean depth of coverage to simulate.
Please note that this represents the mean coverage across all sites. If a
vector is supplied instead, the function assumes that each entry of the
vector is the mean for a different population.}

\item{variance}{an integer that defines the variance of the depth of coverage
across all sites. If a vector is supplied instead, the function assumes
that each entry of the vector is the variance for a different population.}

\item{nSNPs}{an integer representing the number of polymorphic sites per
locus to simulate. This is an optional input but either this or the
\code{genotypes} matrix must be supplied.}

\item{genotypes}{a matrix of simulated genotypes, where each column is a
different SNP and each row is a different individual. This is an optional
input but either this or the \code{nSNPs} must be supplied.}
}
\value{
a matrix with the total coverage per population and per site.
Different rows represent different populations and each column is a
different site.
}
\description{
Simulates the total number of reads, for each polymorphic site of a given
locus using a negative binomial distribution.
}
\details{
The total number of reads is simulated with a negative binomial and according
to a user-defined mean depth of coverage and variance. This function is
intended to work with a matrix of genotypes, simulating the depth of coverage
for each site present in the genotypes. However, it can also be used to
simulate coverage distributions independent of genotypes, by choosing how
many sites should be simulated (with the \code{nSNPs} option).
}
\examples{
# coverage for one population at 10 sites
simReads(mean = 20, variance = 100, nSNPs = 10)

# simulate coverage at one locus with 10 SNPs for two populations:
# the first with 100x and the second with 50x
simReads(mean = c(100, 50), variance = c(250, 150), nSNPs = 10)

}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions_pool.R
\name{remove_by_reads}
\alias{remove_by_reads}
\title{Apply a coverage-based filter over a list}
\usage{
remove_by_reads(nLoci, reads, minimum, maximum, genotypes = NA)
}
\arguments{
\item{nLoci}{an integer that represents how many independent loci were
simulated.}

\item{reads}{a list with the total depth of coverage. Each entry of the list
should be a matrix corresponding to a different locus. Each row of that
matrix should be the coverage of a different population and each column a
different site.}

\item{minimum}{an integer representing the minimum coverage allowed. Sites
where any population has a depth of coverage below this threshold are
removed from the data.}

\item{maximum}{an integer representing the maximum coverage allowed. Sites
where any population has a depth of coverage above this threshold are
removed from the data.}

\item{genotypes}{an optional list input with the genotypes. Each entry of the
list should be a matrix corresponding to a different locus. Each column of
the matrix should be a different site and each row a different individual.}
}
\value{
a list with the total depth of coverage similar to the \code{reads} input
argument but without sites where the coverage was below the \code{minimum} or
above the \code{maximum}. If the genotypes were included, a second list entry
will also be included in the output, containing the genotypes minus the
sites that were removed.
}
\description{
This function removes sites that have a coverage below a \code{minimum} value and
sites with a coverage above a \code{maximum} value. This is done over multiple
loci, assuming that each entry of the \code{reads} list is a different locus. If a
list of genotypes is also supplied, then those same sites are also removed
from each locus of the genotypes.
}
\examples{
set.seed(10)

# simulate coverage for 10 locus
reads <- simulateCoverage(mean = c(25, 25), variance = c(200, 200), nSNPs = 10, nLoci = 10)

# remove sites with coverage below 10x or above 100x
reads <- remove_by_reads(nLoci = 10, reads = reads, minimum = 5, maximum = 100)
# notice that some locus no longer have 10 SNPs - those sites were removed
reads

}

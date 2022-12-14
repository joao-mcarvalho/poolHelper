% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions_pool.R
\name{splitMatrix}
\alias{splitMatrix}
\title{Split matrix of genotypes}
\usage{
splitMatrix(matrix, size)
}
\arguments{
\item{matrix}{is obviously a matrix, ideally one containing genotypes. Each
column of the matrix should be a different site and each row a different
individual.}

\item{size}{a list with one entry per population. Each entry should be a
vector containing the size (in number of diploid individuals) of each pool.
Thus, if a population was sequenced using a single pool, the vector should
contain only one entry. If a population was sequenced using two pools, each
with 10 individuals, this vector should contain two entries and both will
be 10.}
}
\value{
a list with one entry per entry of the size input argument. Each
entry contains the information of the original matrix for the number of
individuals specified by the corresponding entry of the size input
argument.
}
\description{
This function splits a matrix into different list entries. The matrix is
split according to a set of row indexes defined by the \code{size} input.
}
\details{
The \code{size} input is utilized to create the index of the rows that go into the
different list entries. It specifies the size, in terms of number of
individuals, of each population.
}
\examples{
set.seed(10)

# create a random matrix
mymatrix <- matrix(round(runif(n = 50, min = 10, max = 25)), nrow = 10)

# split that matrix assuming 8 individuals in the first population and two in the second
splitMatrix(matrix = mymatrix, size = list(8, 2))

}
\keyword{internal}

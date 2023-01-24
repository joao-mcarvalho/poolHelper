#' Simulate a single population
#'
#' Simulates the evolution of biological sequences for a single population with
#' variable theta values.
#'
#' @param nDip an integer representing the total number of diploid individuals
#'   to simulate. Note that [scrm::scrm()] actually simulates haplotypes, so the
#'   number of simulated haplotypes is double of this.
#' @param nloci is an integer that represents how many independent loci should
#'   be simulated.
#' @param theta a value for the mutation rate assuming theta = 4Nu, where u is
#'   the neutral mutation rate per locus.
#'
#' @return a list with genotypes. Each entry of the list corresponds to a
#'   different locus. For each locus, the genotypes are in a matrix, with each
#'   row representing a different individual and each column a different site.
#'
#' @examples
#' run_scrm(nDip = 100, nloci = 10)
#' run_scrm(nDip = 100, nloci = 10, theta = 5)
#'
#' @export
run_scrm <- function(nDip, nloci, theta = 10) {

  # run scrm to obtain haplotypes
  haplo <- scrm::scrm(paste(nDip*2, nloci, "-t", theta))

  # combine those haplotypes into genotypes
  genotypes <- GetGenotypes(haplotypes = haplo$seg_sites, nDip = nDip)

  # output the genotypes
  genotypes
}


#' Compute allele frequencies from genotypes
#'
#' Computes alternative allele frequencies from genotypes by dividing the total
#' number of alternative alleles by the total number of gene copies.
#'
#' @param nDip an integer representing the total number of diploid individuals
#'   to simulate. Note that [scrm::scrm()] actually simulates haplotypes, so the
#'   number of simulated haplotypes is double of this.
#' @param genotypes a list of simulated genotypes, where each entry is a matrix
#'   corresponding to a different locus. At each matrix, each column is a
#'   different SNP and each row is a different individual.
#'
#' @return a list of allele frequencies. Each entry of the list corresponds to a
#'   different locus.
#'
#' @examples
#' genotypes <- run_scrm(nDip = 10, nloci = 10)
#' Ifreqs(nDip = 10, genotypes)
#'
#' @export
Ifreqs <- function(nDip, genotypes) {

  # perform a sum over all columns to obtain the number of alternative alleles per site
  alternative <- lapply(genotypes, colSums)
  # compute allele frequencies by dividing the number of alternative alleles by the total number of alleles
  ifreqs <- lapply(alternative, function(l) l/(nDip*2))

  # output the allele frequencies computed directly from genotypes
  ifreqs
}


#' Apply a minor allele reads threshold
#'
#' Removes sites where the total number of minor-allele reads is below a certain
#' threshold.
#'
#' If a site has less minor-allele reads than \code{min.minor} across all
#' populations, that site is removed from the data.
#'
#' @param freqs a vector of allele frequencies where each entry corresponds to a
#'   different site.
#' @param alternative a matrix with the number of reads with the alternative
#'   allele. Each row should be a different population and each column a
#'   different site.
#' @param coverage a matrix with the total coverage. Each row should be a
#'   different population and each column a different site.
#' @param minor a matrix with the number of minor-allele reads. Each row should
#'   be a different population and each column a different site.
#' @param min.minor is an integer representing the minimum allowed number of
#'   minor-allele reads. Sites that, across all populations, have less
#'   minor-allele reads than this threshold will be removed from the data.
#'
#' @return a list with three named entries:
#'
#'   \item{freqs}{is a vector with the allele frequencies minus the frequency of
#'   the removed sites.}
#'
#'   \item{alternative}{is a matrix with the number of alternative-allele reads
#'   per site, minus any removed sites.}
#'
#'   \item{coverage}{is a matrix with the depth of coverage minus the coverage
#'   of the removed sites.}
#'
#' @examples
#' # create a vector of allele frequencies
#' freqs <- runif(20)
#'
#' set.seed(10)
#' # create a matrix with the number of reads with the alternative allele
#' alternative <- matrix(sample(x = c(0,5,10), size = 20, replace = TRUE), nrow = 1)
#' # create a matrix with the depth of coverage
#' coverage <- matrix(sample(100:150, size = 20), nrow = 1)
#' # the number of reads with the reference allele is obtained by subtracting
#' # the number of alternative allele reads from the depth of coverage
#' reference <- coverage - alternative
#'
#' # find the minor allele at each site
#' minor <- findMinor(reference = reference, alternative = alternative, coverage = coverage)
#' # keep only the matrix with the minor allele reads
#' minor <- minor[["minor"]]
#'
#' # remove sites where the number of minor-allele reads is below the threshold
#' removeSites(freqs = freqs, alternative = alternative, coverage = coverage,
#' minor = minor, min.minor = 2)
#'
#' @export
removeSites <- function(freqs, alternative, coverage, minor, min.minor) {

  # check if the input is correct - freqs should always be supplied as a vector
  if(!inherits(freqs, "numeric"))
    stop(paste("freqs should be supplied on a numeric vector, with each entry corresponding to a site. Please check"))

  # get the total number of minor allele reads in the data
  tminor <- colSums(minor)
  # find out in which columns the total sum of the reads with the minor allele is below the threshold
  toremove <- tminor < min.minor

  # if there are sites where the sum of the reads with the minor allele is below the threshold
  if(sum(toremove) != 0) {

    # remove those columns from the matrix containing the depth of coverage
    coverage <- coverage[, !toremove, drop = FALSE]
    # remove those columns from the matrix containing the number of reads with the alternative allele
    alternative <- alternative[, !toremove, drop = FALSE]

    # remove those entries from the vector containing the allele frequencies computed directly from genotypes
    freqs <- freqs[!toremove, drop = FALSE]
  }

  # create the output containing the frequencies computed from the genotypes
  # and the number of Pool-seq alternative allele reads and total coverage
  out <- list(freqs = freqs, alternative = alternative, coverage = coverage)
  # output the results of the function
  out
}


#' Compute allele frequencies from pooled sequencing data
#'
#' Computes the frequency of the alternative allele in Pool-seq data and removes
#' any site with too few minor-allele reads from both the pool frequencies and
#' the frequencies computed directly from genotypes.
#'
#' The frequency at a given SNP is calculated according to: `pi = c/r`, where c
#' = number of alternative allele reads and r = total number of observed reads.
#' Additionally, if a site has less minor-allele reads than \code{min.minor}
#' across all populations, that site is removed from the data.
#'
#' @param reference a matrix with the number of reference allele reads. Each row
#'   should be a different population and each column a different site.
#' @param alternative a matrix with the number of alternative allele reads. Each
#'   row should be a different population and each column a different site.
#' @param coverage a matrix with the total coverage. Each row should be a
#'   different population and each column a different site.
#' @param min.minor is an integer representing the minimum allowed number of
#'   minor-allele reads. Sites that, across all populations, have less
#'   minor-allele reads than this threshold will be removed from the data.
#' @param ifreqs a vector of allele frequencies computed directly from the
#'   genotypes where each entry corresponds to a different site.
#'
#'
#' @return a list with two entries. The \code{ifreqs} entry contains the allele
#'   frequencies computed directly from genotypes and \code{pfreqs} the allele
#'   frequencies computed from pooled sequencing data.
#'
#' @examples
#' set.seed(10)
#' # create a vector of allele frequencies
#' freqs <- runif(20)
#' set.seed(10)
#' # create a matrix with the number of reads with the alternative allele
#' alternative <- matrix(sample(x = c(0,5,10), size = 20, replace = TRUE), nrow = 1)
#' # create a matrix with the depth of coverage
#' coverage <- matrix(sample(100:150, size = 20), nrow = 1)
#' # the number of reads with the reference allele is obtained by subtracting
#' # the number of alternative allele reads from the depth of coverage
#' reference <- coverage - alternative
#' # compute allele frequencies from pooled sequencing data
#' Pfreqs(reference = reference, alternative = alternative, coverage = coverage,
#' min.minor = 2, ifreqs = freqs)
#'
#' @export
Pfreqs <- function(reference, alternative, coverage, min.minor, ifreqs) {

  # if the minimum number of minor allele reads is not set to zero
  if(min.minor != 0) {

    # check which of the two simulated alleles (reference or alternative) corresponds to the minor allele
    minor <- findMinor(reference = reference, alternative = alternative, coverage = coverage)
    # keep only the matrix with the number of minor allele reads
    minor <- minor[["minor"]]

    # use the removeSites function to remove sites with less than `min.minor` minor allele reads
    temp <- removeSites(freqs = ifreqs, alternative = alternative, coverage = coverage,
                        minor = minor, min.minor = min.minor)

    # get the number of reads with the alternative allele
    alternative <- temp[["alternative"]]
    # and the total coverage
    coverage <- temp[["coverage"]]

    # get the frequencies computed from genotypes
    ifreqs <- temp[["freqs"]]
  }

  # compute the allele frequencies for Pool-seq data
  freqs <- alternative/coverage
  # create the output in this instance - with the allele frequencies computed from genotypes and from Pool-seq data
  freqs <- list(ifreqs = ifreqs, pfreqs = freqs)

  # output the allele frequencies
  freqs
}


#' Average absolute difference between allele frequencies
#'
#' Calculates the average absolute difference between the allele frequencies
#' computed directly from genotypes and from pooled sequencing data.
#'
#' Different combinations of parameters can be tested to check the effect of the
#' various parameters. The average absolute difference is computed with the
#' \link[Metrics]{mae} function, assuming the frequencies computed directly from
#' the genotypes as the \code{actual} input argument and the frequencies from
#' pooled data as the \code{predicted} input argument.
#'
#' @param nDip an integer representing the total number of diploid individuals
#'   to simulate. Note that [scrm::scrm()] actually simulates haplotypes, so the
#'   number of simulated haplotypes is double of this.
#' @param nloci is an integer that represents how many independent loci should
#'   be simulated.
#' @param pools a list with a vector containing the size (in number of diploid
#'   individuals) of each pool. Thus, if a population was sequenced using a
#'   single pool, the vector should contain only one entry. If a population was
#'   sequenced using two pools, each with 10 individuals, this vector should
#'   contain two entries and both will be 10.
#' @param pError an integer representing the value of the error associated with
#'   DNA pooling. This value is related with the unequal contribution of both
#'   individuals and pools towards the total number of reads observed for a
#'   given population - the higher the value the more unequal are the individual
#'   and pool contributions.
#' @param sError a numeric value with error rate associated with the sequencing
#'   and mapping process. This error rate is assumed to be symmetric:
#'   error(reference -> alternative) = error(alternative -> reference). This
#'   number should be between 0 and 1.
#' @param mCov an integer that defines the mean depth of coverage to simulate.
#'   Please note that this represents the mean coverage across all sites.
#' @param vCov an integer that defines the variance of the depth of coverage
#'   across all sites.
#' @param min.minor is an integer representing the minimum allowed number of
#'   minor-allele reads. Sites that, across all populations, have less
#'   minor-allele reads than this threshold will be removed from the data.
#' @param minimum an optional integer representing the minimum coverage allowed.
#'   Sites where the population has a depth of coverage below this threshold are
#'   removed from the data.
#' @param maximum an optional integer representing the maximum coverage allowed.
#'   Sites where the population has a depth of coverage above this threshold are
#'   removed from the data.
#'
#' @return a data.frame with columns detailing the number of diploid
#'   individuals, the pool error, the number of pools, the number of individuals
#'   per pool, the mean coverage, the variance of the coverage and the average
#'   absolute difference between the frequencies computed from genotypes and
#'   from pooled data.
#'
#' @examples
#' # single population sequenced with a single pool of 100 individuals
#' maePool(nDip = 100, nloci = 10, pools = list(100), pError = 100, sError = 0.01,
#' mCov = 100, vCov = 250, min.minor = 2)
#'
#' # single population sequenced with two pools, each with 50 individuals
#' maePool(nDip = 100, nloci = 10, pools = list(c(50, 50)), pError = 100, sError = 0.01,
#' mCov = 100, vCov = 250, min.minor = 2)
#'
#' # single population sequenced with two pools, each with 50 individuals
#' # removing sites with coverage below 10x or above 180x
#' maePool(nDip = 100, nloci = 10, pools = list(c(50, 50)), pError = 100, sError = 0.01,
#' mCov = 100, vCov = 250, min.minor = 2, minimum = 10, maximum = 180)
#'
#' @export
maePool <- function(nDip, nloci, pools, pError, sError, mCov, vCov, min.minor, minimum = NA, maximum = NA) {

  # run SCRM and obtain genotypes for a single population
  genotypes <- run_scrm(nDip = nDip, nloci = nloci)

  # simulate number of reads
  reads <- simulateCoverage(mean = mCov, variance = vCov, genotypes = genotypes)

  # if the minimum coverage is defined
  if(!is.na(minimum)) {

    # check if the maximum coverage is also defined
    if(is.na(maximum))
      stop("please define the maximum coverage")

    # remove sites with a depth of coverage above or below the defined threshold
    reads <- remove_by_reads(nLoci = nloci, reads, minimum = minimum, maximum = maximum, genotypes = genotypes)

    # get the genotypes - without sites simulated with a coverage below or above the threshold
    genotypes <- lapply(reads, function(locus) locus[[2]])

    # check the dimensions of the matrices with the genotypes
    dimensions <- matrix(unlist(lapply(genotypes, dim)), ncol = 2, byrow = TRUE)
    # we only wish to keep the locus where we have at least one polymorphic site
    tokeep <- dimensions[, 2] != 0
    # remove all loci without polymorphic sites
    genotypes <- genotypes[tokeep]

    # get the reads - without sites simulated with a coverage below or above the threshold
    reads <- lapply(reads, function(locus) locus[[1]])
    # use the same index to remove entries of the reads list that correspond to locus without polymorphic sites
    reads <- reads[tokeep]
    # ensure that each entry is a matrix
    reads <- lapply(reads, function(locus) matrix(locus, nrow = 1))
  }

  # compute allele frequencies directly from the genotypes
  ifreqs <- Ifreqs(nDip = nDip, genotypes = genotypes)

  # simulate individual contribution to the total number of reads
  indContribution <- lapply(1:nloci, function(locus)
    popsReads(list_np = pools, coverage = reads[[locus]], pError = pError))

  # simulate the number of reference reads
  reference <- lapply(1:nloci, function(locus)
    numberReferencePop(genotypes = genotypes[[locus]], indContribution = indContribution[[locus]],
                       size = pools, error = sError))

  # simulate pooled sequencing data
  pool <- poolPops(nPops = 1, nLoci = nloci, indContribution = indContribution, readsReference = reference)

  # compute the allele frequencies obtained with pooled sequencing
  pfreqs <- lapply(1:nloci, function(locus)
    Pfreqs(reference = pool[["reference"]][[locus]], alternative = pool[["alternative"]][[locus]],
           coverage = pool[["total"]][[locus]], min.minor = min.minor, ifreqs = ifreqs[[locus]]))

  # get the allele frequencies computed directly from genotypes after removing sites that did not pass the threshold
  ifreqs <- lapply(pfreqs, `[[`, 1)
  # get the allele frequencies computed from Pool-seq after removing sites that did not pass the threshold
  pfreqs <- lapply(pfreqs, `[[`, 2)

  # compute the mean absolute error for each locus
  abs_error <- lapply(1:nloci, function(locus) Metrics::mae(actual = ifreqs[[locus]], predicted = pfreqs[[locus]]))
  # replace any NaN values with the mean of the remaining values
  abs_error[is.na(abs_error)] <- mean(unlist(abs_error), na.rm = TRUE)

  # get the number of pools used to sequence a population
  nPools <- length(pools[[1]])
  # get the number of individuals per pool
  indsPool <- unique(pools[[1]])

  # create a dataframe with the values for this particular combination of parameters
  out <- data.frame(nDip = nDip, PoolError = pError, nPools = nPools, indsPool = indsPool, mean = mCov, var = vCov,
                    absError = unlist(abs_error))

  # output the results of the function
  out
}


#' Average absolute difference between allele frequencies computed from
#' genotypes and from Pool-seq data
#'
#' Calculates the average absolute difference between the allele frequencies
#' computed directly from genotypes and from pooled sequencing data.
#'
#' The average absolute difference is computed with the \link[Metrics]{mae}
#' function, assuming the frequencies computed directly from the genotypes as
#' the \code{actual} input argument and the frequencies from pooled data as the
#' \code{predicted} input argument.
#'
#' Note that this functions allows for different combinations of parameters.
#' Thus, the effect of different combinations of parameters on the average
#' absolute difference can be tested. For instance, it is possible to check what
#' is the effect of different coverages by including more than one value in the
#' \code{mCov} input argument. This function will run and compute the average
#' absolute difference for all combinations of the \code{nDip}, \code{pError}
#' and \code{mCov} input arguments. This function assumes that a single pool of
#' size \code{nDip} was used to sequence the population.
#'
#' @param nDip is an integer or a vector representing the total number of
#'   diploid individuals to simulate. Note that [scrm::scrm()] actually
#'   simulates haplotypes, so the number of simulated haplotypes is double of
#'   this. If it is a vector, then each vector entry will be simulated
#'   independently. For instance, if \code{nDip = c(100, 200)}, simulations will
#'   be carried out for samples of 100 and 200 individuals.
#' @param nloci is an integer that represents how many independent loci should
#'   be simulated.
#' @param pError an integer or a vector representing the value of the error
#'   associated with DNA pooling. This value is related with the unequal
#'   contribution of both individuals and pools towards the total number of
#'   reads observed for a given population - the higher the value the more
#'   unequal are the individual and pool contributions. If it is a vector, then
#'   each vector entry will be simulated independently.
#' @param sError a numeric value with error rate associated with the sequencing
#'   and mapping process. This error rate is assumed to be symmetric:
#'   error(reference -> alternative) = error(alternative -> reference). This
#'   number should be between 0 and 1.
#' @param mCov an integer or a vector that defines the mean depth of coverage to
#'   simulate. Please note that this represents the mean coverage across all
#'   sites. If it is a vector, then each vector entry will be simulated
#'   independently.
#' @param vCov an integer or a vector that defines the variance of the depth of
#'   coverage across all sites. If the \code{mCov} is a vector, then \code{vCov}
#'   should also be a vector, with each entry corresponding to the variance of
#'   the respective entry in the \code{mCov} vector. Thus, the first entry of
#'   the \code{vCov} vector will be the variance associated with the first entry
#'   of the \code{mCov} vector.
#' @param min.minor is an integer representing the minimum allowed number of
#'   minor-allele reads. Sites that, across all populations, have less
#'   minor-allele reads than this threshold will be removed from the data.
#' @param minimum an optional integer representing the minimum coverage allowed.
#'   Sites where the population has a depth of coverage below this threshold are
#'   removed from the data.
#' @param maximum an optional integer representing the maximum coverage allowed.
#'   Sites where the population has a depth of coverage above this threshold are
#'   removed from the data.
#'
#' @return a data.frame with columns detailing the number of diploid
#'   individuals, the pool error, the number of pools, the number of individuals
#'   per pool, the mean coverage, the variance of the coverage and the average
#'   absolute difference between the frequencies computed from genotypes and
#'   from pooled data.
#'
#' @examples
#' # a simple test with a simple combination of parameters
#' maeFreqs(nDip = 100, nloci = 10, pError = 100, sError = 0.01, mCov = 100, vCov = 200, min.minor = 1)
#'
#' # effect of two different pool error values in conjugation with a fixed coverage and pool size
#' maeFreqs(nDip = 100, nloci = 10, pError = c(100, 200), sError = 0.01,
#' mCov = 100, vCov = 200, min.minor = 1)
#'
#' # effect of two different pool error values in conjugation with a fixed pool size
#' # and two different coverages
#' maeFreqs(nDip = 100, nloci = 10, pError = c(100, 200), sError = 0.01,
#' mCov = c(100, 200), vCov = c(200, 500), min.minor = 1)
#'
#' @export
maeFreqs <- function(nDip, nloci, pError, sError, mCov, vCov, min.minor, minimum = NA, maximum = NA) {

  # create a matrix to save the values of the mean absolute error for the various conditions
  final <- matrix(data = NA, nrow = 1, ncol = 7)
  # add names to the columns of the matrix
  colnames(final) <- c("nDip", "PoolError", "nPools", "indsPool", "mean", "var", "absError")

  # create a matrix containing all combinations of factor variables
  # each row will contain one combination of number of diploids, pool error and mean depth of coverage
  combinations <- expand.grid(nDip, pError, mCov, stringsAsFactors = FALSE)

  # do a loop over all the possible combinations
  for (i in 1:nrow(combinations)) {

    # get the number of diploid individuals for this combination
    dip <- combinations[i, 1]
    # get the pooling error for this combination
    pError <- combinations[i, 2]
    # get the mean depth of coverage for this combination
    meanCov <- combinations[i, 3]
    # get the variance of the depth of coverage associated with that particular mean coverage
    varCov <- vCov[which(mCov == combinations[i, 3])]

    # compute the average absolute difference between the allele frequencies from genotypes and from Pool-seq data
    temp <- maePool(nDip = dip, nloci = nloci, pError = pError, pools = list(dip), sError = sError, mCov = meanCov,
                    vCov = varCov, min.minor = min.minor, minimum = minimum, maximum = maximum)

    # add those values to the dataframe containing all the results
    final <- rbind(final, temp)
  }

  # remove the first row of the final matrix - this row contains NAs
  final <- final[-1 ,]
  # output the final dataframe with the MAE values for the different parameter combinations
  final
}

#' Compute expected heterozygosity per site
#'
#' Computes the expected heterozygosity for a given site.
#'
#' @param geno_site is a vector where each entry contains the genotype for a
#'   single individual, coded as 0,1,2 and using NA for the missing data.
#'
#' @return a numerical value corresponding to the expected heterozygosity at
#'   that site.
#'
#' @keywords internal
#'
#' @export
ExpHet_site <- function(geno_site) {

  # get the number of individuals with data
  ngenecopies <- 2*sum(!is.na(geno_site))

  # get the frequency of the alternative allele
  freq <- sum(geno_site, na.rm = T)/ngenecopies

  # compute the expected heterozygosity
  he <- (ngenecopies/(ngenecopies-1))*2*freq*(1-freq)

  # output the expected heterozygosity
  he
}


#' Compute expected heterozygosity within a population
#'
#' This functions calculates the value of the expected heterozygosity for each
#' SNP.
#'
#' @param Pop_Pi is a matrix or list of allele frequencies. When dealing with a
#'   single locus, this input is a matrix and when dealing with multiple loci it
#'   is a list. Each entry of that list is a matrix representing a different
#'   locus. Each row of that matrix should correspond to a different population
#'   and each column to a different SNP.
#'
#' @return if the input is a single matrix, the output will be a matrix where
#'   each row represents a different population and each column is the expected
#'   heterozygosity of a population at that site. If the input is a list, the
#'   output will also be a list, with each entry corresponding to a different
#'   locus. Each of those entries will be a matrix with different populations in
#'   different rows and the expected heterozygosity of different sites at
#'   different columns.
#'
#' @keywords internal
#'
#' @export
Expected_Het <- function(Pop_Pi) {

  # dealing with a single matrix of population allelic frequencies - a single locus or simulation
  if(class(Pop_Pi)[1] == "matrix") {

    # compute the expected heterozygosity for a site - this code goes across all sites
    het <- apply(Pop_Pi, c(1,2), function (frequency) 2*frequency*(1 - frequency))

  } else { # dealing with more than one locus or simulation

    het <- lapply (Pop_Pi, FUN = function(x) {
      apply(x, c(1,2), function (frequency) 2*frequency*(1 - frequency))})
  }

  # output the expected heterozygosity
  het
}


#' Average absolute difference between expected heterozygosity
#'
#' Calculates the average absolute difference between the expected
#' heterozygosity computed directly from genotypes and from pooled sequencing
#' data.
#'
#' Different combinations of parameters can be tested to check the effect of the
#' various parameters. The average absolute difference is computed with the
#' \link[Metrics]{mae} function, assuming the expected heterozygosity computed
#' directly from the genotypes as the \code{actual} input argument and the
#' expected heterozygosity from pooled data as the \code{predicted} input
#' argument.
#'
#' @param nDip an integer representing the total number of diploid individuals
#'   to simulate. Note that [scrm::scrm()] actually simulates haplotypes, so the
#'   number of simulated haplotypes is double of this.
#' @param nloci is an integer that represents how many independent loci should
#'   be simulated.
#' @param pools a list with a vector containing the size (in number of diploid
#'   individuals) of each pool. Thus, if a population was sequenced using a
#'   single pool, the vector should contain only one entry. If a population was
#'   sequenced using two pools, each with 10 individuals, this vector should
#'   contain two entries and both will be 10.
#' @param pError an integer representing the value of the error associated with
#'   DNA pooling. This value is related with the unequal contribution of both
#'   individuals and pools towards the total number of reads observed for a
#'   given population - the higher the value the more unequal are the individual
#'   and pool contributions.
#' @param sError a numeric value with error rate associated with the sequencing
#'   and mapping process. This error rate is assumed to be symmetric:
#'   error(reference -> alternative) = error(alternative -> reference). This
#'   number should be between 0 and 1.
#' @param mCov an integer that defines the mean depth of coverage to simulate.
#'   Please note that this represents the mean coverage across all sites.
#' @param vCov an integer that defines the variance of the depth of coverage
#'   across all sites.
#' @param min.minor is an integer representing the minimum allowed number of
#'   minor-allele reads. Sites that, across all populations, have less
#'   minor-allele reads than this threshold will be removed from the data.
#' @param minimum an optional integer representing the minimum coverage allowed.
#'   Sites where the population has a depth of coverage below this threshold are
#'   removed from the data.
#' @param maximum an optional integer representing the maximum coverage allowed.
#'   Sites where the population has a depth of coverage above this threshold are
#'   removed from the data.
#'
#' @return a data.frame with columns detailing the number of diploid
#'   individuals, the pool error, the number of pools, the number of individuals
#'   per pool, the mean coverage, the variance of the coverage and the average
#'   absolute difference between the expected heterozygosity computed from
#'   genotypes and from pooled data.
#'
#' @examples
#' # single population sequenced with a single pool of 100 individuals
#' errorHet(nDip = 100, nloci = 10, pools = list(100), pError = 100, sError = 0.01,
#' mCov = 100, vCov = 250, min.minor = 2)
#'
#' # single population sequenced with two pools, each with 50 individuals
#' errorHet(nDip = 100, nloci = 10, pools = list(c(50, 50)), pError = 100, sError = 0.01,
#' mCov = 100, vCov = 250, min.minor = 2)
#'
#' # single population sequenced with two pools, each with 50 individuals
#' # removing sites with coverage below 10x or above 180x
#' errorHet(nDip = 100, nloci = 10, pools = list(c(50, 50)), pError = 100, sError = 0.01,
#' mCov = 100, vCov = 250, min.minor = 2, minimum = 10, maximum = 180)
#'
#' @export
errorHet <- function(nDip, nloci, pools, pError, sError, mCov, vCov, min.minor, minimum = NA, maximum = NA) {

  # run SCRM and obtain genotypes for a single population
  genotypes <- run_scrm(nDip = nDip, nloci = nloci)

  # simulate number of reads
  reads <- simulateCoverage(mean = mCov, variance = vCov, genotypes = genotypes)

  # if the minimum coverage is defined
  if(!is.na(minimum)) {

    # check if the maximum coverage is also defined
    if(is.na(maximum))
      stop("please define the maximum coverage")

    # remove sites with a depth of coverage above or below the defined threshold
    reads <- remove_by_reads(nLoci = nloci, reads, minimum = minimum, maximum = maximum, genotypes = genotypes)

    # get the genotypes - without sites simulated with a coverage below or above the threshold
    genotypes <- lapply(reads, function(locus) locus[[2]])

    # check the dimensions of the matrices with the genotypes
    dimensions <- matrix(unlist(lapply(genotypes, dim)), ncol = 2, byrow = TRUE)
    # we only wish to keep the locus where we have at least one polymorphic site
    tokeep <- dimensions[, 2] != 0
    # remove all loci without polymorphic sites
    genotypes <- genotypes[tokeep]

    # get the reads - without sites simulated with a coverage below or above the threshold
    reads <- lapply(reads, function(locus) locus[[1]])
    # use the same index to remove entries of the reads list that correspond to locus without polymorphic sites
    reads <- reads[tokeep]
    # ensure that each entry is a matrix
    reads <- lapply(reads, function(locus) matrix(locus, nrow = 1))
  }

  # compute the expected heterozygosity directly from genotypes
  indHets <- lapply(genotypes, function(locus) apply(X = locus, MARGIN = 2, FUN = function(site)
    ExpHet_site(geno_site = site)))

  # simulate individual contribution to the total number of reads
  indContribution <- lapply(1:nloci, function(locus)
    popsReads(list_np = pools, coverage = reads[[locus]], pError = pError))

  # simulate the number of reference reads
  reference <- lapply(1:nloci, function(locus)
    numberReferencePop(genotypes = genotypes[[locus]], indContribution = indContribution[[locus]],
                       size = pools, error = sError))

  # simulate pooled sequencing data
  pool <- poolPops(nPops = 1, nLoci = nloci, indContribution = indContribution, readsReference = reference)

  # compute the allele frequencies obtained with pooled sequencing
  pfreqs <- lapply(1:nloci, function(locus)
    Pfreqs(reference = pool[["reference"]][[locus]], alternative = pool[["alternative"]][[locus]],
           coverage = pool[["total"]][[locus]], min.minor = min.minor, ifreqs = indHets[[locus]]))

  # get the expected heterozygosity computed directly from genotypes after removing sites that did not pass the threshold
  indHets <- lapply(pfreqs, `[[`, 1)
  # get the allele frequencies computed from Pool-seq after removing sites that did not pass the threshold
  pfreqs <- lapply(pfreqs, `[[`, 2)

  # compute the mean expected heterozygosity for each population and locus - using the Pool-seq data
  poolHets <- Expected_Het(pfreqs)

  # compute the mean absolute error for each locus
  abs_error <- lapply(1:nloci, function(locus) Metrics::mae(actual = indHets[[locus]], predicted = poolHets[[locus]]))
  # replace any NaN values with the mean of the remaining values
  abs_error[is.na(abs_error)] <- mean(unlist(abs_error), na.rm = TRUE)

  # get the number of pools used to sequence a population
  nPools <- length(pools[[1]])
  # get the number of individuals per pool
  indsPool <- unique(pools[[1]])

  # create a dataframe with the values for this particular combination of parameters
  out <- data.frame(nDip=nDip, PoolError=pError, nPools=nPools, indsPool=indsPool, mean=mCov, var=vCov,
                    absError=unlist(abs_error))

  # output the results of the function
  out
}


#' Average absolute difference between the expected heterozygosity computed from
#' genotypes and from Pool-seq data
#'
#' Calculates the average absolute difference between the expected
#' heterozygosity computed directly from genotypes and from pooled sequencing
#' data.
#'
#' The average absolute difference is computed with the \link[Metrics]{mae}
#' function, assuming the expected heterozygosity computed directly from the
#' genotypes as the \code{actual} input argument and the expected heterozygosity
#' from pooled data as the \code{predicted} input argument.
#'
#' Note that this functions allows for different combinations of parameters.
#' Thus, the effect of different combinations of parameters on the average
#' absolute difference can be tested. For instance, it is possible to check what
#' is the effect of different coverages by including more than one value in the
#' \code{mCov} input argument. This function will run and compute the average
#' absolute difference for all combinations of the \code{nDip}, \code{pError}
#' and \code{mCov} input arguments. This function assumes that a single pool of
#' size \code{nDip} was used to sequence the population.
#'
#' @param nDip is an integer or a vector representing the total number of
#'   diploid individuals to simulate. Note that [scrm::scrm()] actually
#'   simulates haplotypes, so the number of simulated haplotypes is double of
#'   this. If it is a vector, then each vector entry will be simulated
#'   independently. For instance, if \code{nDip = c(100, 200)}, simulations will
#'   be carried out for samples of 100 and 200 individuals.
#' @param nloci is an integer that represents how many independent loci should
#'   be simulated.
#' @param pError an integer or a vector representing the value of the error
#'   associated with DNA pooling. This value is related with the unequal
#'   contribution of both individuals and pools towards the total number of
#'   reads observed for a given population - the higher the value the more
#'   unequal are the individual and pool contributions. If it is a vector, then
#'   each vector entry will be simulated independently.
#' @param sError a numeric value with error rate associated with the sequencing
#'   and mapping process. This error rate is assumed to be symmetric:
#'   error(reference -> alternative) = error(alternative -> reference). This
#'   number should be between 0 and 1.
#' @param mCov an integer or a vector that defines the mean depth of coverage to
#'   simulate. Please note that this represents the mean coverage across all
#'   sites. If it is a vector, then each vector entry will be simulated
#'   independently.
#' @param vCov an integer or a vector that defines the variance of the depth of
#'   coverage across all sites. If the \code{mCov} is a vector, then \code{vCov}
#'   should also be a vector, with each entry corresponding to the variance of
#'   the respective entry in the \code{mCov} vector. Thus, the first entry of
#'   the \code{vCov} vector will be the variance associated with the first entry
#'   of the \code{mCov} vector.
#' @param min.minor is an integer representing the minimum allowed number of
#'   minor-allele reads. Sites that, across all populations, have less
#'   minor-allele reads than this threshold will be removed from the data.
#' @param minimum an optional integer representing the minimum coverage allowed.
#'   Sites where the population has a depth of coverage below this threshold are
#'   removed from the data.
#' @param maximum an optional integer representing the maximum coverage allowed.
#'   Sites where the population has a depth of coverage above this threshold are
#'   removed from the data.
#'
#' @return a data.frame with columns detailing the number of diploid
#'   individuals, the pool error, the number of pools, the number of individuals
#'   per pool, the mean coverage, the variance of the coverage and the average
#'   absolute difference between the expected heterozygosity computed from
#'   genotypes and from pooled data.
#'
#' @examples
#' # a simple test with a simple combination of parameters
#' maeHet(nDip = 100, nloci = 10, pError = 100, sError = 0.01, mCov = 100, vCov = 200, min.minor = 1)
#'
#' # effect of two different pool error values in conjugation with a fixed coverage and pool size
#' maeHet(nDip = 100, nloci = 10, pError = c(100, 200), sError = 0.01,
#' mCov = 100, vCov = 200, min.minor = 1)
#'
#' # effect of two different pool error values in conjugation with a fixed pool size
#' # and two different coverages
#' maeHet(nDip = 100, nloci = 10, pError = c(100, 200), sError = 0.01,
#' mCov = c(100, 200), vCov = c(200, 500), min.minor = 1)
#'
#' @export
maeHet <- function(nDip, nloci, pError, sError, mCov, vCov, min.minor, minimum = NA, maximum = NA) {

  # create a matrix to save the values of the mean absolute error for the various conditions
  final <- matrix(data = NA, nrow = 1, ncol = 7)
  # add names to the columns of the matrix
  colnames(final) <- c("nDip", "PoolError", "nPools", "indsPool", "mean", "var", "absError")

  # create a matrix containing all combinations of factor variables
  # each row will contain one combination of number of diploids, pool error and mean depth of coverage
  combinations <- expand.grid(nDip, pError, mCov, stringsAsFactors = FALSE)

  # do a loop over all the possible combinations
  for (i in 1:nrow(combinations)) {

    # get the number of diploid individuals for this combination
    dip <- combinations[i, 1]
    # get the pooling error for this combination
    pError <- combinations[i, 2]
    # get the mean depth of coverage for this combination
    meanCov <- combinations[i, 3]
    # get the variance of the depth of coverage associated with that particular mean coverage
    varCov <- vCov[which(mCov == combinations[i, 3])]

    # compute the average absolute difference between the allele frequencies from genotypes and from Pool-seq data
    temp <- errorHet(nDip = dip, nloci = nloci, pError = pError, pools = list(dip), sError = sError, mCov = meanCov,
                     vCov = varCov, min.minor = min.minor, minimum = minimum, maximum = maximum)

    # add those values to the dataframe containing all the results
    final <- rbind(final, temp)
  }

  # remove the first row of the final matrix - this row contains NAs
  final <- final[-1 ,]

  # output the final dataframe with the MAE values for the different parameter combinations
  final
}


#' Create invariable sites
#'
#' This function applies a correction for the situations where [scrm::scrm()]
#' does not produce a single polymorphic site for a given locus. In this
#' situation, two artificial sites are created at that locus. All individuals
#' are assumed to be homozygous for the reference allele at those sites.
#'
#' @param haplotypes a list of haplotypes obtained from the simulations done
#'   with [scrm::scrm()]. Each entry of the list is a matrix that corresponds to
#'   a given locus. At each matrix, each column is a different site and each row
#'   is a different haplotype.
#' @param nHap an integer representing the total number of haplotypes simulated.
#'
#' @return a list of haplotypes identical to `haplotypes`, but without empty
#'   loci.
#'
#' @keywords internal
#'
#' @export
haplo.fix <- function(haplotypes, nHap) {

  # get the dimensions of each list entry - organized into a matrix with number of rows in the first column
  # and number of columns in the second
  size <- matrix(unlist(lapply(haplotypes, dim)), ncol = 2, byrow = TRUE)

  # get the index of list entries with zero columns
  index <- size[, 2] == 0

  # replace those entries with a matrix containing a single column of zeros
  haplotypes[index] <- list(matrix(rep(0, times = nHap), ncol = 1))

  # get the dimensions of each list entry - organized in the same manner as before
  size <- matrix(unlist(lapply(haplotypes, dim)), ncol = 2, byrow = TRUE)

  # get the index of list entries with a single column
  index <- size[, 2] == 1

  # if there are any columns with a single entry
  if(sum(index) != 0) {

    # add another column with only zeros to those entries
    haplotypes[index] <- sapply(haplotypes[index], function(x)
      list(cbind(as.matrix(unlist(x), byrow = TRUE), matrix(rep(0, times = nHap)))))
  }

  # output the haplotypes
  haplotypes
}


#' Convert haplotypes to genotypes
#'
#' This function converts haplotypes simulated with [scrm::scrm()] into
#' genotypes by adding the entries on one row with the entries of the subsequent
#' row.
#'
#' @param haplo a matrix of haplotypes obtained from the simulations done with
#'   [scrm::scrm()]. Each column of the matrix is a different site and each row
#'   is a different haplotype.
#'
#' @return a matrix of genotypes with half the rows of the `haplo` matrix. Each
#'   column of this matrix is a different site and each row is a different
#'   genotype.
#'
#' @keywords internal
#'
#' @export
hap2geno <- function(haplo) {

  # add each row of the matrix to the row that is immediately below it
  genotypes <- haplo[seq(1, by = 2, to = nrow(haplo)),] + haplo[seq(2, by = 2, to = nrow(haplo)), , drop = FALSE]
  # output the matrix of genotypes
  genotypes
}


#' Create genotypes from a output with haplotypes
#'
#' This function applies the \code{\link{hap2geno}} function to all entries of a
#' list. Each entry of that list is a different locus simulated with
#' [scrm::scrm()]. Thus, this function converts the haplotypes of all simulated
#' loci into genotypes.
#'
#' @param haplotypes a list of haplotypes obtained from the simulations done
#'   with [scrm::scrm()]. Each entry of the list is a matrix that corresponds to
#'   a given locus. At each matrix, each column is a different site and each row
#'   is a different haplotype.
#' @param nDip an integer representing the total number of diploid individuals
#'   to simulate. Note that this is the total number of diploid individuals and
#'   not the number of individuals per population.
#'
#' @return a list of genotypes. Each entry of the list is a matrix corresponding
#'   to a different locus. locus. At each matrix, each column is a different
#'   site and each row is a different genotype
#'
#' @keywords internal
#'
#' @export
GetGenotypes <- function(haplotypes, nDip) {

  # apply the hap2geno function across all entries of the list - sum haplotypes to get genotypes
  genotypes <- lapply(haplotypes, function(x) hap2geno(x))
  # remove the name (position) of each site - this is something that scrm creates
  genotypes <- lapply(genotypes, function(x) unname(x))
  # output the genotypes
  genotypes
}


#' Simulate coverage at a single locus
#'
#' Simulates the total number of reads, for each polymorphic site of a given
#' locus using a negative binomial distribution.
#'
#' The total number of reads is simulated with a negative binomial and according
#' to a user-defined mean depth of coverage and variance. This function is
#' intended to work with a matrix of genotypes, simulating the depth of coverage
#' for each site present in the genotypes. However, it can also be used to
#' simulate coverage distributions independent of genotypes, by choosing how
#' many sites should be simulated (with the `nSNPs` option).
#'
#' @param mean an integer that defines the mean depth of coverage to simulate.
#'   Please note that this represents the mean coverage across all sites. If a
#'   vector is supplied instead, the function assumes that each entry of the
#'   vector is the mean for a different population.
#' @param variance an integer that defines the variance of the depth of coverage
#'   across all sites. If a vector is supplied instead, the function assumes
#'   that each entry of the vector is the variance for a different population.
#' @param nSNPs an integer representing the number of polymorphic sites per
#'   locus to simulate. This is an optional input but either this or the
#'   `genotypes` matrix must be supplied.
#' @param genotypes a matrix of simulated genotypes, where each column is a
#'   different SNP and each row is a different individual. This is an optional
#'   input but either this or the `nSNPs` must be supplied.
#'
#' @return a matrix with the total coverage per population and per site.
#'   Different rows represent different populations and each column is a
#'   different site.
#'
#' @examples
#' # coverage for one population at 10 sites
#' simReads(mean = 20, variance = 100, nSNPs = 10)
#'
#' # simulate coverage at one locus with 10 SNPs for two populations:
#' # the first with 100x and the second with 50x
#' simReads(mean = c(100, 50), variance = c(250, 150), nSNPs = 10)
#'
#' @export
simReads <- function(mean, variance, nSNPs = NA, genotypes = NA) {

  # check if either the number of SNPs or the genotypes were supplied as input
  if(all(is.na(nSNPs), is.na(genotypes)))
    stop("You should define the number of SNPs to simulate or supply a matrix of genotypes. Please check")

  # check if the variance and mean are reasonable
  if(any(variance - mean > 0) == FALSE)
    stop("Error: variance equal to mean, or variance smaller than mean")

  # if genotypes are supplied as input argument, assume that the number of SNPs
  # is equal to the number of columns of the genotypes matrix
  if(!all(is.na(genotypes)))
    nSNPs <- ncol(genotypes)

  # calculate the parameters for the negative binomial
  pnb <- mean/variance
  rnb <- (mean^2)/(variance - mean)

  # use a negative binomial to draw random values, per site and per population, for the total number of observed reads
  readnumbers <- t(mapply(FUN = function(size, prob) stats::rnbinom(n = nSNPs, size = size, prob = prob), rnb, pnb))

  # if there is only a single SNP - we need to transpose the previous result to get each population in a different row
  if(nSNPs == 1)
    readnumbers <- t(readnumbers)

  # get the output - number of reads per site and per population
  readnumbers
}


#' Simulate total number of reads per site
#'
#' This function simulates the total number of reads, for each polymorphic site
#' using a negative binomial distribution.
#'
#' The total number of reads is simulated with a negative binomial and according
#' to a user-defined mean depth of coverage and variance. This function is
#' intended to work with a list of genotypes, simulating the depth of coverage
#' for each site present in the genotypes. However, it can also be used to
#' simulate coverage distributions independent of genotypes, by choosing how
#' many loci to simulate (with the `nLoci` option) and choosing how many sites
#' per locus should be simulated (with the `nSNPs` option).
#'
#' @param mean an integer that defines the mean depth of coverage to simulate.
#'   Please note that this represents the mean coverage across all sites. If a
#'   vector is supplied instead, the function assumes that each entry of the
#'   vector is the mean for a different population.
#' @param variance an integer that defines the variance of the depth of coverage
#'   across all sites. If a vector is supplied instead, the function assumes
#'   that each entry of the vector is the variance for a different population.
#' @param genotypes a list of simulated genotypes, where each entry is a matrix
#'   corresponding to a different locus. At each matrix, each column is a
#'   different SNP and each row is a different individual. This is an optional
#'   input but either this or the `nSNPs` must be supplied.
#' @param nSNPs an integer representing the number of polymorphic sites per
#'   locus to simulate. This is an optional input but either this or the
#'   `genotypes` list must be supplied.
#' @param nLoci an optional integer that represents how many independent loci
#'   should be simulated.
#'
#' @return a list with the total coverage per population and per site. Each list
#'   entry is a matrix corresponding to a different locus. For each matrix,
#'   different rows represent different populations and each column is a
#'   different site.
#'
#' @examples
#' # simulate 10 loci, each with 10 SNPs for a single population
#' simulateCoverage(mean = 100, variance = 250, nSNPs = 10, nLoci = 10)
#'
#' # simulate 10 loci, each with 10 SNPs for two populations:
#' # the first with 100x and the second with 50x
#' simulateCoverage(mean = c(100, 50), variance = c(250, 150), nSNPs = 10, nLoci = 10)
#'
#' # simulate coverage given a set of genotypes
#' # run scrm and obtain genotypes
#' genotypes <- run_scrm(nDip = 100, nloci = 10)
#' # simulate coverage
#' simulateCoverage(mean = 50, variance = 200, genotypes = genotypes)
#'
#' @export
simulateCoverage <- function(mean, variance, nSNPs = NA, nLoci = NA, genotypes = NA) {

  # check if either the number of SNPs or the genotypes were supplied as input
  if(all(is.na(nSNPs), is.na(genotypes)))
    stop("You should define the number of SNPs to simulate or supply a list of genotypes. Please check")

  # if the genotypes are supplied as input to the function
  if(!all(is.na(genotypes))) {

    # check if the input is correct - genotypes should always be supplied as a list
    if(any(class(genotypes) != "list"))
      stop(paste("genotypes should be supplied on a list format, with each entry corresponding to a locus. Please check"))

    # use a negative binomial to draw random values, per site and per population, for the total number of observed reads
    # this outputs a list where each entry corresponds to a locus
    readnumbers <- lapply(genotypes, FUN = function(geno) simReads(mean, variance, genotypes = geno))

  } else {

    # check if the input is correct - when genotypes are not supplied as input, the number of loci should be defined
    if(is.na(nLoci))
      stop(paste("Please define the number of loci to simulate"))

    # use a negative binomial to draw random values, per site and per population, for the total number of observed reads
    # this outputs a list where each entry corresponds to a locus
    readnumbers <- lapply(1:nLoci, FUN = function(geno) simReads(mean, variance, nSNPs))
  }

  # get the output - number of reads per site and per population
  readnumbers
}


#' Apply a coverage-based filter to a matrix
#'
#' This function removes sites that have a coverage below a `minimum` value and
#' sites with a coverage above a `maximum` value. If a matrix of genotypes is
#' also supplied, then those same sites are also removed from that matrix.
#'
#' @param reads a matrix with the total depth of coverage. Each row of the
#'   matrix should be the coverage of a different population and each column a
#'   different site.
#' @param minimum an integer representing the minimum coverage allowed. Sites
#'   where any population has a depth of coverage below this threshold are
#'   removed from the data.
#' @param maximum an integer representing the maximum coverage allowed. Sites
#'   where any population has a depth of coverage above this threshold are
#'   removed from the data.
#' @param genotypes an optional matrix input with the genotypes. Each column of
#'   the matrix should be a different site and each row a different individual.
#'
#' @return a matrix with the total depth of coverage minus the sites (i.e.
#'   columns) where the coverage for any of the populations was below the
#'   minimum or above the maximum. If genotypes were supplied, then the output
#'   will be a list, with one entry per locus. Each entry will contain the
#'   filtered coverage in the first entry and the genotypes, minus the removed
#'   sites, in the second entry.
#'
#' @examples
#' set.seed(10)
#'
#' # simulate coverage for a single locus - select the first entry to obtain a matrix
#' reads <- simulateCoverage(mean = c(25, 25), variance = c(200, 200), nSNPs = 10, nLoci = 1)[[1]]
#'
#' # check the coverage matrix
#' reads
#'
#' # remove sites with coverage below 10x or above 100x
#' remove_by_reads_matrix(reads = reads, minimum = 10, maximum = 100)
#'
#' @export
remove_by_reads_matrix <- function(reads, minimum, maximum, genotypes = NA) {

  # check which sites, if any, have a coverage below or above the threshold
  toremove <- apply(X = reads, MARGIN = 2, function(col) any(col < minimum | col > maximum))

  # if there are any sites with a coverage below or above the required value - remove those sites from the matrix
  if(length(toremove) != 0)
    reads <- reads[, !toremove, drop = FALSE]

  # set the output if the only input were the number of reads
  output <- reads

  # when genotypes were also supplied as input
  if(!all(is.na(genotypes))) {

    # remove the same sites from the matrix containing the genotypes
    genotypes <- genotypes[, !toremove, drop = FALSE]
    # create the new output, containing both the reads and the genotypes
    output <- list(reads, genotypes)
  }

  # output the number of reads
  output
}


#' Apply a coverage-based filter over a list
#'
#' This function removes sites that have a coverage below a `minimum` value and
#' sites with a coverage above a `maximum` value. This is done over multiple
#' loci, assuming that each entry of the `reads` list is a different locus. If a
#' list of genotypes is also supplied, then those same sites are also removed
#' from each locus of the genotypes.
#'
#' @param nLoci an integer that represents how many independent loci were
#'   simulated.
#' @param reads a list with the total depth of coverage. Each entry of the list
#'   should be a matrix corresponding to a different locus. Each row of that
#'   matrix should be the coverage of a different population and each column a
#'   different site.
#' @param minimum an integer representing the minimum coverage allowed. Sites
#'   where any population has a depth of coverage below this threshold are
#'   removed from the data.
#' @param maximum an integer representing the maximum coverage allowed. Sites
#'   where any population has a depth of coverage above this threshold are
#'   removed from the data.
#' @param genotypes an optional list input with the genotypes. Each entry of the
#'   list should be a matrix corresponding to a different locus. Each column of
#'   the matrix should be a different site and each row a different individual.
#'
#' @return a list with the total depth of coverage similar to the `reads` input
#'   argument but without sites where the coverage was below the `minimum` or
#'   above the `maximum`. If the genotypes were included, a second list entry
#'   will also be included in the output, containing the genotypes minus the
#'   sites that were removed.
#'
#' @examples
#' set.seed(10)
#'
#' # simulate coverage for 10 locus
#' reads <- simulateCoverage(mean = c(25, 25), variance = c(200, 200), nSNPs = 10, nLoci = 10)
#'
#' # remove sites with coverage below 10x or above 100x
#' reads <- remove_by_reads(nLoci = 10, reads = reads, minimum = 5, maximum = 100)
#' # notice that some locus no longer have 10 SNPs - those sites were removed
#' reads
#'
#' @export
remove_by_reads <- function(nLoci, reads, minimum, maximum, genotypes = NA) {

  # check if the input is correct - reads should always be supplied as a list
  if(!inherits(reads, "list"))
    stop(paste("reads should be supplied on a list format, with each entry corresponding to a locus. Please check"))

  # this applies the remove_by_reads_matrix function to all the list entries - i.e. to all the different loci
  # note that this will also remove those sites from the genotypes list - if genotypes are supplied as input
  # if genotypes were not supplied as input
  if(any(is.na(genotypes))) {
    # remove sites with a depth of coverage above or below the defined threshold
    out <- lapply(1:nLoci, function(locus) remove_by_reads_matrix(reads = reads[[locus]], minimum, maximum))

  } else { # if genotypes are supplied as input

    # remove sites with a depth of coverage above or below the defined threshold
    out <- lapply(1:nLoci, function(locus)
      remove_by_reads_matrix(reads = reads[[locus]], minimum, maximum, genotypes = genotypes[[locus]]))
  }

  # output the number of reads - and genotypes if relevant
  out
}


#' Correct alpha value
#'
#' This function corrects the alpha value of the Dirichlet distribution used to
#' simulate the probability of contribution of different pools and individuals.
#'
#' The alpha value corresponds to the vector of shape parameters of the
#' Dirichlet distribution. When the alpha is very small, the random generation
#' of numbers from the Dirichlet distribution produces `NaN`. Thus, this
#' function replaces small values of alpha with a minimum threshold value to
#' avoid that.
#'
#' @param alpha_i is a vector of shape parameters for the Dirichlet
#'   distribution.
#'
#' @return a vector of corrected shape parameters for the Dirichlet
#'   distribution. Very small values are replaced by a minimum value.
#'
#' @keywords internal
#'
#' @export
set_alpha <- function(alpha_i) {

  # define threshold for the minimum alpha value
  min.alpha <- 1e-2

  # if any value of alpha_i is below the minimum allowed value, set alpha to min.alpha value
  alpha_i[alpha_i < min.alpha] <- min.alpha

  # output the alpha_i for the Dirichlet distribution
  alpha_i
}


#' Probability of contribution of each pool
#'
#' This function computes the probability of contribution of each pool towards
#' the total depth of coverage of a single population. If multiple pools where
#' used to sequence a single population, it is possible that some pools
#' contribute more than others.
#'
#' @param nPools an integer indicating how many pools were used to sequence the
#'   population.
#' @param vector_np is a vector where each entry contains the number of diploid
#'   individuals of a given pool. Thus, if a population was sequenced using two
#'   pools, each with 10 individuals, this vector would contain two entries and
#'   both will be 10.
#' @param nSNPs an integer indicating how many SNPs exist in the data.
#' @param pError an integer representing the value of the error associated with
#'   DNA pooling. This value is related with the unequal pool contribution
#'   towards the total number of reads of a population - the higher the value
#'   the more unequal are the pool contributions.
#'
#' @return a matrix with the probabilities of contribution for each pool. Each
#'   row represents a different pool and each column is a different site.
#'
#' @examples
#' # probability of contribution at 8 SNPs for 5 pools, each with 10 individuals
#' poolProbs(nPools = 5, vector_np = rep(10, 5), nSNPs = 8, pError = 50)
#'
#' @export
poolProbs <- function(nPools, vector_np, nSNPs, pError) {

  # Now, we need to calculate the probability of contributing for each individual - in each population
  # This package contains a function to perform random draws from a Dirichlet distribution

  # check if we are dealing with a single population - this function should be used on a single population
  # also check if the input is on the correct format
  if(!inherits(vector_np, "numeric") | length(vector_np) == 1)
    stop(paste("The vector_np input should be a vector. It should also contain more than one entry. Please check"))

  # check if we are dealing with a single population - this function should be used on a single population
  # also check if the input is on the correct format
  if(nPools != length(vector_np))
    stop(paste("The nPools input is", paste(nPools), "and so, the length of the size vector should also be",
               paste(nPools, ".", sep = ""), "Please check"))

  # the total number of individuals in the population (n) can be obtained by adding the individuals in all pools
  n <- sum(vector_np)

  # pooling error is defined in % - change this to a proportion
  pError <- pError/100

  # if we use Dir(rho*np/n), then the alpha_i for Dirichlet can be written as
  numerator <- (n - 1 - (pError^2))*vector_np
  denominator <- n*(pError^2)
  alpha <- numerator/denominator

  # check, and correct if needed, if any alpha_i value is above or below the threshold
  alpha <- set_alpha(alpha_i = alpha)

  # use a Dirichlet distribution to get the probability of contribution for each pool across all sites
  probs <- t(MCMCpack::rdirichlet(n = nSNPs, alpha = alpha))

  # output the probability of contributing for each pool
  probs
}


#' Reads contributed by each pool
#'
#' This function simulates the contribution, in terms of reads, of each pool.
#' The number of reads contributed from all pools is equal to the total coverage
#' of the population.
#'
#' @param nPools an integer indicating how many pools were used to sequence the
#'   population.
#' @param coverage a vector containing the total depth of coverage of the
#'   population. Each entry of the vector represents a different site.
#' @param probs a matrix containing the probability of contribution of each pool
#'   used to sequence the population. This matrix can be obtained with the
#'   `poolProbs` function.
#'
#' @return a matrix with the number of reads contributed by each pool towards
#'   the total coverage of the population. Each row of the matrix is a different
#'   pool and each column a different site.
#'
#' @examples
#' # simulate the probability of contribution of each pool
#' probs <- poolProbs(nPools = 5, vector_np = rep(10, 5), nSNPs = 8, pError = 50)
#'
#' # simulate the number of reads contributed, assuming 10x coverage for each site
#' poolReads(nPools = 5, coverage = rep(10, 8), probs = probs)
#'
#' @export
poolReads <- function(nPools, coverage, probs) {

  # set the output of this function if a given locus has no sites
  if(length(coverage) == 0) {

    # set the contribution to NA
    contribution <- NA

  } else {

    # check if the dimensions make sense
    if(length(coverage) != ncol(probs))
      stop("different number of sites in the coverage and probs input arguments")

    # simulate the contribution of each pool towards the total depth of coverage
    contribution <- vapply(1:length(coverage), FUN = function(i) {
      stats::rmultinom(1, size = coverage[i], prob = probs[, i])
    }, FUN.VALUE = numeric(nPools))
  }

  # output a matrix containing the contribution (in number of reads) for each pool and across all sites
  contribution
}


#' Probability of contribution of each individual
#'
#' This function computes the probability of contribution for each individual of
#' a given pool. Please note that this function works for a single pool and
#' should not be directly applied to situations where multiple pools were used.
#'
#' @param np an integer specifying how many individuals were pooled.
#' @param nSNPs an integer indicating how many SNPs exist in the data.
#' @param pError an integer representing the value of the error associated with
#'   DNA pooling. This value is related with the unequal individual contribution
#'   towards the total number of reads contributed by a single pool - the higher
#'   the value the more unequal are the individual contributions.
#'
#' @return a matrix with the probabilities of contribution for each individual.
#'   Each row represents a different individual and each column is a different
#'   site.
#'
#' @examples
#' # probability of contribution for 10 individuals at 5 sites
#' indProbs(np = 10, nSNPs = 5, pError = 100)
#'
#' @export
indProbs <- function(np, nSNPs, pError) {

  # pooling error is defined in % - change this to a proportion
  pError <- pError/100

  # if we use Dir(rho/np), then the alpha_i for Dirichlet can be written as
  numerator <- (np - 1 - (pError^2))
  denominator <- np*(pError^2)
  alpha <- numerator/denominator

  # check, and correct if needed, if any alpha_i value is above or below the threshold
  alpha <- set_alpha(alpha_i = alpha)

  # use a dirichlet distribution to get the probability of contribution for each individual across all sites
  probs <- t(MCMCpack::rdirichlet(n = nSNPs, alpha = rep(alpha, times = np)))

  # output the probability of contributing for each individual
  probs
}


#' Reads contributed by each individual
#'
#' This function simulates the contribution, in terms of reads, of each
#' individual of a given pool. Please note that this function works for a single
#' pool and should not be directly applied to situations where multiple pools
#' were used.
#'
#' @param np an integer specifying how many individuals were pooled.
#' @param coverage a vector containing the total depth of coverage of a given
#'   pool. Each entry of the vector represents a different site.
#' @param probs a matrix containing the probability of contribution of each
#'   individual. This matrix can be obtained with the `indProbs` function.
#'
#' @return a matrix with the number of reads contributed by each individual
#'   towards the coverage of its pool. Each row of the matrix is a different
#'   individual and each column a different site.
#'
#' @examples
#' # probability of contribution for 10 individuals at 5 sites
#' probs <- indProbs(np = 10, nSNPs = 5, pError = 100)
#'
#' # simulate the number of reads contributed, assuming 10x coverage for each site
#' indReads(np = 10, coverage = rep(10, 5), probs = probs)
#'
#' @export
indReads <- function(np, coverage, probs) {

  # set the output of this function if a given locus has no sites
  if(length(coverage) == 0) {

    # set the contribution to NA
    contribution <- NA

  } else {

    # check if the dimensions make sense
    if(length(coverage) != ncol(probs))
      stop("different number of sites in the coverage and probs input arguments")

    # simulate the contribution of each individual towards the total depth of coverage
    contribution <- vapply(1:length(coverage), FUN = function(i) {
      stats::rmultinom(1, size = coverage[i], prob = probs[, i])
    }, FUN.VALUE = numeric(np))
  }

  # output a matrix containing the contribution (in number of reads) for each individual and across all sites
  contribution
}


#' Compute number of reads for each individual and across all sites
#'
#' This function computes the contribution of each individual towards the total
#' coverage of a given population.
#'
#' If multiple pools were used to sequence a population, this will compute the
#' contribution of each pool and then use that to calculate how many reads does
#' that pool contribute. Next, the probability of contribution of each
#' individual is computed and utilized to calculate the number of reads that
#' each individual contributes towards the total number of reads observed in the
#' corresponding pool.
#'
#' @param vector_np is a vector where each entry contains the number of diploid
#'   individuals of a given pool. Thus, if a population was sequenced using two
#'   pools, each with 10 individuals, this vector would contain two entries and
#'   both will be 10.
#' @param coverage a vector containing the total depth of coverage of the
#'   population. Each entry of the vector represents a different site.
#' @param pError an integer representing the value of the error associated with
#'   DNA pooling. This value is related with the unequal contribution of both
#'   individuals and pools towards the total number of reads observed for a
#'   given population - the higher the value the more unequal are the individual
#'   and pool contributions.
#'
#' @return a matrix with the number of reads contributed by each individual.
#'   Each row of the matrix corresponds to a different individual and each
#'   column to a different site.
#'
#' @examples
#' # simulate number of reads contributed by each individual towards the total population coverage
#' # assuming a coverage of 10x at 5 sites and two pools, each with 5 individuals
#' popReads(vector_np = c(5, 5), coverage = rep(10, 5), pError = 100)
#'
#' @export
popReads <- function(vector_np, coverage, pError) {

  # when the coverage input is a list, convert it to a vector
  if(inherits(coverage, "list"))
    coverage <- unlist(coverage)

  # get the number of diploid individuals
  nDip <- sum(vector_np)

  # get the number of pools used to sequence the population
  nPools <- length(vector_np)

  # get the number of polymorphic sites
  nSNPs <- length(coverage)

  # if no polymorphic sites exist at any given locus - set the output to NA
  if(nSNPs == 0) {
    # set the output
    iReads <- matrix(data = NA, nrow = 1)

  } else {

    # when the population was sequenced using a single pool of individuals
    if (nPools == 1) {

      # calculate the probability of contribution for each individual
      probs <- indProbs(np = nDip, nSNPs = nSNPs, pError = pError)

      # compute the contribution of each individual and across all sites
      # towards the total depth of coverage of the population
      iReads <- indReads(np = nDip, coverage = coverage, probs = probs)

    } else { # When more than one pool of individuals was used to sequence the population

      # calculate the probability of contribution for each pool used to sequence the population
      probs_pool <- poolProbs(nPools = nPools, vector_np = vector_np, nSNPs = nSNPs, pError = pError)

      # compute the number of reads that each pool contributes per site
      # This takes into account the total depth of coverage of the population for a given site
      # and divides that coverage amongst all the pools - considering also the probability of contribution of each pool
      pool_reads <- poolReads(nPools = nPools, coverage, probs = probs_pool)

      # calculate individual probabilities for each pool
      probs_ind <- lapply(vector_np, FUN = function(pool) indProbs(np = pool, nSNPs = nSNPs, pError = pError))

      # Taking into account the coverage of each pool, simulate how many reads does each individual inside a pool
      # contributes towards that number i.e if a pool has 25 reads at a given site:
      # simulate how many of these reads come from the first individual, how many from the second, etc
      iReads <- lapply(1:nPools, function(i)
        indReads(np = vector_np[i], coverage = pool_reads[i, ], probs = probs_ind[[i]]))

      # combine the entries from the various pools to obtain how many reads does a given individual
      # contributes to the total depth of coverage of the population
      iReads <- do.call(rbind, iReads)
    }
  }

  # output a matrix containing the contribution (in number of reads) for each individual and across all sites
  iReads
}


#' Simulate total number of reads for multiple populations
#'
#' Simulates the contribution of each individual towards the total coverage of
#' its population.
#'
#' If multiple pools were used to sequence a population, this will compute the
#' contribution of each pool and then use that to calculate how many reads does
#' that pool contribute. Next, the probability of contribution of each
#' individual is computed and utilized to calculate the number of reads that
#' each individual contributes towards the total number of reads observed in the
#' corresponding pool. These steps will be performed for each population, thus
#' obtaining the number of reads contributed by each individual for each
#' population.
#'
#' @param list_np is a list where each entry corresponds to a different
#'   population. Each entry is a vector and each vector entry contains the
#'   number of diploid individuals of a given pool. Thus, if a population was
#'   sequenced using two pools, each with 10 individuals, this vector would
#'   contain two entries and both will be 10.
#' @param coverage a matrix containing the total depth of coverage of all
#'   populations. Each row corresponds to a different population and each column
#'   to a different site.
#' @param pError an integer representing the value of the error associated with
#'   DNA pooling. This value is related with the unequal contribution of both
#'   individuals and pools towards the total number of reads observed for a
#'   given population - the higher the value the more unequal are the individual
#'   and pool contributions.
#'
#' @return a list with one entry per population. Each entry represents the
#'   number of reads contributed by each individual towards the total coverage
#'   of its population. Different individuals correspond to different rows and
#'   different sites to different columns.
#'
#' @examples
#' # simulate coverage for two populations sequenced at 10x at 5 sites
#' reads <- simulateCoverage(mean = c(10, 10), variance = c(20, 20), nSNPs = 5, nLoci = 1)
#'
#' # simulate the individual contribution towards that coverage
#' # assuming that the first population was sequenced using two pools of 5 individuals
#' # and the second using a single pool with 10 individuals
#' popsReads(list_np = list(c(5, 5), 10), coverage = reads, pError = 5)
#'
#' @export
popsReads <- function(list_np, coverage, pError) {

  # get the number of populations
  nPops <- length(list_np)

  # when the coverage input is a list, convert with to a matrix, with each row corresponding to a population
  if(any(class(coverage) == "list"))
    coverage <- matrix(unlist(coverage), nrow = nPops, byrow = FALSE)

  # Taking into account the depth of coverage for each population, simulate how many reads does each individual contributes
  # towards the total coverage of the population
  indCoverage <- lapply(1:nPops, function(pop) popReads(vector_np = list_np[[pop]], coverage[pop, ], pError))

  # output the contribution (in number of reads) for each individual and across all sites
  indCoverage
}


#' Split matrix of genotypes
#'
#' This function splits a matrix into different list entries. The matrix is
#' split according to a set of row indexes defined by the `size` input.
#'
#' The `size` input is utilized to create the index of the rows that go into the
#' different list entries. It specifies the size, in terms of number of
#' individuals, of each population.
#'
#' @param matrix is obviously a matrix, ideally one containing genotypes. Each
#'   column of the matrix should be a different site and each row a different
#'   individual.
#' @param size a list with one entry per population. Each entry should be a
#'   vector containing the size (in number of diploid individuals) of each pool.
#'   Thus, if a population was sequenced using a single pool, the vector should
#'   contain only one entry. If a population was sequenced using two pools, each
#'   with 10 individuals, this vector should contain two entries and both will
#'   be 10.
#'
#' @return a list with one entry per entry of the size input argument. Each
#'   entry contains the information of the original matrix for the number of
#'   individuals specified by the corresponding entry of the size input
#'   argument.
#'
#' @keywords internal
#'
#' @examples
#' set.seed(10)
#'
#' # create a random matrix
#' mymatrix <- matrix(round(runif(n = 50, min = 10, max = 25)), nrow = 10)
#'
#' # split that matrix assuming 8 individuals in the first population and two in the second
#' splitMatrix(matrix = mymatrix, size = list(8, 2))
#'
#' @export
splitMatrix <- function(matrix, size) {

  # general check to see if the input is correctly supplied
  if(!inherits(size, "list"))
    stop(paste("The size input should be a list. Please check"))

  # get the number of populations - it's the number of entries in the size input
  nPops <- length(size)

  # Perform a cumulative sum - this will create a vector, starting at the number one
  # each subsequent entry on the vector is the index of the last individual of a population
  popsize <- c(0, cumsum(lapply(size, sum)))

  # Use the previous index to create vectors containing the index of all individuals per population.
  # For example, if you have a population starting at index 30 and ending at 50, this will create a vector
  # containing the numbers 30, 31, 32... etc until 50. This creates a list, where each entry is a vector with
  # the indices of a single population
  index <- lapply(1:nPops, function(i) seq(popsize[i]+1, popsize[i+1]))

  # use the vectors containing the index of all individuals of a given population to subset the matrix
  # The sub-setting is done by rows and it creates a list, where each entry contains the information for one population
  output <- lapply(1:nPops, function(pop) matrix[index[[pop]], , drop = FALSE])

  # Output the list containing, in each entry, the information for each population
  output
}


#' Compute the number of reference reads
#'
#' This function takes as input the total depth of coverage and computes how
#' many of those reads are reference allele reads.
#'
#' More precisely, this function computes the number of reference reads per site
#' for one individual, given the genotype of the individual at each site, the
#' total number of reads observed for the individual at that site and an error
#' rate.
#'
#' @param genotype_v is a vector with the genotype of a given individual. Each
#'   entry of the vector should be a different site. Genotypes should be encoded
#'   as 0: reference homozygote, 1: heterozygote and 2: alternative homozygote.
#' @param readCount_v is a vector with the number of reads contributed by the
#'   same given individual. Each entry of that vector should be a different
#'   site.
#' @param error a numeric value with error rate associated with the sequencing
#'   and mapping process. This error rate is assumed to be symmetric:
#'   error(reference -> alternative) = error(alternative -> reference). This
#'   number should be between 0 and 1.
#'
#' @return a vector with the number of reference allele reads. Each entry of the
#'   vector corresponds to a different individual.
#'
#' @examples
#' # number of reference allele reads for three individuals, each with 10x coverage
#' # one individual is homozygote for the reference allele (0), other is heterozygote (1)
#' # and the last is homozygote for the alternative allele (2)
#' getNumReadsR_vector(genotype_v = c(0,1,2), readCount_v = c(10, 10, 10), error = 0.01)
#'
#' @export
getNumReadsR_vector <- function(genotype_v, readCount_v, error) {

  # ensure that both vectors have the same length
  # you should have one genotype and one value of depth of coverage per site
  if(length(genotype_v) != length(readCount_v))
    stop(paste("The lengths of the genotype_v and the readCount_v inputs should be the same. Please check"))

  # initialize the number of A reads as a vector of 0 with same size as Genotype_v
  numAreads <- numeric(length(genotype_v))

  # Get the entries of genotype for hom reference
  hom_anc <- genotype_v == 0
  # Get the entries of heterozygotes (1)
  het <- genotype_v == 1
  # Get the entries for hom alternative (2)
  hom_der <- genotype_v == 2

  # Call the rbinom in a vector way, using size as the vector
  # homozygote reference
  numAreads[hom_anc] <- stats::rbinom(n = sum(hom_anc), size = readCount_v[hom_anc], prob = 1-error)
  # hets
  numAreads[het] <- stats::rbinom(n = sum(het), size = readCount_v[het], prob = 0.5)
  # homozygote alternative
  numAreads[hom_der] <- stats::rbinom(n = sum(hom_der), size = readCount_v[hom_der], prob = error)
  # return the number of reads
  numAreads
}


#' Compute the number of reference reads over a matrix
#'
#' This function works over all the rows and columns of a matrix and computes
#' the number of reads containing the reference allele at each site and for each
#' individual.
#'
#' @param genotypes is a matrix of genotypes. Each column of the matrix should
#'   be a different site and each row a different individual. Genotypes should
#'   be encoded as 0: reference homozygote, 1: heterozygote and 2: alternative
#'   homozygote.
#' @param indContribution is a matrix of individual contributions. Each row of
#'   that matrix is a different individual and each column is a different site.
#'   Thus, each entry of the matrix should contain the number of reads
#'   contributed by that individual at that particular site.
#' @param error a numeric value with error rate associated with the sequencing
#'   and mapping process. This error rate is assumed to be symmetric:
#'   error(reference -> alternative) = error(alternative -> reference). This
#'   number should be between 0 and 1.
#'
#' @return a matrix with the number of reference allele reads contributed by
#'   each individual. Each row of the matrix represents a different individual
#'   and each column is a different site.
#'
#' @examples
#' # probability of contribution for 10 individuals at 5 sites
#' probs <- indProbs(np = 10, nSNPs = 5, pError = 5)
#'
#' # simulate the number of reads contributed, assuming 20 coverage for each site
#' indContribution <- indReads(np = 10, coverage = rep(20, 5), probs = probs)
#'
#' # set seed and create a random matrix of genotypes
#' set.seed(10)
#' genotypes <- matrix(rpois(50, 0.5), nrow = 10)
#'
#' # simulate the number of reads with the reference allele
#' computeReference(genotypes = genotypes, indContribution = indContribution, error = 0.01)
#'
#' @export
computeReference <- function(genotypes, indContribution, error) {

  # ensure that both matrices have the same dimension
  # you should have one genotype and one value of depth of coverage per site - both matrices should have the same dimension
  if(identical(dim(genotypes), dim(indContribution)) == FALSE)
    stop(paste("The dimensions of the genotypes and indContribution matrices should be the same. Please check"))

  # the function getNumReadsR_vector gets as input vectors
  # hence, I can apply it to each column of genotypes and individual contribution
  tempReference <- sapply(1:ncol(genotypes), function(i) {
    getNumReadsR_vector(genotype_v = genotypes[,i], readCount_v = indContribution[,i], error = error)})

  # check that indContribution-tmp does not give negative values
  if(sum( (indContribution - tempReference) < 0 ) > 0 ) {
    stop("Error in creating number of reads!")
  }

  # output the number of reference reads
  tempReference
}


#' Compute the number of reference reads at multiple loci
#'
#' This function computes the number of reference reads over multiple loci and
#' for a single population.
#'
#' Note that this function will also work on a single locus, provided that the
#' input is in the list format.
#'
#' @param genotypes is a list, where each entry corresponds to a different
#'   locus. Each entry should be a matrix containing the genotypes (coded as 0,
#'   1 or 2). Each column of that matrix should be a different site and each row
#'   a different individual.
#' @param indContribution either a list or a matrix (that the function will
#'   convert to a list). Each list entry should be a matrix of individual
#'   contributions.  Each row of that matrix is a different individual and each
#'   column is a different site. Thus, each entry of the matrix should contain
#'   the number of reads contributed by that individual at that particular site.
#' @param error a numeric value with error rate associated with the sequencing
#'   and mapping process. This error rate is assumed to be symmetric:
#'   error(reference -> alternative) = error(alternative -> reference). This
#'   number should be between 0 and 1.
#'
#' @return a list with one entry per locus. Each of those entries is a matrix
#'   with the number of reference allele reads contributed by each individual.
#'   Each matrix row represents a different individual and each column is a
#'   different site.
#'
#' @keywords internal
#'
#' @export
numberReference <- function(genotypes, indContribution, error) {

  if(!inherits(genotypes, "list"))
    stop(paste("The genotypes input should be on a list format. Please check"))

  if(any(class(indContribution) != "list"))
    indContribution <- list(indContribution)

  if(length(genotypes) != length(indContribution))
    stop(paste("The genotypes and indContribution lists should have the same number of entries. Please check"))

  # get the number of loci
  nLoci <- length(genotypes)

  # When dealing with a single population - all the genotypes in the matrix correspond to a single population
  readsReference <- lapply(1:nLoci, function(locus) computeReference(genotypes[[locus]], indContribution[[locus]], error))

  # output the number of reads with the reference allele per individual and per site
  readsReference
}


#' Compute the number of reference reads for multiple populations
#'
#' This function computes the number of reference reads over a single locus for
#' multiple populations.
#'
#' Note that this function will not work as intended if the input consists of
#' multiple loci.
#'
#' @param genotypes either a list with a single entry (one locus) or a matrix
#'   (that the function will convert to a list) containing the genotypes (coded
#'   as 0, 1 or 2). Each column of that matrix should be a different site and
#'   each row a different individual.
#' @param indContribution a list where each entry contains the information for a
#'   single population. Each entry should be a matrix, with as many rows as the
#'   number of individuals of that population. Each row contains the number of
#'   contributed reads for a given individual and across all sites.
#' @param size a list with one entry per population. Each entry should be a
#'   vector containing the size (in number of diploid individuals) of each pool.
#'   Thus, if a population was sequenced using a single pool, the vector should
#'   contain only one entry. If a population was sequenced using two pools, each
#'   with 10 individuals, this vector should contain two entries and both will
#'   be 10.
#' @param error a numeric value with error rate associated with the sequencing
#'   and mapping process. This error rate is assumed to be symmetric:
#'   error(reference -> alternative) = error(alternative -> reference). This
#'   number should be between 0 and 1.
#'
#' @return a list with one entry per population. Each entry contains the number
#'   of reference allele reads for the individuals of that population and for
#'   that locus. Different individuals are in different rows and each columns
#'   represents a different site.
#'
#' @examples
#' # simulate coverage at 5 SNPs for two populations, assuming 20x mean coverage
#' reads <- simulateCoverage(mean = c(20, 20), variance = c(100, 100), nSNPs = 5, nLoci = 1)
#'
#' # simulate the number of reads contributed by each individual
#' # for each population there are two pools, each with 5 individuals
#' indContribution <- popsReads(list_np = rep(list(rep(5, 2)), 2), coverage = reads, pError = 5)
#'
#' # set seed and create a random matrix of genotypes for the 20 individuals - 10 per population
#' set.seed(10)
#' genotypes <- matrix(rpois(100, 0.5), nrow = 20)
#'
#' # simulate the number of reference reads for the two populations
#' numberReferencePop(genotypes = genotypes, indContribution = indContribution,
#' size = rep(list(rep(5, 2)), 2), error = 0.01)
#'
#' @export
numberReferencePop <- function(genotypes, indContribution, size, error) {

  # check if the indContribution input is the right format
  if(any(class(indContribution) != "list"))
    stop(paste("The indContribution input should be on a list format. Please check"))

  # get the number of populations
  nPops <- length(size)

  # set the output when a given locus has no SNPs
  if(any(is.na(unlist(indContribution)))) {
    # set the output to NA
    reference <- list(matrix(data = NA, nrow = nPops), matrix(data = NA, nrow = nPops))

  } else {

    # when the genotypes input is not a list, convert with to a list
    if(any(class(genotypes) != "list"))
      genotypes <- list(genotypes)

    # check if the indContribution list contains one list entry per population
    if(length(indContribution) != nPops)
      stop(paste("The indContribution input should have one list entry per population. Please check"))

    # use the splitMatrix function to split the genotypes into different list entries for each population
    tempGeno <- sapply(genotypes, FUN = function(geno) splitMatrix(geno, size))

    # compute the number of reads with the reference allele, for each individual and across all sites and populations
    reference <- sapply(1:nPops, function(pop)
      numberReference(genotypes = list(tempGeno[[pop]]), indContribution = indContribution[[pop]], error))
  }

  # output the number of reference reads
  reference
}


#' Create Pooled DNA sequencing data for multiple populations
#'
#' This function combines the information for each individual of each population
#' into information at the population level.
#'
#' In other words, the information of all individuals in a given population is
#' combined into a single population value and this is done for the various
#' populations. In this situation, each entry of the `indContribution` and
#' `readsReference` lists should contain one entry per population - being, in
#' essence, a list within a list. Please note that this function is intended to
#' work for multiple populations and should not be used with a single
#' population.
#'
#' @param nPops An integer representing the total number of populations in the
#'   dataset.
#' @param nLoci An integer that represents the total number of independent loci
#'   in the dataset.
#' @param indContribution  Either a list or a matrix (when dealing with a single
#'   locus).
#' @param readsReference A list, where each entry contains the information for a
#'   single locus. Each list entry should then have one separate entry per
#'   population. Each of these entries should be a matrix, with each row
#'   corresponding to a single individual and each column a different site.
#'   Thus, each entry of the matrix contains the number of observed reads with
#'   the reference allele for that individual at a given site. The output of the
#'   `numberReference` or `numberReferencePop` functions should be the input
#'   here.
#'
#' @return a list with three names entries
#'
#'   \item{reference}{a list with one entry per locus. Each entry is a matrix
#'   with the number of reference allele reads for each population. Each column
#'   represents a different site and each row a different population.}
#'
#'   \item{alternative}{a list with one entry per locus. Each entry is a matrix
#'   with the number of alternative allele reads for each population. Each
#'   column represents a different site and each row a different population.}
#'
#'   \item{total}{a list with one entry per locus. Each entry is a matrix with
#'   the coverage of each population. Each column represents a different site
#'   and each row a different population.}
#'
#' @examples
#' # simulate coverage at 5 SNPs for two populations, assuming 20x mean coverage
#' reads <- simulateCoverage(mean = c(20, 20), variance = c(100, 100), nSNPs = 5, nLoci = 1)
#'
#' # simulate the number of reads contributed by each individual
#' # for each population there are two pools, each with 5 individuals
#' indContribution <- popsReads(list_np = rep(list(rep(5, 2)), 2), coverage = reads, pError = 5)
#'
#' # set seed and create a random matrix of genotypes for the 20 individuals - 10 per population
#' set.seed(10)
#' genotypes <- matrix(rpois(100, 0.5), nrow = 20)
#'
#' # simulate the number of reference reads for the two populations
#' readsReference <- numberReferencePop(genotypes = genotypes, indContribution = indContribution,
#' size = rep(list(rep(5, 2)), 2), error = 0.01)
#'
#' # create Pooled DNA sequencing data for these two populations and for a single locus
#' poolPops(nPops = 2, nLoci = 1, indContribution = indContribution, readsReference = readsReference)
#'
#' @export
poolPops <- function(nPops, nLoci, indContribution, readsReference) {

  # when dealing with a single locus it's possible that the indContribution input is not on a list format
  if(nLoci == 1) {
    indContribution <- list(indContribution); readsReference <- list(readsReference)
  }

  # number of reads with the alternative allele is simply the total number of reads per individual minus the number of reads
  # with the reference allele for that individual
  readsAlternative <- lapply(1:nLoci, function(locus)
    mapply(function(individual, reference) FUN = individual - reference, SIMPLIFY = FALSE,
           individual = indContribution[[locus]], reference = readsReference[[locus]]))

  # now, since each entry (each locus) has independent entries for each population, we can simple perform colSums across the
  # various entries. This will sum the number of reads (with the reference, the alternative and the total number of reads)
  # across all individuals of a population
  # in this way, you get the total number of reads with each allele for each population
  referencePool <- lapply(readsReference, function(reference) matrix(t(sapply(reference, colSums)), nrow = nPops))
  alternativePool <- lapply(readsAlternative, function(alternative) matrix(t(sapply(alternative, colSums)), nrow = nPops))
  totalPool <- lapply(indContribution, function(total) matrix(t(sapply(total, colSums)), nrow = nPops))

  # combine the information about the different read types into a list, create names for each entry and output that list
  final_list <- list(referencePool, alternativePool, totalPool)
  names(final_list)=c("reference", "alternative", "total")
  return(final_list)
}


#' Define major and minor alleles
#'
#' This function checks which of the two simulated alleles (reference or
#' alternative) corresponds to the minor allele. This function can also be used
#' to remove sites according to a minor-allele reads threshold.
#'
#' More precisely, this function counts the number of reads with the reference
#' or alternative allele at each site and then sets the minor allele as the
#' least frequent of the two. This is done across all populations and so the
#' major and minor alleles are defined at a global level. Then if the
#' `min.minor` input is not NA, sites where the number of minor allele reads,
#' across all populations, is below the user-defined threshold are removed.
#'
#' @param reference is a matrix of reference allele reads. Each row of the
#'   matrix should be a different population and each column a different site.
#'   Thus, each entry of the matrix contains the number of observed reads with
#'   the reference allele for that population at a given site.
#' @param alternative is a matrix of alternative allele reads. Each row of the
#'   matrix should be a different population and each column a different site.
#'   Thus, each entry of the matrix contains the number of observed reads with
#'   the alternative allele for that population at a given site.
#' @param coverage is a matrix of total coverage. Each row of the matrix should
#'   be a different population and each column a different site. Thus, each
#'   entry of the matrix contains the total number of observed reads for that
#'   population at a given site.
#' @param min.minor is an integer representing the minimum allowed number of
#'   minor-allele reads. Sites that, across all populations, have less
#'   minor-allele reads than this threshold will be removed from the data.
#'
#' @return a list with three names entries
#'
#'   \item{major}{a list with one entry per locus. Each entry is a matrix with
#'   the number of major allele reads for each population. Each column
#'   represents a different site and each row a different population.}
#'
#'   \item{minor}{a list with one entry per locus. Each entry is a matrix with
#'   the number of minor allele reads for each population. Each column
#'   represents a different site and each row a different population.}
#'
#'   \item{total}{a list with one entry per locus. Each entry is a matrix with
#'   the coverage of each population. Each column represents a different site
#'   and each row a different population.}
#'
#' @examples
#' # simulate coverage at 5 SNPs for two populations, assuming 20x mean coverage
#' reads <- simulateCoverage(mean = c(20, 20), variance = c(100, 100), nSNPs = 5, nLoci = 1)
#'
#' # simulate the number of reads contributed by each individual
#' # for each population there are two pools, each with 5 individuals
#' indContribution <- popsReads(list_np = rep(list(rep(5, 2)), 2), coverage = reads, pError = 5)
#'
#' # set seed and create a random matrix of genotypes for the 20 individuals - 10 per population
#' set.seed(10)
#' genotypes <- matrix(rpois(100, 0.5), nrow = 20)
#'
#' # simulate the number of reference reads for the two populations
#' readsReference <- numberReferencePop(genotypes = genotypes, indContribution = indContribution,
#' size = rep(list(rep(5, 2)), 2), error = 0.01)
#'
#' # create Pooled DNA sequencing data for these two populations and for a single locus
#' pools <- poolPops(nPops = 2, nLoci = 1, indContribution = indContribution,
#' readsReference = readsReference)
#'
#' # define the major and minor alleles for this Pool-seq data
#' # we have to select the first entry of the pools list because this function works for matrices
#' findMinor(reference = pools$reference[[1]], alternative = pools$alternative[[1]],
#' coverage = pools$total[[1]])
#'
#' @export
findMinor <- function(reference, alternative, coverage, min.minor = NA) {

  # set the output for the situations where there is no SNP at the locus
  # check for NAs in one of the matrices
  if(any(is.na(unlist(reference)))) {
    # set the output to NA
    out <- list(reference = NA, alternative = NA, total = NA)

  } else {

    # create two temporary matrices - one containing the reads with the reference allele for each population and across all sites
    # and another with the reads with the alternative allele
    Ranc <- reference; Rder <- alternative

    # perform an evaluation to check if the total number or reads for each site
    # is bigger in the matrix containing the "reference" allele reads or in the matrix containing the "alternative" allele
    eval <- colSums(Ranc) < colSums(Rder)

    # at the sites where the number of reads for the "alternative" allele is bigger than the number of reads for the "reference" allele
    # replace those columns with the rows from the matrix of the "alternative" allele
    reference[, eval] <- Rder[, eval]
    # and then replace those same columns in the matrix with the reads of the "alternative" allele
    alternative[, eval] <- Ranc[, eval]

    # output the matrices with the various types of read numbers
    out <- list(major = reference, minor = alternative, total = coverage)
  }

  # output the matrices with the various types of read numbers
  out
}


#' Calculate population frequency at each SNP
#'
#' The frequency at a given SNP is calculated according to: `pi = c/r`, where c
#' = number of minor-allele reads and r = total number of observed reads.
#'
#' This function takes as input a list that contains the number of reads with
#' the minor allele and the number of total reads per population at a given
#' site. The names of the respective elements of the list should be minor and
#' total. It works with lists containing just one set of minor and total reads,
#' corresponding to a single locus, and with lists where each entry contains a
#' different set of minor and total number of reads, corresponding to different
#' loci.
#'
#' @param listPool a list containing the "minor" element, representing the
#'   number of reads with the minor-allele and the "total" element that contains
#'   information about the total number of Reads. The list should also contain a
#'   "major" entry with the information about reads containing the major-allele.
#'   The output of the `poolPops` function should be used as input here.
#' @param nLoci an integer that represents the total number of independent loci
#'   in the dataset.
#'
#' @return a list with two named entries
#'
#'   \item{pi}{a list with the allele frequencies of each population. Each list
#'   entry is a matrix, corresponding to a different locus. Each row of a matrix
#'   corresponds to a different population and each column to a different site.}
#'
#'   \item{pool}{a list with three different entries: major, minor and total.
#'   This list is similar to the one obtained with the \code{\link{findMinor}}
#'   function.}
#'
#' @examples
#' # simulate coverage at 5 SNPs for two populations, assuming 20x mean coverage
#' reads <- simulateCoverage(mean = c(20, 20), variance = c(100, 100), nSNPs = 5, nLoci = 1)
#'
#' # simulate the number of reads contributed by each individual
#' # for each population there are two pools, each with 5 individuals
#' indContribution <- popsReads(list_np = rep(list(rep(5, 2)), 2), coverage = reads, pError = 5)
#'
#' # set seed and create a random matrix of genotypes for the 20 individuals - 10 per population
#' set.seed(10)
#' genotypes <- matrix(rpois(100, 0.5), nrow = 20)
#'
#' # simulate the number of reference reads for the two populations
#' readsReference <- numberReferencePop(genotypes = genotypes, indContribution = indContribution,
#' size = rep(list(rep(5, 2)), 2), error = 0.01)
#'
#' # create Pooled DNA sequencing data for these two populations and for a single locus
#' pools <- poolPops(nPops = 2, nLoci = 1, indContribution = indContribution,
#' readsReference = readsReference)
#'
#' # define the major and minor alleles for this pool-seq data
#' # note that we have to select the first entry of the pools list
#' # because this function works for matrices
#' pools <- findMinor(reference = pools$reference[[1]], alternative = pools$alternative[[1]],
#' coverage = pools$total[[1]])
#'
#' # calculate population frequency at each SNP of this locus
#' calculatePi(listPool = pools, nLoci = 1)
#'
#' @export
calculatePi <- function(listPool, nLoci) {

  # the input should be the list containing information about reads per pool that was created by the previous function
  # a failsafe to detect if you're using the right list
  if(("minor" %in% names(listPool) | "total" %in% names(listPool)) == FALSE)
    stop(paste("Using an incorrect list"))

  # it is possible that each entry of the list is a single matrix - particularly when dealing with a single locus
  # if this is the case, then this step will convert those entries into lists
  if(any(sapply(listPool, class) == "matrix") == TRUE)
    listPool <- lapply(listPool, list)

  # by doing that transformation, we can use an lapply, whether we have one locus or multiple loci
  # divide the number of reads with the minor allele by the total number of reads - to obtain allelic frequencies
  Site_Pi <- lapply(1:nLoci, function(locus) listPool[["minor"]][[locus]]/listPool[["total"]][[locus]])

  # remove sites without information from the various categories of reads
  listPool[["minor"]] <- lapply(1:nLoci, function(locus)
    listPool[["minor"]][[locus]][, colSums(is.na(Site_Pi[[locus]])) == 0, drop = FALSE])

  listPool[["major"]] <- lapply(1:nLoci, function(locus)
    listPool[["major"]][[locus]][, colSums(is.na(Site_Pi[[locus]])) == 0, drop = FALSE])

  listPool[["total"]] <- lapply(1:nLoci, function(locus)
    listPool[["total"]][[locus]][, colSums(is.na(Site_Pi[[locus]])) == 0, drop = FALSE])

  # remove sites without information from the matrix containing the allelic frequencies
  Site_Pi <- lapply(1:nLoci, function(locus)
    Site_Pi[[locus]][, colSums(is.na(Site_Pi[[locus]])) == 0, drop = FALSE])

  # the output is a list that contains the information about the site frequencies
  # but also information about the total number of reads and the number of reads with the major/minor alleles
  list(pi = Site_Pi, pool = listPool)
}


#' Randomly select the required number of loci from the pooled sequencing data
#'
#' This function removes loci without polymorphic sites and then randomly
#' selects the required number of loci from a larger dataset.
#'
#' If extra loci were simulated to try to have sufficient loci to keep the
#' required number of loci after filtering, then this function is used to remove
#' extra loci. This is done by randomly selecting the required number of loci
#' from the full contingent of extra simulated loci.
#'
#' @param nSims is an integer representing how many types of simulations were
#'   performed. The possible types of simulations include loci without barriers
#'   against migration between divergent ecotypes, loci without migration from
#'   the C towards the W ecotype, loci without migration from the W towards the
#'   C ecotypes and loci where no migration occurs between divergent ecotypes.
#' @param pool a list containing the "minor" element, representing the number of
#'   reads with the minor-allele and the "total" element that contains
#'   information about the total number of Reads. The list should also contain a
#'   "major" entry with the information about reads containing the major-allele.
#' @param target is a vector with the required number of loci per category. If
#'   extra loci were simulated, this vector informs how many loci of each
#'   simulation type should be randomly selected.
#'
#' @return a list with three names entries
#'
#'   \item{major}{a list with one entry per locus. Each entry is a matrix with
#'   the number of major allele reads for each population. Each column
#'   represents a different site and each row a different population.}
#'
#'   \item{minor}{a list with one entry per locus. Each entry is a matrix with
#'   the number of minor allele reads for each population. Each column
#'   represents a different site and each row a different population.}
#'
#'   \item{total}{a list with one entry per locus. Each entry is a matrix with
#'   the coverage of each population. Each column represents a different site
#'   and each row a different population.}
#'
#' @keywords internal
#'
#' @export
forcePool <- function(nSims, pool, target) {

  # check the number of columns of the matrices with the number of minor-allele reads
  dimensions <- lapply(1:nSims, function(sim) sapply(pool[[sim]][["minor"]], ncol))

  # we only wish to keep loci where we have at least one polymorphic site
  # with this we obtain a vector that marks the loci without polymorphic sites with a FALSE
  tokeep <- lapply(dimensions, function(sim) sim != 0)
  # remove those loci from the data
  pool <- lapply(1:nSims, function(sim) lapply(pool[[sim]], function(i) i[tokeep[[sim]]]))

  # create a random vector of TRUE/FALSE to keep only the required number of loci per category
  tokeep <- lapply(1:nSims, function(sim)
    sample(x = c(rep(TRUE, target[sim]), rep(FALSE, length(pool[[sim]][["minor"]]) - target[sim])),
           size = length(pool[[sim]][["minor"]]), replace = F))
  # keep only the randomly selected loci
  pool <- lapply(1:nSims, function(sim) lapply(pool[[sim]], function(i) i[tokeep[[sim]]]))

  # convert the pool list back to the previous format:
  # one entry for major allele, one for minor allele and a final one for total coverage
  pool <- list(major = unlist(lapply(pool, function(x) x[["major"]]), recursive = FALSE),
               minor = unlist(lapply(pool, function(x) x[["minor"]]), recursive = FALSE),
               total = unlist(lapply(pool, function(x) x[["total"]]), recursive = FALSE))

  # output the pooled sequencing data
  pool
}

#' Simulate a single population
#'
#' Simulates the evolution of biological sequences for a single population with
#' variable theta values.
#'
#' @param nDip an integer representing the total number of diploid individuals
#'   to simulate. Note that scrm actually simulates haplotypes, so the number of
#'   simulated haplotypes is double of this.
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


#' Change to minor allele frequency
#'
#' Ensures that the allele frequencies computed directly from genotypes
#' correspond to the minor allele frequencies at each site.
#'
#' @param freqs a vector of allele frequencies where each entry corresponds to a
#'   different site.
#'
#' @return a vector with the same length as \code{freqs} with the allele
#'   frequencies above 0.5 replaced by the minor allele frequency.
#'
#' @examples
#' set.seed(10)
#' freqs <- runif(20)
#' changeFreqs(freqs)
#'
#' @export
changeFreqs <- function(freqs) {

  # check to see which sites have an allele frequency above 0.5
  toreplace <- freqs > 0.5
  # replace those sites by 1 - frequency to ensure that we have the minor allele frequency for each site
  freqs[toreplace] <- 1 - freqs[toreplace]

  # output the minor allele frequency computed directly from genotypes
  freqs
}


#' Compute allele frequencies from genotypes
#'
#' Computes alternative allele frequencies from genotypes by dividing the total
#' number of alternative alleles by the total number of gene copies.
#'
#' @param nDip an integer representing the total number of diploid individuals
#'   to simulate. Note that scrm actually simulates haplotypes, so the number of
#'   simulated haplotypes is double of this.
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

  # ensure that the allele frequencies correspond to the minor allele frequencies
  # by changing the frequencies above 0.5
  ifreqs <- lapply(ifreqs, FUN = function(locus) changeFreqs(locus))

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
#' @param minor a matrix with the number of minor-allele reads. Each row should
#'   be a different population and each column a different site.
#' @param coverage a matrix with the total coverage. Each row should be a
#'   different population and each column a different site.
#' @param min.minor is an integer representing the minimum allowed number of
#'   minor-allele reads. Sites that, across all populations, have less
#'   minor-allele reads than this threshold will be removed from the data.
#'
#' @return a list with three named entries. \code{freqs} contains the allele
#'   frequencies minus the frequency of any removed site. \code{minor} contains
#'   the number of minor-allele reads minus any removed site and the
#'   \code{coverage} entry contains the total coverage minus any removed site.
#'
#' @examples
#' freqs <- runif(20)
#' set.seed(10)
#' minor <- matrix(sample(x = c(0,5,10), size = 20, replace = TRUE), nrow = 1)
#' coverage <- matrix(sample(100:150, size = 20), nrow = 1)
#' removeSites(freqs = freqs, minor = minor, coverage, min.minor = 2)
#'
#' @export
removeSites <- function(freqs, minor, coverage, min.minor) {

  # get the total number of minor allele reads in the data
  tminor <- colSums(minor)
  # find out in which columns the total sum of the reads with the minor allele is below the threshold
  toremove <- tminor < min.minor

  # if there are sites where the sum of the reads with the minor allele is below the threshold
  if(sum(toremove) != 0) {

    # remove those columns from the matrix containing the depth of coverage
    coverage <- coverage[, !toremove, drop = FALSE]
    # remove those columns from the matrix containing the number of reads with the minor allele
    minor <- minor[, !toremove, drop = FALSE]

    # remove those entries from the vector containing the allele frequencies computed directly from genotypes
    freqs <- freqs[!toremove, drop = FALSE]
  }

  # create the output containing the frequencies computed from the genotypes
  # and the number of Pool-seq minor allele reads and total coverage
  out <- list(freqs = freqs, minor = minor, coverage = coverage)
  # output the results of the function
  out
}


#' Compute allele frequencies from pooled sequencing data
#'
#' Computes the minor-allele frequency from pooled data and removes any site
#' with too few minor-allele reads from both the pool frequencies and
#' frequencies computed directly from genotypes.
#'
#' The frequency at a given SNP is calculated according to: π = c/r, where c =
#' number of minor-allele reads and r = total number of observed reads.
#' Additionally, if a site has less minor-allele reads than \code{min.minor}
#' across all populations, that site is removed from the data.
#'
#' @param minor a matrix with the number of minor-allele reads. Each row should
#'   be a different population and each column a different site.
#' @param coverage a matrix with the total coverage. Each row should be a
#'   different population and each column a different site.
#' @param min.minor is an integer representing the minimum allowed number of
#'   minor-allele reads. Sites that, across all populations, have less
#'   minor-allele reads than this threshold will be removed from the data.
#' @param ifreqs a vector of allele frequencies computed directly from the
#'   genotypes where each entry corresponds to a different site.
#'
#' @return a list with two entries. The \code{ifreqs} entry contains the allele
#'   frequencies computed directly from genotypes and \code{pfreqs} the allele
#'   frequencies computed from pooled sequencing data.
#'
#' @examples
#' set.seed(10)
#' freqs <- runif(20)
#' freqs <- changeFreqs(freqs)
#' set.seed(10)
#' minor <- matrix(sample(x = c(0,5,10), size = 20, replace = TRUE), nrow = 1)
#' coverage <- matrix(sample(100:150, size = 20), nrow = 1)
#' Pfreqs(minor = minor, coverage = coverage, min.minor = 2, ifreqs = freqs)
#'
#' @export
Pfreqs <- function(minor, coverage, min.minor, ifreqs) {

  # use the removeSites function to remove sites with less than `min.minor` minor allele reads
  temp <- removeSites(freqs = ifreqs, minor = minor, coverage = coverage, min.minor = min.minor)

  # compute the allele frequencies for Pool-seq data
  freqs <- temp[["minor"]]/temp[["coverage"]]
  # create the output in this instance - with the allele frequencies computed from genotypes and from Pool-seq data
  freqs <- list(ifreqs = temp[["freqs"]], pfreqs = freqs)

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
#'   to simulate. Note that scrm actually simulates haplotypes, so the number of
#'   simulated haplotypes is double of this.
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
#' @export
maePool <- function(nDip, nloci, pools, pError, sError, mCov, vCov, min.minor) {

  # run SCRM and obtain genotypes for a single population
  genotypes <- run_scrm(nDip = nDip, nloci = nloci)

  # compute allele frequencies directly from the genotypes
  ifreqs <- Ifreqs(nDip = nDip, genotypes = genotypes)

  # simulate number of reads
  reads <- simulateCoverage(mean = mCov, variance = vCov, genotypes = genotypes)

  # simulate individual contribution to the total number of reads
  indContribution <- lapply(1:nloci, function(locus)
    popsReads(list_np = pools, coverage = reads[[locus]], pError = pError))

  # simulate the number of reference reads
  reference <- lapply(1:nloci, function(locus)
    numberReferencePop(genotypes = genotypes[[locus]], indContribution = indContribution[[locus]], size = pools, error = sError))

  # simulate pooled sequencing data
  pool <- poolPops(nPops = 1, nLoci = nloci, indContribution = indContribution, readsReference = reference)

  # use an lapply to ensure that the ancestral allele of the simulations is also the major allele - do this for each locus
  pool <- lapply(1:nloci, function(locus)
    minorPool(reference = pool[["reference"]][[locus]], alternative = pool[["alternative"]][[locus]],
              coverage = pool[["total"]][[locus]]))

  # convert the pool list back to the previous format
  # one entry for ancestral matrices, one for derived and a final one for total matrices
  pool <- list(major = lapply(pool, function(locus) locus[["major"]]),
               minor = lapply(pool, function(locus) locus[["minor"]]),
               total = lapply(pool, function(locus) locus[["total"]]))

  # compute the allele frequencies obtained with pooled sequencing
  pfreqs <- lapply(1:nloci, function(locus)
    Pfreqs(minor = pool[["minor"]][[locus]], coverage = pool[["total"]][[locus]], min.minor = min.minor, ifreqs = ifreqs[[locus]]))

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
#'   diploid individuals to simulate. Note that scrm actually simulates
#'   haplotypes, so the number of simulated haplotypes is double of this. If it
#'   is a vector, then each vector entry will be simulated independently. For
#'   instance, if \code{nDip = c(100, 200)}, simulations will be carried out for
#'   samples of 100 and 200 individuals.
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
maeFreqs <- function(nDip, nloci, pError, sError, mCov, vCov, min.minor) {

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
                    vCov = varCov, min.minor = min.minor)

    # add those values to the dataframe containing all the results
    final <- rbind(final, temp)
  }

  # remove the first row of the final matrix - this row contains NAs
  final <- final[-1 ,]
  # output the final dataframe with the MAE values for the different parameter combinations
  final
}


#' Draw parameters from the priors
#'
#' This function creates a named vector of parameters that can be used as input
#' in the command line of the scrm package. Please note that this function needs
#' to be adjusted if you wish to test the effect of different prior
#' distributions.
#'
#' @param Nref The minimum and maximum value of the uniform distribution for the
#'   effective population size of the reference population (Nref).
#' @param ratio The minimum and maximum value of the distribution from which the
#'   relative size of the present-day and ancestral populations are drawn. The
#'   size of these populations is set as a ratio of the size of the Nref
#'   population. All of these ratios are drawn from a log10 uniform
#'   distribution.
#' @param split The minimum and maximum values, at the 4Nref scale, of the
#'   uniform distribution from which the values of the times of the split events
#'   are draw. Both the time of the recent split event and the distance between
#'   the two split events are drawn from this distribution.
#' @param pool The minimum and maximum values of the uniform distribution from
#'   which the value of the error associated with DNA pooling is drawn. More
#'   specifically, this value is related with the unequal individual
#'   contribution to the pool.
#' @param seq The minimum and maximum values of the uniform distribution from
#'   which the value of the error associated with DNA sequencing is drawn. This
#'   parameter should be supplied as a decimal number between zero and one.
#' @param CW The minimum and maximum value of the uniform distribution from
#'   which the migration rate between the two divergent ecotypes inhabiting the
#'   same location is drawn. We consider that this parameter is drawn on a m
#'   scale. This is the migration rate from ecotype C to ecotype W.
#' @param WC The minimum and maximum value of the uniform distribution from
#'   which the migration rate between the two divergent ecotypes inhabiting the
#'   same location is drawn. We consider that this parameter is drawn on a m
#'   scale. This is the migration rate from ecotype W to ecotype C.
#' @param CC The minimum and maximum value of the uniform distribution from
#'   which the migration rate between similar ecotypes inhabiting different
#'   locations is drawn. We consider that this parameter is drawn on a m scale.
#'   This is the migration between the two C ecotypes at two different
#'   locations.
#' @param WW The minimum and maximum value of the uniform distribution from
#'   which the migration rate between similar ecotypes inhabiting different
#'   locations is drawn. We consider that this parameter is drawn on a m scale.
#'   This is the migration between the two W ecotypes at two different
#'   locations.
#' @param ANC The minimum and maximum value of the uniform distribution from
#'   which the migration rate between the two ancestral populations is drawn. We
#'   consider that this parameter is drawn on a m scale.
#' @param bT The minimum and maximum values of the distribution from which the
#'   proportion of the simulated loci where no migration occurs between
#'   divergent ecotypes is drawn. The maximum value should not be higher than
#'   one.
#' @param bCW The minimum and maximum values of the distribution from which the
#'   proportion of the simulated loci where no migration occurs from the C
#'   ecotype towards the W ecotype is drawn. The maximum value should not be
#'   higher than one.
#' @param bWC The minimum and maximum values of the distribution from which the
#'   proportion of the simulated loci where no migration occurs from the W
#'   ecotype towards the C ecotype is drawn. The maximum value should not be
#'   higher than one.
#' @param model Either "2pops", "Single" or "Parallel" indicating for which
#'   model should parameters be drawn.
#' @param digits An optional integer indicating the number of decimal places to
#'   use when rounding certain parameters. The default is five.
#'
#' @return a vector with one named entry per relevant parameter. Each entry is
#'   the sampled value from the prior for that particular parameter.
#'
#' @examples
#' # for a model with two populations
#' createParams(Nref = c(25000, 25000), ratio = c(0.1, 3), pool = c(5, 250), seq = c(0.0001, 0.001),
#' split = c(0, 3), CW = c(1e-13, 1e-3), WC = c(1e-13, 1e-3), bT = c(0, 0.2), model = "2pops")
#'
#' # for a single origin scenario
#' createParams(Nref = c(25000, 25000), ratio = c(0.1, 3), pool = c(5, 250), seq = c(0.0001, 0.001),
#' split = c(0, 3), CW = c(1e-13, 1e-3), WC = c(1e-13, 1e-3), CC =  c(1e-13, 1e-3),
#' WW = c(1e-13, 1e-3), ANC = c(1e-13, 1e-3), bT = c(0, 0.2), bCW = c(0, 0.5),
#' bWC = c(0, 0.5), model = "Single")
#'
#' @export
createParams <- function(Nref, ratio, split, pool, seq, CW, WC, CC = NA, WW = NA, ANC = NA, bT, bCW = NA, bWC = NA,
                         model, digits = 5) {

  # check if the input is correct - the selected model should be one of the following
  if(model %in% c("2pops", "Single", "Parallel") == FALSE)
    stop("The selected model should be either 2pops, Single or Parallel. Please check")

  # draw a value for the population size of the most ancestral population - it's also the Nref
  Nrf <- stats::runif(n = 1, min = Nref[1], max = Nref[2])
  # Draw the values for the extant and ancestral population sizes
  # Values are drawn as a ratio of the Ne size - a value of 2 means that the population is twice the size of the Ne
  N1 <- round(stats::runif(n = 1, min = log10(ratio[1]), max = log10(ratio[2])), digits = 5)
  N2 <- round(stats::runif(n = 1, min = log10(ratio[1]), max = log10(ratio[2])), digits = 5)
  # set the values in the natural scale
  N1 <- 10^N1; N2 <- 10^N2

  # Draw the split times from a uniform distribution - times are drawn on a 4Nref scale
  Split <- stats::runif(n = 1, min = split[1], max = split[2])
  # Draw the errors for the pooling parameter
  Pool_Error <- stats::runif(n = 1, min = pool[1], max = pool[2])
  # and the sequencing error
  Error <- stats::runif(n = 1, min = seq[1], max = seq[2])

  # draw the migration rate (at the m scale) - this corresponds to the migration between different ecotypes in the same site
  # the migration rate from crab to wave is drawn as mCW - for the first site
  mCW1 <- stats::runif(n = 1, min = CW[1], max = CW[2])
  # and the migration rate from wave to crab is drawn as mWC - for the first site
  mWC1 <- stats::runif(n = 1, min = WC[1], max = WC[2])

  # proportion of the genome without migration - total barrier
  total <- stats::rbeta(n = 1, shape1 = 1, shape2 = 10)
  # replace values below the minimum threshold with the minimum
  total[total < 1e-2] <- bT[1]
  # and values above the maximum threshold with the maximum
  total[total > bT[2]] <- bT[2]

  # stop the function if we are working with the two-population model
  if(model == "2pops") {
    # assume that the proportion of the genome with unrestricted migration is 1 - proportion without migration
    pMig <- 1 - total
    # create the parameters vector for this particular model
    parameters <- c(Nrf, N1, N2, Split, Pool_Error, Error, mCW1, mWC1, pMig, total)
    # add names to the entries of the vector
    names(parameters) <- c("Nref", "N1", "N2", "Split", "PoolError", "SeqError", "mCW", "mWC", "pM", "pNO")
    # stop the function and output the parameters vector
    stop(return(parameters))
  }

  # draw the values for the remaining present-day population sizes
  N3 <- round(stats::runif(n = 1, min = log10(ratio[1]), max = log10(ratio[2])), digits = 5)
  N4 <- round(stats::runif(n = 1, min = log10(ratio[1]), max = log10(ratio[2])), digits = 5)
  # set the values in the natural scale
  N3 <- 10^N3; N4 <- 10^N4

  # draw the values for the ancestral population sizes
  # values are drawn as a ratio of the Ne size - a value of 2 means that the population is twice the size of the Ne
  NA1 <- round(stats::runif(n = 1, min = log10(ratio[1]), max = log10(ratio[2])), digits = 5)
  NA2 <- round(stats::runif(n = 1, min = log10(ratio[1]), max = log10(ratio[2])), digits = 5)
  # set the values in the natural scale
  NA1 <- 10^NA1; NA2 <- 10^NA2

  # draw the second split time from a uniform distribution - times are drawn on a 4Ne scale
  Dsplit <- stats::runif(n = 1, min = split[1], max = split[2])

  # draw the migration rate (at the m scale) - this corresponds to the migration between different ecotypes in the same site
  # the migration rate from crab to wave is drawn as mCW - for the second site
  mCW2 <- stats::runif(n = 1, min = CW[1], max = CW[2])
  # and the migration rate from wave to crab is drawn as mWC - for the second site
  mWC2 <- stats::runif(n = 1, min = WC[1], max = WC[2])

  # if required, draw additional migration rates
  if(!any(is.na(CC)))
    mCC <- stats::runif(n = 1, min = CC[1], max = CC[2]) # between crab ecotypes in different locations
  else
    mCC <- NA  # set mCC to NA

  # if required, draw additional migration rates
  if(!any(is.na(WW)))
    mWW <- stats::runif(n = 1, min = WW[1], max = WW[2]) # between wave ecotypes in different locations
  else
    mWW <- NA  # set mWW to NA

  # if required, draw additional migration rates
  if(!any(is.na(ANC)))
    mAA <- stats::runif(n = 1, min = ANC[1], max = ANC[2])  # migration rates between the ancestral populations
  else
    mAA <- NA # set mAA to NA

  # if required, draw the proportion of the genome without migration
  if(!any(is.na(bCW))) {
    # from the crab to the wave ecotype
    pCW <- stats::rbeta(n = 1, shape1 = 1, shape2 = 10)
    # replace values below the minimum threshold with the minimum
    pCW[pCW < 1e-2] <- bCW[1]
    # and values above the maximum threshold with the maximum
    pCW[pCW > bCW[2]] <- bCW[2]

  } else {

    # set pCW to NA
    pCW <- NA
  }

  # if required, draw the proportion of the genome without migration
  if(!any(is.na(bWC))) {
    # from the wave to the crab ecotype
    pWC <- stats::rbeta(n = 1, shape1 = 1, shape2 = 10)
    # replace values below the minimum threshold with the minimum
    pWC[pWC < 1e-2] <- bWC[1]
    # and values above the maximum threshold with the maximum
    pWC[pWC > bWC[2]] <- bWC[2]

  } else {

    # set pWC to NA
    pWC <- NA
  }

  # assume that the proportion of the genome with unrestricted migration is 1 - proportion without migration
  pMig <- 1 - sum(c(total, pCW, pWC), na.rm = TRUE)

  # create the parameters vector for the single and parallel models
  parameters <- c(Nrf, N1, N2, N3, N4, NA1, NA2, Split, Dsplit, Pool_Error, Error, mCW1, mCW2, mWC1, mWC2, mCC, mWW, mAA,
                  pMig, pCW, pWC, total)

  # add names to the entries of the vector
  names(parameters) <- c("Nref", "N1", "N2", "N3", "N4", "NA1", "NA2", "Split", "Dsplit", "PoolError", "SeqError", "mCW1", "mCW2",
                         "mWC1", "mWC2", "mCC", "mWW", "mAA", "pM", "pCW", "pWC", "pNO")

  # remove any parameters that are set as NA
  parameters <- parameters[!is.na(parameters)]

  # output the parameters vector
  parameters
}


#' Create SCRM command line for a model with two populations
#'
#' This function creates a command line tailored for an isolation with migration
#' model with two populations. The command line can then be fed to the scrm
#' package to run the model.
#'
#' @param parameters A vector where each entry corresponds to a different
#'   parameter, e.g. one entry is the size of the reference population, another
#'   is the time of recent split, etc. Please note that this functions depends
#'   on the ordering of the parameters in the vector and thus, it should only be
#'   used with a vector created with the `createParams` function.
#' @param nSites An integer representing the number of base pairs that each
#'   locus should have.
#' @param nLoci An integer that represents how many independent loci should be
#'   simulated.
#' @param nDip An integer representing the total number of diploid individuals
#'   to simulate. Note that scrm actually simulates haplotypes, so the number of
#'   simulated haplotypes is double of this. Also note that this is the total
#'   number of diploid individuals and this function will distribute the
#'   individuals equally by the two populations.
#' @param mutrate A number representing the mutation rate assumed for the
#'   simulations.
#' @param extra is a logical value indicating whether the required number of
#'   loci should be enforced. The default is FALSE but, if set to TRUE, then
#'   additional loci will be simulated. These additional loci are simulated to
#'   try to have sufficient loci to keep the required number of loci after
#'   filtering.
#'
#' @return a character vector with two entries. The first entry is the scrm
#'   command line for the loci without any barriers against migration, while the
#'   second entry is the scrm command line for the loci without migration
#'   between divergent ecotypes.
#'
#' @examples
#' # create a vector with parameter values for a two populations model
#' params <- createParams(Nref = c(25000, 25000), ratio = c(0.1, 3), pool = c(5, 250),
#' seq = c(0.0001, 0.001), split = c(0, 3), CW = c(1e-13, 1e-3), WC = c(1e-13, 1e-3),
#' bT = c(0, 0.2), model = "2pops")
#'
#' # create the command line for the scrm package
#' cmd2pops(parameters = params, nSites = 2000, nLoci = 100, nDip = 100, mutrate = 2e-8)
#'
#' @export
cmd2pops <- function(parameters, nSites, nLoci, nDip, mutrate, extra = FALSE) {

  # this function is intended to be used with a two-population model
  nPops <- 2

  # Read the vector with the parameters and assign each parameter to the correct command name
  Ne <- parameters[1]
  # set the relative size of each population
  N1 <- parameters[2]; N2 <- parameters[3]

  # get the proportion of loci with migration - no barriers against migration between the different ecotypes
  pM <- parameters[9]
  # and the proportion of loci without any migration - total barrier against migration between the different ecotypes
  pNO <- parameters[10]

  # get the migration rates - between different ecotypes at the same site - this is the m value on the M = 4N0m formula
  # from the crab ecotype to the wave ecotype
  mCW <-  parameters[7]
  # from the wave ecotype to the crab ecotype
  mWC <-  parameters[8]

  # m <i> <j> <M>: Set the migration rate from population j to population i to M (looking forward in time)
  # and REMEMBER that M = 4N0m
  # set the migration rates to the scale of Nref - between different ecotypes at the same site
  mig_CW <- 4*Ne*mCW # -m 2 1 mig_CW

  # the migration from wave to crab is parametrized as a ratio of the migration from crab to wave
  # so we need to multiply the migration from crab to wave by this ratio
  mig_WC <- 4*Ne*mWC # -m 1 2 mig_WC

  # get the time of the split event
  split <- round(parameters[4], digits = 3)

  # Compute the value of theta
  mutrate_locus <- nSites*mutrate
  theta <- 4*Ne*mutrate_locus
  # Create a vector with information about how many haplotypes are sampled from each population
  n <- c(rep(nDip/nPops, times = nPops))

  # use a multinomial distribution to get the number of loci simulated under each category
  lociTotal <- as.vector(stats::rmultinom(n = 1, size = nLoci, c(pM, pNO)))

  # if extra is TRUE, then more loci than required per category will be simulated
  if(extra == TRUE) {
    # save the required number of simulated loci per category
    targetLoci <- lociTotal
    # simulate more 25 loci per category
    lociTotal <- lociTotal + 25
  }

  # cheat code: pop1 - crab in site 1; pop2 - wave in site 1; pop3 - crab in site 2; pop4 - wave in site 2
  # create command line with no barriers to migration
  # set the basic elements for scrm - nhap: total number of haplotypes that are simulated at each locus and
  # nrep: the number of independent loci that will be produced
  with.mig <- paste(paste(nDip*2, collapse = " "), lociTotal[1], "-t",theta, "-I", paste(nPops), paste(n*2, collapse = " "),
                    # set the size of the present day populations - n <i> <n> Set the size of population i to n*N0.
                    "-n 1", N1, "-n 2", N2,
                    # m <i> <j> <M>: Set the migration rate from population j to population i to M (looking forward in time)
                    # set the migration rate between different ecotypes inhabiting the same location
                    "-m 2 1", mig_CW, "-m 1 2", mig_WC,
                    # now, set the migration rate right before the split event to zero by using the switch:
                    # -eM <t> <M>: assume a symmetric migration rate of M/(npop-1) at time t.
                    "-eM", split, "0",
                    # add a split event -ej <t> <j> <i> in population i that creates population j (forwards in time)
                    # finally, set the size of the ancestral pop equal to the size of the reference population with:
                    # -eN <t> <n>: set the size of all populations to n*N0 at time t.
                    "-ej", split, "2 1 -eN", split, 1)

  # create a command line for the loci without any migration (between the different ecotypes)
  # set the basic elements for scrm - nhap: total number of haplotypes that are simulated at each locus and
  # nrep: the number of independent loci that will be produced
  without.mig <- paste(paste(nDip*2, collapse = " "), lociTotal[2], "-t",theta, "-I", paste(nPops), paste(n*2, collapse = " "), "0",
                       # set the size of the present day populations - n <i> <n> Set the size of population i to n*N0.
                       "-n 1", N1, "-n 2", N2,
                       # add a split event -ej <t> <j> <i> in population i that creates population j (forwards in time)
                       # finally, set the size of the ancestral pop equal to the size of the reference population with:
                       # -eN <t> <n>: set the size of all populations to n*N0 at time t.
                       "-ej", split, "2 1 -eN", split, 1)


  # combine the two different types of commands
  cmd_2pops <- c(with.mig, without.mig)

  # if extra is equal to TRUE, then we simulated more loci than required
  if(extra == TRUE)
    # include the required number of loci per category in the output
    cmd_2pops <- list(commands = cmd_2pops, targetLoci = targetLoci)

  # output the command line for the two population model
  cmd_2pops
}

#' Create SCRM command line for a single origin scenario
#'
#' This function creates a command line tailored for a scenario of single origin
#' to explain ecotype formation. The command line can then be fed to the scrm
#' package to run the model.
#'
#' For convenience, imagine we have two divergent ecotypes, named C and W. This
#' model assumes that the first population corresponds to the C ecotype at the
#' first location, the second population to the C ecotype in the first location,
#' the third population to the W ecotype in the second location and the fourth
#' population to the W ecotype in the second location.
#'
#' @param parameters A vector where each entry corresponds to a different
#'   parameter, e.g. one entry is the size of the reference population, another
#'   is the time of recent split, etc. Please note that this functions depends
#'   on the ordering of the parameters in the vector and thus, it should only be
#'   used with a vector created with the `createParams` function.
#' @param nSites An integer representing the number of base pairs that each
#'   locus should have.
#' @param nLoci An integer that represents how many independent loci should be
#'   simulated.
#' @param nDip An integer representing the total number of diploid individuals
#'   to simulate. Note that scrm actually simulates haplotypes, so the number of
#'   simulated haplotypes is double of this. Also note that this is the total
#'   number of diploid individuals and this function will distribute the
#'   individuals equally by the two populations.
#' @param mutrate A number representing the mutation rate assumed for the
#'   simulations.
#' @param extra is a logical value indicating whether the required number of
#'   loci should be enforced. The default is FALSE but, if set to TRUE, then
#'   additional loci will be simulated. These additional loci are simulated to
#'   try to have sufficient loci to keep the required number of loci after
#'   filtering.
#'
#' @return a character vector with four entries. The first entry is the scrm
#'   command line for the loci without any barriers against migration. The
#'   second entry is the command line for the loci without migration from the C
#'   towards the W ecotype. The third entry is command line for the loci without
#'   migration from the W towards the C ecotype and the last entry is the scrm
#'   command line for the loci without migration between divergent ecotypes.
#'
#' @examples
#' # create a vector with parameter values for the single origin scenario
#' params <- createParams(Nref = c(25000, 25000), ratio = c(0.1, 3), pool = c(5, 250),
#' seq = c(0.0001, 0.001), split = c(0, 3), CW = c(1e-13, 1e-3), WC = c(1e-13, 1e-3),
#' CC =  c(1e-13, 1e-3), WW = c(1e-13, 1e-3), ANC = c(1e-13, 1e-3), bT = c(0, 0.2),
#' bCW = c(0, 0.5), bWC = c(0, 0.5), model = "Single")
#'
#' # create the command line for the scrm package
#' cmdSingle(parameters = params, nSites = 2000, nLoci = 100, nDip = 400, mutrate = 2-8)
#'
#' @export
cmdSingle <- function(parameters, nSites, nLoci, nDip, mutrate, extra = FALSE) {

  # this function is intended to be used with a four-population model
  nPops <- 4

  # read the vector with the parameters and assign each parameter to the correct variable name
  Ne <- parameters[1]
  # set the relative size of each population - for the extant populations
  N1 <- parameters[2]; N2 <- parameters[3]; N3 <- parameters[4]; N4 <- parameters[5]
  # and the ancient populations
  NA1 <- parameters[6]; NA2 <- parameters[7]

  # get the proportion of loci with migration - no barriers against migration between the different ecotypes
  pM <- parameters[19]
  # get the proportion of loci without migration - from the crab to the wave ecotype at the same location
  pCW <- parameters[20]
  # get the proportion of loci without migration - from the wave to the crab ecotype at the same location
  pWC <- parameters[21]
  # and the proportion of loci without any migration - total barrier against migration between the different ecotypes
  pNO <- parameters[22]

  # get the migration rates - between different ecotypes at the same site - this is the m value on the M = 4N0m formula
  # from the crab ecotype to the wave ecotype - at the first site
  mCW1 <- parameters[12]
  # from the crab ecotype to the wave ecotype - at the second site
  mCW2 <- parameters[13]
  # from the wave ecotype to the crab ecotype - at the fist site
  mWC1 <- parameters[14]
  # from the wave ecotype to the crab ecotype - at the second site
  mWC2 <- parameters[15]

  # between crab populations inhabiting different locations
  mCC <- parameters[16]
  # between wave populations inhabiting different locations
  mWW <- parameters[17]
  # and between the two ancestral populations
  mAA <- parameters[18]

  # m <i> <j> <M>: Set the migration rate from population j to population i to M (looking forward in time)
  # REMEMBER that M = 4N0m
  # set the migration rates to the scale of Nref - between different ecotypes at the same site
  # from the crab ecotype to the wave ecotype - at the first site
  mig_CW1 <- 4*Ne*mCW1 # -m 3 1 mig_CW
  # from the crab ecotype to the wave ecotype - at the second site
  mig_CW2 <- 4*Ne*mCW2 # -m 4 2 mig_CW
  # from the wave ecotype to the crab ecotype - at the first site
  mig_WC1 <- 4*Ne*mWC1 # -m 1 3 mig_WC
  # from the wave ecotype to the crab ecotype - at the second site
  mig_WC2 <- 4*Ne*mWC2 # -m 2 4 mig_WC

  # between crab populations at different locations
  mig_CC <- 4*Ne*mCC
  # between wave populations at different locations
  mig_WW <- 4*Ne*mWW
  # and between the two ancestral populations
  mig_AA <- 4*Ne*mAA

  # get the time of the recent split event
  Rsplit <- round(parameters[8], digits = 3)
  # and the ancient split event
  Asplit <- round(parameters[9], digits = 3)
  # the actual time of the ancient split event is obtained by doing:
  Asplit <- Rsplit + Asplit

  # Compute the value of theta
  mutrate_locus <- nSites*mutrate
  theta <- 4*Ne*mutrate_locus
  # Create a vector with information about how many haplotypes are sampled from each population
  n <- c(rep(nDip/nPops, times = nPops))

  # create variables to ensure that the changes in migration rates occur at different times than the split time
  # a variable to inform when does migration start between ancestral populations
  tmAA <- Rsplit + 0.0001
  # create a variable to inform when does migration stop before the recent split
  # if Rsplit is zero, we can not subtract something from it
  if(Rsplit != 0) {
    # when Rsplit is not zero, set the end of the migration before the split
    tmRS <- Rsplit - 0.0001
  } else {
    # when Rsplit is zero, set the end of the migration at the split
    tmRS <- Rsplit
  }

  # create also a variable to inform when does migration start between ancestral populations
  tmAA <- Rsplit + 0.0001

  # use a multinomial distribution to get the number of loci simulated under each category
  lociTotal <- as.vector(stats::rmultinom(n = 1, size = nLoci, c(pM, pCW, pWC, pNO)))

  # if extra is TRUE, then more loci than required per category will be simulated
  if(extra == TRUE) {
    # save the required number of simulated loci per category
    targetLoci <- lociTotal
    # simulate more 25 loci per category
    lociTotal <- lociTotal + 25
  }

  # cheat code: pop1 - crab in site 1; pop2 - crab in site 2; pop3 - wave in site 1; pop4 - wave in site 2

  # create command line with no barriers to migration
  # set the basic elements for scrm - nhap: total number of haplotypes that are simulated at each locus and
  # nrep: the number of independent loci that will be produced
  with.mig <- paste(paste(nDip*2, collapse = " "), lociTotal[1], "-t",theta, "-I", paste(nPops), paste(n*2, collapse = " "),
                    # set the size of the present day populations - n <i> <n> Set the size of population i to n*N0.
                    "-n 1", N1, "-n 2", N2, "-n 3", N3, "-n 4", N4,
                    # m <i> <j> <M>: Set the migration rate from population j to population i to M (looking forward in time)
                    # set the migration rate between different ecotypes inhabiting the same location
                    "-m 3 1", mig_CW1, "-m 4 2", mig_CW2, "-m 1 3", mig_WC1, "-m 2 4", mig_WC2,
                    # set the migration rate between the same ecotypes inhabiting different locations
                    "-m 2 1", mig_CC, "-m 1 2", mig_CC, "-m 4 3", mig_WW, "-m 3 4", mig_WW,
                    # set the migration, between all populations, to zero - immediately before the recent split
                    "-eM", tmRS, "0",
                    # add a split event -ej <t> <j> <i> in population i that creates population j (forwards in time)
                    # looking backwards in time, it moves all lines from population j into population i at time t
                    # Migration rates into population j are set to 0 for the time further back into the past
                    "-ej", Rsplit, "1 2 -ej", Rsplit, "4 3",
                    # set the size of the ancient populations
                    # -en <t> <i> <n>: Set the size of population i to n*N0 at time t.
                    "-en", Rsplit, "2", NA1, "-en", Rsplit, "3", NA2,
                    # set the migration rate between the ancestral populations
                    "-em", tmAA, "2 3", mig_AA, "-em", tmAA, "3 2", mig_AA,
                    # set the migration, between all populations, to zero - immediately at the ancient split
                    "-eM", Asplit, "0",
                    # add a split event - this event creates the two ancestral populations
                    "-ej", Asplit, "2 3",
                    # finally, set the size of the most ancestral pop equal to the size of the reference population with:
                    # -eN <t> <n>: set the size of all populations to n*N0 at time t.
                    "-eN", Asplit, 1)

  # create command line with no migration from the crab to the wave ecotype within the same location
  no.mCW <- paste(paste(nDip*2, collapse = " "), lociTotal[2], "-t",theta, "-I", paste(nPops), paste(n*2, collapse = " "),
                  # set the size of the present day populations - n <i> <n> Set the size of population i to n*N0
                  "-n 1", N1, "-n 2", N2, "-n 3", N3, "-n 4", N4,
                  # set the migration rate from wave to crab (looking forward in time) inhabiting the same location
                  "-m 3 1 0 -m 4 2 0 -m 1 3", mig_WC1, "-m 2 4", mig_WC2,
                  # set the migration rate between the same ecotypes inhabiting different locations
                  "-m 2 1", mig_CC, "-m 1 2", mig_CC, "-m 4 3", mig_WW, "-m 3 4", mig_WW,
                  # set the migration, between all populations, to zero - immediately before the recent split
                  "-eM", tmRS, "0",
                  # add a split event -ej <t> <j> <i> in population i that creates population j (forwards in time)
                  "-ej", Rsplit, "1 2 -ej", Rsplit, "4 3",
                  # set the size of the ancient populations
                  "-en", Rsplit, "2", NA1, "-en", Rsplit, "3", NA2,
                  # set the migration rate between the ancestral populations
                  "-em", tmAA, "2 3", mig_AA,
                  # set the migration, between all populations, to zero - immediately at the ancient split
                  "-eM", Asplit, "0",
                  # add a split event - this event creates the two ancestral populations
                  "-ej", Asplit, "2 3",
                  # finally, set the size of the most ancestral pop equal to the size of the reference population with:
                  "-eN", Asplit, 1)

  # create command line with no migration from the wave to the crab ecotype within the same location
  no.mWC <- paste(paste(nDip*2, collapse = " "), lociTotal[3], "-t",theta, "-I", paste(nPops), paste(n*2, collapse = " "),
                  # set the size of the present day populations
                  "-n 1", N1, "-n 2", N2, "-n 3", N3, "-n 4", N4,
                  # set the migration rate from crab to wave (looking forward in time) inhabiting the same location
                  "-m 3 1", mig_CW1, "-m 4 2", mig_CW1, "-m 1 3 0 -m 2 4 0",
                  # set the migration rate between the same ecotypes inhabiting different locations
                  "-m 2 1", mig_CC, "-m 1 2", mig_CC, "-m 4 3", mig_WW, "-m 3 4", mig_WW,
                  # set the migration, between all populations, to zero - immediately before the recent split
                  "-eM", tmRS, "0",
                  # add a split event -ej <t> <j> <i> in population i that creates population j (forwards in time)
                  "-ej", Rsplit, "1 2 -ej", Rsplit, "4 3",
                  # set the size of the ancient populations
                  "-en", Rsplit, "2", NA1, "-en", Rsplit, "3", NA2,
                  # set the migration rate between the ancestral populations
                  "-em", tmAA, "3 2", mig_AA,
                  # set the migration, between all populations, to zero - immediately at the ancient split
                  "-eM", Asplit, "0",
                  # add a split event - this event creates the two ancestral populations
                  "-ej", Asplit, "2 3",
                  # finally, set the size of the most ancestral pop equal to the size of the reference population with:
                  "-eN", Asplit, 1)

  # create a command line for the loci without any migration (between the different ecotypes)
  # set the basic elements for scrm - nhap: total number of haplotypes that are simulated at each locus and
  # nrep: the number of independent loci that will be produced
  without.mig <- paste(paste(nDip*2, collapse = " "), lociTotal[4], "-t",theta, "-I", paste(nPops), paste(n*2, collapse = " "), "0",
                       # set the size of the present day populations - n <i> <n> Set the size of population i to n*N0.
                       "-n 1", N1, "-n 2", N2, "-n 3", N3, "-n 4", N4,
                       # set the migration rate between the same ecotypes inhabiting different locations
                       "-m 2 1", mig_CC, "-m 1 2", mig_CC, "-m 4 3", mig_WW, "-m 3 4", mig_WW,
                       # set the migration, between all populations, to zero - immediately before the recent split
                       "-eM", tmRS, "0",
                       # add a split event -ej <t> <j> <i> in population i that creates population j (forwards in time)
                       "-ej", Rsplit, "1 2 -ej", Rsplit, "4 3",
                       # set the size of the ancient populations
                       "-en", Rsplit, "2", NA1, "-en", Rsplit, "3", NA2,
                       # set the migration, between all populations, to zero - immediately at the ancient split
                       "-eM", Asplit, "0",
                       # add a split event - this event creates the two ancestral populations
                       "-ej", Asplit, "2 3",
                       # finally, set the size of the most ancestral pop equal to the size of the reference population with:
                       "-eN", Asplit, 1)

  # combine the two different types of commands
  cmdSingle <- c(with.mig, no.mCW, no.mWC, without.mig)

  # if extra is equal to TRUE, then we simulated more loci than required
  if(extra == TRUE)
    # include the required number of loci per category in the output
    cmdSingle <- list(commands = cmdSingle, targetLoci = targetLoci)

  # output the command line for the single origin model
  cmdSingle
}


#' Create SCRM command line for a parallel origin scenario
#'
#' This function creates a command line tailored for a scenario of parallel
#' origin to explain ecotype formation. The command line can then be fed to the
#' scrm package to run the model.
#'
#' For convenience, imagine we have two divergent ecotypes, named C and W. This
#' model assumes that the first population corresponds to the C ecotype at the
#' first location, the second population to the W ecotype in the first location,
#' the third population to the C ecotype in the second location and the fourth
#' population to the W ecotype in the second location.
#'
#' @param parameters A vector where each entry corresponds to a different
#'   parameter, e.g. one entry is the size of the reference population, another
#'   is the time of recent split, etc. Please note that this functions depends
#'   on the ordering of the parameters in the vector and thus, it should only be
#'   used with a vector created with the `createParams` function.
#' @param nSites An integer representing the number of base pairs that each
#'   locus should have.
#' @param nLoci An integer that represents how many independent loci should be
#'   simulated.
#' @param nDip An integer representing the total number of diploid individuals
#'   to simulate. Note that scrm actually simulates haplotypes, so the number of
#'   simulated haplotypes is double of this. Also note that this is the total
#'   number of diploid individuals and this function will distribute the
#'   individuals equally by the two populations.
#' @param mutrate A number representing the mutation rate assumed for the
#'   simulations.
#' @param extra is a logical value indicating whether the required number of
#'   loci should be enforced. The default is FALSE but, if set to TRUE, then
#'   additional loci will be simulated. These additional loci are simulated to
#'   try to have sufficient loci to keep the required number of loci after
#'   filtering.
#'
#' @return a character vector with four entries. The first entry is the scrm
#'   command line for the loci without any barriers against migration. The
#'   second entry is the command line for the loci without migration from the C
#'   towards the W ecotype. The third entry is command line for the loci without
#'   migration from the W towards the C ecotype and the last entry is the scrm
#'   command line for the loci without migration between divergent ecotypes.
#'
#' @examples
#' # create a vector with parameter values for the parallel origin scenario
#' params <- createParams(Nref = c(25000, 25000), ratio = c(0.1, 3), pool = c(5, 250),
#' seq = c(0.0001, 0.001), split = c(0, 3), CW = c(1e-13, 1e-3), WC = c(1e-13, 1e-3),
#' CC =  c(1e-13, 1e-3), WW = c(1e-13, 1e-3), ANC = c(1e-13, 1e-3), bT = c(0, 0.2),
#' bCW = c(0, 0.5), bWC = c(0, 0.5), model = "Parallel")
#'
#' # create the command line for the scrm package
#' cmdParallel(parameters = params, nSites = 2000, nLoci = 100, nDip = 400, mutrate = 2-8)
#'
#' @export
cmdParallel <- function(parameters, nSites, nLoci, nDip, mutrate, extra = FALSE) {

  # this function is intended to be used with a four-population model
  nPops <- 4

  # Read the vector with the parameters and assign each parameter to the correct command name
  Ne <- parameters[1]
  # set the relative size of each population - for the extant populations
  N1 <- parameters[2]; N2 <- parameters[3]; N3 <- parameters[4]; N4 <- parameters[5]
  # and the ancient populations
  NA1 <- parameters[6]; NA2 <- parameters[7]

  # get the proportion of loci with migration - no barriers against migration between the different ecotypes
  pM <- parameters[19]
  # get the proportion of loci without migration - from the crab to the wave ecotype at the same location
  pCW <- parameters[20]
  # get the proportion of loci without migration - from the wave to the crab ecotype at the same location
  pWC <- parameters[21]
  # and the proportion of loci without any migration - total barrier against migration between the different ecotypes
  pNO <- parameters[22]

  # get the migration rates - between different ecotypes at the same site - this is the m value on the M = 4N0m formula
  # from the crab ecotype to the wave ecotype - at the first site
  mCW1 <- parameters[12]
  # from the crab ecotype to the wave ecotype - at the second site
  mCW2 <- parameters[13]
  # from the wave ecotype to the crab ecotype - at the fist site
  mWC1 <- parameters[14]
  # from the wave ecotype to the crab ecotype - at the second site
  mWC2 <- parameters[15]

  # between crab populations inhabiting different locations
  mCC <- parameters[16]
  # between wave populations inhabiting different locations
  mWW <- parameters[17]
  # and between the two ancestral populations
  mAA <- parameters[18]

  # m <i> <j> <M>: Set the migration rate from population j to population i to M (looking forward in time)
  # and REMEMBER that M = 4N0m
  # set the migration rates to the scale of Nref - between different ecotypes at the same site
  # from the crab ecotype to the wave ecotype - at the first site
  mig_CW1 <- 4*Ne*mCW1 # -m 2 1 mig_CW
  # from the crab ecotype to the wave ecotype - at the second site
  mig_CW2 <- 4*Ne*mCW2 # -m 4 3 mig_CW
  # from the wave ecotype to the crab ecotype - at the first site
  mig_WC1 <- 4*Ne*mWC1 # -m 1 2 mig_WC
  # from the wave ecotype to the crab ecotype - at the second site
  mig_WC2 <- 4*Ne*mWC2 # -m 3 4 mig_WC

  # between crab populations at different locations
  mig_CC <- 4*Ne*mCC
  # between wave populations at different locations
  mig_WW <- 4*Ne*mWW
  # and between the two ancestral populations
  mig_AA <- 4*Ne*mAA

  # get the time of the recent split event
  Rsplit <- round(parameters[8], digits = 3)
  # and the ancient split event
  Asplit <- round(parameters[9], digits = 3)
  # the actual time of the ancient split event is obtained by doing:
  Asplit <- Rsplit + Asplit

  # Compute the value of theta
  mutrate_locus <- nSites*mutrate
  theta <- 4*Ne*mutrate_locus
  # Create a vector with information about how many haplotypes are sampled from each population
  n <- c(rep(nDip/nPops, times = nPops))

  # create a variable to be used for a minor correction related to the split time
  # this will ensure that the changes in migration rates occur at different times
  # a variable to inform when does migration start between ancestral populations
  tmAA <- Rsplit + 0.0001
  # if Rsplit is zero, we can not subtract something from it
  if(Rsplit != 0) {
    # when Rsplit is not zero, set the end of the migration before the split
    tmRS <- Rsplit - 0.0001
  } else {
    # when Rsplit is zero, set the end of the migration at the split
    tmRS <- Rsplit
  }

  # use a multinomial distribution to get the number of loci simulated under each category
  lociTotal <- as.vector(stats::rmultinom(n = 1, size = nLoci, c(pM, pCW, pWC, pNO)))

  # if extra is TRUE, then more loci than required per category will be simulated
  if(extra == TRUE) {
    # save the required number of simulated loci per category
    targetLoci <- lociTotal
    # simulate more 25 loci per category
    lociTotal <- lociTotal + 25
  }

  # cheat code: pop1 - crab in site 1; pop2 - wave in site 1; pop3 - crab in site 2; pop4 - wave in site 2

  # create command line with no barriers to migration
  # set the basic elements for scrm - nhap: total number of haplotypes that are simulated at each locus and
  # nrep: the number of independent loci that will be produced
  with.mig <- paste(paste(nDip*2, collapse = " "), lociTotal[1], "-t",theta, "-I", paste(nPops), paste(n*2, collapse = " "),
                    # set the size of the present day populations - n <i> <n> Set the size of population i to n*N0.
                    "-n 1", N1, "-n 2", N2, "-n 3", N3, "-n 4", N4,
                    # m <i> <j> <M>: Set the migration rate from population j to population i to M (looking forward in time)
                    # set the migration rate between different ecotypes inhabiting the same location
                    "-m 2 1", mig_CW1, "-m 4 3", mig_CW2, "-m 1 2", mig_WC1, "-m 3 4", mig_WC2,
                    # set the migration rate between the same ecotypes inhabiting different locations
                    "-m 3 1", mig_CC, "-m 1 3", mig_CC, "-m 4 2", mig_WW, "-m 2 4", mig_WW,
                    # set the migration, between all populations, to zero - immediately before the recent split
                    "-eM", tmRS, "0",
                    # add a split event -ej <t> <j> <i> in population i that creates population j (forwards in time)
                    "-ej", Rsplit, "1 2 -ej", Rsplit, "4 3",
                    # set the size of the ancient populations
                    # -en <t> <i> <n>: Set the size of population i to n*N0 at time t.
                    "-en", Rsplit, "2", NA1, "-en", Rsplit, "3", NA2,
                    # set the migration rate between the ancestral populations
                    "-em", tmAA, "2 3", mig_AA, "-em", tmAA, "3 2", mig_AA,
                    # set the migration, between all populations, to zero - immediately at the ancient split
                    "-eM", Asplit, "0",
                    # add a split event - this event creates the two ancestral populations
                    "-ej", Asplit, "2 3",
                    # finally, set the size of the most ancestral pop equal to the size of the reference population with:
                    # -eN <t> <n>: set the size of all populations to n*N0 at time t.
                    "-eN", Asplit, 1)

  # create command line with no migration from the crab to the wave ecotype within the same location
  no.mCW <- paste(paste(nDip*2, collapse = " "), lociTotal[2], "-t",theta, "-I", paste(nPops), paste(n*2, collapse = " "),
                  # set the size of the present day populations - n <i> <n> Set the size of population i to n*N0
                  "-n 1", N1, "-n 2", N2, "-n 3", N3, "-n 4", N4,
                  # set the migration rate from wave to crab (looking forward in time) inhabiting the same location
                  "-m 2 1 0 -m 4 3 0 -m 1 2", mig_WC1, "-m 3 4", mig_WC2,
                  # set the migration rate between the same ecotypes inhabiting different locations
                  "-m 3 1", mig_CC, "-m 1 3", mig_CC, "-m 4 2", mig_WW, "-m 2 4", mig_WW,
                  # set the migration, between all populations, to zero - immediately before the recent split
                  "-eM", tmRS, "0",
                  # add a split event -ej <t> <j> <i> in population i that creates population j (forwards in time)
                  "-ej", Rsplit, "1 2 -ej", Rsplit, "4 3",
                  # set the size of the ancient populations
                  "-en", Rsplit, "2", NA1, "-en", Rsplit, "3", NA2,
                  # set the migration rate between the ancestral populations
                  "-em", tmAA, "2 3", mig_AA, "-em", tmAA, "3 2", mig_AA,
                  # set the migration, between all populations, to zero - immediately at the ancient split
                  "-eM", Asplit, "0",
                  # add a split event - this event creates the two ancestral populations
                  "-ej", Asplit, "2 3",
                  # finally, set the size of the most ancestral pop equal to the size of the reference population with:
                  "-eN", Asplit, 1)

  # create command line with no migration from the wave to the crab ecotype within the same location
  no.mWC <- paste(paste(nDip*2, collapse = " "), lociTotal[3], "-t",theta, "-I", paste(nPops), paste(n*2, collapse = " "),
                  # set the size of the present day populations
                  "-n 1", N1, "-n 2", N2, "-n 3", N3, "-n 4", N4,
                  # set the migration rate from crab to wave (looking forward in time) inhabiting the same location
                  "-m 2 1", mig_CW1, "-m 4 3", mig_CW2, "-m 1 2 0 -m 3 4 0",
                  # set the migration rate between the same ecotypes inhabiting different locations
                  "-m 3 1", mig_CC, "-m 1 3", mig_CC, "-m 4 2", mig_WW, "-m 2 4", mig_WW,
                  # set the migration, between all populations, to zero - immediately before the recent split
                  "-eM", tmRS, "0",
                  # add a split event -ej <t> <j> <i> in population i that creates population j (forwards in time)
                  "-ej", Rsplit, "1 2 -ej", Rsplit, "4 3",
                  # set the size of the ancient populations
                  "-en", Rsplit, "2", NA1, "-en", Rsplit, "3", NA2,
                  # set the migration rate between the ancestral populations
                  "-em", tmAA, "2 3", mig_AA, "-em", tmAA, "3 2", mig_AA,
                  # set the migration, between all populations, to zero - immediately at the ancient split
                  "-eM", Asplit, "0",
                  # add a split event - this event creates the two ancestral populations
                  "-ej", Asplit, "2 3",
                  # finally, set the size of the most ancestral pop equal to the size of the reference population with:
                  "-eN", Asplit, 1)

  # create a command line for the loci without any migration (between the different ecotypes)
  without.mig <- paste(paste(nDip*2, collapse = " "), lociTotal[4], "-t",theta, "-I", paste(nPops), paste(n*2, collapse = " "), "0",
                       # set the size of the present day populations - n <i> <n> Set the size of population i to n*N0.
                       "-n 1", N1, "-n 2", N2, "-n 3", N3, "-n 4", N4,
                       # set the migration rate between the same ecotypes inhabiting different locations
                       "-m 3 1", mig_CC, "-m 1 3", mig_CC, "-m 4 2", mig_WW, "-m 2 4", mig_WW,
                       # set the migration, between all populations, to zero - immediately before the recent split
                       "-eM", tmRS, "0",
                       # add a split event -ej <t> <j> <i> in population i that creates population j (forwards in time)
                       "-ej", Rsplit, "1 2 -ej", Rsplit, "4 3",
                       # set the size of the ancient populations
                       "-en", Rsplit, "2", NA1, "-en", Rsplit, "3", NA2,
                       # set the migration rate between the ancestral populations
                       "-em", tmAA, "2 3", mig_AA, "-em", tmAA, "3 2", mig_AA,
                       # set the migration, between all populations, to zero - immediately at the ancient split
                       "-eM", Asplit, "0",
                       # add a split event - this event creates the two ancestral populations
                       "-ej", Asplit, "2 3",
                       # finally, set the size of the most ancestral pop equal to the size of the reference population with:
                       "-eN", Asplit, 1)

  # combine the four different types of commands
  cmdParallel <- c(with.mig, no.mCW, no.mWC, without.mig)

  # if extra is equal to TRUE, then we simulated more loci than required
  if(extra == TRUE)
    # include the required number of loci per category in the output
    cmdParallel <- list(commands = cmdParallel, targetLoci = targetLoci)

  # output the command line for the parallel origin model
  cmdParallel
}


#' Create invariable sites on the scrm output
#'
#' This function applies a correction for the situations where scrm does not
#' produce a single polymorphic site for a given locus. In this situation, two
#' artificial sites are created at that locus. All individuals are assumed to be
#' homozygous for the reference allele at those sites.
#'
#' @param haplotypes a list of haplotypes obtained from the simulations done
#'   with scrm. Each entry of the list is a matrix that corresponds to a given
#'   locus. At each matrix, each column is a different site and each row is a
#'   different haplotype.
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


#' Organize scrm output
#'
#' This function is utilized to sort out the scrm output. The order of the
#' populations changes accordingly to the model used (i.e. single or parallel
#' origin). Running this function will re-organize the output produced by scrm,
#' so that the populations are in the same order in both models.
#'
#' @param seg_sites a matrix of segregating sites as produced by scrm. Each
#'   column of the matrix is a different site and each row is a different
#'   haplotype.
#' @param nHap an integer representing the total number of haplotypes simulated.
#' @param nPops an integer, representing the total number of populations of the
#'   simulated model.
#'
#' @return a matrix of segregating sites, similar to `seg_sites` but with the
#'   populations organized so that the order is always the same, regardless of
#'   the model used.
#'
#' @keywords internal
#'
#' @export
organizeSCRM <- function(seg_sites, nHap, nPops) {

  # get the number of haplotypes simulated by population
  haPop <- nHap/nPops
  # create a vector with the index representing the beginning of each population
  beginPop <- seq(from = 1, to = nHap, by = haPop)

  # remove the name (position) of each site - this is something that scrm creates
  seg_sites <- unname(seg_sites)

  # in the single model, we need to switch the order of the second and third population
  # get the haplotypes corresponding to each population
  pop2 <- seg_sites[(beginPop[2]):(beginPop[3]-1), ]
  pop3 <- seg_sites[(beginPop[3]):(beginPop[4]-1), ]

  # re-organize the matrix of haplotypes with the populations in the correct order
  seg_sites[(beginPop[2]):(beginPop[3]-1), ] <- pop3
  seg_sites[(beginPop[3]):(beginPop[4]-1), ] <- pop2

  # output the matrix of haplotypes with the populations in the correct order
  seg_sites
}


#' Convert Haplotypes to Genotypes
#'
#' This function converts haplotypes simulated with scrm into genotypes by
#' adding the entries on one row with the entries of the subsequent row.
#'
#' @param haplo a matrix of haplotypes obtained from the simulations done with
#'   scrm. Each column of the matrix is a different site and each row is a
#'   different haplotype.
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


#' Create Genotypes from the scrm output
#'
#' This function applies the \code{\link{hap2geno}} function to all entries of a
#' list. Each entry of that list is a different locus simulated with scrm. Thus,
#' this function converts the haplotypes of all simulated loci into genotypes.
#'
#' @param haplotypes a list of haplotypes obtained from the simulations done
#'   with scrm. Each entry of the list is a matrix that corresponds to a given
#'   locus. At each matrix, each column is a different site and each row is a
#'   different haplotype.
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


#' Run scrm and obtain genotypes
#'
#' This function will run the scrm package, according to the command line
#' supplied as input. It will also combine haplotypes into genotypes and
#' re-organize the output if the simulations were performed under a single
#' origin scenario. This is to ensure that the output of the four-population
#' models will always follow the same order: the two divergent ecotypes in the
#' first location, followed by the two divergent ecotypes in the second
#' location.
#'
#' @param commands A character string containing the commands for the scrm
#'   package. This string can be created using the `cmd2pops`, the `cmdSingle`
#'   or the `cmdParallel` functions.
#' @param nDip An integer representing the total number of diploid individuals
#'   to simulate. Note that scrm actually simulates haplotypes, so the number of
#'   simulated haplotypes is double of this.
#' @param nPops An integer that informs of how many populations exist on the
#'   model you are trying to run.
#' @param model Either "2pops", "Single" or "Parallel" indicating which model
#'   should be simulated.
#'
#' @return a list with the simulated genotypes. Each entry is a different locus
#'   and, for each locus, different rows represent different individuals and
#'   each column is a different site.
#'
#' @examples
#' # create a vector with parameter values for a two populations model
#' params <- createParams(Nref = c(25000, 25000), ratio = c(0.1, 3), pool = c(5, 250),
#' seq = c(0.0001, 0.001), split = c(0, 3), CW = c(1e-13, 1e-3), WC = c(1e-13, 1e-3),
#' bT = c(0, 0.2), model = "2pops")
#'
#' # create the command line for the scrm package
#' cmds <- cmd2pops(parameters = params, nSites = 2000, nLoci = 100, nDip = 100, mutrate = 2e-8)
#'
#' # run scrm and obtain genotypes
#' runSCRM(commands = cmds, nDip = 100, nPops = 2, model = "2pops")
#'
#' @export
runSCRM <- function(commands, nDip, nPops, model) {

  # check if the input is correct - the selected model should be one of the following
  if(model %in% c("2pops", "Single", "Parallel") == FALSE)
    stop(paste("The selected model should be either 2pops, Single or Parallel. Please check"))

  # binding the variable locally to the function
  temp1 <- NULL

  # run the scrm package for each set of commands - with and without migration
  simulation <- lapply(commands, FUN = function(x) scrm::scrm(x))

  # extract the information from each simulation and store it on a temporary matrix
  for (i in 1:length(simulation)) {

    # create the temporary matrix for each simulation
    assign(paste("temp", i, sep = ""), simulation[[i]][["seg_sites"]])
  }

  if(length(simulation) != 1) {

    # combine all simulations into a matrix of haplotypes
    haplotypes <- append(temp1, unlist(mget(paste0("temp", 2:length(simulation))), recursive = FALSE, use.names = FALSE))

  } else {

    # if only set of simulations was performed, only one set of haplotypes exist
    haplotypes <- temp1
  }

  # get the total number of haplotypes
  nHap <- nDip*2

  # apply a correction for the situations where scrm does not produce a single polymorphic site
  # first check the dimensions of each list entry
  size <- matrix(unlist(lapply(haplotypes, dim)), ncol = 2, byrow = TRUE)

  # if one entry has no columns, i.e. no sites, then add columns containing only zeros to that entry
  if(any(size[, 2] == 0)) {

    # add two columns containing zeros to that locus
    haplotypes <- haplo.fix(haplotypes = haplotypes, nHap = nHap)
  }

  # re-organize output for the single model
  if (model == "Single") {

    # change the order of the populations in the single origin model
    # so that the order is always: ecotype C and W in the first location and ecotype C and W in the second location
    haplotypes <- lapply(haplotypes, function(segSites) organizeSCRM(segSites, nHap, nPops))
  }

  # convert the haplotypes to genotypes
  genotypes <- GetGenotypes(haplotypes, nDip = nDip)

  # output the genotypes
  genotypes
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
#' @param nSNPs An integer representing the number of polymorphic sites per
#'   locus to simulate. This is an optional input but either this or the
#'   `genotypes` list must be supplied.
#' @param nLoci An optional integer that represents how many independent loci
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
#' # create a vector with parameter values for a two populations model
#' params <- createParams(Nref = c(25000, 25000), ratio = c(0.1, 3), pool = c(5, 250),
#' seq = c(0.0001, 0.001), split = c(0, 3), CW = c(1e-13, 1e-3), WC = c(1e-13, 1e-3),
#' bT = c(0, 0.2), model = "2pops")
#' # create the command line for the scrm package
#' cmds <- cmd2pops(parameters = params, nSites = 2000, nLoci = 100, nDip = 100, mutrate = 2e-8)
#' # run scrm and obtain genotypes
#' genotypes <- runSCRM(commands = cmds, nDip = 100, nPops = 2, model = "2pops")
#' # simulate coverage
#' simulateCoverage(mean = c(100, 50), variance = c(250, 150), genotypes = genotypes)
#'
#' @export
simulateCoverage <- function(mean, variance, nSNPs = NA, nLoci = NA, genotypes = NA) {

  # check if either the number of SNPs or the genotypes were supplied as input
  if(all(is.na(nSNPs), is.na(genotypes)))
    stop("You should define the number of SNPs to simulate or supply a list of genotypes. Please check")

  # check if the variance and mean are reasonable
  if(any(variance - mean > 0) == FALSE)
    stop("Error: variance equal to mean, or variance smaller than mean.")

  # calculate the parameters for the negative binomial
  pnb <- mean/variance
  rnb <- (mean^2)/(variance - mean)

  # if the genotypes are supplied as input to the function
  if(!all(is.na(genotypes))) {

    # check if the input is correct - genotypes should always be supplied as a list
    if(any(class(genotypes) != "list"))
      stop(paste("Genotypes should be supplied on a list format, with each entry corresponding to a locus. Please check"))

    # use a negative binomial to draw random values, per site and per population, for the total number of observed reads
    # this outputs a list where each entry corresponds to a locus
    readnumbers <- lapply(genotypes, FUN = function(geno)
      t(mapply(FUN = function(size, prob) stats::rnbinom(n = ncol(geno), size = size, prob = prob), rnb, pnb)))

  } else {

    # check if the input is correct - when genotypes are not supplied as input, the number of loci should be defined
    if(is.na(nLoci))
      stop(paste("Please define the number of loci to simulate"))

    # use a negative binomial to draw random values, per site and per population, for the total number of observed reads
    # this outputs a list where each entry corresponds to a locus
    readnumbers <- lapply(1:nLoci, FUN = function(geno)
      t(mapply(FUN = function(size, prob) stats::rnbinom(n = nSNPs, size = size, prob = prob), rnb, pnb)))
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
  if(class(reads) != "list")
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
#' of numbers from the Dirichlet distribution produces NaN. Thus, this function
#' replaces small values of alpha with a minimum threshold value to avoid that.
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
  if(class(vector_np) != "numeric" | length(vector_np) == 1)
    stop(paste("The vector_np input should be a vector. It should also contain more than one entry. Please check"))

  # check if we are dealing with a single population - this function should be used on a single population
  # also check if the input is on the correct format
  if(nPools != length(vector_np))
    stop(paste("The nPools input is", paste(nPools), "and so, the length of the size vector should also be",
               paste(nPools, ".", sep = ""), "Please check"))

  # Silence warnings - when using only two pools, the rdirichlet functions prints a warning about reducing
  # to a beta function. This ensures that the warning is not printed
  options(warn = -1)

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

  # set the warning level back to normal
  options(warn = 0)

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

  # This package contains a function to perform random draws from a Dirichlet distribution

  # Silence warnings - when using only two individuals, the rdirichlet functions prints a warning about reducing
  # to a beta function. This ensures that the warning is not printed
  options(warn = -1)

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

  # set the warning level back to normal
  options(warn = 0)

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
  if(class(coverage) == "list")
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
  if(class(size) != "list")
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

  if(class(genotypes) != "list")
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
#' alternative) corresponds to the major allele. This function can also be used
#' to remove sites according to a minor-allele reads threshold.
#'
#' More precisely, this function counts the number of reads with the reference
#' or alternative allele at each site and then sets the major allele as the most
#' frequent of the two. This is done across all populations and so the major and
#' minor alleles are defined at a global level. Then if the `min.minor` input is
#' not NA, sites where the number of minor allele reads, across all populations,
#' are below the user-defined threshold are removed.
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
#' # define the major and minor alleles for this ool-seq data
#' # we have to select the first entry of the pools list because this function works for matrices
#' minorPool(reference = pools$reference[[1]], alternative = pools$alternative[[1]],
#' coverage = pools$total[[1]])
#'
#' @export
minorPool <- function(reference, alternative, coverage, min.minor = NA) {

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

    # if the min.minor input is not NA - then it should be an integer
    # representing the minimum number of reads with the minor allele that we should observe across all populations
    if (is.na(min.minor) == FALSE) {
      # now we need to find the total coverage - across all populations - of the minor frequency allele
      # since we switched the columns in the previous step, the number of reads with the minor allele - the less frequent allele
      # are stored in the alternative matrix
      minor <- colSums(alternative)
      # find out in which columns the total sum of the reads with the minor allele is below the threshold
      toremove <- minor < min.minor

      # if there are sites where the sum of the reads with the minor allele is below the threshold
      if (length(toremove) != 0) {
        # remove those columns from the matrix containing the depth of coverage
        coverage <- coverage[, !toremove, drop = FALSE]
        # remove those columns from the matrix containing the number of reads with the reference allele
        reference <- reference[, !toremove, drop = FALSE]
        # remove those columns from the matrix containing the number of reads with the alternative allele
        alternative <- alternative[, !toremove, drop = FALSE]
      }
    }

    # output the matrices with the various types of read numbers
    out <- list(major = reference, minor = alternative, total = coverage)
  }

  # output the matrices with the various types of read numbers
  out
}


#' Calculate population frequency at each SNP
#'
#' The frequency at a given SNP is calculated according to: π = c/r, where c =
#' number of minor-allele reads and r = total number of observed reads.
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
#' @return a list with two names entries
#'
#'   \item{pi}{a list with the allele frequencies of each population. Each list
#'   entry is a matrix, corresponding to a different locus. Each row of a matrix
#'   corresponds to a different population and each column to a different site.}
#'
#'   \item{pool}{a list with three different entries: major, minor and total.
#'   This list is similar to the one obtained with the \code{\link{minorPool}}
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
#' pools <- minorPool(reference = pools$reference[[1]], alternative = pools$alternative[[1]],
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
  if(any(lapply(listPool, class) == "matrix") == TRUE)
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


#' Force the simulations to contain the required number of loci
#'
#' This function attempts to force the required number of loci after the
#' filtering steps are performed.
#'
#' This is done by simulating extra loci for each of the different types of
#' simulations performed. The possible types of simulations include loci without
#' barriers against migration between divergent ecotypes, loci without migration
#' from the C towards the W ecotype, loci without migration from the W towards
#' the C ecotypes and loci where no migration occurs between divergent ecotypes.
#' Using this function, more loci than required are simulated for each of those
#' types of simulations.
#'
#' Then, a coverage-based filter is applied to the data, followed by a filter
#' based on a required number of minor-allele reads per site. Those filters
#' remove some loci from the data. The extra simulated loci should allow us to
#' keep the required number of loci per type of simulation even after filtering.
#'
#'
#' @param model a character, either 2pops", "Single" or "Parallel" indicating
#'   which model should be simulated.
#' @param parameters a vector of parameters used to create the command line for
#'   the scrm package. Each entry of the vector is a different parameter. Note
#'   that each vector entry should be named with the name of the corresponding
#'   parameter. The output of the `CreateParameters` function is the intended
#'   input.
#' @param nSites is an integer that specifies how many base pairs should scrm
#'   simulate, i.e. how many sites per locus to simulate.
#' @param nLoci an integer that represents how many independent loci should be
#'   simulated.
#' @param nDip an integer representing the total number of diploid individuals
#'   to simulate. Note that scrm actually simulates haplotypes, so the number of
#'   simulated haplotypes is double of this. Also note that this is the total
#'   number of diploid individuals and this function will distribute the
#'   individuals equally by the simulated populations.
#' @param mutrate an integer representing the mutation rate assumed for the
#'   simulations.
#' @param mean an integer or a vector defining the mean value of the negative
#'   binomial distribution from which different number of reads are drawn. It
#'   represents the mean coverage across all sites. If a vector is supplied, the
#'   function assumes that each entry of the vector is the mean for a different
#'   population.
#' @param variance an integer or a vector defining the variance of the negative
#'   binomial distribution from which different number of reads are drawn. It
#'   represents the variance of the total coverage across all sites. If a vector
#'   is supplied, the function assumes that each entry of the vector is the
#'   variance for a different population.
#' @param minimum an integer representing the minimum coverage allowed. Sites
#'   where any population has a depth of coverage below this threshold are
#'   removed from the data.
#' @param maximum an integer representing the maximum coverage allowed. Sites
#'   where any population has a depth of coverage above this threshold are
#'   removed from the data.
#' @param size a list with one entry per population. Each entry should be a
#'   vector containing the size (in number of diploid individuals) of each pool.
#'   Thus, if a population was sequenced using a single pool, the vector should
#'   contain only one entry. If a population was sequenced using two pools, each
#'   with 10 individuals, this vector should contain two entries and both will
#'   be 10.
#' @param min.minor is an integer representing the minimum allowed number of
#'   minor-allele reads. Sites that, across all populations, have less
#'   minor-allele reads than this threshold will be removed from the data.
#'
#' @return a list with two names entries
#'
#'   \item{pool}{a list with three different entries: major, minor and total.
#'   This list is obtained by running the \code{\link{forcePool}} function.}
#'
#'   \item{nPoly}{a numeric value indicating the mean number of polymorphic
#'   sites across all simulated locus.}
#'
#' @examples
#' # create a vector of parameters for a model with two populations
#' parameters <- createParams(Nref = c(25000, 25000), ratio = c(0.1, 3), pool = c(5, 250),
#' seq = c(0.0001, 0.001), split = c(0, 3), CW = c(1e-13, 1e-3), WC = c(1e-13, 1e-3), bT = c(0, 0.2),
#' model = "2pops")
#'
#' # simulate a two populations model, forcing the number of loci (100 loci)
#' forceLocus(model = "2pops", parameters, nSites = 2000, nLoci = 100, nDip = 200, mutrate = 2e-8,
#' mean = c(100, 80), variance = c(200, 180), minimum = 10, maximum = 150,
#' size = rep(list(rep(50, 2)), 2), min.minor = 1)
#'
#' @export
forceLocus <- function(model, parameters, nSites, nLoci, nDip, mutrate, mean, variance, minimum, maximum, size, min.minor) {

  # create the command line to run the scrm package
  # the command line varies according to the selected model
  if(model == "2pops") {

    # create the command line for a star shaped model
    commands <- cmd2pops(parameters, nSites, nLoci, nDip, mutrate, extra = TRUE)
    # set the number of populations
    nPops <- 2

  } else if (model == "Parallel") {

    # create the command line for the parallel origin model
    commands <- cmdParallel(parameters, nSites, nLoci, nDip, mutrate, extra = TRUE)
    # set the number of populations
    nPops <- 4

  } else if (model == "Single") {

    # create the command line for the single origin model
    commands <- cmdSingle(parameters, nSites, nLoci, nDip, mutrate, extra = TRUE)
    # set the number of populations
    nPops <- 4

  } else {

    # if a correct model is not supplied as input for the function - stop and warn
    stop(paste("model should be 2pops, Parallel or Single. Please check!"))
  }

  # get the required number of loci per category - this is the number of loci we want in the end
  target <- commands[["targetLoci"]]
  # get the command line for scrm
  commands <- commands[["commands"]]

  # the length of the commands object indicates how many different categories of simulations we are performing
  # the possible categories are: loci with no barriers to migration, loci with no migration from one ecotype to the other
  # and loci without any migration between the different ecotypes
  nSims <- length(commands)

  # run the scrm package and obtain the genotypes
  genotypes <- lapply(commands, function(sim) runSCRM(commands = sim, nDip, nPops, model))

  # get the mean number of polymorphic sites
  nPoly <- mean(unlist(sapply(genotypes, function(sim) sapply(sim, ncol))))

  # simulate total number of reads per site
  reads <- lapply(genotypes, function(sim) simulateCoverage(mean, variance, genotypes = sim))

  # remove sites with a depth of coverage above or below the defined threshold
  reads <- lapply(1:nSims, function(sim)
    remove_by_reads(nLoci = length(reads[[sim]]), reads[[sim]], minimum, maximum, genotypes = genotypes[[sim]]))

  # get the genotypes - without sites simulated with a coverage below or above the threshold
  genotypes <- lapply(reads, FUN = function(sim) lapply(sim, "[[", 2))
  # get the reads - without sites simulated with a coverage below or above the threshold
  reads <- lapply(reads, FUN = function(sim) lapply(sim, "[[", 1))

  # check the dimensions of the matrices with the genotypes
  # it is possible that some loci do not have any polymorphic site after the coverage filter
  dimensions <- lapply(genotypes, function(sim) sapply(sim, ncol))
  # keep only those loci where we have at least one polymorphic site
  tokeep <- lapply(dimensions, function(sim) sim != 0)
  # remove loci without polymorphic sites from the genotypes
  genotypes <- lapply(1:nSims, function(sim) genotypes[[sim]][tokeep[[sim]]])
  # remove loci without polymorphic sites from the matrices with the coverage
  reads <- lapply(1:nSims, function(sim) reads[[sim]][tokeep[[sim]]])

  # simulate individual contribution to the total number of reads
  indContribution <- lapply(1:nSims, FUN = function(sim) lapply(reads[[sim]], function(locus)
    popsReads(list_np = size, coverage = locus, pError = parameters["PoolError"])))

  # simulate the number of reference reads
  reference <- lapply(1:nSims, FUN = function(sim) lapply(1:length(genotypes[[sim]]), function(locus)
    numberReferencePop(genotypes = genotypes[[sim]][[locus]], indContribution = indContribution[[sim]][[locus]],
                       size = size, error = parameters["SeqError"])))

  # simulate pooled sequencing data
  pool <- lapply(1:nSims, function(sim)
    poolPops(nPops, nLoci=length(indContribution[[sim]]), indContribution=indContribution[[sim]], readsReference=reference[[sim]]))

  # define major and minor alleles
  pool <- lapply(1:nSims, function(sim) lapply(1:length(pool[[sim]][["total"]]), function(locus)
    minorPool(reference = pool[[sim]][["reference"]][[locus]], alternative = pool[[sim]][["alternative"]][[locus]],
              coverage = pool[[sim]][["total"]][[locus]], min.minor)))

  # convert the pool list back to the previous format:
  # one entry for major allele, one for minor allele and a final one for total coverage
  pool <- lapply(1:nSims, function(sim)
    list(major = lapply(pool[[sim]], function(locus) locus[["major"]]), minor = lapply(pool[[sim]], function(locus) locus[["minor"]]),
         total = lapply(pool[[sim]], function(locus) locus[["total"]])))

  # remove loci without polymorphic sites and randomly select the required number of loci per category
  pool <- forcePool(nSims, pool = pool, target)

  # output the pooled sequencing data and the number of polymorphic sites prior to filtering
  list(pool = pool, nPoly = nPoly)
}

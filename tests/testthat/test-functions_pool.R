test_that("number of simulated loci is correct", {
  expect_length(run_scrm(nDip = 100, nloci = 10), n = 10)
})

test_that("individual frequencies are correct", {
  expect_equal(Ifreqs(nDip = 10, genotypes = list(matrix(data = 1, nrow = 10, ncol = 5))), list(rep(0.5, 5)))
})

test_that("an error occurs when freqs is not a numeric vector", {
  expect_error(removeSites(freqs = list(runif(4)), alternative = matrix(rep(10, 4), nrow = 1),
                           coverage = matrix(rep(100, 4), nrow = 1), minor = matrix(rep(10, 4), nrow = 1), min.minor = 5))
})

test_that("column names are correct", {
  expect_equal(object = colnames(maePool(nDip = 10, nloci = 5, pools = list(c(5, 5)), pError = 100, sError = 0.01,
                                         mCov = 100, vCov = 250, min.minor = 2, minimum = 10, maximum = 200)),
               expected = c("nDip", "PoolError", "nPools", "indsPool", "mean", "var", "absError"))
})

test_that("number of loci is correct", {
  expect_equal(object = nrow(maePool(nDip = 10, nloci = 5, pools = list(c(5, 5)), pError = 100, sError = 0.01,
                                         mCov = 100, vCov = 250, min.minor = 2, minimum = 10, maximum = 200)),
               expected = 5)
})

test_that("maximum coverage should be defined", {
  expect_error(maePool(nDip = 10, nloci = 5, pools = list(c(5, 5)), pError = 100, sError = 0.01, mCov = 100, vCov = 250,
                       min.minor = 2, minimum = 10))
})

test_that("number of combinations is correct", {
  expect_equal(object = nrow(maeFreqs(nDip = 50, nloci = 10, pError = c(5, 100), sError = 0.01, mCov = c(50, 100),
                                      vCov = c(200, 500), min.minor = 1)),
               expected = nrow(expand.grid(50, c(5, 100), c(50, 100), stringsAsFactors = FALSE))*10)
})

test_that("number of loci is correct", {
  expect_equal(object = nrow(errorHet(nDip = 10, nloci = 5, pools = list(c(5, 5)), pError = 100, sError = 0.01,
                                      mCov = 100, vCov = 250, min.minor = 2, minimum = 10, maximum = 200)),
               expected = 5)
})

test_that("number of combinations is correct", {
  expect_equal(object = nrow(maeHet(nDip = 50, nloci = 10, pError = c(5, 100), sError = 0.01, mCov = c(50, 100),
                                    vCov = c(200, 500), min.minor = 1)),
               expected = nrow(expand.grid(50, c(5, 100), c(50, 100), stringsAsFactors = FALSE))*10)
})

test_that("no loci are removed", {
  expect_length(object = haplo.fix(list(matrix(sample(0:1, size = 20, replace = T), nrow = 5),
                                        matrix(sample(0:1, size = 20, replace = T), nrow = 5),
                                        matrix(nrow = 5, ncol = 0)), nHap = 5), n = 3)
})

test_that("variance is larger than mean", {
  expect_error(simReads(mean = 20, variance = 5, nSNPs = 10))
})

test_that("number of simulated SNPs is correct", {
  expect_length(object = unlist(simulateCoverage(mean = 50, variance = 100, nSNPs = 10, nLoci = 1)), n = 10)
})

test_that("check that number of loci is required", {
  expect_error(simulateCoverage(mean = 50, variance = 100, nSNPs = 10))
})

test_that("sites are correctly removed from a matrix", {
  expect_equal(object = remove_by_reads_matrix(matrix(c(100, 10, 110, 120, 125, 130, 140, 150, 180, 250), nrow = 2),
                                               minimum = 20, maximum = 200),
               expected = matrix(c(110, 120, 125, 130, 140, 150), nrow = 2))
})

test_that("alpha is corrected", {
  expect_equal(object = set_alpha(alpha_i = 0.00001), expected = 0.01)
})

test_that("pool probablities sum to one", {
  expect_equal(object = colSums(poolProbs(nPools = 5, vector_np = rep(10, 5), nSNPs = 5, pError = 50)),
               expected = rep(x = 1, 5))
})

test_that("individual probablities sum to one", {
  expect_equal(object = colSums(indProbs(np = 10, nSNPs = 5, pError = 100)), expected = rep(x = 1, 5))
})

test_that("input must be correctly named", {
  expect_error(object = calculatePi(listPool = list(100, 15, 115), nLoci = 1))
})

test_that("output must be named", {
  expect_named(object = calculatePi(listPool = list(major = matrix(100, 2, 5), minor = matrix(15, 2, 5),
                                                    total = matrix(115, 2, 5)), nLoci = 1))
})

test_that("number of retained loci is correct", {

  # create pool for one type of simulation
  pool1 <- list(major = rep(list(matrix(100, 2, 5)), 10), minor = rep(list(matrix(15, 2, 5)), 10),
                total = rep(list(matrix(115, 2, 5)), 10))

  # create pool for another type of simulation
  pool2 <- list(major = rep(list(matrix(100, 2, 5)), 10), minor = rep(list(matrix(15, 2, 5)), 10),
                total = rep(list(matrix(115, 2, 5)), 10))

  # combine both into a single list
  pool <- list(pool1, pool2)

  expect_length(object = forcePool(nSims = 2, pool = pool, target = c(2, 2))[["total"]], n = 4)
})

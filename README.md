# poolHelper
Simulates Pooled Sequencing Genetic Data

This is the source code for the poolHelper R package. The main goal of this package is to provide tools to help researchers design their pooled sequencing studies.

The package provides functions to simulate pooled sequencing (Pool-seq) data under a variety of conditions. Users can define the average coverage, the number of individuals in the pool, the number of pools used and the Pool-seq error. The poolHelper package simulates the allele frequencies obtained with Pool-seq under different combinations of those parameters. These allele frequencies are then compared with the allele frequencies computed directly from the genotypes in the sample and the average absolute difference between both sets of allele frequencies is calculated. 

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7520303.svg)](https://doi.org/10.5281/zenodo.7520303)

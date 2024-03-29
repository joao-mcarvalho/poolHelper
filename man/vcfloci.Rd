% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions_pool.R
\name{vcfloci}
\alias{vcfloci}
\title{Create vcf string for all SNPs in multiple loci}
\usage{
vcfloci(reference, alternative, total)
}
\arguments{
\item{reference}{is a list where each entry corresponds to a different locus.
Each list entry is a vector with the number of reads with the reference
allele. Each entry of the vector corresponds to a different SNP.}

\item{alternative}{is a list where each entry corresponds to a different
locus. Each list entry is a vector with the number of reads with the
alternative allele. Each entry of the vector corresponds to a different
SNP.}

\item{total}{is a list where each entry corresponds to a different locus.
Each list entry is a vector with the total number of reads observed at each
SNP. Each entry of the vector corresponds to a different SNP.}
}
\value{
is a list where each entry corresponds to a different locus. Each
entry of the list is a character vector with as many entries as the number
of SNPs in the locus. Each entry of this character vector contains the
information for a single SNP coded as R,A:DP.
}
\description{
Creates a string with the information for all SNPs across multiple loci. The
information is coded as R,A:DP. R is the number of reads of the reference
allele, A is the number of reads of the alternative allele and DP is the
total depth of coverage. Each entry of the character string corresponds to a
different SNP and each entry of the list to a different locus.
}
\keyword{internal}

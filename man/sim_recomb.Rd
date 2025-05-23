% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/main.R
\name{sim_recomb}
\alias{sim_recomb}
\title{Simulate segments in a recombining genotype}
\usage{
sim_recomb(
  k,
  rho = 7.4e-07,
  set_props,
  chrom_sizes,
  zero_based_positioning = T
)
}
\arguments{
\item{k}{the number of serial meiosis events.}

\item{rho}{the recombination rate (per-site, per-meiosis). By default uses
the value 7.4e-7 from Miles et al. (2016).}

\item{set_props}{proportions of each of the members of the sample set.
Dictates the probability they are chosen as an ancestor.}

\item{chrom_sizes}{lengths of each chromosome, default taken from `get_pf3d7_chrom_sizes()`
by default.}
}
\description{
Produces a tibble of segments along the genome produced by
  meiosis of several distinct ancestors.
}
\details{
`k` sets the number of serial meiosis events. The number of
  ancestors chosen from the sample set is equal to 2^k. Although the same
  ancestor may be chosen multiple times by chance, it is not allowed for all
  2^k parents to come from the same ancestor, i.e. there is an element of
  sampling without replacement to avoid this.
}
\references{
Miles A, Iqbal Z, Vauterin P, Pearson R, Campino S, Theron M, Gould K, Mead D, Drury E, O'Brien J, et al. (2016).
*Indels, structural variation, and recombination drive genomic diversity in Plasmodium falciparum*.
Genome Research, 26(9), 1288–1299.
\doi{10.1101/gr.203711.115}
}

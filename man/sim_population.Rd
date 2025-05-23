% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/main.R
\name{sim_population}
\alias{sim_population}
\title{Population simulator}
\usage{
sim_population(
  input_samples,
  n_samples_out,
  pop_alpha,
  coi_r,
  coi_p,
  k_s,
  max_coi = 100,
  rho = 7.4e-07,
  chrom_sizes = get_pf3d7_chrom_sizes()
)
}
\arguments{
\item{input_samples}{a character list of samples, this will be the ancestral genotypes for the simulated population}

\item{n_samples_out}{the number of samples to generate}

\item{pop_alpha}{the alpha parameter to set for the dirichlet function for the population proportions, lower alphas generate more closely related final populations}

\item{coi_r, coi_p}{the r and p parameters to be given to the zero truncated negative binomial distribution COI generator `recombuddy::rztnbinom()` see `?rztnbinom` for more details}

\item{k_s}{the s parameter to be given to the type 1 geometric distribution random generator function `rgeom()` to select for serial meiosis}

\item{max_coi}{the maximum allowable COI}

\item{rho}{the recombination rate}

\item{chrom_sizes}{a named vector of lengths of chromosome lengths}
}
\value{
a list the simulated population and look up tables of population proportions, sample indexes, and chromosome sizes
}
\description{
Simulate a number of samples of a new population from a previous ancestral population
}
\examples{
# simulate pop_alpha 9 (~10\% between sample relatedness), coi_r = 0.25, coi_p = 0.7 (COI mean of 1.256, ~80.4 proportion will be monoclonal), k_s = 0.5 (50\% of genotypes will be recombinant)
pop1 = sim_population(paste0("sample", seq(0,100,1)), 5, pop_alpha = 9, coi_r = 0.25, coi_p = 0.7, k_s = 0.5)
}

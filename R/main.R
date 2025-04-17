
#------------------------------------------------
# Declare global variable bindings to prevent R CMD check notes caused by the
# use of non-standard evaluation (NSE) in packages like dplyr and ggplot2. These
# variables are referenced in pipelines and tidyverse code but are not
# explicitly defined within function scopes.
if (getRversion() >= "2.15.1") {
  utils::globalVariables(c(
    "pos", "WSAF", "REF", "ALT", "PLAF", "chrom",
    "index", "start", "end", "genotype", "prop",
    "chrom_start", "chrom_col", "group"
  ))
}

#------------------------------------------------
#' @title Get chromosome lengths for P.falciparum 3D7
#'
#' @description Utility function that returns a vector of chromosome lengths for
#'   all 14 chromosomes of the 3D7 genome (version 2020-09-01)
#'
#' @return returns a named vector of the lengths of P.falciparum 3D7 chromosomes
#'
#' @export

get_pf3d7_chrom_sizes <- function() {
  c(Pf3D7_01_v3 = 640851, Pf3D7_02_v3 = 947102, Pf3D7_03_v3 = 1067971, Pf3D7_04_v3 = 1200490, Pf3D7_05_v3 = 1343557,
    Pf3D7_06_v3 = 1418242, Pf3D7_07_v3 = 1445207, Pf3D7_08_v3 = 1472805, Pf3D7_09_v3 = 1541735, Pf3D7_10_v3 = 1687656,
    Pf3D7_11_v3 = 2038340, Pf3D7_12_v3 = 2271494, Pf3D7_13_v3 = 2925236, Pf3D7_14_v3 = 3291936)
}

#------------------------------------------------
#' @title Get chromosome lengths for P.vivax P01
#'
#' @description Utility function that returns a vector of chromosome lengths for
#'   all 14 chromosomes of the P01 genome (version 2020-09-01)
#'
#' @return returns a named vector of the lengths of P.vivax P01 chromosomes
#'
#' @export

get_pvp01_chrom_sizes <- function() {
  c(PvP01_01_v2 = 1021664, PvP01_02_v2 = 956327, PvP01_03_v2 = 896704, PvP01_04_v2 = 1012024, PvP01_05_v2 = 1524814,
    PvP01_06_v2 = 1042791, PvP01_07_v2 = 1652210, PvP01_08_v2 = 1761288, PvP01_09_v2 = 2237066, PvP01_10_v2 = 1548844,
    PvP01_11_v2 = 2131221, PvP01_12_v2 = 3182763, PvP01_13_v2 = 2093556, PvP01_14_v2 = 3153402)
}

#------------------------------------------------
#' @title Draw from symmetric Dirichlet distribution
#'
#' @description Utility function for drawing from a symmetric Dirichlet
#'   distribution, where all concentration parameters are equal to \code{alpha / n}.
#'
#' @param n the number of dimensions ofthe Dirichlet distribution.
#' @param alpha the unscaled concentration parameter. The final concentration
#'   parameter is \code{alpha / n}
#'
#' @importFrom MCMCpack rdirichlet
#' @export

rdirichlet_single <- function(n, alpha = 1.0) {
  MCMCpack::rdirichlet(1, rep(alpha / n, n)) |>
    as.vector()
}

#------------------------------------------------
#' @title Simulate segments in a non-recombining genotype
#'
#' @description Produces a tibble of segments along the genome, assuming the
#'   genotype is copied over exactly from a single ancestor (no meiosis).
#'
#' @param index the index of the ancestor in the sample set.
#' @param chrom_sizes a vector of chromosome sizes.
#' @param zero_based_positioning whether to make the output positions zero based
#'
#' @importFrom tibble tibble
#' @export

sim_nonrecomb <- function(index, chrom_sizes, zero_based_positioning = T) {
  tibble(chrom = names(chrom_sizes), start = ifelse(zero_based_positioning, 0, 1), end = chrom_sizes, index = index)
}

#------------------------------------------------
#' @title Simulate segments in a recombining genotype
#'
#' @description Produces a tibble of segments along the genome produced by
#'   meiosis of several distinct ancestors.
#'
#' @details `k` sets the number of serial meiosis events. The number of
#'   ancestors chosen from the sample set is equal to 2^k. Although the same
#'   ancestor may be chosen multiple times by chance, it is not allowed for all
#'   2^k parents to come from the same ancestor, i.e. there is an element of
#'   sampling without replacement to avoid this.
#'
#' @param k the number of serial meiosis events.
#' @param rho the recombination rate (per-site, per-meiosis). By default uses
#'   the value 7.4e-7 from Miles et al. (2016).
#' @param set_props proportions of each of the members of the sample set.
#'   Dictates the probability they are chosen as an ancestor.
#' @param chrom_sizes lengths of each chromosome, default taken from `get_pf3d7_chrom_sizes()`
#'   by default.
#'
#' @references
#' Miles A, Iqbal Z, Vauterin P, Pearson R, Campino S, Theron M, Gould K, Mead D, Drury E, O'Brien J, et al. (2016).
#' *Indels, structural variation, and recombination drive genomic diversity in Plasmodium falciparum*.
#' Genome Research, 26(9), 1288–1299.
#' \doi{10.1101/gr.203711.115}
#'
#' @importFrom stats rgeom
#' @importFrom tibble tibble
#' @import dplyr
#' @export

sim_recomb <- function(k, rho = 7.4e-7, set_props, chrom_sizes, zero_based_positioning = T) {

  # we want to draw 2^k parents, but they cannot be all identical. Therefore,
  # draw (2^k - 1) with replacement and then draw the last one to be distinct
  n_set <- length(set_props)
  n_parents <- 2^k
  parents <- sample(x = n_set, size = n_parents, replace = TRUE, prob = set_props)
  if (all(parents == parents[1])) {
    parents[1] <- sample(x = (1:n_set)[-parents[1]], size = 1, prob = set_props[-parents[1]])
  }

  ret_list <- list()
  for (chrom in seq_along(chrom_sizes)) {

    # draw positions of recombination breakpoints
    p <- 1 - exp(-k*rho)
    x0_vec <- x1_vec <- NULL
    x0 <- x1 <- 0
    while (x1 < chrom_sizes[chrom]) {
      x0 <- x1 + 1
      x1 <- x0 + rgeom(1, prob = p)
      x1 <- min(x1, chrom_sizes[chrom])
      x0_vec <- c(x0_vec, x0)
      x1_vec <- c(x1_vec, x1)
    }
    n_segments <- length(x0_vec)

    # draw parent index for each segment
    w_vec <- rep(NA, n_segments)
    w_vec[1] <- sample(n_parents, 1)
    if (n_segments > 1) {
      for (i in 2:n_segments) {
        w_vec[i] <- sample((1:n_parents)[-w_vec[i-1]], 1)
      }
    }
    parent_vec <- parents[w_vec]

    # create raw table, before merging segments
    tab_raw <- tibble(start = x0_vec,
                      end = x1_vec,
                      index = parent_vec)

    # merge adjacent segments with same index
    tab_merge <- tab_raw |>
      mutate(group = cumsum(index != lag(index, default = first(index)))) |>
      group_by(group) |>
      summarise(
        start = min(start) - ifelse(zero_based_positioning, 1, 0),
        end = max(end),
        index = index[1]
      ) |>
      select(-group) |>
      mutate(chrom = names(chrom_sizes)[chrom], .before = 0)

    ret_list[[chrom]] <- tab_merge
  }

  ret_list |>
    bind_rows()
}

#------------------------------------------------
#' @title Simulate a single sample
#'
#' @description Simulates one or more haploid genotypes within a single sample
#'   by drawing ancestors from a sample set and applying recombination.
#'
#' @inheritParams sim_recomb
#'
#' @return A list with one element per haploid genome. Within each element is
#'   another list with three elements; the first two specify whether the
#'   genotype was copied over without recombination, and if so, which ancestor
#'   was it copied from. The final element `segments` contains the richest
#'   information - a tibble giving segments along the genome and which memeber
#'   of the sample set is ancestral to each segment.
#'
#' @references
#' Miles A, Iqbal Z, Vauterin P, Pearson R, Campino S, Theron M, Gould K, Mead D, Drury E, O'Brien J, et al. (2016).
#' *Indels, structural variation, and recombination drive genomic diversity in Plasmodium falciparum*.
#' Genome Research, 26(9), 1288–1299.
#' \doi{10.1101/gr.203711.115}
#'
#' @importFrom tibble tibble
#' @export

sim_sample <- function(k, rho = 7.4e-7, set_props, chrom_sizes = get_pf3d7_chrom_sizes()) {

  # get basic measures
  COI <- length(k)
  n_nonrecomb <- sum(k == 0)
  n_set <- length(set_props)
  n_chrom <- length(chrom_sizes)

  # draw samples without replacement for non-recombinants
  index_nonrecomb <- rep(NA, COI)
  if (n_nonrecomb > n_set) {
    stop(sprintf("Cannot generate %s distinct non-recombinant genotypes from sample set of size %s", n_nonrecomb, n_set))
  }
  if (n_nonrecomb > 0) {
    index_nonrecomb[which(k == 0)] <- sample(x = n_set, size = n_nonrecomb, prob = set_props)
  }

  ret <- list()
  ret$k = k
  ret$COI = COI
  ret$rho = rho
  ret$genotypes = list()
  for (i in 1:COI) {
    ret$genotypes[[i]] <- list()
    ret$genotypes[[i]]$is_nonrecomb <- (k[i] == 0)
    ret$genotypes[[i]]$index_nonrecomb <- index_nonrecomb[i]

    if (k[i] == 0) {
      ret$genotypes[[i]]$segments <- sim_nonrecomb(index = index_nonrecomb[i], chrom_sizes = chrom_sizes)
    } else {
      ret$genotypes[[i]]$segments <- sim_recomb(k = k[i], rho = rho, set_props = set_props, chrom_sizes = chrom_sizes)
    }
  }
  return(ret)
}

#------------------------------------------------
#' @title Collect all `segment` tables together
#'
#' @description For a given sample, extracts the `segments` element of each
#'   haploid genotype into a single tibble.
#'
#' @param sim a single simulated sample, as produced by `sim_sample()`.
#'
#' @import dplyr
#' @export

get_all_genotypes <- function(sim) {

  mapply(function(i) {
    sim$genotypes[[i]]$segments |>
      mutate(genotype = i,
             .before = 0)
  }, seq_along(sim$genotypes), SIMPLIFY = FALSE) |>
    bind_rows()
}

#------------------------------------------------
#' @title Block visualization of genome segments
#'
#' @description Takes a table of genomic segments (see `?get_all_genotypes`) and
#'   produces a simple plot showing chromosomes as blocks where the colour of a
#'   segment indicates which member of the sample set is ancestral. This
#'   provides an easy way to view blocks of shared ancestry.
#'
#' @param df_genotypes a data.frame (or tibble) of genomic segments, as produced
#'   by `get_all_genotypes()`.
#'
#' @importFrom grDevices grey
#' @import ggplot2
#' @export

plot_genotypes <- function(df_genotypes) {

  df_genotypes |>
    mutate(chrom = factor(chrom)) |>
    ggplot() + theme_void() +
    geom_rect(aes(xmin = start, xmax = end, ymin = as.numeric(chrom) - 0.4, ymax = as.numeric(chrom) + 0.4,
                  fill = as.factor(index)), color = grey(0)) +
    facet_wrap(~sprintf("Genotype_%s", genotype)) +
    scale_fill_discrete(name = "Ancestral Index") +
    scale_y_continuous(breaks = 1:n_distinct(df_genotypes$chrom),
                       labels = levels(factor(df_genotypes$chrom)))
}

#------------------------------------------------
#' @title Sample loci at random along the genome
#'
#' @description Creates a data.frame (tibble) of marker positions (loci) by
#'   sampling sites uniformly along the genome. Can be useful when using a
#'   simulated sample set, rather than a true reference panel.
#'
#' @param n_loci number of markers.
#' @inheritParams sim_recomb
#'
#' @importFrom stats rmultinom
#' @importFrom tibble tibble
#' @import dplyr
#'
#' @export

draw_loci <- function(n_loci = 1e3,
                      chrom_sizes = get_pf3d7_chrom_sizes()) {

  # assume loci of interest are distributed evenly around the genome. Therefore,
  # draw the number of loci per chromosome based on chromosome lengths
  n_chrom <- length(chrom_sizes)
  n_loci_per_chrom <- rmultinom(1, size = n_loci, prob = chrom_sizes)[,1]

  # draw positions per chromosome
  mapply(function(n, i) {
    if (n == 0) {
      return(NULL)
    }
    tibble(chrom = names(chrom_sizes)[i],
           pos = sort(sample(chrom_sizes[i], size = n, replace = FALSE)))
  }, n_loci_per_chrom, 1:n_chrom, SIMPLIFY = FALSE) |>
    bind_rows()
}

#------------------------------------------------
#' @title Sample population-level allele frequencies (PLAFs) at given loci
#'
#' @description For each locus in a data.frame, draws a PLAF from a Beta
#'   distribution with specified shape parameters. These parameters should be
#'   chosen to represent the allele frequency spectrum of *loci of interest*,
#'   which may differ from the frequency spectrum of background mutations (i.e.
#'   tending towards intermediate frequencies).
#'
#' @param df_loci a data.frame (tibble) of locus positions, as produced by
#'   `draw_loci()`.
#' @param beta_shape1,beta_shape2 shape parameters of the Beta distribution on
#'   PLAFs.
#'
#' @importFrom stats rbeta
#' @import dplyr
#' @export

draw_PLAF <- function(df_loci,
                      beta_shape1 = 0.5,
                      beta_shape2 = 0.5) {
  df_loci |>
    mutate(PLAF = rbeta(n = length(pos), shape1 = beta_shape1, shape2 = beta_shape2))
}

#------------------------------------------------
#' @title Draw a sample set given known population-level allele frequencies (PLAFs)
#'
#' @description Takes a data.frame (tibble) of known PLAFs, as produced by
#'   `draw_PLAF()`, and simulates a sample set with alleles drawn from these
#'   frequencies. The PLAF is taken to be the frequency of the REF allele.
#'
#' @param df_PLAF a data.frame (tibble) of PLAFs, as produced by `draw_PLAF()`.
#' @param n_set the number of members of the sample set.
#'
#' @importFrom stats runif
#' @import dplyr
#' @import tidyr
#'
#' @export

draw_sample_set_WSAF <- function(df_PLAF, n_set) {
  df_PLAF |>
    expand_grid(index = 1:n_set) |>
    mutate(WSAF = as.numeric(runif(length(PLAF)) < PLAF)) |>
    select(index, chrom, pos, WSAF) |>
    arrange(index, chrom, pos)
}

#------------------------------------------------
#' @title Sample genotype proportions
#'
#' @description Although a sample may contain multiple haploid genotypes, these
#'   are not usually present in equal proportions. This function draws genotype
#'   proportions from a symmetric Dirichlet distribution, giving control over
#'   the level of skew.
#'
#' @param COI the complexity of infection of the sample.
#' @param alpha the concentration parameter of the symmetric Dirichlet
#'   distribution on genotype proportions.
#'
#' @importFrom tibble tibble
#' @export

draw_genotype_props <- function(COI, alpha = 10.0) {
  tibble(genotype = 1:COI,
         prop = rdirichlet_single(n = COI, alpha = alpha))
}

#------------------------------------------------
#' @title Calculates within-sample allele frequencies (WSAFs) of a sample
#'
#' @description Combines information over multiple data.frames to calculate the
#'   overall within-sample allele frequency (WSAF) of the sample. This is the
#'   "convolution" of all the haploid genotypes in the sample.
#'
#' @details Takes a tibble containing segments along the genome of all haploid
#'   genotypes (as produced by `get_all_genotypes()`) and combines this with a
#'   tibble of the observed alleles in a sample set (as produced by
#'   `draw_sample_set_WSAF()`) to calculate the overall WSAF. Haploid genotypes
#'   are combined in proportions defined by another tibble, as produced by
#'   `draw_genotype_props()`.
#'
#' @param df_all_genotypes a data.frame (tibble) containing segments along the
#'   genome of all haploid genotypes, as produced by `get_all_genotypes()`.
#' @param df_sample_set a data.frame (tibble) of the observed alleles in a
#'   sample set, as produced by `draw_sample_set_WSAF()`.
#' @param df_genotype_props a data.frame (tibble) of genomic proportions, as
#'   produced by `draw_genotype_props()`.
#'
#' @import dplyr
#' @import ggplot2
#'
#' @export

get_WSAF <- function(df_all_genotypes, df_sample_set, df_genotype_props) {

  # merge genotypes table with sample set info
  df_all_genotypes |>
    left_join(df_sample_set, relationship = "many-to-many",
              by = join_by(chrom, index)) |>
    filter(pos >= start) |>
    filter(pos <= end) |>
    left_join(df_genotype_props, by = join_by(genotype)) |>
    group_by(chrom, pos) |>
    summarise(WSAF = sum(WSAF * prop), .groups = "drop")
}

#------------------------------------------------
#' @title Simulate read counts given known within-sample allele frequencies (WSAFs).
#'
#' @description Uses a simple statistical model to draw simulated read counts.
#'   The expected value is given by the WSAF at a locus, the distribution is
#'   beta-binomial with a given denominator size (read depth) and
#'   over-dispersion parameter.
#'
#' @details draws from `rbbinom()` in the `extraDistr` package using the
#'   parameters:
#'   \itemize{
#'     \item `alpha = WSAF*(1 - overdisp) / overdisp`
#'     \item `beta = (1 - WSAF)*(1 - overdisp) / overdisp)`
#'   }
#'
#' @param df_WSAF a data.frame (tibble) of WSAFs, as produced by `get_WSAF()`.
#' @param depth the *size* argument in the beta-binomial.
#' @param overdisp the over-dispersion of the beta-binomial distribution.
#'
#' @importFrom extraDistr rbbinom
#' @import dplyr
#'
#' @export

draw_read_counts <- function(df_WSAF, depth = 100, overdisp = 0.01) {

  df_WSAF |>
    mutate(alpha = WSAF*(1 - overdisp) / overdisp,
           beta = (1 - WSAF)*(1 - overdisp) / overdisp) |>
    rowwise() |>
    mutate(REF = extraDistr::rbbinom(1, size = depth, alpha = alpha, beta = beta),
           ALT = depth - REF) |>
    select(chrom, pos, REF, ALT)
}

#------------------------------------------------
#' @title Plots draws from the read counts model
#'
#' @description Takes values produced by `draw_read_counts()` and produces a
#'   simple plot showing *estimated* WSAFs along the genome. This simple plot
#'   often contains patterns that demonstrate clear blocks of relatedness.
#'
#' @param df_counts a data.frame (tibble) of read counts, as produced by
#'   `draw_read_counts()`.
#' @inheritParams sim_recomb
#'
#' @import dplyr
#' @import ggplot2
#'
#' @export

plot_read_counts <- function(df_counts, chrom_sizes = get_pf3d7_chrom_sizes()) {

  df_counts |>
    left_join(data.frame(chrom = names(chrom_sizes),
                         chrom_start = cumsum(chrom_sizes) - chrom_sizes + 1),
              by = join_by(chrom)) |>
    mutate(chrom = factor(chrom)) |>
    mutate(pos = pos + chrom_start - 1,
           chrom_col = as.factor((as.numeric(chrom) %% 2) + 1)) |>
    ggplot() + theme_bw() +
    geom_point(aes(x = pos, y = REF / (REF + ALT), col = chrom_col)) +
    scale_color_manual(values = c("black", "firebrick1")) +
    guides(col = "none") +
    ylab("Within-sample allele frequency (REF)") + xlab("Position on contiguous genome")
}

#------------------------------------------------
#' @title Zero-truncated negative binomial distribution
#'
#' @description Implements the zero-truncated negative binomial distribution, in
#'   which there is zero chance of seeing `n=0`.
#'
#' @details The negative binomial distribution has several alternative
#'   parameterizations. We follow the same parameterization implemented in
#'   `dnbinom()`, with the small caveat that defining via `mu` is not available,
#'   only via `size` and `prob`. Note that the mean of the truncated
#'   distribution is different from the mean of the un-truncated distribution,
#'   now being given by `size*(1 - prob) / (prob*(1 - prob^size))`.
#'
#' @param x,size,prob,log these parameters are passed directly to `dnbinom()`.
#'   See `?dnbinom` for further details of these parameters.
#'
#' @importFrom stats dnbinom qnbinom pnbinom
#' @export

dztnbinom <- function(x, size, prob, log = FALSE) {
  s <- dnbinom(x = 0, size = size, prob = prob, log = FALSE)
  ret <- dnbinom(x = x, size = size, prob = prob, log = FALSE) / (1 - s)
  ret[x == 0] <- 0
  if (log) {
    ret <- log(ret)
  }
  return(ret)
}


#' @title Zero-truncated negative binomial distribution random number generator
#'
#' @param n the number of random numbers to generate
#' @param size passed directly to `dnbinom()`.See `?dnbinom` for further details of these parameters.
#' @param prob passed directly to `dnbinom()`.See `?dnbinom` for further details of these parameters.
#'
#' @returns `n` random numbers from a Zero-truncated negative binomial distribution
#' @importFrom stats dnbinom qnbinom pnbinom
#' @export
#'
#' @examples rztnbinom(10, 0.25, 0.7)
rztnbinom <- function(n, size, prob) {
  prob0 <- dnbinom(x = 0, size = size, prob = prob, log = FALSE)
  q_thresh <- prob0 + (1 - prob0)*(1 - 1e-4)
  n_max <- qnbinom(q_thresh, size = size, prob = prob)
  cumul_prob <- pnbinom(0:n_max, size = size, prob = prob) - prob0
  X <- runif(n, min = 0, max = max(cumul_prob))
  ret <- as.numeric(cut(X, breaks = cumul_prob, include.lowest = TRUE))
  return(ret)
}



#' @title Generate a vector of recombinant parameters
#'
#' @description
#' Will generate a vector of Serial meiosis parameters, the length will be the COI of the sample
#'
#'
#' @param coi_r,coi_p the r and p parameters to be given to the zero truncated negative binomial distribution COI generator `recombuddy::rztnbinom()` see `?rztnbinom` for more details
#' @param k_s the s parameter to be given to the type 1 geometric distribution random generator function `rgeom()` to select for serial meiosis
#' @param max_coi the maximum allowable COI
#'
#' @returns a vector of ks, to be given to `recombuddy::sim_sample()`
#' @export
#' @importFrom stats rgeom
#' @examples generate_serial_meiosis_k(0.25, 0.7, 0.5)
generate_serial_meiosis_k <-function(coi_r, coi_p, k_s, max_coi = 100){
  # randomly sample from a zero truncated negative binomial distribution
  sample_coi = rztnbinom(1, coi_r, coi_p)
  while(sample_coi > max_coi){
    sample_coi = rztnbinom(1, coi_r, coi_p)
  }
  # randomly sample the serial meiosis from a type 1 geometric distribution
  sample_k = rgeom(sample_coi, k_s)
  return (sample_k)
}

#' @title Population simulator
#'
#' @description
#' Simulate a number of samples of a new population from a previous ancestral population
#'
#'
#' @param input_samples a character list of samples, this will be the ancestral genotypes for the simulated population
#' @param n_samples_out the number of samples to generate
#' @param pop_alpha the alpha parameter to set for the dirichlet function for the population proportions, lower alphas generate more closely related final populations
#' @param coi_r,coi_p the r and p parameters to be given to the zero truncated negative binomial distribution COI generator `recombuddy::rztnbinom()` see `?rztnbinom` for more details
#' @param k_s the s parameter to be given to the type 1 geometric distribution random generator function `rgeom()` to select for serial meiosis
#' @param max_coi the maximum allowable COI
#' @param rho the recombination rate
#' @param chrom_sizes a named vector of lengths of chromosome lengths
#'
#' @returns a list the simulated population and look up tables of population proportions, sample indexes, and chromosome sizes
#' @export
#' @import dplyr
#' @importFrom tibble tibble
#' @examples
#' # simulate pop_alpha 9 (~10% between sample relatedness), coi_r = 0.25, coi_p = 0.7 (COI mean of 1.256, ~80.4 proportion will be monoclonal), k_s = 0.5 (50% of genotypes will be recombinant)
#' pop1 = sim_population(paste0("sample", seq(0,100,1)), 5, pop_alpha = 9, coi_r = 0.25, coi_p = 0.7, k_s = 0.5)
sim_population <- function(input_samples, n_samples_out, pop_alpha, coi_r, coi_p, k_s, max_coi = 100, rho = 7.4e-7, chrom_sizes = get_pf3d7_chrom_sizes()){
  ret = list()
  # generate a index key tibble for samples
  input_samples_df = tibble(ancestral_genotype = input_samples) |> mutate(index = row_number())
  ret[["ancestral_indexes"]] = input_samples_df
  # set proportions
  set_props <- rdirichlet_single(nrow(input_samples_df), alpha = pop_alpha)
  ret[["set_props"]] = set_props
  ret[["chrom_sizes"]] = chrom_sizes
  # simulate samples
  ret[["simulated_samples"]] = list()
  for(samp in 1:n_samples_out){
    ret[["simulated_samples"]][[samp]] = sim_sample(k = generate_serial_meiosis_k(coi_r = coi_r, coi_p = coi_p, k_s = k_s, max_coi = max_coi),
                                                    rho = rho,
                                                    set_props = set_props,
                                                    chrom_sizes = chrom_sizes)
  }
  return (ret)
}



#' @title Validate panel locations
#'
#' @description
#' Run validator with rules to make sure that there 4 columns and the following types chrom (character), start (numeric), end(numeric), target(character), will error out if fails to validate
#'
#' @param panel_locs a data frame with 4 columns, chrom (character), start (numeric), end(numeric), target(character)
#'
#' @returns nothing
#' @export
#' @importFrom validate validator confront summary
#' @import dplyr
#' @examples
#' panel <- tibble(chrom = c("Pf3D7_04_v3", "Pf3D7_05_v3"), start = c(748173, 958071), end = c(748361, 958206), target = c("dhfr-1", "dhps-1"))
#' validate_panel_locs_df(panel)
validate_panel_locs_df <- function(panel_locs){
  # validate columns of the intersecting panel_locs
  rules <- validate::validator(
    is.character(chrom),
    is.numeric(start),
    is.numeric(end),
    is.character(target),
    ! is.na(chrom),
    ! is.na(start),
    ! is.na(end),
    ! is.na(target)
  )
  fails <- validate::confront(panel_locs, rules) |>
    validate::summary() |>
    dplyr::filter(error)
  warns = c()
  if (nrow(fails) > 0) {
    stop(paste0(
      "panel_locs failed one or more validation checks, check if required columns are present and that they pass these checks: \n",
      paste0(fails$expression, collapse = "\n")
    ) )
  }
}

#' @title Intersect target locations for a panel with the recombined segments
#' @param segs the segments from one of the genotypes generated from `sim_sample`
#' @param panel_locs a data frame with 4 columns, chrom, start, end, target, chrom values must match the recombined segment chroms
#'
#' @returns a table with the sample indexes from which each panel intersects this, this is the ancestral genotype for this target in the panel
#' @export
#' @import dplyr
#' @importFrom tibble tibble
#' @examples
#' samp1 <- sim_sample(k = c(0, 2, 4), rho = 7.4e-7, set_props = rdirichlet_single(10, alpha = 9))
#' panel <- tibble(chrom = c("Pf3D7_04_v3", "Pf3D7_05_v3"), start = c(748173, 958071), end = c(748361, 958206), target = c("dhfr-1", "dhps-1"))
#' intersect_panel_with_segs(panel, samp1$genotypes[[3]]$segments)
intersect_panel_with_segs <-function(panel_locs, segs){
  validate_panel_locs_df(panel_locs)
  panel_with_sample_index = panel_locs |>
    mutate(index = NA)
  # fill in index with first overlapping segment, if more than one segment overlaps aka recombination happened within target than it will be first come first serve
  for(panel_row in 1:nrow(panel_with_sample_index)){
    for(seg_row in 1:nrow(segs)){
      if(panel_with_sample_index$chrom[panel_row] == segs$chrom[seg_row] &&
         ( (panel_with_sample_index$start[panel_row] >= segs$start[seg_row] && panel_with_sample_index$start[panel_row] < segs$end[seg_row]) ||
           (panel_with_sample_index$end[panel_row] > segs$start[seg_row] && panel_with_sample_index$end[panel_row] <= segs$end[seg_row]))
      ){
        panel_with_sample_index$index[panel_row] = segs$index[seg_row]
        break
      }
    }
  }
  return(panel_with_sample_index)
}



#' @title Intersect a panel with the all genotypes from a simulated sample
#'
#' @param panel_locs a data frame with 4 columns, chrom, start, end, target, chrom values must match the recombined segment chroms
#' @param simulated_sample a simulated generated from `recombuddy::sim_sample`
#'
#' @returns a table with the genotype and sample index from the simulated sample corresponding to the panel of intersect
#' @export
#'
#' @examples
#' samp1 <- sim_sample(k = c(0, 2, 4), rho = 7.4e-7, set_props = rdirichlet_single(10, alpha = 9))
#' panel <- tibble(chrom = c("Pf3D7_04_v3", "Pf3D7_05_v3"), start = c(748173, 958071), end = c(748361, 958206), target = c("dhfr-1", "dhps-1"))
#' intersect_panel_with_sample_genotypes(panel, samp1)
intersect_panel_with_sample_genotypes <-function(panel_locs, simulated_sample){
  validate_panel_locs_df(panel_locs)
  genos = get_all_genotypes(simulated_sample)
  genotypes <- tibble(genotype = unique(sort(genos$genotype)))
  # Create the expanded data frame
  panel_with_sample_index <- crossing(panel_locs, genotypes)|>
    mutate(index = NA) |>
    arrange(genotype)
  # fill in index with first overlapping segment, if more than one segment overlaps aka recombination happened within target than it will be first come first serve
  for(panel_row in 1:nrow(panel_with_sample_index)){
    for(geno_row in 1:nrow(genos)){
      if(
        panel_with_sample_index$genotype[panel_row] == genos$genotype[geno_row] &&
        panel_with_sample_index$chrom[panel_row] == genos$chrom[geno_row] &&
         ( (panel_with_sample_index$start[panel_row] >= genos$start[geno_row] && panel_with_sample_index$start[panel_row] < genos$end[geno_row]) ||
           (panel_with_sample_index$end[panel_row] > genos$start[geno_row] && panel_with_sample_index$end[panel_row] <= genos$end[geno_row]))
      ){
        panel_with_sample_index$index[panel_row] = genos$index[geno_row]
        break
      }
    }
  }
  return(panel_with_sample_index)
}


#' @title Intersect locations with populations
#'
#' @description
#' Intersect a table of genomic locations with simulated samples from a simulated population to get a table of final genotypes joined with what ancestral genotype they came from
#'
#' @param panel_locs a data frame with 4 columns, chrom (character), start (numeric), end(numeric), target(character)
#' @param simulated_population a population generated with
#'
#' @returns a table with chrom,start,end,target,within_sample_genotype(a within sample genotype),ancestral_index, simulated_sample (a name for the simulate sample, will just be a number in order of simulated population),ancestral_genotype(original genotype for target)
#' @export
#' @import dplyr
#' @importFrom tibble tibble
#'
#'
#' @examples
#' panel <- tibble(chrom = c("Pf3D7_04_v3", "Pf3D7_05_v3"), start = c(748173, 958071), end = c(748361, 958206), target = c("dhfr-1", "dhps-1"))
#' # simulate pop_alpha 9 (~10% between sample relatedness), coi_r = 0.25, coi_p = 0.7 (COI mean of 1.256, ~80.4 proportion will be monoclonal), k_s = 0.5 (50% of genotypes will be recombinant)
#' pop1 = sim_population(paste0("sample", seq(0,100,1)), 5, pop_alpha = 9, coi_r = 0.25, coi_p = 0.7, k_s = 0.5)
#' intersect_panel_with_simulated_population(panel, pop1)
intersect_panel_with_simulated_population <-function(panel_locs, simulated_population){
  validate_panel_locs_df(panel_locs)
  all_sample_genotypes = tibble()
  for(samp in seq_along(simulated_population[["simulated_samples"]])){
    panel_with_genos_samp_indexes = intersect_panel_with_sample_genotypes(panel_locs, simulated_population[["simulated_samples"]][[samp]]) |>
      mutate(simulated_sample = samp)
    all_sample_genotypes = bind_rows(
      all_sample_genotypes,
      panel_with_genos_samp_indexes
    )
  }
  return(all_sample_genotypes |> left_join(simulated_population[["ancestral_indexes"]]) |>  dplyr::rename(within_sample_genotype = genotype, ancestral_index = index))
}


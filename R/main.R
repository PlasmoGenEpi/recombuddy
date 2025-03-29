
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
#'   all 14 chromosomes of the 3D7 genome.
#'
#' @export

get_chrom_sizes <- function() {
  c(chr1 = 640851, chr2 = 947102, chr3 = 1067971, chr4 = 1200490, chr5 = 1343557,
    chr6 = 1418242, chr7 = 1445207, chr8 = 1472805, chr9 = 1541735, chr10 = 1687656,
    chr11 = 2038340, chr12 = 2271494, chr13 = 2925236, chr14 = 3291936)
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
#'
#' @importFrom tibble tibble
#' @export

sim_nonrecomb <- function(index, chrom_sizes) {
  tibble(chrom = seq_along(chrom_sizes), start = 1, end = chrom_sizes, index = index)
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
#' @param chrom_sizes lengths of each chromosome, taken from `get_chrom_sizes()`
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

sim_recomb <- function(k, rho = 7.4e-7, set_props, chrom_sizes) {

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
        start = min(start),
        end = max(end),
        index = index[1]
      ) |>
      select(-group) |>
      mutate(chrom = chrom, .before = 0)

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

sim_sample <- function(k, rho = 7.4e-7, set_props, chrom_sizes = get_chrom_sizes()) {

  # get basic measures
  MOI <- length(k)
  n_nonrecomb <- sum(k == 0)
  n_set <- length(set_props)
  n_chrom <- length(chrom_sizes)

  # draw samples without replacement for non-recombinants
  index_nonrecomb <- rep(NA, MOI)
  if (n_nonrecomb > n_set) {
    stop(sprintf("Cannot generate %s distinct non-recombinant genotypes from sample set of size %s", n_nonrecomb, n_set))
  }
  if (n_nonrecomb > 0) {
    index_nonrecomb[which(k == 0)] <- sample(x = n_set, size = n_nonrecomb, prob = set_props)
  }

  ret <- list()
  for (i in 1:MOI) {
    ret[[i]] <- list()
    ret[[i]]$is_nonrecomb <- (k[i] == 0)
    ret[[i]]$index_nonrecomb <- index_nonrecomb[i]

    if (k[i] == 0) {
      ret[[i]]$segments <- sim_nonrecomb(index = index_nonrecomb[i], chrom_sizes = chrom_sizes)
    } else {
      ret[[i]]$segments <- sim_recomb(k = k[i], rho = rho, set_props = set_props, chrom_sizes = chrom_sizes)
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
    sim[[i]]$segments |>
      mutate(genotype = i,
             .before = 0)
  }, seq_along(sim), SIMPLIFY = FALSE) |>
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
    ggplot() + theme_void() +
    geom_rect(aes(xmin = start, xmax = end, ymin = chrom - 0.9, ymax = chrom - 0.1,
                  fill = as.factor(index)), color = grey(0)) +
    facet_wrap(~sprintf("Genotype_%s", genotype)) +
    scale_fill_discrete(name = "Parent Index")
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
                      chrom_sizes = get_chrom_sizes()) {

  # assume loci of interest are distributed evenly around the genome. Therefore,
  # draw the number of loci per chromosome based on chromosome lengths
  n_chrom <- length(chrom_sizes)
  n_loci_per_chrom <- rmultinom(1, size = n_loci, prob = chrom_sizes)[,1]

  # draw positions per chromosome
  mapply(function(n, i) {
    if (n == 0) {
      return(NULL)
    }
    tibble(chrom = i,
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
#' @param MOI the multiplicity of infection of the sample.
#' @param alpha the concentration parameter of the symmetric Dirichlet
#'   distribution on genotype proportions.
#'
#' @importFrom tibble tibble
#' @export

draw_genotype_props <- function(MOI, alpha) {
  tibble(genotype = 1:MOI,
         prop = rdirichlet_single(n = MOI, alpha = 10.0))
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

plot_read_counts <- function(df_counts, chrom_sizes = get_chrom_sizes()) {

  df_counts |>
    left_join(data.frame(chrom = seq_along(chrom_sizes),
                         chrom_start = cumsum(chrom_sizes) - chrom_sizes + 1),
              by = join_by(chrom)) |>
    mutate(pos = pos + chrom_start - 1,
           chrom_col = as.factor((chrom %% 2) + 1)) |>
    ggplot() + theme_bw() +
    geom_point(aes(x = pos, y = REF / (REF + ALT), col = chrom_col)) +
    scale_color_manual(values = c("black", "firebrick1")) +
    guides(col = "none") +
    ylab("Within-sample allele frequency (REF)") + xlab("Position on contiguous genome")
}



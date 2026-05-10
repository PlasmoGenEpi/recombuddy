#------------------------------------------------
#' @title Collect all `segment` tables together for all samples in a population
#'
#' @description For a population of samples, extracts the `segments` element of each
#'   haploid genotype per sample into a single population tibble.
#'
#' @param sim_pop a simulated population, as produced by `sim_population()` (see `?sim_population()`).
#'
#' @return
#' Returns a tibble with the following columns:
#' \itemize{
#'   \item `sample_id`
#'   \item `genotype`
#'   \item `chrom`
#'   \item `start`
#'   \item `end`
#'   \item `index`
#' }
#'
#' @details
#' This is a wrapper around `get_all_genotypes()` that applies the same logic
#' across all simulated samples in a population and appends a `sample_id` column
#' to identify which sample each genotype belongs to.
#'
#' @examples
#' # simulate a population of 5 samples
#' pop1 <- sim_population(
#'   paste0("sample", seq(0,100,1)),
#'   5,
#'   pop_alpha = 9,
#'   coi_r = 0.25,
#'   coi_p = 0.7,
#'   k_s = 0.5
#' )
#'
#' df_population_genotypes <- get_population_genotypes(pop1)
#' df_population_genotypes
#' @import dplyr
#' @export
get_population_genotypes <- function(sim_pop) {

  mapply(function(i) {
    get_all_genotypes(sim_pop$simulated_samples[[i]]) |>
      mutate(sample_id = i,
             .before = 0)
  }, seq_along(sim_pop$simulated_samples), SIMPLIFY = FALSE) |>
    bind_rows()
}

#------------------------------------------------
#' @title Add segment lengths to a segments table
#'
#' @description
#' Internal helper function that appends a `segment_length` column to a table of
#' genotype ancestry segments.
#'
#' @param df_genotypes a data.frame (or tibble) of genomic segments containing columns
#'   `start` and `end`, such as output from `get_all_genotypes()` (see `?get_all_genotypes`),
#'   or `get_population_genotypes()`.
#'
#' @return
#' Returns the input data.frame with a numeric `segment_length` column appended.
#' If `segment_length` is already present, the input is returned unchanged.
#'
#' @details
#' This function is intended for internal use by summary and plotting helpers
#' that operate on genomic segment data frames produced by
#' `get_all_genotypes()` or `get_population_genotypes()`.
#'
#' @keywords internal
#' @import dplyr
add_segment_lengths <- function(df_genotypes) {
  if (!"segment_length" %in% names(df_genotypes)) {
    df_genotypes <- df_genotypes |>
      mutate(segment_length = end - start)
  }
  df_genotypes
}

#------------------------------------------------
#' @title Get ancestry composition
#'
#' @description
#' Calculates the overall genomic proportion inherited from each ancestral index.
#'
#' @param df_genotypes a data.frame (or tibble) of ancestry segments, such as that produced
#'   by `get_all_genotypes()` for individual samples (see `?get_all_genotypes()`),
#'   or by `get_population_genotypes()` for simulated populations (see `?get_population_genotypes()`).
#'
#' @return
#' Returns a tibble with the following columns:
#' \itemize{
#'   \item `sample_id`
#'   \item `genotype`
#'   \item `index`
#'   \item `ancestry_total_bp`
#'   \item `ancestry_prop`
#' }
#'
#' @details
#' The input is first passed through an internal helper that appends a
#' `segment_length` column if one is not already present. Ancestry proportions
#' are then calculated for each haploid genotype within sample(s). The output of
#' `get_ancestry_composition` can then serve as input for the plotting helper
#' function `plot_ancestry_composition`.
#'
#' @examples
#' n_set <- 10
#' set_props <- rdirichlet_single(n_set, alpha = 9)
#' samp1 <- sim_sample(k = c(0, 1, 2), set_props = set_props)
#' df_all_genotypes <- get_all_genotypes(samp1)
#'
#' # For a single simulated sample
#' df_ancestry_sample <- get_ancestry_composition(df_all_genotypes)
#' df_ancestry_sample
#'
#' # for a simulated population
#' # df_ancestry_population <- get_ancestry_composition(df_population_genotypes)
#' # df_ancestry_population
#'
#' @import dplyr
#' @export
get_ancestry_composition <- function(df_genotypes) {

  # adds a sample_id=1 to keep plotting functionality grouped by sample_id
  if (!"sample_id" %in% names(df_genotypes)) {
    df_genotypes <- df_genotypes |>
      mutate(sample_id = 1) |>
      relocate(sample_id, .before = genotype)
  }

  # appends segment_lengths if not present already.
  if (!"segment_length" %in% names(df_genotypes)) {
    df_genotypes <- df_genotypes |>
      add_segment_lengths()
  }

  df_ancestry <- df_genotypes |>
    group_by(sample_id, genotype, index) |>
    summarise(ancestry_total_bp = sum(segment_length),
              .groups = "drop") |>
    group_by(sample_id, genotype) |>
    mutate(ancestry_prop = ancestry_total_bp / sum(ancestry_total_bp)) |>
    ungroup()

  df_ancestry
}

#------------------------------------------------
#' @title Plot ancestry composition
#'
#' @description
#' Produces a horizontal stacked bar plot showing the ancestry composition of
#' simulated genotypes. The length of each colored bar segment
#' corresponds to the genomic proportion inherited from a given ancestral index.
#' Separate panels are shown for each sample.
#'
#' @param df_ancestry a data.frame (or tibble) of ancestry composition, as
#'   produced by `get_ancestry_composition()` (see `?get_ancestry_composition()`).
#' @param by the hierarchical level at which to display ancestry composition.
#'   Must be either `"genotype"` or `"sample"`. The default is `"genotype"`.
#'
#' @return
#' Returns a `ggplot2` object.
#'
#' @details
#' Requires columns `sample_id`, `genotype`, `index` and `ancestry_prop`.
#' If `by = "genotype"`, genotype-level ancestry proportions
#' are plotted and separated into panels by sample. If `by = "sample"`,
#' ancestry composition is first summed across genotypes within each sample and
#' then plotted at the sample level.
#'
#' @examples
#' n_set <- 10
#' set_props <- rdirichlet_single(n_set, alpha = 9)
#' samp1 <- sim_sample(k = c(0, 1, 2), set_props = set_props)
#' df_all_genotypes <- get_all_genotypes(samp1)

#' # genotype-level ancestry composition
#' df_all_genotypes |>
#'   get_ancestry_composition() |>
#'   plot_ancestry_composition(by = "genotype")
#'
#' # sample-level ancestry composition
#' # df_population_genotypes |>
#' #  get_ancestry_composition() |>
#' #  plot_ancestry_composition(by = "sample")
#'
#' @import dplyr
#' @import ggplot2
#' @export
plot_ancestry_composition <- function(df_ancestry, by = c("genotype", "sample")) {

  by <- match.arg(by)

  if(by == "genotype"){
    df_plot <- df_ancestry |>
      mutate(sample_label = factor(paste0("Sample ", sample_id),
                                   levels = paste0("Sample ", sort(unique(as.numeric(sample_id))))),
             genotype_label = factor(paste0("Genotype ", genotype)))

    p <-  df_plot |>
      ggplot(aes(x = ancestry_prop, y = genotype_label, fill = factor(index))) +
        geom_col(width = 0.8, color = grey(0)) +
        scale_x_continuous(breaks = scales::pretty_breaks(n=3),
                           labels = scales::percent_format(),
                           expand = c(0,0)) +
        scale_fill_viridis_d(name = "Ancestral Index", option = "turbo") +
        facet_wrap(~sample_label, scales = "free_y") +
        labs(x = "Ancestry proportion",
             y = "") +
        theme_bw() +
        theme(axis.text.y = element_text(size = 6),
              strip.text = element_text(size = 6))
  }

  if(by == "sample"){

    sample_labels <- df_ancestry |>
      group_by(sample_id) |>
      summarise(n_genotypes = n_distinct(genotype),
                .groups = "drop") |>
      mutate(sample_label = paste0("Sample ", sample_id, " (COI=", n_genotypes, ")"))

    sample_levels <- sample_labels |>
      arrange(sample_id) |>
      pull(sample_label)

    df_plot <- df_ancestry |>
      group_by(sample_id, index) |>
      summarise(index_total_bp = sum(ancestry_total_bp),
                .groups = "drop") |>
      left_join(sample_labels, by = "sample_id") |>
      group_by(sample_id) |>
      mutate(sample_ancestry_prop = index_total_bp / sum(index_total_bp),
             sample_label = factor(sample_label, levels = sample_levels)) |>
      ungroup()

    p <- df_plot |>
      ggplot(aes(x = sample_ancestry_prop, y = sample_label, fill = factor(index))) +
        geom_col(width = 0.8, color = grey(0)) +
        scale_x_continuous(breaks = scales::pretty_breaks(n=5),
                           labels = scales::percent_format(),
                           expand = c(0,0)) +
        scale_fill_viridis_d(name = "Ancestral Index", option = "turbo") +
        labs(x = "Ancestry proportion",
             y = "") +
        theme_minimal()
  }

  p
}

#------------------------------------------------
#' @title Plot ancestry segment lengths
#'
#' @description
#' Produces a boxplot of ancestry segment lengths from a data frame of simulated ancestry
#' segments.
#'
#' @param df_genotypes a data.frame (or tibble) of ancestry segments, such as that produced
#'   by `get_all_genotypes()` for individual samples (see `?get_all_genotypes()`),
#'   or by `get_population_genotypes()` for simulated populations (see `?get_population_genotypes()`).
#'
#' @return
#' Returns a `ggplot2` object.
#'
#' @details
#' The input is first passed through an internal helper that appends a
#' `segment_length` column if one is not already present. Note that boxplots for
#' non-recombinant genotypes reflect chromosome sizes since there is only
#' one ancestral index for the entire genome.
#'
#' @examples
#' n_set <- 10
#' set_props <- rdirichlet_single(n_set, alpha = 9)
#' samp1 <- sim_sample(k = c(0, 1, 2), set_props = set_props)
#' df_all_genotypes <- get_all_genotypes(samp1)
#' df_all_genotypes |> plot_segment_lengths()
#'
#' @import dplyr
#' @import ggplot2
#' @export
plot_segment_lengths <- function(df_genotypes) {

  # adds a sample_id=1 to keep plotting functionality grouped by sample_id
  if (!"sample_id" %in% names(df_genotypes)) {
    df_genotypes <- df_genotypes |>
      mutate(sample_id = 1) |>
      relocate(sample_id, .before = genotype)
  }

  # appends segment_lengths if not present already
  if (!"segment_length" %in% names(df_genotypes)) {
    df_genotypes <- df_genotypes |>
      add_segment_lengths()
  }

  df_plot <- df_genotypes |>
      mutate(
        sample_label = factor(
          paste0("Sample ", sample_id),
          levels = paste0("Sample ", sort(unique(as.numeric(sample_id))))
        ),
        index = factor(index),
        genotype = factor(genotype)
      )

    p <- df_plot |>
      ggplot(aes(x = genotype,
                 y = segment_length,
                 group = interaction(genotype, index),
                 color = index,
                 fill = index)) +
      geom_boxplot(
        alpha = 0.1,
        position = position_dodge(width = 0.8),
        outlier.shape = NA
      ) +
      geom_jitter(
        position = position_jitterdodge(
          jitter.width = 0.15,
          dodge.width = 0.8
        ),
        alpha = 0.6
      ) +
      scale_fill_viridis_d(name = "Ancestral Index", option = "turbo") +
      scale_color_viridis_d(name = "Ancestral Index", option = "turbo") +
      scale_y_continuous(breaks = scales::pretty_breaks(5)) +
      labs(
        x = "Genotype",
        y = "Segment length (bp)"
      ) +
      facet_wrap(~sample_label, scales = "free_x") +
      theme_bw() +
      theme(axis.text.y = element_text(size = 6),
            strip.text = element_text(size = 6))

    p
}

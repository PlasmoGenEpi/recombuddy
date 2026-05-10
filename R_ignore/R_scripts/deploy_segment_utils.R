# deploy_segment_utils.R
#
# Author: Shazia Ruybal-Pesántez
# Date: 2026-05-04
#
# Purpose:
# Showing the usage of the new utility functions in the recombuddy package.
# This is a first pass and will likely need updates.
#
# I have built the functions to handle single simulated samples and also
# simulated populations (but as sample N grows, the plots get uglier!). So likely
# will need to edit or subset large plots.
#
# ------------------------------------------------------------------

# Single sample simulation ------------------------------------------------
# define the number of samples in the initial sample set, and in the final
# population that we are producing
n_set <- 10
n_pop <- 20

# define concentration parameter of the Dirichlet distribution
f <- 0.1
alpha <- (1 - f) / f

# draw proportions of each member of the sample set
set_props <- rdirichlet_single(n_set, alpha = alpha)

# Simulate a single sample
samp1 <- sim_sample(k = c(0, 1, 2), rho = 6e-7, set_props = set_props)

# get all genotypes into a single table
df_all_genotypes <- get_all_genotypes(samp1)

# plot genotypes with recombuddy function, NOTE: I edited the color pal
plot_genotypes(df_all_genotypes)

# New functionality for easily visualizing a genomic summary of ancestral composition:
df_ancestry_samp1 <- df_all_genotypes |> get_ancestry_composition()
df_ancestry_samp1

## By genotype
df_ancestry_samp1 |> plot_ancestry_composition(by = "genotype")

## For the whole sample (collapsing all genotypes). Given this is ground truth and we have the full mosaic genome
## I think this is OK - but of course this would be adjusted with weights for actual sample proportions of each genotype
df_ancestry_samp1 |> plot_ancestry_composition(by = "sample")

# New functionality for easily visualizing. NOTE: I don't love this, but as I was playing
# with simulating populations for the MAGIC grant writing I was interested to see segment block length.
# This will be interested when thinking of "f2 generation" imported infections that are now local.

df_all_genotypes |> plot_segment_lengths()


# Population simulation --------------------------------------------------
# Example 1:
# generate off a set of 100 samples, generate 10 new samples
# pop_alpha 9 (~10% between sample relatedness)
# coi_r = 0.25, coi_p = 0.7 (COI mean of 1.256, ~80.4% proportion will be monoclonal)
# k_s = 0.5 (50% of genotypes will be recombinant)
pop1 <- sim_population(paste0("sample", seq(0,100,1)), 10, pop_alpha = 9, coi_r = 0.25, coi_p = 0.7, k_s = 0.5)

# New functionality for getting genotype table for a population (with multiple samples), basically a wrapper
# for get_all_genotypes()
df_pop1_genotypes <- pop1 |> get_population_genotypes()

# Plotting ancestry composition (the same function works for both single samples and populations of samples)
df_pop1_ancestry <- df_pop1_genotypes |>
  get_ancestry_composition()

df_pop1_ancestry

## By genotype
df_pop1_ancestry |>
  plot_ancestry_composition(by = "genotype") + labs(title = "Population 1")

## By sample
df_pop1_ancestry |>
  plot_ancestry_composition(by = "sample") + labs(title = "Population 1")

# Plotting segment length (the same function works for both single samples and populations of samples)
df_pop1_genotypes |>
  plot_segment_lengths() + labs(title = "Population 1")

# Example 2:
# generate off a set of 100 samples, generate 50 new samples
# pop_alpha 9 (~10% between sample relatedness)
# coi_r = 0.25, coi_p = 0.7 (COI mean of 1.256, ~80.4% proportion will be monoclonal)
# k_s = 0.5 (50% of genotypes will be recombinant)
pop2 <- sim_population(paste0("sample", seq(0,100,1)), 50, pop_alpha = 9, coi_r = 0.25, coi_p = 0.7, k_s = 0.5)

df_pop2_genotypes <- pop2 |> get_population_genotypes()

# Ancestry composition
df_pop2_ancestry <- df_pop2_genotypes |>
  get_ancestry_composition()
df_pop2_ancestry

## By genotype
df_pop2_ancestry |>
  plot_ancestry_composition("genotype") + labs(title = "Population 2")

## By sample
df_pop2_ancestry |>
  plot_ancestry_composition(by = "sample") + labs(title = "Population 2")

# Segment lengths
df_pop2_genotypes |>
  plot_segment_lengths() + labs(title = "Population 2")

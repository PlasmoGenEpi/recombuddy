# deploy_recombuddy.R
#
# Author: Bob Verity
# Date: 2025-03-20
#
# Purpose:
# A sandbox to demonstrate the use of the recombuddy package.
#
# ------------------------------------------------------------------

# define the number of samples in the initial sample set, and in the final
# population that we are producing
n_set <- 10
n_pop <- 20

# define concentration parameter of the Dirichlet distribution
f <- 0.1
alpha <- (1 - f) / f

# draw proportions of each member of the sample set
set_props <- rdirichlet_single(n_set, alpha = alpha)

# simulate a single sample
sim <- sim_sample(k = c(0, 1, 0), rho = 6e-7, set_props = set_props)

# get all genotypes into a single table
df_all_genotypes <- get_all_genotypes(sim)

plot_genotypes(df_all_genotypes)



# draw positions of loci along the genome
df_loci <- draw_loci(n_loci = 2e3)

# draw population-level allele frequencies (PLAFs) at specified loci
df_PLAF <- draw_PLAF(df_loci = df_loci)

# draw alleles for each member of sample set
df_sample_set <- draw_sample_set_WSAF(df_PLAF = df_PLAF, n_set = n_set)

# draw proportions for each genotype
df_genotype_props <- draw_genotype_props(MOI = MOI, alpha = 10)

# calculate true within-sample allele frequencies (WSAFs) for the entire sample
df_WSAF <- get_WSAF(df_all_genotypes = df_all_genotypes,
                    df_sample_set = df_sample_set,
                    df_genotype_props = df_genotype_props)

# apply a simple read count model
df_counts <- draw_read_counts(df_WSAF = df_WSAF)

# plot read counts
plot_read_counts(df_counts = df_counts)

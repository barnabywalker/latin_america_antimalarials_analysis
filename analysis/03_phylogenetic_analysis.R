#' Phylogenetic analysis using techniques commonly used
#' in community phylogenetics.
#' 
#' This analysis calculates measures of overall relatedness
#' (mean pairwise distance and mean nearest taxon distance),
#' phylogenetic signal (Pagel's lambda, Blomberg's K, and the D-statistic),
#' and over/under representation at nodes in the tree (nodesig).

# libraries ----

library(here)   # handle data
library(dplyr)  # manipulate data
library(readr)  # read and write text files
library(tidyr)  # reshape data
library(purrr)  # map functions
library(tibble) # table formatting
library(ape)    # handle phylogenies
library(PhyloMeasures)  # fast calculation of mpd and mntd
library(ggtree)  # tidy handling of phylogenetic trees
# accessed but not loaded:
# phytools - calculate phylogenetic signal (lambda and K)
# caper - calculate phylogenetic signal (D)

source(here("R/helper_functions.R"))

# load data and tree ----
d <- read_csv(here("output/antimalarial_genus_counts.csv"))
tree <- read.tree(here("output/genus_level_tree.txt"))

# set some constants ----

angiosperms <- c("Magnoliopsida",
                 "Liliopsida")
# re-format some of the data ----

angiosperm_data <-
  d %>%
  filter(class %in% angiosperms, 
         genus %in% tree$tip.label)

# get the number of antimalarials in each genus
observed_counts <- 
  angiosperm_data %>%
  select(genus, n=n_antimalarial)

# get just the overall number of species in each genus
expected_counts <- 
  angiosperm_data %>%
  select(genus, n=n_species)

# remove tips from the tree that aren't in our data  
angiosperm_tree <- subset_tree(tree, angiosperm_data$genus)

# calculate MPD and MNTD for count data ----

observed_mat <-
  observed_counts %>%
  pivot_wider(names_from=genus, values_from=n) %>%
  as.matrix()

expected_mat <-
  expected_counts %>%
  pivot_wider(names_from=genus, values_from=n) %>%
  as.matrix()

# calculate the phylogenetic distance matrix
dis <- cophenetic(angiosperm_tree)

# calculate the observed metrics
mpd_obs <- mpd(observed_mat, dis, abundance.weighted=TRUE, self.distances=TRUE)
mntd_obs <- mntd(observed_mat, dis, abundance.weighted=TRUE)

# sampling for the null models:
# need to draw samples same size as antimalarial species pool from the total phylogenetic tree
# but each genus has a probability of being in the sample proportional to the total number
# of species in that genus

# total number of species to calculate the sampling probabilities
total_species <- sum(expected_counts$n)
# total number of antimalarial species as the number of species to sample
total_antimalarial <- sum(observed_counts$n)

# probability of a species being in each genus
sample_weights <- expected_mat / total_species
names(sample_weights) <- colnames(expected_mat)

# draw null models
null_mpd <- 
  null_distribution(sample_weights, total_antimalarial, dis, mpd, n_repeats=1000)

null_mntd <- 
  null_distribution(sample_weights, total_antimalarial, dis, mntd, n_repeats=1000)

# put all the results in a single table
phylo_diversity_results <-
  tibble(mpd=null_mpd,
         mntd=null_mntd) %>%
  gather(measure, value) %>%
  group_by(measure) %>%
  summarise(null_mean=mean(value),
            null_std=sd(value),
            .groups="drop") %>%
  # calling the standardised metrics "ses" as a placeholder for NTI(MNTD) and NRI(MPD)
  mutate(obs=c(mntd_obs, mpd_obs),
         ses=-1*(obs - null_mean) / null_std)

# save results table
write_csv(phylo_diversity_results, here("output/phylo_diversity_results.csv"))

# calculate MPD and MNTD for binary data ----

observed_bin <- observed_mat > 0

# don't need to do our own sampling because this is already implemented
binary_phylo_results <-
  tribble(
  ~measure, ~obs, ~ses,
  "mntd", mntd.query(angiosperm_tree, observed_bin), 
    mntd.query(angiosperm_tree, observed_bin, standardize=TRUE)*-1,
  "mpd", mpd.query(angiosperm_tree, observed_bin), 
    mpd.query(angiosperm_tree, observed_bin, standardize=TRUE)*-1
  )

write_csv(binary_phylo_results, here("output/phylo_diversity_binary_results.csv"))

# calculate phylogenetic signal ----

# first try treating the proportion of antimalarial species as a continuous trait
proportions <- angiosperm_data$n_antimalarial / angiosperm_data$n_species
names(proportions) <- angiosperm_data$genus

# calculate Blomberg's K and Pagel's lambda
signal_K <- phytools::phylosig(angiosperm_tree, proportions, method="K", test=TRUE)
signal_lambda <- phytools::phylosig(angiosperm_tree, proportions, method="lambda", test=TRUE)

# then use an arcsine transformation on the trait to make it actually continuous
arcsin_proportions <- asin(sqrt(proportions))

signal_K_transformed <- phytools::phylosig(angiosperm_tree, arcsin_proportions, method="K", test=TRUE)
signal_lambda_transformed <- phytools::phylosig(angiosperm_tree, arcsin_proportions, method="lambda", test=TRUE)

# then treat it as a binary trait
observed_bin <-
  observed_counts %>%
  mutate(antimalarial=n > 0) %>%
  as.data.frame()

# need to replace the node names with numbers so lots don't just have a blank name
new_tree <- angiosperm_tree
node_numbers <-
  angiosperm_tree$node.label %>% 
  enframe(name="node", value="label")

new_tree$node.label <- as.character(node_numbers$node)

# and then calculate D statistic
signal_D <- caper::phylo.d(observed_bin, new_tree, names.col=genus, binvar=antimalarial)

# put all measures together in table. D has two p-values
signal_table <-
  tribble(
    ~metric, ~transformation, ~value, ~pvalue_random, ~pvalue_brownian,
    "K", "none", signal_K$K, signal_K$P, NA_real_,
    "lambda", "none", signal_lambda$lambda, signal_lambda$P, NA_real_,
    "K", "arcsine", signal_K_transformed$K, signal_K_transformed$P, NA_real_,
    "lambda", "arcsine", signal_lambda_transformed$lambda, signal_lambda_transformed$P, NA_real_,
    "D", "none", signal_D$DEstimate, signal_D$Pval0, signal_D$Pval1
  )

# save the table
write_csv(signal_table, here("output/phylogenetic_signal_results.csv"))

# calculate "hot nodes" ----

# calculate richness at nodes, as per nodesig
node_richness <- calculate_node_richness(angiosperm_tree, observed_counts, 
                                         expected_counts, n_samples=1000)

# save whole output
write_csv(node_richness, here("output/node_richness_results_all.csv"))


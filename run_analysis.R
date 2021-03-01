#' Script to run the whole analysis for the paper:
#'   "Plants used traditionally as antimalarials in Latin America: 
#'    mining the Tree of Life for potential new medicines"
library(here)

# prepare the names and phylogeny ----
source(here("analysis/01_prepare_names.R"))

# calculate effect size compared to binomial distributions ----
source(here("analysis/02_binomial_analysis.R"))

# calculate phylogenetic diversity, signal, and over-representation ----
source(here("analysis/03_phylogenetic_analysis.R"))

# make figures and table data ----
source(here("analysis/04_make_figures.R"))


#' Use binomial tests to identify family and genera that are
#' over or under represented in terms of anti-malarial species.
#' 
#' These tests calculate the effect size of the difference from 
#' a null distribution - the distribution you would expect if all
#' species had an equal probability of being anti-malarial, regardless
#' of their genus.
#' 
#' Calculates the data underlying funnel plots for the manuscript.

# libraries ----

library(here)   # handle file paths
library(dplyr)  # manipulate data
library(readr)  # read and write text files
library(tidyr)  # reshape data

# load in aggregated medicinal data ----
d <- read_csv(here("output/antimalarial_genus_counts.csv"))
head(d)

# set some constants ----
seedplants <- c("Magnoliopsida",
                "Liliopsida",
                "Gnetopsida",
                "Pinopsida",
                "Cycadopsida")

# genus-wise effect sizes ---- 

genus_effects <-
  d %>%
  filter(class %in% seedplants) %>%
  mutate(mean=sum(n_antimalarial) / sum(n_species),
         std_err=sqrt(mean * (1 - mean) / n_species),
         eff=((n_antimalarial/n_species) - mean) / std_err) %>%
  arrange(desc(eff))

write_csv(genus_effects, here("output/binomial_effects_genus.csv"))  

# family-wise effect sizes ----

family_effects <-
  d %>%
  filter(class %in% seedplants) %>%
  group_by(family) %>%
  summarise(n_antimalarial=sum(n_antimalarial), 
            n_species=sum(n_species),
            .groups="drop") %>%
  mutate(mean=sum(n_antimalarial) / sum(n_species),
         std_err=sqrt(mean * (1 - mean) / n_species),
         eff=((n_antimalarial/n_species) - mean) / std_err) %>%
  arrange(desc(eff))

write_csv(family_effects, here("output/binomial_effects_family.csv"))  

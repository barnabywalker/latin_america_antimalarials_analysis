#' Plot figures and collate data to go with the manuscript.
#' 
#' Figures:
#'   1. (A) Genus- and (B) family-level funnel plots.
#'   2. Phylogenetic tree displaying anti-malarial and "hot" nodes.
#'   3. Subtrees showing anti-malarial and "hot" nodes for 
#'      (A) Rubiaceae and (B) Apocynaceae
#'   4. Subtrees showing anti-malarial and "hot" nodes for 
#'      (A) Menispermaceae and (B) Simaroubaceae
#'
#' Tables:
#'   1. Something putting everything together.

# TODO: add the final labels to funnels

# libraries ----
library(here)
library(readr)
library(writexl)
library(dplyr)
library(tidyr)
library(ape)
library(ggplot2)
library(ggtree)
library(ggrepel)
library(patchwork)

source(here("R/helper_functions.R"))
source(here("R/plotting_functions.R"))

# plotting constants ----
point_colours <- c("Under"="#bb2a34", 
                   "Over"="#2870b1", 
                   "Expected"="#a9a9a9")

funnel_linestyles <- c("Mean"=1,
                       "95% CI"=4, 
                       "99% CI"=3)

angiosperms <- c("Magnoliopsida",
                 "Liliopsida")

# load data ---
genus_effects <- read_csv(here("output/binomial_effects_genus.csv"))
family_effects <- read_csv(here("output/binomial_effects_family.csv"))

node_richness <- read_csv(here("output/node_richness_results_all.csv"))

phylo_diversity <- read_csv(here("output/phylo_diversity_results.csv"))
phylo_diversity_binary <- read_csv(here("output/phylo_diversity_binary_results.csv"))
phylo_signal <- read_csv(here("output/phylogenetic_signal_results.csv"))

tree <- read.tree(here("output/genus_level_tree.txt"))

# Fig 1. funnel plots ----

# genus-level
genus_labels <- c(
  "Iriartella",
  "Anadenanthera",
  "Geissospermum",
  "Aspidosperma",
  "Euterpe",
  "Ladenbergia",
  "Cyathula",
  "Cilosemina",
  "Uncaria",
  "Oreocallis",
  "Iriartea",
  "Solanum",
  "Piper",
  "Croton",
  "Aristolochia",
  "Momocordia",
  "Cinchona",
  "Stenna",
  "Stachytarpheta"
)

genus_rep <- 
  genus_effects %>%
  mutate(representation=case_when(eff <= -3 ~ "Under",
                                  eff >= 3 ~ "Over",
                                  TRUE ~ "Expected"),
         label=ifelse(genus %in% genus_labels, genus, NA_character_)) %>%
  mutate(representation=factor(representation, levels=c("Under", "Expected", "Over")))

genus_funnel <- make_funnel_data(genus_rep$n_antimalarial, genus_rep$n_species)

genus_funnel_plot <-
  plot_funnel(genus_rep, genus_funnel, label_italics=TRUE) +
  labs(x="Number of species", y="Proportion with antimalarial uses") +
  scale_colour_manual(values=point_colours,
                      name="Representation") +
  scale_linetype_manual(values=funnel_linestyles,
                        name="Expected distribution",
                        guide=guide_legend(override.aes=list(size=c(0.5,0.5,0.5)))) +
  scale_x_log10(labels=scales::label_comma()) +
  scale_y_continuous(limits=c(-0.05, 1)) +
  guides(colour=guide_legend(title.position="top"),
         linetype=guide_legend(title.position="top")) +
  theme_funnel()

ggsave(here("figures/funnel_plot_genus.svg"), plot=genus_funnel_plot, width=8, height=4.5)

# family-level
family_labels <- c(
  "Rubiaceae",
  "Apocynaceae",
  "Meliaceae",
  "Simaroubaceae",
  "Asteraceae",
  "Bignoniaceae",
  "Menispermaceae",
  "Fabaceae",
  "Solanceae",
  "Arecaceae",
  "Euphorbiaceae"
)

family_rep <- 
  family_effects %>%
  mutate(representation=case_when(eff <= -3 ~ "Under",
                                  eff >= 3 ~ "Over",
                                  TRUE ~ "Expected"),
         label=ifelse(family %in% family_labels, family, NA_character_)) %>%
  mutate(representation=factor(representation, levels=c("Under", "Expected", "Over"))) 
  
family_funnel <- make_funnel_data(family_rep$n_antimalarial, family_rep$n_species)

family_funnel_plot <-
  plot_funnel(family_rep, family_funnel) +
  labs(x="Number of species", y="Proportion with antimalarial uses") +
  scale_colour_manual(values=point_colours, 
                      name="Representation") +
  scale_linetype_manual(values=funnel_linestyles,
                        name="Expected distribution") +
  scale_x_log10(labels=scales::label_comma(accuracy=1)) +
  scale_y_continuous(limits=c(-0.02, 0.2)) +
  guides(colour=guide_legend(title.position="top"),
         linetype=guide_legend(title.position="top")) +
  theme_funnel()

ggsave(here("figures/funnel_plot_family.svg"), plot=family_funnel_plot, width=8, height=4.5)

# join together
funnel_figure <- 
  (genus_funnel_plot / family_funnel_plot) + 
  plot_layout(guides="collect") +
  plot_annotation(tag_levels="A")

# save then edit labels in editing software
ggsave(here("figures/funnel_plot_figure.svg"),
       plot=funnel_figure,
       height=7,
       width=8)

# Fig 2. overall tree ----

# get a tree of just the angiosperms
angiosperm_genera <-
  genus_effects %>%
  filter(class %in% angiosperms) %>%
  pull(genus)

angiosperm_tree <- subset_tree(tree, angiosperm_genera)

node_reps <-
  node_richness %>%
  select(node, rank_richness, obs_richness) %>%
  mutate(overrepresented=rank_richness >= 0.975,
         antimalarial_not_hot=rank_richness < 0.975 & obs_richness > 0)

hot_node_tree <-
  ggtree(angiosperm_tree, mapping=aes(colour=overrepresented, alpha=overrepresented), layout="circular") %<+% 
  node_reps +
  geom_cladelabel(8680, "Magnoliidae", barsize=1, colour="grey30") +
  geom_cladelabel(8220, "Asparagales", barsize=1, colour="grey50") +
  geom_cladelabel(4804, "Malpighiales", hjust=1, barsize=1, colour="grey50") +
  geom_cladelabel(5191, "Rosales", hjust=1, barsize=1, colour="grey30") +
  geom_cladelabel(5280, "Fabales", hjust=1, barsize=1, colour="grey50") +
  geom_cladelabel(5580, "Brassicales", hjust=1, barsize=1, colour="grey30") +
  geom_cladelabel(5847, "Myrtales", hjust=1, barsize=1, colour="grey50") +
  geom_cladelabel(6364, "Gentianales", barsize=1, colour="grey50") +
  geom_cladelabel(6066, "Lamiales", barsize=1, colour="grey50") +
  geom_cladelabel(6711, "Asterales", barsize=1, colour="grey30") +
  geom_cladelabel(7305, "Ericales", hjust=1, barsize=1, colour="grey50") +
  geom_cladelabel(7448, "Caryophyllales", hjust=1, barsize=1, colour="grey30") +
  geom_cladelabel(7912, "Poales", barsize=1, colour="grey30") +
  geom_cladelabel(6621, "Solanales", barsize=1, colour="grey30") +
  scale_colour_manual(values=c(`TRUE`="#ff0d57", `FALSE`="#a9a9a9")) +
  scale_alpha_manual(values=c(`TRUE`=1, `FALSE`=0)) +
  guides(alpha=FALSE, colour=FALSE)

just_antimalarial_tree <-
  ggtree(angiosperm_tree, mapping=aes(colour=antimalarial_not_hot, alpha=antimalarial_not_hot), layout="circular") %<+% 
  node_reps +
  geom_cladelabel(8680, "Magnoliidae", barsize=1, colour="grey30") +
  geom_cladelabel(8220, "Asparagales", barsize=1, colour="grey50") +
  geom_cladelabel(4804, "Malpighiales", hjust=1, barsize=1, colour="grey50") +
  geom_cladelabel(5191, "Rosales", hjust=1, barsize=1, colour="grey30") +
  geom_cladelabel(5280, "Fabales", hjust=1, barsize=1, colour="grey50") +
  geom_cladelabel(5580, "Brassicales", hjust=1, barsize=1, colour="grey30") +
  geom_cladelabel(5847, "Myrtales", hjust=1, barsize=1, colour="grey50") +
  geom_cladelabel(6364, "Gentianales", barsize=1, colour="grey50") +
  geom_cladelabel(6066, "Lamiales", barsize=1, colour="grey50") +
  geom_cladelabel(6711, "Asterales", barsize=1, colour="grey30") +
  geom_cladelabel(7305, "Ericales", hjust=1, barsize=1, colour="grey50") +
  geom_cladelabel(7448, "Caryophyllales", hjust=1, barsize=1, colour="grey30") +
  geom_cladelabel(7912, "Poales", barsize=1, colour="grey30") +
  geom_cladelabel(6621, "Solanales", barsize=1, colour="grey30") +
  scale_colour_manual(values=c(`TRUE`="#1E88E5", `FALSE`="#a9a9a9")) +
  scale_alpha_manual(values=c(`TRUE`=1, `FALSE`=0)) +
  guides(alpha=FALSE, colour=FALSE)

background_tree <-
  ggtree(angiosperm_tree, colour="#a9a9a9", layout="circular") +
  geom_cladelabel(8680, "Magnoliidae", barsize=1, colour="grey30") +
  geom_cladelabel(8220, "Asparagales", barsize=1, colour="grey50") +
  geom_cladelabel(4804, "Malpighiales", hjust=1, barsize=1, colour="grey50") +
  geom_cladelabel(5191, "Rosales", hjust=1, barsize=1, colour="grey30") +
  geom_cladelabel(5280, "Fabales", hjust=1, barsize=1, colour="grey50") +
  geom_cladelabel(5580, "Brassicales", hjust=1, barsize=1, colour="grey30") +
  geom_cladelabel(5847, "Myrtales", hjust=1, barsize=1, colour="grey50") +
  geom_cladelabel(6364, "Gentianales", barsize=1, colour="grey50") +
  geom_cladelabel(6066, "Lamiales", barsize=1, colour="grey50") +
  geom_cladelabel(6711, "Asterales", barsize=1, colour="grey30") +
  geom_cladelabel(7305, "Ericales", hjust=1, barsize=1, colour="grey50") +
  geom_cladelabel(7448, "Caryophyllales", hjust=1, barsize=1, colour="grey30") +
  geom_cladelabel(7912, "Poales", barsize=1, colour="grey30") +
  geom_cladelabel(6621, "Solanales", barsize=1, colour="grey30")

# save trees, then use editing software to overlay all three
ggsave(here("figures/all_tree.svg"), plot=background_tree, height=6, width=6)
ggsave(here("figures/hot_node_tree.svg"), plot=hot_node_tree, height=6, width=6)
ggsave(here("figures/just_antimalarial_tree.svg"), plot=just_antimalarial_tree, height=6, width=6)

# Fig 3. rubiaceae and apocynaceae subtrees ----

# rubiaceae subtree
node <- 6365

rube_data <- make_subtree_data(node_reps, angiosperm_tree, node)

rube_subtree <- subset_tree(angiosperm_tree, rube_data$label)

rube_plot <-
  ggtree(rube_subtree, 
         layout="circular",
         mapping=aes(colour=flag)) %<+%
  rube_data +
  scale_colour_manual(values=c(antimalarial="#1E88E5", nothing="#a9a9a9", hot="#ff0d57"), name="") +
  geom_tiplab2(size=2, show.legend=FALSE) +
  theme(legend.position="bottom")

ggsave(here("figures/rubiaceae_tree.svg"), rube_plot)

# apocynaceae subtree
node <- 6525
apoc_data <- make_subtree_data(node_reps, angiosperm_tree, node)

apoc_subtree <- subset_tree(angiosperm_tree, apoc_data$label)

apoc_plot <-
  ggtree(apoc_subtree, 
         layout="circular",
         mapping=aes(colour=flag)) %<+%
  apoc_data +
  scale_colour_manual(values=c(antimalarial="#1E88E5", nothing="#a9a9a9", hot="#ff0d57"), name="") +
  geom_tiplab2(size=2, show.legend=FALSE) +
  theme(legend.position="bottom")

ggsave(here("figures/apocynaceae_tree.svg"), apoc_plot)

# join them
rube_apoc_plot <-
  (rube_plot / apoc_plot) + 
  plot_annotation(tag_levels="A") + 
  plot_layout(guides="collect") & 
  theme(legend.position = "bottom")

ggsave(here("figures/figure5.svg"), rube_apoc_plot, height=10, width=8)

# Fig 4. ranunculales and sapindales subtrees ----

# ranunc subtree
node <- 7766

ranunc_data <- make_subtree_data(node_reps, angiosperm_tree, node)

ranunc_subtree <- subset_tree(angiosperm_tree, ranunc_data$label)

ranunc_plot <-
  ggtree(ranunc_subtree, 
         layout="circular",
         mapping=aes(colour=flag)) %<+%
  ranunc_data +
  scale_colour_manual(values=c(antimalarial="#1E88E5", nothing="#a9a9a9", hot="#ff0d57"), name="") +
  geom_tiplab2(size=2, show.legend=FALSE) +
  theme(legend.position="bottom")

ggsave(here("figures/ranunculales_tree.svg"), ranunc_plot)

# sapindales subtree
node <- 5779
sapind_data <- make_subtree_data(node_reps, angiosperm_tree, node)

sapind_subtree <- subset_tree(angiosperm_tree, sapind_data$label)

sapind_plot <-
  ggtree(sapind_subtree, 
         layout="circular",
         mapping=aes(colour=flag)) %<+%
  sapind_data +
  scale_colour_manual(values=c(antimalarial="#1E88E5", nothing="#a9a9a9", hot="#ff0d57"), name="") +
  geom_tiplab2(size=2, show.legend=FALSE) +
  theme(legend.position="bottom")

ggsave(here("figures/sapindales_tree.svg"), sapind_plot)

# join them
ranunc_sapind_plot <-
  (ranunc_plot / sapind_plot) + 
  plot_annotation(tag_levels="A") + 
  plot_layout(guides="collect") & 
  theme(legend.position = "bottom")

ggsave(here("figures/figure6.svg"), ranunc_sapind_plot, height=10, width=8)

# Table 1. Compilation of results ----

# comparison of binomial effects to node richness
tip_richness <-
  node_richness %>%
  filter(! node %in% as_tibble(angiosperm_tree)$parent) %>%
  select(label, ses_richness, rank_richness)

genus_comparison <-
  genus_effects %>%
  select(-class) %>%
  inner_join(
    tip_richness,
    by=c("genus"="label")
  ) %>%
  mutate(hot_binomial=eff >= 3,
         hot_nodesig=rank_richness >= 0.975) %>%
  arrange(desc(eff)) %>%
  select("total species"="n_species",
         "antimalarial species"="n_antimalarial",
         "binomial effect size"="eff",
         "nodesig effect size"="ses_richness",
         "nodesig p-value"="rank_richness",
         "hot node (binomial)"="hot_binomial",
         "hot node (nodesig)"="hot_nodesig")

# clean up diversity results
diversity <-
  phylo_diversity %>%
  select(measure, obs, ses) %>%
  mutate(variable="count") %>%
  bind_rows(
    phylo_diversity_binary %>%
      mutate(variable="binary")
  ) %>%
  select("Measure"="measure", "Variable"="variable", 
         "Observed"=obs, "Standardised effect"=ses)

# write everything to an excel file
write_xlsx(
  list(
    "binomial effects (genus)"=genus_effects,
    "binomial effects (family)"=family_effects,
    "node richness (nodesig)"=node_richness,
    "phylo diversity measures"=diversity,
    "phylo signal"=phylo_signal
  ),
  here("output/analysis_results.xlsx")
)

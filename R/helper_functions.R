# data request functions ----

#' Utility to check if a taxon name is in a lookup of IDs.
check_id <- function(name, id) {
  if (name %in% names(missing_names)) {
    id <- missing_names[[name]]
  }
  
  id
}

#' Utility to get accepted name information for an IPNI ID.
get_accepted_info <- function(ipni_id) {
  
  # empty table to return if nothing comes back
  info <- tibble(
    status=NA_character_,
    accepted_id=NA_character_,
    accepted_name=NA_character_,
    accepted_authors=NA_character_,
    accepted_rank=NA_character_,
    accepted_species=NA_character_,
    accepted_species_id=NA_character_
  )
  
  if (is.na(ipni_id)) {
    return(list(info))
  }
  
  # get WCVP entry for ID
  r <- kewr::lookup_wcvp(ipni_id)
  
  info$status <- r$status
  
  # extract info about accepted name
  if (info$status %in% c("accepted", "unplaced")) {
    info$accepted_id <- r$id
    info$accepted_name <- r$name
    info$accepted_authors <- r$author
    info$accepted_rank <- r$rank
  } else {
    info$accepted_id <- r$accepted$id
    info$accepted_name <- r$accepted$name
    info$accepted_authors <- r$accepted$authors
    info$accepted_rank <- r$accepted$rank
  }
  
  # handle infraspecifics by getting the parent name
  if (info$accepted_rank %in% c("Variety", "Subspecies")) {
    r <- kewr::lookup_wcvp(info$accepted_id)
    info$accepted_species <- r$parent$name
    info$accepted_species_id <- r$parent$id
  } else if (info$accepted_rank == "Species") {
    info$accepted_species <- info$accepted_name
    info$accepted_species_id <- info$accepted_id
  }
  
  # wrap in a list to make it interface with dplyr's rowwise nicely
  list(info)
}

#' Check if a species is native to a study region by looking up
#' the distribution in POWO.
check_native <- function(ipni_id, study_region) {
  if (is.na(ipni_id)) {
    return(NA)
  }
  Sys.sleep(0.2)
  r <- kewr::lookup_powo(ipni_id, distribution=TRUE)
  
  native_dist <- r$distribution$natives
  in_study_region <- purrr::map_lgl(native_dist, ~.x$tdwgCode %in% study_region)
  
  any(in_study_region)
}

#' Get a checklist of accepted species in a WGSPRD region.
get_region_species <- function(region, limit=1000) {
  query <- list(distribution=region)
  filters <- c("accepted", "species")
  
  r <- search_powo(query, filters=filters, limit=1)
  total_species <- r$total
  n_pages <- ceiling(total_species / limit)
  
  pb <- progress_bar$new(
    format=paste(region, "[:bar] [:current]/[:total]"),
    total=n_pages
  )
  
  pb$tick()
  results <- list(search_powo(query, filters=filters, limit=limit))
  for (i in 1:(n_pages-1)) {
    pb$tick()
    Sys.sleep(0.2)
    results[[i+1]] <- request_next(results[[i]])
  }
  
  purrr::map_dfr(results, format)
}

# phylogenetics functions ----

#' Drop all nodes from a phylogenetic tree except those in a list.
subset_tree <- function(tree, node_list) {
  drop_list <- tree$tip.label[! tree$tip.label %in% node_list]
  
  drop.tip(tree, drop_list)
}

#' Calculate the mean pairwise distance for a community.
#' 
#' This implementation accepts abundances in the sample,
#' as well as different abundance weights. 
#' 
#' Includes optional calculation between species at tips, to include
#' pairwise distance for species in the same genus/family
#' when using a genus-level or family-level phylogeny.
mpd <- function(samp, dis, abundance.weighted=FALSE, self.distances=FALSE) {
  
  tips_in_sample <- samp[, samp > 0, drop=FALSE]
  tip_names <- colnames(tips_in_sample)
  
  sample_distances <- dis[tip_names, tip_names]
  
  if (abundance.weighted) {
    sample_weights <- t(as.matrix(tips_in_sample)) %*% as.matrix(tips_in_sample)
    if (self.distances) {
      diag(sample_weights) <- choose(as.matrix(tips_in_sample), 2)
    } 
    
    mpd <- weighted.mean(sample_distances, sample_weights)
  } else {
    mpd <- mean(sample_distances[lower.tri(sample_distances)])
  }
  
  return(mpd)
}

#' Calculate the mean nearest taxon distance for a community.
#' 
#' This implementation accepts abundances in the sample,
#' as well as different abundance weights.
mntd <- function(samp, dis, abundance.weighted = TRUE) {
  
  tips_in_sample <- samp[, samp > 0, drop=FALSE]
  tip_names <- colnames(tips_in_sample)
  
  sample_distances <- dis[tip_names, tip_names]
  
  if (abundance.weighted) {
    diag(sample_distances) <- apply(tips_in_sample, 2, function(x) ifelse(x == 1, NA, x))
    mntds <- apply(sample_distances, 2, min, na.rm=TRUE)
    mntd <- weighted.mean(mntds, tips_in_sample)
  } else {
    diag(sample_distances) <- NA
    mntd <- mean(apply(sample_distances, 2, min, na.rm=TRUE))
  }
  
  return(mntd)
}

#' Sample random taxa from a community with an optional probability weighting.
sample_taxa <- function(taxa, sample_size=NULL, probabilites=NULL) {
  if (is.null(sample_size)) {
    sample_size <- length(taxa)
  }
  
  sampled_taxa <- sample(taxa, size=sample_size, prob=probabilites, replace=TRUE)
  
  tibble(taxa_group=sampled_taxa) %>%
    count(taxa_group) %>%
    spread(taxa_group, n) %>%
    as.data.frame()
}

#' Generate a null distribution of for a particular measure.
null_distribution <- function(class_probabilities, sample_count, distance_matrix, summary_function, n_repeats=1000) {
  class_names <- names(class_probabilities)
  purrr::map_dbl(1:n_repeats,
          ~summary_function(sample_taxa(class_names, sample_count, class_probabilities), distance_matrix, abundance.weighted=TRUE))
}

#' Draw sampled counts from an expected distribution.
#' 
#' Draws samples of the same total size as the observed distribution,
#' with the probability of each class drawn according to the expected distribution.
sample_counts <- function(observed_counts, expected_counts, n_samples) {
  total_observed_count <- sum(observed_counts$n)
  
  class_probabilities <- 
    expected_counts %>% 
    mutate(p = n / sum(n)) %>% 
    pull(p)
  
  class_names <- expected_counts$genus
  
  n_samples %>%
    purrr::rerun(sample_taxa(class_names, total_observed_count, class_probabilities)) %>%
    bind_rows() %>%
    mutate(sample = 1:nrow(.)) %>%
    gather(genus, n, -sample) %>%
    replace_na(list(n=0))
}

#' Tally the richness of species at each node in a phylogenetic tree.
tally_richness <- function(richness, node_number, tree_df) {
  node_name <- tree_df[tree_df$node == node_number,]$label
  
  descendants <- tidytree::offspring(tree_df, node_number)
  descendant_tips <- filter(descendants, ! node %in% tree_df$parent)
  
  # handling tallying the richness for tip nodes
  if (nrow(descendant_tips) == 0 & node_name %in% richness$genus) {
    descendant_tips <- filter(tree_df, node == node_number)
  }
  
  descendant_richness <- filter(richness, genus %in% descendant_tips$label)
  
  # handling tallying richness with more than one sample
  if (nrow(descendant_richness) > nrow(descendant_tips)) {
    descendant_richness <- group_by(descendant_richness, sample)
  }
  
  descendant_richness %>%
    summarise(n=sum(n),
              .groups="drop") %>%
    mutate(label=node_name, node=node_number)
}

#' Calculate the standardised residual of an observed value from simulated values.
#' 
#' Calculated as the difference of the observed value from the mean simulated value,
#' divided by the standard deviation of the simulated values.
standardise_residual <- function(observed_values, simulated_values, var) {
  var <- enquo(var)
  
  simulated_values %>%
    group_by(node) %>%
    summarise(mu_sim=mean(!!var), sd_sim=sd(!!var),
              .groups="drop") %>%
    right_join(observed_values, by="node") %>%
    rename(obs = !!var) %>%
    mutate(standardised_residual=(obs - mu_sim) / sd_sim)
}

#' Calculate the proportional rank of observed values compared to simulated values.
#' 
#' The proportional rank can be used as a measure of the probability of observing
#' a particular value.
calculate_proportional_rank <- function(observed_values, simulated_values, var) {
  var <- enquo(var)
  
  observed_values %>%
    mutate(sample=0) %>%
    bind_rows(simulated_values) %>%
    group_by(node) %>%
    mutate(proportional_rank=rank(!!var) / n()) %>%
    filter(sample == 0) %>%
    select(label, node, proportional_rank)
}

#' Carry out an analysis of the richness at nodes in the tree,
#' like nodesig.
calculate_node_richness <- function(tree, observed_counts, expected_counts, n_samples=10) {
 
  sampled_counts <- sample_counts(observed_counts, expected_counts, n_samples)
  
  tree_df <- as_tibble(tree)
  
  nodes <- tree_df
  
  observed_richness <- purrr::map_dfr(nodes$node, ~tally_richness(observed_counts, .x, tree_df))
  simulated_richness <- purrr::map_dfr(nodes$node, ~tally_richness(sampled_counts, .x, tree_df))
  
  richness_sr <- standardise_residual(observed_richness, simulated_richness, n)
  richness_pr <- calculate_proportional_rank(observed_richness, simulated_richness, n)
  
  richness_sr %>%
    rename(ses_richness=standardised_residual,
           obs_richness=obs,
           mu_sim_richness=mu_sim,
           sd_sim_richness=sd_sim) %>%
    full_join(richness_pr, by=c("label", "node")) %>%
    rename(rank_richness=proportional_rank)
}







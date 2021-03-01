# TODO: finish docs

#' Calculate a funnel representing the expected distribution
#' of anti-malarial species in groups, assuming a binomial 
#' distribution where each species has an equal probability of
#' being anti-malarial. Also returns the 95 % and 99 % confidence 
#' intervals.
make_funnel_data <- function(successes, totals) {

  expectation <- sum(successes) / sum(totals)
  n <- seq(from=min(totals), to=max(totals), by=1)
  std_err <- sqrt(expectation * (1 - expectation) / n)

  funnel <- tibble(
    n_species=n,
    std_err=std_err,
    expected=expectation,
    over_95=expectation + 1.96*std_err,
    over_99=expectation + 3*std_err,
    under_95=expectation - 1.96*std_err,
    under_99=expectation - 3*std_err
  )
  
  funnel %>%
    pivot_longer(expected:under_99, names_to="line", values_to="value") %>%
    mutate(value=ifelse(value < 0, NA_real_, value)) %>%
    mutate(representation=stringr::str_extract(line, "^[a-z]+"),
           confidence=stringr::str_extract(line, "(?<=_)\\d+")) %>%
    mutate(confidence=case_when(stringr::str_detect(confidence, "\\d+") ~ glue::glue("{confidence}% CI"),
                                is.na(confidence) ~ "Mean",
                                TRUE ~ confidence)) %>%
    mutate(confidence=factor(confidence,
                             levels=c("Mean", "95% CI", "99% CI")))
}

#' Utility to plot the proportion of species in each group, overlayed with
#' a funnel representing a binomial distribution with equal probability for all
#' species.
plot_funnel <- function(d, funnel, label_italics=FALSE) {
  if (label_italics) {
    face <- "italic"
  } else {
    face <- "plain"
  }

  d %>%
    mutate(p=n_antimalarial / n_species) %>%
    ggplot(mapping=aes(x=n_species, y=p)) +
    geom_point(mapping=aes(colour=representation), alpha=0.5) +
    geom_line(data=funnel, 
              mapping=aes(y=value, group=paste(representation, confidence), linetype=confidence), 
              size=1) +
    geom_text_repel(mapping=aes(colour=representation, label=label, y=p),
                    force=10, fontface=face, size=3, show.legend=FALSE)
}

theme_funnel <- function() {
  theme_bw() +
  theme(
    panel.border=element_blank(),
    panel.grid.major.x=element_blank(),
    panel.grid.minor.x=element_blank(),
    panel.grid.minor.y=element_blank(),
    panel.grid.major.y=element_line(linetype=3, colour="grey60"),
    legend.position="right",
        legend.key.width=grid::unit(2, "lines"))
}

make_funnel <- function(.data, success_var, total_var, colour_var=NULL, label_var=NULL) {
  success_var <- enquo(success_var)
  total_var <- enquo(total_var)
  colour_var <- enquo(colour_var)
  label_var <- enquo(label_var)
  
  successes <- pull(.data, !! success_var)
  totals <- pull(.data, !! total_var)
  
  funnel <- tibble(expected=sum(successes) / sum(totals),
                   !! total_var :=seq(min(totals), max(totals), 1))
  
  funnel <- mutate(funnel, std_err=sqrt(expected * (1 - expected) / !! total_var))
  
  funnel <- 
    funnel %>%
    mutate(over_95=expected + 1.96*std_err,
           over_99=expected + 3*std_err,
           under_95=expected - 1.96*std_err,
           under_99=expected - 3*std_err) %>%
    gather(line, value, -n_species, -std_err) %>%
    mutate(value=ifelse(value < 0, NA_real_, value)) %>%
    mutate(representation=stringr::str_extract(line, "^[a-z]+"),
           confidence=stringr::str_extract(line, "(?<=_)\\d+")) %>%
    mutate(confidence=case_when(stringr::str_detect(confidence, "\\d+") ~ glue::glue("{confidence}% CI"),
                                is.na(confidence) ~ "Mean",
                                TRUE ~ confidence)) %>%
    mutate(confidence=factor(confidence,
                             levels=c("Mean", "95% CI", "99% CI")))
  
  if (rlang::quo_is_null(colour_var)) {
    p <- 
      .data %>%
      mutate(p= !! success_var / !! total_var) %>%
      ggplot(mapping=aes(x=!! total_var, y=p)) +
      geom_point(color="#a9a9a9")  
  } else {
    p <-
      .data %>%
      mutate(p= !! success_var / !! total_var) %>%
      ggplot(mapping=aes(x=!! total_var, y=p)) +
      geom_point(mapping=aes(colour=!! colour_var), alpha=0.5)
  }
  
  
  p <-
    p +
    geom_line(data=funnel, mapping=aes(y=value, group=paste(representation, confidence), linetype=confidence), size=1)
  
  if (! rlang::quo_is_null(label_var)) {
    p <-
      p +
      geom_text_repel(mapping=aes(colour=!! colour_var, label=!! label_var, y=p), force=10, show.legend=FALSE)
  }
  
  p
}

#' Utility to construct a table of nodes in a subtree,
#' along with node richness-based representation data.
make_subtree_data <- function(data, tree, node) {
  node_tbl <-
    tree %>%
    tidytree::as_tibble() %>%
    left_join(data, by=c("node"="node")) %>%
    tidytree::offspring(node, self_include=TRUE) %>%
    select(-parent) %>%
    tibble::rowid_to_column(var="node_id")
  
  subtree <- subset_tree(tree, node_tbl$label)
  
  subtree_data <- 
    node_tbl %>%
    left_join(
      subtree %>%
        tidytree::as_tibble() %>%
        tibble::rowid_to_column(var="node_id"),
      by="node_id",
      suffix=c("_old", "")
    ) %>%
    select(node, label, overrepresented, obs_richness) %>%
    mutate(flag=case_when(overrepresented ~ "hot",
                          obs_richness > 0 ~ "antimalarial",
                          TRUE ~ "nothing"))
}

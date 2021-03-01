#' Prepare a list of taxon names for the analysis.
#' This involves:
#'   1. Clean the name list by removing special characters
#'      and any taxa only determined as "aff." or "cf.".
#'   2. Match the cleaned names to accepted species names in Kew's
#'      World Checklist of Vascular Plants and remove any not native
#'      to study area (Central and South America, excluding the Southern Cone).
#'   3. Download a checklist of accepted species in the study
#'      area.
#'   4. Count the number of total species and anti-malarial species in
#'      each genus.
#'   5. Trim the Smith and Brown phylogeny to match out list of genera.


# libraries ----
library(here)      # handle file paths
library(dplyr)     # manipulate data
library(readr)     # read text files
library(jsonlite)  # handle json format
library(tidyr)     # reshape data
library(purrr)     # map functions
library(ape)       # handle phylogenies
library(rgbif)     # get higher taxonomy from GBIF backbone
library(stringr)   # handle string data
library(kewr)      # request kew data
library(tibble)    # get data into nice tables
library(progress)  # make nice progress bars

source(here("R/helper_functions.R"))

# load taxon list ----

antimalarial_path <- here("data/antimalarial_taxa.csv")

antimalarial_taxa <- read_csv(antimalarial_path)

# define study region ----

# WGSRPD level 2 codes
l2_codes <- c(
  79, # Mexico
  80, # Central America
  81, # Caribbean
  82, # Northern South America
  83, # Western South America
  84  # Brazil
)

# WGSRPD level 3 codes
study_region <- c(
  # mexico
  "MXC", "MXE", "MXG", "MXI", "MXN", "MXS", "MXT",
  # central america
  "BLZ", "COS", "CPI", "ELS", "GUA", "HON", "NIC", "PAN",
  # caribbean
  "ARU", "BAH", "BER", "CAY", "CUB", "DOM", "HAI", "JAM", "LEE",
  "NLA", "PUE", "SWC", "TCI", "TRT", "VNA", "WIN",
  # northern south america
  "FRG", "GUY", "SUR", "VEN",
  # western south america
  "BOL", "CLM", "ECU", "GAL", "PER",
  # brazil
  "BZC", "BZE", "BZL", "BZN", "BZS"
)

# 1. clean taxon list ----

# rename columns and remove unnecessary ones
antimalarial_taxa <-
  antimalarial_taxa %>%
  rename(species=Species, family=Family) %>%
  select(species, family)

# remove cf. or aff. species
antimalarial_taxa <-
  antimalarial_taxa %>%
  filter(! str_detect(species, "(cf\\.|aff\\.)"))

# remove special characters
antimalarial_taxa <-
  antimalarial_taxa %>%
  mutate(species=str_remove_all(species, "(\\*|ø)")) %>%
  mutate(species=str_trim(species))

# remove duplicate taxa
antimalarial_taxa <-
  antimalarial_taxa %>%
  distinct(species, .keep_all=TRUE)

# 2. match names against WCVP ----
names_to_match <-
  antimalarial_taxa %>%
  pull(species)

full_matches <- match_knms(names_to_match)
full_matches

full_matches <- tidy(full_matches)

# show the taxa with no matches
unmatched <- filter(full_matches, ! matched)
unmatched

# resolve unmatched names using a manual matching file
missing_names <- read_json(here("data/name_match_missing.json"))

matched_names <-
  full_matches %>%
  mutate(ipni_id=map2_chr(submitted, ipni_id, check_id))

# show the taxa with multiple matches
multiple_matches <-
  full_matches %>%
  add_count(submitted) %>%
  filter(n > 1)

multiple_matches

# resolve multiple matches with a manual matching file
match_resolutions <- read_json(here("data/name_match_multiples.json"))

matched_names <-
  matched_names %>%
  filter((! submitted %in% names(match_resolutions)) | ipni_id %in% match_resolutions)

# manually fix a couple matches
match_correction <- read_json(here("data/manual_match_correction.json"))

matched_names <-
  matched_names %>%
  mutate(ipni_id=recode(ipni_id, !!! match_correction),
         matched_record=ifelse(ipni_id %in% names(match_correction),
                               NA_character_,
                               matched_record))

# get accepted name info for each match
accepted_names <-
  matched_names %>%
  nest_by(submitted, match_id=ipni_id) %>%
  mutate(info=get_accepted_info(match_id)) %>%
  select(-data) %>%
  unnest(cols=c(info)) %>%
  ungroup()

# flag hybrids
cleaned_names <-
  accepted_names %>%
  mutate(hybrid=str_detect(accepted_name, "×"))

# add indicator if in study region
cleaned_names <-
  cleaned_names %>%
  rowwise() %>%
  mutate(in_study_region=check_native(accepted_id, study_region)) %>%
  ungroup()

# extract genus from accepted name
cleaned_names <-
  cleaned_names %>%
  mutate(accepted_genus=str_extract(accepted_name, "^[A-Z][a-z\\-]+"))

# remove any citruses that have slipped through
cleaned_names <-
  cleaned_names %>%
  mutate(in_study_region=ifelse(accepted_genus == "Citrus", 
                                FALSE,
                                in_study_region))

# save for posterity
cleaned_names %>%
  mutate(removal_reason=case_when(hybrid ~ "hybrid",
                                  ! in_study_region ~ "not native to SA",
                                  is.na(match_id) ~ "not matched to accepted name",
                                  TRUE ~ NA_character_)) %>%
  mutate(included_in_analysis=is.na(removal_reason)) %>%
  write_csv(here("output/antimalarials_matched.csv"))

# 3. get full checklist of species ----

# get list of species for each region
mexico_checklist <- get_region_species("Mexico", limit=2500)
caribbean_checklist <- get_region_species("Caribbean", limit=2500)
brazil_checklist <- get_region_species("Brazil", limit=2500)
central_am_checklist <- get_region_species("Central America", limit=5000)
north_sa_checklist <- get_region_species("Northern South America", limit=5000)
west_sa_checklist <- get_region_species("Western South America", limit=5000)

# join lists
checklist <- bind_rows(
  mexico_checklist,
  caribbean_checklist,
  brazil_checklist,
  central_am_checklist,
  north_sa_checklist,
  west_sa_checklist
)

# deduplicate by species
checklist <- distinct(checklist, fqId, .keep_all=TRUE)

# remove hybrids and clean up names a bit
checklist <- 
  checklist %>%
  filter(! str_detect(name, "×")) %>%
  mutate(ipni_id=str_extract(fqId, "[\\d\\-]+$"),
         genus=str_extract(name, "^[A-Z][a-z\\-]+(?= )")) %>%
  select(ipni_id, name, author, family, genus) %>%
  # correct some bad family data
  mutate(family=case_when(genus == "Prockia" ~ "Salicaceae",
                          genus == "Calyptridium" ~ "Montiaceae",
                          TRUE ~ family))

# check that species in the checklist are native
# WARNING: takes a long time - checks POWO for each species in checklist
checklist <-
  checklist %>%
  rowwise() %>%
  mutate(native=check_native(ipni_id, study_region)) %>% 
  ungroup()

# remove species that aren't listed as native to our study region
checklist <-
  checklist %>%
  # no citruses should be native
  filter(genus != "Citrus") %>%
  filter(native) %>%
  select(-native)

# get class for each family, to remove non-angiosperms later
checklist <-
  checklist %>%
  nest_by(family) %>%
  mutate(class=name_backbone(family)$class) %>%
  unnest(cols=c(data)) %>%
  ungroup()

# append a known unplaced correction
checklist <-
  checklist %>%
  bind_rows(
    tibble(family="Gentianaceae", 
           ipni_id="42702-2",
           name="Calolisianthus pendulus",
           author="(Mart.) Gilg",
           genus="Calolisianthus",
           class="Magnoliopsida")
  )

write_csv(checklist, here("output/south_america_checklist.csv"))

# 4. generate counts of medicinal and overall species ----

overall_counts <-
  checklist %>%
  count(class, family, genus, name="n_species")

antimalarial_counts <-
  cleaned_names %>%
  filter(! hybrid, in_study_region, ! is.na(match_id)) %>%
  distinct(accepted_genus, accepted_species, .keep_all=TRUE) %>%
  group_by(accepted_genus) %>%
  filter(! is.na(accepted_species) | n() == 1) %>%
  ungroup() %>%
  count(accepted_genus, name="n_antimalarial")

# join two counts together
aggregated_data <-
  overall_counts %>%
  left_join(antimalarial_counts,
            by=c("genus"="accepted_genus")) %>%
  replace_na(list(n_antimalarial=0))

write_csv(aggregated_data, here("output/antimalarial_genus_counts.csv"))

# 5. trim the phylogeny ----

# choose species list
species_list <-
  checklist %>%
  filter(class %in% c("Magnoliopsida", "Liliopsida",
                      "Gnetopsida", "Pinopsida", 
                      "Cycadopsida"))

# load smith and brown tree
tree <- read.tree(here("data/ALLMB.txt"))
tree

# match labels to species in southern america
label_matches <-
  tree %>%
  "$"("tip.label") %>%
  enframe(name=NULL, value="label") %>%
  mutate(taxon_name=str_replace_all(label, "_", " ")) %>%
  mutate(taxon_name=str_replace(taxon_name, "LeucothoÃ«", "Leucothoe")) %>%
  mutate(taxon_name=str_replace(taxon_name, "HierochloÃ«", "Hierochloe")) %>%
  right_join(
    species_list, 
    by=c("taxon_name"="name")
  )

# find all records for genera that have no label matches
unmatched <- 
  label_matches %>%
  group_by(genus) %>%
  filter(all(is.na(label))) %>%
  ungroup() %>% 
  arrange(genus)

# find synonyms for those records
unmatched_synonyms <-
  unmatched %>%
  filter(! is.na(ipni_id)) %>%
  select(taxon_name, ipni_id) %>%
  rowwise() %>%
  mutate(info=list(lookup_wcvp(ipni_id)))
  
unmatched_synonyms <-
  unmatched_synonyms %>%
  mutate(synonyms=tidy(info)$synonyms) %>%
  select(-info) %>%
  unnest(cols=c(synonyms)) %>%
  select(taxon_name, ipni_id, synonym_id=id, 
         synonym_name=name)
  
# match labels in the tree to those synonyms
synonym_matches <- 
  tree %>%
  "$"("tip.label") %>%
  enframe(name=NULL, value="label") %>%
  mutate(taxon_name=str_replace_all(label, "_", " ")) %>%
  mutate(taxon_name=str_replace(taxon_name, "LeucothoÃ«", "Leucothoe")) %>%
  mutate(taxon_name=str_replace(taxon_name, "HierochloÃ«", "Hierochloe")) %>%
  inner_join(unmatched_synonyms, by=c("taxon_name"="synonym_name"), suffix=c("", "_accepted")) %>%
  distinct(taxon_name_accepted, .keep_all=TRUE) %>%
  left_join(species_list, by=c("taxon_name_accepted"="name", "ipni_id")) %>%
  select(-taxon_name, -synonym_id) %>%
  rename(taxon_name=taxon_name_accepted)

# find all records for genera that have no label matches again
unmatched_second <- 
  label_matches %>%
  bind_rows(synonym_matches) %>%
  group_by(genus) %>%
  filter(all(is.na(label))) %>%
  ungroup() %>% 
  arrange(genus)

# make a regular expression to find any label starting with the missing genera
regex <-
  unmatched_second %>%
  pull(genus) %>%
  unique() %>%
  paste(collapse="|")

regex <- paste0("^(", regex, ") [a-z\\-]+")

# find any tips belonging to the missing genera
candidate_tips <-
  tree %>%
  "$"("tip.label") %>%
  enframe(name=NULL, value="label") %>%
  mutate(taxon_name=str_replace_all(label, "_", " ")) %>%
  mutate(taxon_name=str_replace(taxon_name, "LeucothoÃ«", "Leucothoe")) %>%
  mutate(taxon_name=str_replace(taxon_name, "HierochloÃ«", "Hierochloe")) %>%
  filter(str_detect(taxon_name, regex)) %>%
  mutate(genus=str_extract(taxon_name, "^[A-Za-z\\-]+(?= )"))

# join all the matched labels together and deduplicate
matched_labels <-
  label_matches %>%
  filter(! is.na(label)) %>%
  bind_rows(synonym_matches) %>%
  select(label, genus) %>%
  bind_rows(candidate_tips) %>%
  distinct(genus, .keep_all=TRUE)

# drop all tips we haven't found matches for
genus_tree <- subset_tree(tree, matched_labels$label)

# replace all tip labels with genus names
label_map <- matched_labels$genus
names(label_map) <- matched_labels$label
  
genus_tree$tip.label <- map_chr(genus_tree$tip.label, ~label_map[.x])

genus_tree

write.tree(genus_tree, here("output/genus_level_tree.txt"))

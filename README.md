# Plants used traditionally as antimalarials in Latin America: mining the Tree of Life for potential new medicines

This repository contains analysis scripts for the paper ***Plants used traditionally as antimalarials in Latin America: mining the Tree of Life for potential new medicines***.

Some scripts take a long time to run, so the analysis outputs have been preserved in the *outputs* folder.

## Running the analysis

To run the analysis, download or clone this repository and run the script `run_analysis.R`.

This will run all of the analysis scripts stored in the *analysis* folder. These scripts are:

1. `01_prepare_names.R`
    * takes the names of taxa identified as being used as antimalarials
    * matches them to accepted species names in the [World Checklist of Vascular Plants (WCVP)](https://wcvp.science.kew.org/)
    * removes any not listed as native to Latin America in [Plants of the World Online (POWO)](http://powo.science.kew.org/)
    * downloads a checklist of all accepted species in Latin America
    * tallies the number of antimalarial and overall species in each genus
    * trims the inclusive seedplant phylogeny from [Smith and Brown](https://bsapubs.onlinelibrary.wiley.com/doi/full/10.1002/ajb2.1019) to a genus level tree that covers these genera

2. `02_binomial_analysis.R`
    * compares the proportion of antimalarial species in each genus to what would be expected from a binomial distribution if all species had the same probability of being reported as antimalarial
    * repeats this for taxonomic families

3. `03_phylogenetic_analysis.R`
    * calculates measures of phylogenetic diversity (NRI and NTI) both taking into account the number of species reported as antimalarial in each genus and for binary presence/absence of a reported antimalarial in each genus
    * calculates measures of phylogenetic signal treating the proportion of antimalarial species in each genus as a continuous trait and as a binary presence/absence
    * calculates a measure of the over- or under-representation of antimalarial species at each node of the phylogenetic tree, as performed by [`nodesig` in *Phylocom*](https://github.com/phylocom/phylocom/)

4. `04_make_figures.R`
    * makes figures used in the paper
    * makes an Excel file containing the results of the analysis

## Data

There are two sources of data needed for this analysis, both stored in the *data* folder:

1. `antimalarial_taxa.csv` - the list of taxa that have been reported to be used as antimalarials in Latin America. This is different from the taxa list supplied as *Table A.1* in the Supplementary Materials of the paper. *Table A.1* was updated after the name matching in `01_prepare_names.R` to follow the taxonomy of the WCVP and includes some unpublished records that have been excluded from this analysis. The mapping of names from this file to the names in *Table A.1* can be found in `output/antimalarials_matched.csv`.

2. `ALLMB.txt` - version 0.1 of the broadly inclusive seedplant phylogeney constructed by [Smith and Brown](https://bsapubs.onlinelibrary.wiley.com/doi/full/10.1002/ajb2.1019). This is the tree that they constructed from GenBank sequence data and Open Tree of Life taxonomic data, using a backbone from [Magallon et al](https://nph.onlinelibrary.wiley.com/doi/10.1111/nph.13264). The tree used here, and possibly updated versions, can be [downloaded here](https://github.com/FePhyFoFum/big_seed_plant_trees).

The *data* folder also contains some `json` files. These were put together during the name preparation process to resolve:

* names that matched to multiple records in WCVP (`name_match_multiples.json`)
* names that did not match to anything in WCVP (`name_match_missing.json`)
* names that matched to records in WCVP but needed correcting (`manual_match_correction.json`)

## Requirements

This analysis was carried out in R version 4.0.3, which you can download [here](https://www.r-project.org/).

The analysis uses these packages:

* [here](https://here.r-lib.org/) (version 0.1)
* [tidyverse](https://www.tidyverse.org/) (version 1.3.0)
* [writexl](https://docs.ropensci.org/writexl/) (version 1.3.1)
* [jsonlite](https://github.com/jeroen/jsonlite) (version 1.7.1)
* [ape](http://ape-package.ird.fr/) (version 5.4.1)
* [rgbif](https://github.com/ropensci/rgbif) (version 3.3.0)
* [kewr](https://barnabywalker.github.io/kewr/) (version 0.4.0)
* [progress](https://github.com/r-lib/progress) (version 1.2.2)
* [PhyloMeasures](https://cran.r-project.org/web/packages/PhyloMeasures/index.html) (version 2.1)
* [phytools](https://github.com/liamrevell/phytools) (version 0.7.70)
* [caper](https://cran.r-project.org/package=caper) (version 1.0.1)
* [ggtree](https://guangchuangyu.github.io/software/ggtree/) (version 2.4.1)
* [ggrepel](https://ggrepel.slowkow.com/index.html) (version 0.8.2)
* [patchwork](https://patchwork.data-imaginist.com/) (version 1.1.0)

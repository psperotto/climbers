rm(list=ls())
library(ape)
source("neotropical_climbers_functions.R") 

# Taxizing tree
sb_tree <- read.tree("Data/GBMB.tre")
sb_tree_taxized_tips <- gbif.taxize(sb_tree$tip.label)
sb_tree$tip.label <- sb_tree_taxized_tips
saveRDS(sb_tree, file="Data/taxized_GBMB.Rdata")

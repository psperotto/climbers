# Taxize names in trees for analyses
# rm(list=ls())
setwd("~/Desktop/WCVP_special_issue/Patricia_Climbers/climbers")

source("00_utility_functions.R")
# Taxize names in climber dataset
# Taxizing tree
sb_tree <- read.tree("Data/ALLOTB.tre")
sb_tree_taxized_tips <- gbif.taxize(sb_tree$tip.label)
sb_tree$tip.label <- sb_tree_taxized_tips
saveRDS(sb_tree, file="Data/taxized_ALLOTB.Rdata")

# Fixing previous taxize
climbers <- read.csv("Data/climber_database.csv", stringsAsFactors = F)
climbers$taxized_names <- fix.names.taxize(climbers$taxized_names)
write.csv(climbers, file="Data/climber_database.csv", row.names=F)

sb_tree <- readRDS("Data/taxized_ALLOTB.Rdata")
sb_tree$tip.label <- fix.names.taxize(sb_tree$tip.label)
saveRDS(sb_tree, file="Data/taxized_ALLOTB.Rdata")

sb_tree_gbmb <- readRDS("Data/taxized_GBMB.Rdata")
sb_tree_gbmb$tip.label <- fix.names.taxize(sb_tree_gbmb$tip.label)
saveRDS(sb_tree_gbmb, file="Data/taxized_GBMB.Rdata")





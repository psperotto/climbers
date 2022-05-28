# rm(list=ls())
setwd("~/Desktop/WCVP_special_issue/Patricia_Climbers/climbers")

library(maptools)
library(data.table)
data("wrld_simpl")

source("00_utility_functions.R")
source("/Users/thaisvasconcelos/Desktop/WCVP_special_issue/WCVPtools/WCVPtools_functions.R")


# Hypothesis: Groups with different climbing mechanisms have different support diameter requirements, 
# and thus the distribution of lineages with different climbing mechanism is influenced by the 
# surrounding vegetation type

climber_points <- as.data.frame(fread("gbif_climbers/climbers_cleaned_points.csv"))
climber_mech <- read.csv("Data/climber_database.csv")
colnames(climber_mech) <- c("family","genus","species","CM","taxized_names")
merged_table <- merge(climber_mech, climber_points, by.x = "taxized_names", by.y = "scientificName")

# Directory to save preliminary datasets:
climate_data.dir <- "./climate_data"
# Directory where climate layers are:
climate_layers.dir <- "./climate_layers"

# 1. Thinning occurence data first
thinned_points <- Thinning(merged_table, species="taxized_names", lat = "decimalLatitude", lon="decimalLongitude", n = 3)
# write.csv(thinned_points, file=paste0(climate_data.dir, "/climbers_thinned_points.csv"), row.names=F)
# thinned_points <- read.csv(paste0(climate_data.dir, "/climbers_thinned_points.csv"))

# 2. Getting summary statistics of climatic variables for each species
allpoints <- ClimateFromPoint_custom(thinned_points, species="species",lon="lon", lat="lat", layerdir = climate_layers.dir)
write.csv(allpoints, file=paste0(climate_data.dir, "/climbers_allpoints.csv"), row.names=F)
summstats <- GetClimateSummStats_custom(allpoints, type="raw")
write.csv(summstats, file=paste0(climate_data.dir, "/climbers_summstats_raw.csv"), row.names=F)  
summstats <- GetClimateSummStats_custom(allpoints, type="transformed")
write.csv(summstats, file=paste0(climate_data.dir, "/climbers_summstats.csv"), row.names=F)

# Different ways in which we can test this:
# (1) Phyloanovas
# Higher AI will represent higher potential for vegetation growth (i.e. closed canopy environments)
summstats <- read.csv(paste0(climate_data.dir, "/climbers_summstats.csv"))

ai_table <- summstats[,c("species","mean_aridity")]
ai_table$mean_aridity <- as.numeric(ai_table$mean_aridity)
merged_table <- merge(climber_mech, ai_table, by.x="taxized_names", by.y="species")

boxplot(exp(merged_table$mean_aridity) ~ merged_table$CM)

# Load tree for PCMs
library(ape)
library(phytools)
tree <- readRDS("Data/taxized_GBMB.Rdata")

tree$tip.label
tree$tip.label <- unname(tree$tip.label)



subset_merged_table <- subset(merged_table, merged_table$taxized_names%in%tree$tip.label)
subset_merged_table <- subset(subset_merged_table, !duplicated(subset_merged_table$taxized_names))
pruned_tree <- keep.tip(tree, intersect(subset_merged_table$taxized_names, tree$tip.label))
pruned_tree$tip.label <- paste0(unlist(lapply(strsplit(pruned_tree$tip.label, " "), "[[", 1)), " ", unlist(lapply(strsplit(pruned_tree$tip.label, " "), "[[", 2)))
subset_merged_table$taxized_names <- paste0(unlist(lapply(strsplit(subset_merged_table$taxized_names, " "), "[[", 1)), " ", unlist(lapply(strsplit(subset_merged_table$taxized_names, " "), "[[", 2)))

# Figure
boxplot(exp(subset_merged_table$mean_aridity) ~ subset_merged_table$CM)

aridity <- subset_merged_table$mean_aridity
mechanisms <- subset_merged_table$CM
names(aridity) <- names(mechanisms) <- subset_merged_table$taxized_names

class(mechanisms)
mechanisms[which(mechanisms==1)] <- "a"
mechanisms[which(mechanisms==2)] <- "b"
mechanisms[which(mechanisms==3)] <- "c"
mechanisms[which(mechanisms==4)] <- "d"
mechanisms[which(mechanisms==5)] <- "e"
mechanisms[which(mechanisms==6)] <- "f"
mechanisms[which(mechanisms==7)] <- "g"
mechanisms[which(mechanisms==8)] <- "h"

phyanova <- phylANOVA(pruned_tree, mechanisms, aridity)
# 



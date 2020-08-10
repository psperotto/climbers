# Diversification of Neotropical climbers across space and time

library(data.table)
library(monographaR)
library(tidyverse)
library(raster)
library(CoordinateCleaner)
library(ape)
source("~/Desktop/climbers/neotropical_climbers_functions.R")

# First get grid cell values from filtered GBIF data 
####
#### from TV local
# gbif_dir <- "~/Desktop/full_gbif_quest/" 
# gbif_files <- list.files(paste0(gbif_dir, "z_filtered_gbif"), ".csv")
# full_list <- read.full.gbif(gbif_files, gbif_dir)
# write.csv(full_list, file=paste0(gbif_dir, "full_tracheophyte_filtered_gbif.csv"))

# load table back
# setwd(gbif_dir)
# full_list <- fread(paste0(gbif_dir, "full_tracheophyte_filtered_gbif.csv"))
# full_list <- full_list[,-1]
# save.gbif.neotropics(full_list)

##############
climbers_dir <- "~/Desktop/climbers" # Repo 
setwd(climbers_dir)
full_list <- fread(file.choose()) ## Procurar "neotropics_tracheophyte_filtered_gbif.csv" no computador. NÃO COLOCAR ESSE ARQUIVO NO REPO
full_list <- full_list[,-1]

full_map <- run.mapDiversity.neotropics(full_list)

pdf("all_angio_div.pdf", height=8, width=5)
plot(full_map)
dev.off()

# load table with climbers

climbers <- read.csv("climber_database.csv", stringsAsFactors = F)

# Make map of grid cell diversity for each of the 8 mechanisms
for(mechanism_index in 1:8) {
  # Select mechanism
  mech <- climbers$Species[climbers$CM==mechanism_index]
  # Taxize names
  mech_taxized <- gbif.taxize(mech)
  # Get points from full list
  mech_points <- subset(full_list, full_list$species %in% as.character(mech_taxized))
  # Get raster of grid cells 
  mech_map <- run.mapDiversity.neotropics(mech_points, filename=paste0("Mechanism_", mechanism_index))
  # Plot pdf
  pdf(paste0("Mechanism_", mechanism_index,".pdf"), height=8, width=5)
  plot(mech_map)
  dev.off()
}

# Load raster back
full_map <- readRDS("full_neotropical_diversity.Rdata")
mechanism1 <- readRDS("Mechanism_1.Rdata")


# corrected for PW?
data = mech_points[mech_points$species %in% phy$tip.label,]
phy <- phy.list(input.dir=climbers_dir, names="GBMB", search.for=".taxized.tre")[[1]] 
phy$tip.label <- gsub("_"," ",phy$tip.label)
phy <- keep.tip(phy, which(phy$tip.label %in% unique(data$species)))
mapDiversity.pw(data, phy)


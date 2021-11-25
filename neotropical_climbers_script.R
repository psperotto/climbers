# Diversification of Neotropical climbers across space and time
# setwd("~/Desktop/Colabs/Patricia_Climbers/climbers")
# setwd("C:/Users/patri/Google Drive/Papers/Diversifica??o/climbers") # Repo 

library(data.table)
library(monographaR)
library(tidyverse)
library(raster)
library(CoordinateCleaner)
library(ape)
library(maptools)
data("wrld_simpl")

source("neotropical_climbers_functions.R") 
      # source das funcoes q vamos usar

# First get grid cell values from filtered GBIF data 
#### from TV local
 gbif_dir <- paste0(getwd(), "/full_gbif_quest")
 gbif_files <- list.files(paste0(gbif_dir, "/z_filtered_gbif"), ".csv")
 full_list <- read.full.gbif(gbif_files, gbif_dir)
 write.csv(full_list, file=paste0(gbif_dir, "full_tracheophyte_filtered_gbif.csv"))

# load table back
 full_list <- fread(paste0(gbif_dir, "/full_tracheophyte_filtered_gbif.csv"))
 full_list <- full_list[,-1]
 save.gbif.neotropics(full_list)

##############
# Directory for descriptive results
descriptive_dir <- paste0(getwd(), "/1_descriptive")
   
full_list <- fread(paste0(gbif_dir, "/neotropics_tracheophyte_filtered_gbif.csv"))
full_list <- full_list[,-1]
full_map <- run.mapDiversity.neotropics(full_list, filename="full_neotropical_diversity", dir=descriptive_dir) # function modified from monographaR

# Plot
pdf(paste0(descriptive_dir,"/full_neotropical_diversity.pdf"), height=8, width=5)
plot(full_map)
dev.off()

# Now load table climbers
climbers <- read.csv("Data/climber_database.csv", stringsAsFactors = F)

# Make map of grid cell diversity for each of the 8 mechanisms
for(mechanism_index in 1:8) {
  # Select mechanism
  mech <- climbers$Species[climbers$CM==mechanism_index]
  # Taxize names
  mech_taxized <- gbif.taxize(mech)
  # Get points from full list
  mech_points <- subset(full_list, full_list$species %in% as.character(mech_taxized))
  saveRDS(mech_points, file=paste0("Data/","Mechanism_", mechanism_index,"_taxized_points.Rdata"))
  # Get raster of grid cells 
  mech_map <- run.mapDiversity.neotropics(mech_points, filename=paste0("Mechanism_", mechanism_index), dir=descriptive_dir)
  # Plot pdf
  pdf(paste0(descriptive_dir,"/Mechanism_", mechanism_index,".pdf"), height=8, width=5)
  plot(mech_map)
  dev.off()
}

###########
# Load raster back
rasters <- list.files(descriptive_dir, ".Rdata")
rasters <- rasters[grep("Mechanism", rasters)]
rasters.names <- sub(".Rdata", "", rasters)
rasters <- lapply(paste0(descriptive_dir, "/", rasters), readRDS)
names(rasters) <- rasters.names
  
full_map <- readRDS("1_descriptive/full_neotropical_diversity.Rdata")

template.map <- full_map
template.map[!is.na(template.map[])] <- 0

tmp.raster.list <- list()
for (i in 1:length(rasters)) {
  r1 <- rasters[[i]]
  r1 <- raster::resample(r1, template.map)
  r1[is.na(r1)] <- 0
  tmp.raster.list[[i]] <- raster::mask(r1, template.map)
  print(i)
}
sprichness_climbers_map <- raster::calc(raster::stack(tmp.raster.list), sum)
sprichness_climbers_map <- raster::crop(sprichness_climbers_map, raster::extent(full_map))

#plot(sprichness_climbers_map)
#plot(wrld_simpl, add=T)

plot.res(sprichness_climbers_map, full_map, dir=getwd(), pal.name = "RdBu", output = "residuals_sprich_pw")

### residuals 

pdf("res_mechanism1.pdf", height=5, width=8)
plot.res(mechanism1, full_map, dir=getwd(), pal.name = "RdBu", output = "residuals_sprich_pw")
dev.off()

pdf("res_mechanism2.pdf", height=5, width=8)
plot.res(mechanism2, full_map, dir=getwd(), pal.name = "RdBu", output = "residuals_sprich_pw")
dev.off()

pdf("res_mechanism3.pdf", height=5, width=8)
plot.res(mechanism3, full_map, dir=getwd(), pal.name = "RdBu", output = "residuals_sprich_pw")
dev.off()

pdf("res_mechanism4.pdf", height=5, width=8)
plot.res(mechanism4, full_map, dir=getwd(), pal.name = "RdBu", output = "residuals_sprich_pw")
dev.off()

pdf("res_mechanism1.pdf", height=5, width=8)
plot.res(mechanism1, full_map, dir=getwd(), pal.name = "RdBu", output = "residuals_sprich_pw")
dev.off()

pdf("res_mechanism5.pdf", height=5, width=8)
plot.res(mechanism5, full_map, dir=getwd(), pal.name = "RdBu", output = "residuals_sprich_pw")
dev.off()

pdf("res_mechanism6.pdf", height=5, width=8)
plot.res(mechanism6, full_map, dir=getwd(), pal.name = "RdBu", output = "residuals_sprich_pw")
dev.off()

pdf("res_mechanism7.pdf", height=5, width=8)
plot.res(mechanism7, full_map, dir=getwd(), pal.name = "RdBu", output = "residuals_sprich_pw")
dev.off()

pdf("res_mechanism8.pdf", height=5, width=8)
plot.res(mechanism8, full_map, dir=getwd(), pal.name = "RdBu", output = "residuals_sprich_pw")
dev.off()
# adicionar coisas




# PD and phyloANOVA



# Carregar a arvore
phy <- phy.list(input.dir=climbers_dir, names="GBMB", search.for=".taxized.tre")[[1]] 
# Substitui "_" por " "
phy$tip.label <- gsub("_"," ",phy$tip.label)
# lista de pontos com os nomes que tem na arvore e arvore com tips que tem na lista de pontos
   ## mas tenho que fazer isso ent?o pra cada um dos mechs?

mech_list <- list()
for(mechanism_index in 1:8) {
  # Select mechanism
  mech <- climbers$Species[climbers$CM==mechanism_index]
  # Taxize names
  mech_taxized <- gbif.taxize(mech)
  # Get points from full list
  mech_points <- subset(full_list, full_list$species %in% as.character(mech_taxized))
  ##
  data <- mech_points[mech_points$species %in% phy$tip.label,] 
  phy_temp <- keep.tip(phy, which(phy$tip.label %in% unique(data$species)))
  mech_list[[mechanism_index]] <- mapDiversity.pw(data, phy_temp)
}
  # IT WORKED! *-*
  # mas os mapas estao sem o delineamento dos paises, tenho q ver como coloca

saveRDS(mech_list,file="mech_list.Rdata")

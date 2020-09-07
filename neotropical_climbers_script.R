# Diversification of Neotropical climbers across space and time

library(data.table)
library(monographaR)
library(tidyverse)
library(raster)
library(CoordinateCleaner)
library(ape)


source("neotropical_climbers_functions.R") 
      # source das fun?oes q vamos usar

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
climbers_dir <- "C:/Users/patri/Google Drive/Papers/Diversificação/climbers" # Repo 
setwd(climbers_dir)
full_list <- fread(file.choose()) 
## Procurar "neotropics_tracheophyte_filtered_gbif.csv" no computador. 
## N?O COLOCAR ESSE ARQUIVO NO REPO
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

###########
# Load raster back
rasters <- list.files(climbers_dir, ".Rdata")
rasters.names <- sub(".Rdata", "", rasters)
rasters <- lapply(paste0(climbers_dir, "/", rasters), readRDS)
names(rasters) <- rasters.names
  
full_map <- readRDS("full_neotropical_diversity.Rdata")
mechanism1 <- readRDS("Mechanism_1.Rdata")
mechanism2 <- readRDS("Mechanism_2.Rdata")
mechanism3 <- readRDS("Mechanism_3.Rdata")
mechanism4 <- readRDS("Mechanism_4.Rdata")
mechanism5 <- readRDS("Mechanism_5.Rdata")
mechanism6 <- readRDS("Mechanism_6.Rdata")
mechanism7 <- readRDS("Mechanism_7.Rdata")
mechanism8 <- readRDS("Mechanism_8.Rdata")

### residuals (n?o sei fazer loops entao copiei tudo varias vezes haha)
    ## tudo bem, plotamos os residuals. mas como calcular estatisticamente pra cada gridcell onde
    ## tem mais/menos spp do que o esperado? e a escala est? em qu??
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

# corrected for PW
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

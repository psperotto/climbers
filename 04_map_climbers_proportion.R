setwd("~/Desktop/WCVP_special_issue/Patricia_Climbers/climbers")

library(data.table)
library(monographaR)
library(maptools)
library(raster)
data("wrld_simpl")

all_angios <- fread("full_gbif_quest/neotropics_tracheophyte_filtered_gbif.csv")
all_angios <- all_angios[,2:4]

wcvp_climbers <- readRDS("wcvp_nt_climbers_final.Rdata")
all_climbers <- paste(wcvp_climbers$taxon_name, wcvp_climbers$taxon_authors, sep=" ")

climber_points <- subset(all_angios, all_angios$species %in% all_climbers)
#non_climber_points <- subset(all_angios, !all_angios$species %in% all_climbers)

climber_points<-as.data.frame(climber_points)
all_angios <- as.data.frame(all_angios)

climber_diversity_raster <- mapDiversity(climber_points)
neotropics_diversity_raster <- mapDiversity(all_angios)

neotropics_diversity_raster[which(neotropics_diversity_raster[]<10)] <- NA
climber_diversity_raster[which(is.na(neotropics_diversity_raster[]))] <- NA

proportion_map <- climber_diversity_raster / neotropics_diversity_raster

pal <- hcl.colors(200, palette = "Plasma", alpha = 0.7)
plot(proportion_map, col=pal, zlim=c(0,0.25))
plot(climber_diversity_raster, col=pal)
plot(wrld_simpl, add=T)


# setwd("~/Desktop/WCVP_special_issue/Patricia_Climbers/climbers")
# rm(list=ls())
require(ape);
require(RColorBrewer);
require(sp); 
require(ggplot2);
library(rgdal);
library(data.table);
library(corpcor)
library(plyr)
library(automap)
library(reshape2)

climber_points <- read.csv("climbers_points.csv")
climber_mech <- read.csv("Data/climber_database.csv")
colnames(climber_mech) <- c("family","genus","species","CM","taxized_names")

merged_table <- merge(climber_mech, climber_points, by.x = "taxized_names", by.y = "species")

###################################
# TEST - Bubble map per TWGD
###################################
organize.bubble.plot <- function(trait_table, reference_table, all_vars, twgd_data) {
  tmp_reference_table <- subset(reference_table, reference_table$gbif_name %in% unique(trait_table$taxized_names))
  wcvp_subset <- subset(all_vars, all_vars$taxon_name %in% tmp_reference_table$wcvp_name)
  wcvp_subset <- subset(wcvp_subset, wcvp_subset$introduced==0)
  wcvp_subset <- subset(wcvp_subset, wcvp_subset$extinct==0)
  wcvp_subset <- subset(wcvp_subset, wcvp_subset$location_doubtful==0)
  
  focal_areas <- unique(wcvp_subset$area_code_l3)
  results <- matrix(nrow=0, ncol=5)
  for(i in 1:length(focal_areas)) {
    one_area <- focal_areas[i]
    one_subset <- subset(wcvp_subset, wcvp_subset$area_code_l3==one_area)
    sp_rich <- length(unique(one_subset$taxon_name))
    family_rich <- length(unique(one_subset$family))
    area_plus_buffer <- twgd_data[which(as.character(twgd_data$LEVEL3_COD) %in% one_area),]
    centroids <- rgeos::gCentroid(area_plus_buffer, byid=TRUE)
    lon <- extent(centroids)[1]
    lat <- extent(centroids)[3]
    results <- rbind(results, cbind(sp_rich, family_rich, one_area, lon, lat))
  }
  results <- as.data.frame(results)
  results$sp_rich <- as.numeric(results$sp_rich)
  results$family_rich <- as.numeric(results$family_rich)
  results$lon <- as.numeric(results$lon)
  results$lat <- as.numeric(results$lat)
  return(results)
}

library(data.table)
library(maptools)
library(raster)
library(sp)
library(rgeos)
library(rworldmap)
data("wrld_simpl")

#-----------------------------
# If local
# setwd("~/Desktop/WCVP_special_issue/Eve_MyrtalesPAFTOL/myrtales")
source("/Users/thaisvasconcelos/Desktop/WCVP_special_issue/WCVPtools/WCVPtools_functions.R")
dist_sample <- read.table("../../wcvp_names_and_distribution_special_edition_2022/wcvp_distribution.txt", sep="|", header=TRUE, quote = "", fill=TRUE, encoding = "UTF-8")
names_sample <- read.table("../../wcvp_names_and_distribution_special_edition_2022/wcvp_names.txt", sep="|", header=TRUE, quote = "", fill=TRUE, encoding = "UTF-8")
#-----------------------------
# If labcomputer
# setwd("~/myrtales")
source("../WCVPtools/WCVPtools_functions.R")
dist_sample <- read.table("../life_history_houwie/wcvp_names_and_distribution_special_edition_2022/wcvp_distribution.txt", sep="|", header=TRUE, quote = "", fill=TRUE, encoding = "UTF-8")
names_sample <- read.table("../life_history_houwie/wcvp_names_and_distribution_special_edition_2022/wcvp_names.txt", sep="|", header=TRUE, quote = "", fill=TRUE, encoding = "UTF-8")
#-----------------------------
# Merge them in one big table
all_vars <- merge(dist_sample, names_sample, by="plant_name_id")

# reference table for taxized names
#-----------------------------
# If local
reference_table <- list.files("/Users/thaisvasconcelos/Desktop/WCVP_special_issue/WCVPtools/taxized_reference_tables", full.names = T)
reference_table <- do.call(rbind, lapply(reference_table, read.csv))
#-----------------------------
# If labcomputer
reference_table <- list.files("../WCVPtools/taxized_reference_tables", full.names = T)
reference_table <- do.call(rbind, lapply(reference_table, read.csv))

# If local
path="/Users/thaisvasconcelos/Desktop/WCVP_special_issue/WCVPtools/wgsrpd-master/level3/level3.shp"
#-----------------------------
# If labcomputer
path="../WCVPtools/wgsrpd-master/level3/level3.shp"
#-----------------------------
twgd_data <- suppressWarnings(maptools::readShapeSpatial(path))

organized_table_for_plot1 <- organize.bubble.plot(merged_table, reference_table, all_vars, twgd_data)
# filter neotropics
organized_table_for_plot1 <- subset(organized_table_for_plot1, organized_table_for_plot1$lon < -25 & organized_table_for_plot1$lon > -90)
organized_table_for_plot1 <- subset(organized_table_for_plot1, organized_table_for_plot1$lat < 23 & organized_table_for_plot1$lat > -30)

mech1 <- subset(merged_table, merged_table$CM==1)
organized_table_for_plot1 <- organize.bubble.plot(mech1, reference_table, all_vars, twgd_data)
# filter neotropics
organized_table_for_plot1 <- subset(organized_table_for_plot1, organized_table_for_plot1$lon < -25 & organized_table_for_plot1$lon > -90)
organized_table_for_plot1 <- subset(organized_table_for_plot1, organized_table_for_plot1$lat < 23 & organized_table_for_plot1$lat > -30)

nrow(organized_table_for_plot1)


mybreaks <- c(1, 5, 10, 15, 20, 25)

library(ggplot2)
library(maps)
library(ggthemes)
library(viridis)

world <- ggplot() +
  borders("world", colour = "gray85", fill = "gray80") +
  theme_map() 

# Build the map
map <- world +
  geom_point(data = organized_table_for_plot1, aes(x=lon, y=lat, size=family_rich, color=sp_rich), shape=20, stroke=FALSE) +
  scale_size_continuous(name="number of families", range=c(1,30), breaks=mybreaks) +
  #scale_alpha_continuous(trans="log", range=c(0.5, 0.5)) +
  scale_color_viridis(option="magma", name="species richness", alpha=0.5) 

organized_table_for_plot1$family_rich * 100

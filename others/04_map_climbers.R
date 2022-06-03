setwd("~/Desktop/WCVP_special_issue/Patricia_Climbers/climbers")

library(maptools)
library(data.table)
library(raster)
data("wrld_simpl")

source("00_utility_functions.R")
source("/Users/thaisvasconcelos/Desktop/WCVP_special_issue/WCVPtools/WCVPtools_functions.R")
source("/Users/thaisvasconcelos/Desktop/WCVP_special_issue/WCVPtools/GetRanges.R")
source("/Users/thaisvasconcelos/Desktop/WCVP_special_issue/WCVPtools/GetTraitDistribution.R")

# Some smart way of plotting them? 







#
climber_points <- as.data.frame(fread("gbif_climbers/climbers_cleaned_points.csv"))
climber_mech <- read.csv("Data/climber_database.csv")
colnames(climber_mech) <- c("family","genus","species","CM","taxized_names")
merged_table <- merge(climber_mech, climber_points, by.x = "taxized_names", by.y = "scientificName")




# PILOT APRIL 26
mechanism_5 <- subset(merged_table, merged_table$CM==5)
climber_distribution <- GetRanges(mechanism_5, species="taxized_names", res=10)
climber_distribution[which(unname(unlist(lapply(climber_distribution, is.character))))] <- NULL
x <- GetSpRichness(climber_distribution)
saveRDS(x, file="mech5_pilot.Rdata")



library(picante)

#-----------------------------
dist_sample <- read.table("../../wcvp_names_and_distribution_special_edition_2022/wcvp_distribution.txt", sep="|", header=TRUE, quote = "", fill=TRUE, encoding = "UTF-8")
names_sample <- read.table("../../wcvp_names_and_distribution_special_edition_2022/wcvp_names.txt", sep="|", header=TRUE, quote = "", fill=TRUE, encoding = "UTF-8")

#-----------------------------
# Merge them in one big table
all_vars <- merge(dist_sample, names_sample, by="plant_name_id")

# reference table for taxized names
#-----------------------------
# If local
reference_table <- list.files("/Users/thaisvasconcelos/Desktop/WCVP_special_issue/WCVPtools/taxized_reference_tables", full.names = T)
reference_table <- do.call(rbind, lapply(reference_table, read.csv))

# If local
path="/Users/thaisvasconcelos/Desktop/WCVP_special_issue/WCVPtools/wgsrpd-master/level3/level3.shp"

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
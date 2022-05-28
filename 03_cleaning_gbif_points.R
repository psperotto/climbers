# rm(list=ls())
setwd("~/Desktop/WCVP_special_issue/Patricia_Climbers/climbers")
source("00_utility_functions.R")

#-----------------------------
# If local
source("/Users/thaisvasconcelos/Desktop/WCVP_special_issue/WCVPtools/WCVPtools_functions.R")
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
reference_table$gbif_name <- fix.names.taxize(reference_table$gbif_name)

# Reading gbif file
gbif_data <- fread("gbif_climbers/0306501-210914110416597.csv") # load the table you downloaded from GBIF

# Looking at the WCVP table and TDWG to clean GBIF points
#-----------------------------
# If local
path="/Users/thaisvasconcelos/Desktop/WCVP_special_issue/WCVPtools/wgsrpd-master/level3/level3.shp"

twgd_data <- suppressWarnings(maptools::readShapeSpatial(path))
#issues_to_remove <- read.csv("gbif_issues_to_remove.csv")
#cultivated <- read.csv("cultivated_species.csv")

#------ Cleaning steps:
cleaned_points <- gbif_data
# WCVP filtering
subset_reference_table <- subset(reference_table, reference_table$gbif_name %in% unique(cleaned_points$scientificName))
subset_all_vars <- subset(all_vars, all_vars$taxon_name %in% subset_reference_table$wcvp_name)
cleaned_points <- FilterWCVP_genus(cleaned_points, subset_all_vars, twgd_data)
# 0 points removed
if(nrow(subset_reference_table)>0){
    cleaned_points <- FilterWCVP(cleaned_points, all_vars, subset_reference_table, twgd_data) # This will filter the GBIF points acording to WCVP for species
}
# Restricting to Neotropics  
cleaned_points <- subset(cleaned_points, cleaned_points$decimalLongitude < -25 & cleaned_points$decimalLongitude > -90)
cleaned_points <- subset(cleaned_points, cleaned_points$decimalLatitude < 23 & cleaned_points$decimalLatitude > -30)
# Cleaning common problems:
cleaned_points <- RemoveCentroids(cleaned_points, lon="decimalLongitude", lat="decimalLatitude")
cleaned_points <- RemoveDuplicates(cleaned_points, lon="decimalLongitude", lat="decimalLatitude")
cleaned_points <- RemoveOutliers(cleaned_points, species="scientificName", lon="decimalLongitude", lat="decimalLatitude")
cleaned_points <- RemoveZeros(cleaned_points, lon="decimalLongitude", lat="decimalLatitude")
cleaned_points <- RemoveSeaPoints(cleaned_points, lon="decimalLongitude", lat="decimalLatitude")
write.csv(cleaned_points, file=paste0("gbif_climbers/climbers_cleaned_points.csv"), row.names=F)

#------------------------
all_cleaned_points_files <- list.files("gbif_climbers", full.names = T)
all_cleaned_points_files <- all_cleaned_points_files[grep("cleaned", all_cleaned_points_files)]
points_cleaned <- read.csv(all_cleaned_points_files)

# Plotting to inspect distributions
  all_spp <- unique(points_cleaned$scientificName)
  all_spp <- subset(all_spp, all_spp!="")
  pdf(paste0("climbers_points.pdf"))
  for(genus_index in 1:length(all_spp)){
    tmp_subset <- as.data.frame(points_cleaned[points_cleaned$genus==genera[genus_index],])
    coord <- tmp_subset[,c("decimalLongitude","decimalLatitude")]
    coordinates(coord) <- ~ decimalLongitude + decimalLatitude
    plot(wrld_simpl)
    plot(coord, col="red", add=T)
    title(genera[genus_index])
    print(genus_index)
  }
  dev.off()


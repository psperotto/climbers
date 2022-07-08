
# rm(list=ls())
library(dplyr)
# function to filter dubious ids from gbif data according to what the wcvp says about their distribution
FilterWCVP <- function(points, all_vars, reference_table, twgd_data, lon="decimalLongitude", lat="decimalLatitude") {
  npoints_start <- nrow(points)
  tmp_points = as.data.frame(points)
  colnames(tmp_points)[colnames(tmp_points)==lon] <- "x"
  colnames(tmp_points)[colnames(tmp_points)==lat] <- "y"
  tmp_points = subset(tmp_points, !is.na(tmp_points$x))
  tmp_points = subset(tmp_points, !is.na(tmp_points$y))
  # Load shape files and make sure they have the same name as the WCVP column with the TDWG areas
  # twgd_data <- suppressWarnings(maptools::readShapeSpatial(path))
  dubiousGBIF_ids <- c()
  for(species_index in 1:nrow(reference_table)) {
    gbif_subset <- subset(tmp_points, tmp_points$scientificName == reference_table$gbif_name[species_index])
    if(nrow(gbif_subset)!=0) {
      wcvp_subset <- subset(all_vars, all_vars$taxon_name == reference_table$wcvp_name[species_index])
      wcvp_subset <- subset(wcvp_subset, wcvp_subset$introduced==0)
      wcvp_subset <- subset(wcvp_subset, wcvp_subset$extinct==0)
      wcvp_subset <- subset(wcvp_subset, wcvp_subset$location_doubtful==0)
      occ_areas <- wcvp_subset$area_code_l3
      area_plus_buffer <- twgd_data[which(as.character(twgd_data$LEVEL3_COD) %in% occ_areas),]
      if(nrow(area_plus_buffer)>0) {
        coords <- gbif_subset[,c("x","y")]
        sp::coordinates(coords) <- ~ x + y
        answer <- which(is.na(sp::over(coords, area_plus_buffer)[,3]))
        if(length(answer) != 0) {
          dubiousGBIF_ids <- c(dubiousGBIF_ids, as.character(gbif_subset$gbifID[answer]))
        }
      }
    }
    cat(species_index, "\r")
  }
  cleaned_points <- subset(points, !as.character(points$gbifID) %in% dubiousGBIF_ids)
  npoints_end <- nrow(cleaned_points)
  print(paste0(npoints_start - npoints_end, " points removed."))
  return(cleaned_points)
}

# function to deal with a bug in how taxize treats authors names
fix.names.taxize <- function(focal_species_trees) {
  for(name_index in 1:length(focal_species_trees)){
    one_tmp_string <- focal_species_trees[name_index]
    if(any(grepl("[()]", one_tmp_string))){
      splitted_names <- strsplit(one_tmp_string," ")[[1]]
      begin_author <- which(grepl("[()]", splitted_names))[1]
      species_name <- paste0(splitted_names[1:(begin_author-1)], collapse=" ")
      author <- splitted_names[begin_author:length(splitted_names)]
      old_authors <- author[grep("[()]", author)]
      end_first_half <- floor(length(old_authors)/2)
      before <- old_authors[1:end_first_half]
      after <- old_authors[(end_first_half+1):(length(old_authors))]
      if(paste(before,collapse = " ") == paste(after, collapse=" ")) {
        author <- paste(author[1:(length(author)/2)], collapse=" ")
        focal_species_trees[name_index] <- paste0(species_name, " ", author, collapse=" ")
      } else {
        author <- paste(author, collapse=" ")
        focal_species_trees[name_index] <- paste0(species_name, " ", author, collapse=" ")
      }
    }
  }
  return(focal_species_trees)
}


#setwd("~/Desktop/WCVP_special_issue/Patricia_Climbers/climbers")

#-----------------------------
# load wcvp data
dist_sample <- read.table("C:/Users/patri/Desktop/temp/wcvp_names_and_distribution_special_edition_2022/wcvp_distribution.txt", sep="|", header=TRUE, quote = "", fill=TRUE, encoding = "UTF-8")
names_sample <- read.table("C:/Users/patri/Desktop/temp/wcvp_names_and_distribution_special_edition_2022/wcvp_names.txt", sep="|", header=TRUE, quote = "", fill=TRUE, encoding = "UTF-8")
#-----------------------------
# Merge them in one big table
all_vars <- merge(dist_sample, names_sample, by="plant_name_id")

# reference table for taxized names
#-----------------------------
reference_table <- read.csv("taxized_reference_table.csv")
length(unique(reference_table$wcvp_name))
length(unique(reference_table$gbif_name))
##################################
##################################
# A ideia aqui e' ver quantas especies a gente perde se usar o filtro do wcvp no arquivo com os
# pontos que a gente tinha baixado antes
library(data.table)
# Reading gbif file
gbif_data <- fread("C:/Users/patri/Desktop/temp/neotropics_tracheophyte_filtered_gbif.csv") # load  GBIF
# subgbif<-gbif_data[1:100000,]
 length(unique(gbif_data$scientificName))
# [1] 82614 # n species nos neotropicos
# [1] 4340 -> n spp no subgbif

gbif_data <- as.data.frame(gbif_data)[,c(2:4)]
colnames(gbif_data)[1] <- "scientificName" 

# A way to simplify names in table and trees so that species names match again
simplify.names.taxize <- function(names) {
  results <- c()
  for(name_index in 1:length(names)){
    one_tmp_string <- names[name_index]
    splitted_names <- strsplit(one_tmp_string," ")[[1]]
    genus <- splitted_names[1]
    epiphet <- splitted_names[2]
    if(is.na(epiphet)) {
      full_name <- "tip_to_drop" # indet species
    } else {
      if(any(grepl("indet_sp",splitted_names))) {
        full_name <- "tip_to_drop" # indet species
      } else {
        if(stringr::str_detect(epiphet,"[[:upper:]]")) {
          full_name <- "tip_to_drop" # indet species
        } else {
          if(length(splitted_names) > 2) {
            complement <- splitted_names[3:length(splitted_names)]
            if(grepl("[()]", complement[1])) {
              full_name <- paste(c(genus, epiphet), collapse = " ")
            } else {
              if(stringr::str_detect(complement[1],"[[:upper:]]")) {
                full_name <- paste(c(genus, epiphet), collapse = " ")
              } else {
                complement <- subset(complement, !stringr::str_detect(complement,"[[:upper:]]"))
                complement <- subset(complement, !grepl(paste(c("[()]","&","([0-9]+).*$","^ex$"), collapse="|"), complement))
                if(length(complement)==0){
                  full_name <- paste(c(genus, epiphet), collapse = " ")
                } else {
                  full_name <- paste(c(genus, epiphet, complement), collapse = " ")
                }
              }
            } 
          }
        }
      }
    }
    results[name_index] <- full_name
  }
  return(results)
}
library(stringr)
reference_table2<-reference_table
cleaned_points3<-cleaned_points

reference_table2$gbif_name<-simplify.names.taxize(reference_table2$gbif_name)
reference_table2<-subset(reference_table2, gbif_name!="tip_to_drop")

cleaned_points3$scientificName<-simplify.names.taxize(cleaned_points3$scientificName)
cleaned_points3<-subset(cleaned_points3, scientificName!="tip_to_drop")

# Looking at the WCVP table and TDWG to clean GBIF points
#-----------------------------
# If local
# path="/Users/thaisvasconcelos/Desktop/WCVP_special_issue/WCVPtools/wgsrpd-master/level3/level3.shp"
path="C:/Users/patri/Desktop/temp/maps/wgsrpd-master/level3/level3.shp"

twgd_data <- suppressWarnings(maptools::readShapeSpatial(path))
#issues_to_remove <- read.csv("gbif_issues_to_remove.csv")
#cultivated <- read.csv("cultivated_species.csv")

#------ Cleaning steps:
cleaned_points <- gbif_data
# WCVP filtering
subset_reference_table <- subset(reference_table, reference_table$gbif_name %in% unique(cleaned_points$scientificName))
subset_reference_table2 <- subset(reference_table2, reference_table2$gbif_name %in% unique(cleaned_points3$scientificName))
subset_all_vars <- subset(all_vars, all_vars$taxon_name %in% subset_reference_table$wcvp_name)
# cleaned_points <- FilterWCVP(cleaned_points, all_vars, subset_reference_table, twgd_data,  lon="lon", lat="lat") # This will filter the GBIF points acording to WCVP for species
# function(points, all_vars, reference_table, twgd_data, lon="decimalLongitude", lat="decimalLatitude")

npoints_start <- nrow(cleaned_points)
tmp_points = as.data.frame(cleaned_points)
colnames(tmp_points)[colnames(tmp_points)=="lon"] <- "x"
colnames(tmp_points)[colnames(tmp_points)=="lat"] <- "y"
tmp_points = subset(tmp_points, !is.na(tmp_points$x))
tmp_points = subset(tmp_points, !is.na(tmp_points$y))
# Load shape files and make sure they have the same name as the WCVP column with the TDWG areas
# twgd_data <- suppressWarnings(maptools::readShapeSpatial(path))
cleaned_points<-vector("list",nrow(subset_reference_table))
excluded_spp<-c(length(nrow(subset_reference_table)))
for(species_index in 1:nrow(subset_reference_table)) {
  gbif_subset <- subset(tmp_points, tmp_points$scientificName == subset_reference_table$gbif_name[species_index])
  if(nrow(gbif_subset)!=0) {
    wcvp_subset <- subset(all_vars, all_vars$taxon_name == subset_reference_table$wcvp_name[species_index]) # aqui usa o wcvp names
    wcvp_subset <- subset(wcvp_subset, wcvp_subset$introduced==0)
    wcvp_subset <- subset(wcvp_subset, wcvp_subset$extinct==0)
    wcvp_subset <- subset(wcvp_subset, wcvp_subset$location_doubtful==0)
    occ_areas <- wcvp_subset$area_code_l3
    area_plus_buffer <- twgd_data[which(as.character(twgd_data$LEVEL3_COD) %in% occ_areas),]
    if(nrow(area_plus_buffer)>0) {
      coords <- gbif_subset[,c("x","y")]
      sp::coordinates(coords) <- ~ x + y
      answer <- which(is.na(sp::over(coords, area_plus_buffer)[,3])) # coord (x) são os pontos, area_plus_buffer (y) é a área de ocorrencia segundo o POWO
      if(length(answer) != 0) {
        cleaned_points[[species_index]]<-gbif_subset[-answer,]
        #dubiousGBIF_ids <- c(dubiousGBIF_ids, as.character(gbif_subset$gbifID[answer]))
      }
      if(length(answer) == 0) {
        cleaned_points[[species_index]]<-gbif_subset
      }
    }
  }
  
  cat(species_index, "\r")
}

cleaned_points<-do.call("rbind",cleaned_points)

#cleaned_points <- subset(cleaned_points, !as.character(cleaned_points$gbifID) %in% dubiousGBIF_ids)
#npoints_end <- nrow(cleaned_points)
#print(paste0(npoints_start - npoints_end, " points removed."))
#return(cleaned_points)

length(unique(cleaned_points$scientificName))
length(unique(gbif_data$scientificName))
length(which(unique(cleaned_points$scientificName) %in% reference_table$gbif_name))
which(cleaned_points$scientificName == "Anthurium pahumense Cerón & Croat")
which(reference_table$gbif_name == "Anthurium pahumense Cerón & Croat")

# quantas species no final?
# [1] 2593 com o subset


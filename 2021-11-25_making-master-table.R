# Making master table 
rm(list=ls())
library(data.table)
source("neotropical_climbers_functions.R") 
# setwd("~/climbers")
# setwd("~/Desktop/Colabs/Patricia_Climbers/climbers")

# Load all Neotropical stuff 
gbif_dir <- paste0(getwd(), "/full_gbif_quest")
full_list <- fread(paste0(gbif_dir, "/neotropics_tracheophyte_filtered_gbif.csv"))
full_list <- full_list[,-1]

# Now load table climbers
climbers <- read.csv("Data/climber_database.csv", stringsAsFactors = F)
taxized_climbers <- gbif.taxize(climbers$Species)
climbers$taxized_names <- taxized_climbers
write.csv(climbers, file="Data/climber_database.csv")

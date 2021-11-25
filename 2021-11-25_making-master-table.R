# Making master table 
rm(list=ls())
library(data.table)
source("neotropical_climbers_functions.R") 
# setwd("~/climbers")
# setwd("~/Desktop/Colabs/Patricia_Climbers/climbers")

# Taxize climbers accordinf to GBIF
climbers <- read.csv("Data/climber_database.csv", stringsAsFactors = F)
taxized_climbers <- gbif.taxize(climbers$Species)
climbers$taxized_names <- taxized_climbers
# write.csv(climbers, file="Data/climber_database.csv", row.names = F)

# Load back
climbers <- read.csv("Data/climber_database.csv", stringsAsFactors = F)
climbers$taxized_names[which(duplicated(climbers$taxized_names))]
#[1] "Aphelandra tomentosa Lindau"                                                 
#[2] "Ditassa conceptionis Fontella"                                               
#[3] "Matelea purpureolineata Woodson"                                             
#[4] "Mikania pterocaula Sch.Bip. ex Klatt"                                        
#[5] "Adenocalymma pseudopatulum (A.H.Gentry) L.G.Lohmann L.G.Lohmann (A.H.Gentry)"
#[6] "Tontelea cylindrocarpa (A.C.Sm.) A.C.Sm. A.C.Sm. (A.C.Sm.)"                  
#[7] "Tontelea martiana (Peyr.) A.C.Sm. A.C.Sm. (Peyr.)"                           
#[8] "Cayaponia tubulosa Cogn."                                                    
#[9] "Elateriopsis oerstedii Pittier"                                              
#[10] "Dioscorea perenensis R.Knuth"                                                
#[11] "Niedenzuella castanea (Cuatrec.) W.R.Anderson W.R.Anderson (Cuatrec.)"       
#[12] "Passiflora ernestii Harms"                                                   
#[13] "Passiflora hahnii (E.Fourn.) Mast. Mast. (E.Fourn.)"                         
#[14] "Paullinia bilobulata Radlk."                                                 
#[15] "Paullinia eriocarpa Triana & Planch."                                        
#[16] "Serjania pyramidata Radlk."                                                  
#[17] "Serjania regnellii Schltdl."   

#---------------------
# Making master table
# Load all Neotropical stuff 
gbif_dir <- paste0(getwd(), "/full_gbif_quest")
full_list <- fread(paste0(gbif_dir, "/neotropics_tracheophyte_filtered_gbif.csv"))
full_list <- full_list[,-1]

# length(unique(full_list$species))
# [1] 82614 species in the dataset

# length(which(unique(full_list$species) %in% climbers$taxized_names))
# [1] 4485 species in the climbers dataset have good points

all_neotropical_species <- unique(full_list$species)
climbers_w_points <- intersect(all_neotropical_species, climbers$taxized_names)
full_list$mechanism <- NA

for(i in 1:length(climbers_w_points)) {
  one_climber <- climbers_w_points[i]
  mechanism <- paste0("Mechanism_", climbers$CM[which(climbers$taxized_names==one_climber)])
  full_list$mechanism[which(full_list$species==climbers_w_points[i])] <- mechanism
  print(i)
}

full_list$mechanism[which(is.na(full_list$mechanism))] <- "not_a_climber"
write.csv(full_list, file="full_gbif_quest/master_table.csv", row.names = F)

#table(full_list$mechanism)





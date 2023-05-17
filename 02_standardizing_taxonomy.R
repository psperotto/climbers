# Script for standardizing taxonomy using RBG Kew's World Checklist of Vascular Plants

######################################################################
##### PART 1: matching the names of our initial database to WCVP #####

## Reading initial data
wcvp_names <- read.table("C:/Users/patri/Desktop/temp/wcvp_names_and_distribution_special_edition_2022/wcvp_names.txt", sep="|", header=T, quote = "", fill=TRUE, encoding = "UTF-8")
climbers_temp <- read.csv("C:/Users/patri/Desktop/temp/database.csv", stringsAsFactors = F) # This is our initial database based on Pedro Acevedo's list of neotropical climbing plants.
# "CM" stands for Climbing Mechanism and each number refers to a different CM following Sperotto et al. 2020
# 1 = "Twiner"; 2 = "Tendrils"; 3 = "Scrambler"; 4 = "Root-climber; 5 = "Prehensile branches";
# 6 = "Prehensile petioles"; 7 = "Hooks or grapnels"; 8 = "Twining inflorescences"

## Filtering WCVP for names that are in our initial database
library(dplyr)
x <- wcvp_names %>% filter(taxon_name %in% climbers_temp$Species)

## Filtering WCVP to contain only the accepted plant name ids of 'x'
acc_climbers_temp <- wcvp_names %>% filter(plant_name_id %in% x$accepted_plant_name_id)
length(unique(acc_climbers_temp$taxon_name)) # 8959 accepted names 

## Considering that in taxonomic lists there might be misspellings of names, we searched for names that were in our initial database but not in WCVP in order to optimize the number of species considered in our analyses
table(climbers_temp$Species %in% wcvp_names$taxon_name) 
# FALSE  TRUE 
#   241  9581 -> 241 names were not in WCVP
a<-which(!climbers_temp$Species %in% wcvp_names$taxon_name)
b<-climbers_temp[a,] 
# write.csv(b, file="notinwcvp.csv") -> each name in this file was manually searched and corrected in the "database.csv" file, which was saved as "database_corrected.csv"

## Reading the corrected database
climbers<-read.csv("C:/Users/patri/Desktop/temp/database_corrected.csv", stringsAsFactors = F)
table(climbers$Species %in% wcvp_names$taxon_name) 
# FALSE  TRUE 
#   45  9764 -> 45 names still not found in WCVP

## Excluding names not in WCVP
c<-which(!climbers$Species %in% wcvp_names$taxon_name)
climbers<-climbers[-c,] 
head(climbers) 
# write.csv(climbers, file ="climbers.csv")

# "climbers.csv" is the list with the corrected species names and their respective genera, families and climbing mechanisms ("CM", see line 9 of this script). 
# All species names here can be found in the WCVP, but are not necessarily accepted species. Taxon status will be addressed in Part 2 of this script (below).

##############################################################################################################
##### PART 2: assigning the correct climbing mechanisms to only accepted species names according to WCVP #####

## Filtering WCVP to contain only the species names in 'climbers'
library(dplyr)
y <- wcvp_names %>% filter(taxon_name %in% climbers$Species)

## Adding a climbing mechanism (i.e., "CM") column to the filtered WCVP and assigning the climbing mechanisms to each name
y[,"CM"] <- NA
for (k in 1:length(climbers$Species)) {
  i <- which(y$taxon_name == climbers$Species[k])
  y$CM[i] <- climbers$CM[k]
}
# Reminder:
# "CM" stands for Climbing Mechanism and each number refers to a different CM following Sperotto et al. 2020
# 1 = "Twiner"; 2 = "Tendrils"; 3 = "Scrambler"; 4 = "Root-climber; 5 = "Prehensile branches";
# 6 = "Prehensile petioles"; 7 = "Hooks or grapnels"; 8 = "Twining inflorescences"

## Checking if all accepted plant name IDs of 'y' are in the plant name IDs of the full WCVP
which(!y$accepted_plant_name_id %in% wcvp_names$plant_name_id) 
print(y$taxon_status[which(!y$accepted_plant_name_id %in% wcvp_names$plant_name_id)]) 
# all of these are Unplaced taxa

## Excluding Unplaced taxa
y<-y[-which(!y$accepted_plant_name_id %in% wcvp_names$plant_name_id),] # excluding unplaced taxa

## Filtering WCVP to contain only names corresponding to the accepted plant name IDs of 'y'
climbers2 <- wcvp_names %>% filter(plant_name_id %in% y$accepted_plant_name_id)

## Checking if all accepted plant IDs of 'y' correspond to only one CM
y_acc<-unique(y$accepted_plant_name_id)
for (k in 1:length(y_acc)) {
  i <- which(y$accepted_plant_name_id == y_acc[k])
  if (length(unique(y$CM[i])) > 1) {
    print(y_acc[k])
  }
} # The printed IDs are accepted plant IDs that, due to certain names in our initial database being synonymous, ended up being assigned more than one climbing mechanism 

## Accepted plant id of species that are dubious regarding their climbing mechanism (you can use 'subset(wcvp_names, accepted_plant_name_id == "insert ID here") to check the names)
b<-c("300711-wcs", "1118364-az", "1079222-az", "1123981-az", "1016262-az", "733369-az", "732736-az", 
     "732993-az", "289086-wcs", "520470-wcs", "457445-az", "913748-az", "366371-az", "506018-az", 
     "602525-az", "602531-az", "602541-az", "118041-wcs", "118054-wcs", "1032600-az")

## Excluding IDs of species with more than one climbing mechanism from 'wcvp_climber's and 'y'
for (i in 1:length(b)){
  climbers2<-climbers2 %>% filter(plant_name_id !=b[i])
}
for (i in 1:length(b)){
  y<-y %>% filter(accepted_plant_name_id !=b[i])
}

## Assigning the correct climbing mechanism to each accepted species
library(beepr) # not necessary, just fun! :)
climbers2[,"CM"]<-NA
y_acc2<-y_acc[which(!y_acc %in% b)]
for (k in 1:length(y_acc2)) {
  i <- which(y$accepted_plant_name_id == y_acc2[k])
  for (p in 1:length(climbers2$plant_name_id)){
    if (climbers2$plant_name_id[p] == y_acc2[k]) {
      climbers2$CM[p]<-y$CM[i][[1]]
    }
  }
}
beep("mario")

## Creating the the final database
library(xlsx)
APGIV <- read.xlsx("C:/Users/patri/Desktop/temp/APGIV.xlsx", sheetIndex = 1)
head(APGIV) # List of angiosperm families, orders and clades following APG IV

## Leaving only accepted species in 'climbers2'
climbers2<-subset(climbers2, taxon_rank=="Species")
climbers2<-subset(climbers2, taxon_status=="Accepted")

## Checking if every family in 'climbers2' is in APG IV and correcting if they're not 
which(!climbers2$family %in% APGIV$family)
# [1] 7811 8969 -> THAIS: tem alguma coisa mto estranha que acontece aqui, se tu procura a linha 8969 diz que é uma Violaceae, sendo que Violaceae ta no WCVP e no APGIV

climbers2$family[7811] # Polypodiaceae: it's a fern, must be excluded!
climbers2<-climbers2[-7811,]

which(!climbers2$family %in% APGIV$family)
# [1] 8968 -> THAIS: aí depois de excluir a linha 7811, se checar de novo diz que é a linha 8968 que não ta na APGIV, o que é verdade! 
climbers2$family[8968] # Adoxaceae in APG IV = Viburnaceae in WCVP
climbers2$family[8968] <- c("Adoxaceae")

## Checking again
which(!climbers2$family %in% APGIV$family) # 'integer(0)', all good!

## Selecting only the columns that are of interest for us and
climbers_final<-climbers2[,c("plant_name_id", "taxon_status", "family","genus","lifeform_description","taxon_name","taxon_authors","CM")]

##############################################################################
##### PART 3: finding which species the WCVP already considered climbers #####
## Adding a column saying if WCVP consider the species a climbers
wcvp_climbers <- subset(climbers_final, grepl(paste0(c("Climb","Liana","Scrambl"), collapse="|"), climbers_final$lifeform_description))
climbers_final[,"is_climber_WCVP"]<-NA
for (k in 1:length(wcvp_climbers$taxon_name)) {
  i <- which(climbers_final$taxon_name == wcvp_climbers$taxon_name[k])
  climbers_final$is_climber_WCVP[i] <- "yes"
}

## Saving the final complete database
saveRDS(climbers_final, file="climbers_final.Rdata")
write.csv(climbers_final, file="climbers_final.csv")

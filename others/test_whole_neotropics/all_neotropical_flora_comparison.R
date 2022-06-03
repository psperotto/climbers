
# Let's work with the whole flora:
library(data.table)
library(monographaR)

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




# get AI for all:
library(maptools)
library(data.table)
data("wrld_simpl")

source("00_utility_functions.R")
source("/Users/thaisvasconcelos/Desktop/WCVP_special_issue/WCVPtools/WCVPtools_functions.R")

# Directory to save preliminary datasets:
climate_data.dir <- "./climate_data"
# Directory where climate layers are:
climate_layers.dir <- "./climate_layers"

# 1. Thinning occurence data first
thinned_points <- Thinning(all_angios, species="species", lat = "lat", lon="lon", n = 3)
write.csv(thinned_points, file="full_gbif_quest/neotropics_tracheophyte_filtered_gbif_thinned.csv", row.names=F)
# 2. Getting summary statistics of climatic variables for each species
allpoints <- ClimateFromPoint_custom(thinned_points, species="species",lon="lon", lat="lat", layerdir = climate_layers.dir)
write.csv(allpoints, file=paste0(climate_data.dir, "/neotropics_allpoints.csv"), row.names=F)
summstats <- GetClimateSummStats_custom(allpoints, type="raw")
write.csv(summstats, file=paste0(climate_data.dir, "/neotropics_summstats_raw.csv"), row.names=F)  
summstats <- GetClimateSummStats_custom(allpoints, type="transformed")
write.csv(summstats, file=paste0(climate_data.dir, "/neotropics_summstats.csv"), row.names=F)

# Different ways in which we can test this:
# (1) Phyloanovas
# Higher AI will represent higher potential for vegetation growth (i.e. closed canopy environments)

# Hypothesis: Groups with different climbing mechanisms have different support diameter requirements, 
# and thus the distribution of lineages with different climbing mechanism is influenced by the 
# surrounding vegetation type

ai<-summstats[[4]]
ai <- as.data.frame(ai)


climber_points <- subset(ai, ai$species %in% all_climbers)
non_climber_points <- subset(ai, !ai$species %in% all_climbers)

climber_points$climber <- "climber"
non_climber_points$climber <- "not_a_climber"

merged_table_for_plot  <- rbind(climber_points, non_climber_points)

hist(as.numeric(climber_points$mean_aridity), breaks=50)
hist(as.numeric(non_climber_points$mean_aridity), breaks=50)

library(ggplot2)

merged_table_for_plot %>% 
  ggplot(aes(x=climber,y=mean_aridity, fill=climber)) +
  geom_boxplot() + geom_jitter(width=0.1,alpha=0.2) 

colnames(merged_table_for_plot)

set.seed(1)
Ixos=as.numeric(climber_points$mean_aridity)
Primadur=as.numeric(non_climber_points$mean_aridity)


# First distribution
hist(log(Ixos), breaks=30, col=rgb(1,0,0,0.5), xlab="height", 
     ylab="nbr of plants", main="distribution of height of 2 durum wheat varieties" )

# Second with add=T to plot on top
hist(log(Primadur), breaks=30, xlim=c(0,11), col=rgb(0,0,1,0.5), add=T)


# Add legend
legend("topright", legend=c("Ixos","Primadur"), col=c(rgb(1,0,0,0.5), 
                                                      rgb(0,0,1,0.5)), pt.cex=2, pch=15 )




climber_mech <- read.csv("Data/climber_database.csv")
colnames(climber_mech) <- c("family","genus","species","CM","taxized_names")
merged_table <- merge(climber_mech, climber_points, by.x = "taxized_names", by.y = "scientificName")


ai_table <- summstats[[4]]
ai_table <- as.data.frame(ai_table)
ai_table$mean_aridity <- as.numeric(ai_table$mean_aridity)
merged_table <- merge(climber_mech, ai_table, by.x="taxized_names", by.y="species")

boxplot(merged_table$mean_aridity ~ merged_table$CM)

# Load tree for PCMs
library(ape)
library(phytools)
tree <- readRDS("Data/taxized_GBMB.Rdata")

tree$tip.label

tree$tip.label <- unname(tree$tip.label)

pruned_tree <- keep.tip(tree, grep("Coprosma",tree$tip.label))
plot(pruned_tree, cex=0.5)
axisPhylo()

subset_merged_table <- subset(merged_table, merged_table$taxized_names%in%tree$tip.label)
subset_merged_table <- subset(subset_merged_table, !duplicated(subset_merged_table$taxized_names))
pruned_tree <- keep.tip(tree, intersect(subset_merged_table$taxized_names, tree$tip.label))
pruned_tree$tip.label <- paste0(unlist(lapply(strsplit(pruned_tree$tip.label, " "), "[[", 1)), " ", unlist(lapply(strsplit(pruned_tree$tip.label, " "), "[[", 2)))
subset_merged_table$taxized_names <- paste0(unlist(lapply(strsplit(subset_merged_table$taxized_names, " "), "[[", 1)), " ", unlist(lapply(strsplit(subset_merged_table$taxized_names, " "), "[[", 2)))


aridity <- subset_merged_table$mean_aridity
mechanisms <- subset_merged_table$CM
names(aridity) <- names(mechanisms) <- subset_merged_table$taxized_names

class(mechanisms)
mechanisms[which(mechanisms==1)] <- "a"
mechanisms[which(mechanisms==2)] <- "b"
mechanisms[which(mechanisms==3)] <- "c"
mechanisms[which(mechanisms==4)] <- "d"
mechanisms[which(mechanisms==5)] <- "e"
mechanisms[which(mechanisms==6)] <- "f"
mechanisms[which(mechanisms==7)] <- "g"
mechanisms[which(mechanisms==8)] <- "h"

phyanova <- phylANOVA(pruned_tree, mechanisms, aridity)


# 


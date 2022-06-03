# Load points and extracting variables
library(phytools)

#-------------------------------
# Getting climate data
#-------------------------------
# Load cleaned points back:
data_dir <- paste0(getwd(),"/Data")
all_taxized_points_files <- list.files(data_dir, "taxized_points.Rdata")
all_taxized_points <- lapply(paste0(data_dir, "/",all_taxized_points_files), readRDS)
names(all_taxized_points) <- gsub("_taxized_points.Rdata","", all_taxized_points_files)

all_cleaned_points = all_taxized_points
# Directory to save preliminary datasets:
div_space_dir <- "./3_div_in_space"
# Directory where climate layers are:
climate_layers.dir <- "./climate_layers"

{; for(family_index in 1:length(all_cleaned_points)) {
  # 1. Thinning occurence data first
  thinned_points <- Thinning(all_cleaned_points[[family_index]], species="species", lat = "lat", lon="lon", n = 1)
  # 2. Getting summary statistics of climatic variables for each species
  allpoints <- ClimateFromPoint_custom(thinned_points, species="species",lon="lon", lat="lat", layerdir = climate_layers.dir)
  write.csv(allpoints, file=paste0(data_dir, "/", names(all_cleaned_points)[family_index], "_allpoints.csv"))
  summstats <- GetClimateSummStats_custom(allpoints, type="raw")
  write.csv(summstats, file=paste0(data_dir, "/", names(all_cleaned_points)[family_index], "_summstats_raw.csv"))
  summstats <- GetClimateSummStats_custom(allpoints, type="transformed")
  write.csv(summstats, file=paste0(data_dir, "/", names(all_cleaned_points)[family_index], "_summstats.csv"))
}
  beepr::beep("complete"); } 

#-------------------------------
# Getting organized table for phyloANOVA
#-------------------------------
data_dir <- paste0(getwd(),"/Data")
summstats_files <- list.files(data_dir, "summstats.csv")
summstats <- lapply(paste0(data_dir, "/", summstats_files), read.csv)
names(summstats) <- gsub("_summstats.csv","", summstats_files)

all_cleaned_points <- list()
{; for(family_index in 1:length(summstats)) {
  group <- names(summstats)[family_index]
  group_summstats <- summstats[[grep(group, names(summstats))]]
  group_summstats$Mechanism <- group
  cleaned_table <- group_summstats[,c("species","mean_temp","se_temp","within_sp_var_temp","mean_prec","se_prec","within_sp_var_prec",
                                   "mean_pet","se_pet","within_sp_var_pet","mean_aridity","se_aridity", "within_sp_var_aridity","Mechanism")]
  
  if(any(is.na(cleaned_table$mean_prec))) { 
    cleaned_table <- cleaned_table[-which(is.na(cleaned_table$mean_prec)),]
  }
  if(any(cleaned_table$mean_prec==0)) {
    cleaned_table <- cleaned_table[-which(cleaned_table$mean_prec == 0),]
  }
  all_cleaned_points[[family_index]] <- cleaned_table
  write.csv(cleaned_table, file=paste0(data_dir,"/",group, "_niche.csv"))
}
  beepr::beep("complete"); } 

all_cleaned_points <- do.call(rbind, all_cleaned_points)
write.csv(all_cleaned_points, file=paste0(data_dir,"/all_mechanisms_niche.csv"))
#---------------------------------------------
#---------------------------------------------
#---------------------------------------------
# Now let's do the same for all neotropical points 
full_list <- fread(paste0(gbif_dir, "/neotropics_tracheophyte_filtered_gbif.csv"))

# Directory to save preliminary datasets:
data_dir <- paste0(getwd(),"/Data")
# Directory where climate layers are:
climate_layers.dir <- "./climate_layers"

# 1. Thinning occurence data first
thinned_points <- Thinning(full_list, species="species", lat = "lat", lon="lon", n = 1)
# 2. Getting summary statistics of climatic variables for each species
allpoints <- ClimateFromPoint_custom(thinned_points, species="species",lon="lon", lat="lat", layerdir = climate_layers.dir)
write.csv(allpoints, file=paste0(data_dir, "/neotropics_allpoints.csv"))
summstats <- GetClimateSummStats_custom(allpoints, type="raw")
write.csv(summstats, file=paste0(data_dir, "/neotropics_summstats_raw.csv"))
summstats <- GetClimateSummStats_custom(allpoints, type="transformed")
write.csv(summstats, file=paste0(data_dir, "/neotropics_summstats.csv"))

beepr::beep("complete")


#---------------------------------------------
#---------------------------------------------
#---------------------------------------------
# phylANOVA
# Load tree with taxized tips
phy <- phy.list(input.dir=getwd(), names="GBMB", search.for=".taxized.tre")[[1]] 
phy$tip.label <- gsub("_"," ",phy$tip.label)

# Make datasets match:
all_cleaned_points <- subset(all_cleaned_points, !is.na(all_cleaned_points$mean_aridity))
all_cleaned_points <- subset(all_cleaned_points, !is.na(all_cleaned_points$mean_temp))
all_cleaned_points <- subset(all_cleaned_points, !is.na(all_cleaned_points$mean_prec))
all_cleaned_points <- subset(all_cleaned_points, !is.na(all_cleaned_points$mean_pet))

phy_pruned <- keep.tip(phy, which(phy$tip.label %in% all_cleaned_points$species))
all_cleaned_points_pruned <- subset(all_cleaned_points, all_cleaned_points$species %in% phy_pruned$tip.label)

trait1 <- all_cleaned_points$mean_aridity
mechanisms <- all_cleaned_points$Mechanism
names(trait1) <- names(mechanisms) <- all_cleaned_points$species
phy_pruned <- force.ultrametric(phy_pruned)


anova_results <- phylANOVA(phy_pruned, mechanisms, trait1, nsim=1000)

boxplot(all_cleaned_points$mean_aridity~all_cleaned_points$Mechanism)



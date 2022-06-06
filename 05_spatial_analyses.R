# rm(list=ls())
setwd("~/Desktop/WCVP_special_issue/Patricia_Climbers/climbers")

library(maptools)
library(data.table)
data("wrld_simpl")

source("00_utility_functions.R")
source("/Users/thaisvasconcelos/Desktop/WCVP_special_issue/WCVPtools/WCVPtools_functions.R")


# Hypothesis: Groups with different climbing mechanisms have different support diameter requirements, 
# and thus the distribution of lineages with different climbing mechanism is influenced by the 
# surrounding vegetation type

climber_points <- as.data.frame(fread("gbif_climbers/climbers_cleaned_points.csv"))
climber_mech <- readRDS("wcvp_nt_climbers_final.Rdata")
#climber_mech <- read.csv("Data/climber_database.csv")
climber_mech$taxon_and_author <- paste0(climber_mech$taxon_name, " ", climber_mech$taxon_authors)
#colnames(climber_mech) <- c("family","genus","species","CM","taxized_names")
merged_table <- merge(climber_mech, climber_points, by.x = "taxon_and_author", by.y = "scientificName")

# Directory to save preliminary datasets:
climate_data.dir <- "./climate_data"
# Directory where climate layers are:
climate_layers.dir <- "./climate_layers"

# 1. Thinning occurence data first
thinned_points <- Thinning(merged_table, species="taxon_and_author", lat = "decimalLatitude", lon="decimalLongitude", n = 3)
# write.csv(thinned_points, file=paste0(climate_data.dir, "/climbers_thinned_points.csv"), row.names=F)
# thinned_points <- read.csv(paste0(climate_data.dir, "/climbers_thinned_points.csv"))

# 2. Getting summary statistics of climatic variables for each species
allpoints <- ClimateFromPoint_custom(thinned_points, species="species",lon="lon", lat="lat", layerdir = climate_layers.dir)
write.csv(allpoints, file=paste0(climate_data.dir, "/climbers_allpoints.csv"), row.names=F)
summstats <- GetClimateSummStats_custom(allpoints, type="raw")
write.csv(summstats, file=paste0(climate_data.dir, "/climbers_summstats_raw.csv"), row.names=F)  
summstats <- GetClimateSummStats_custom(allpoints, type="transformed")
write.csv(summstats, file=paste0(climate_data.dir, "/climbers_summstats.csv"), row.names=F)

# Different ways in which we can test this:
# (1) Phyloanovas
# Higher AI will represent higher potential for vegetation growth (i.e. closed canopy environments)
summstats <- read.csv(paste0(climate_data.dir, "/climbers_summstats.csv"))

ai_table <- summstats[,c("species","mean_aridity")]
ai_table$mean_aridity <- as.numeric(ai_table$mean_aridity)
merged_table <- merge(climber_mech, ai_table, by.x="taxon_and_author", by.y="species")

#boxplot(exp(merged_table$mean_aridity) ~ merged_table$CM)

# Load tree for PCMs
library(ape)
library(phytools)
tree <- readRDS("trees/taxized_GBMB.Rdata")
tree$tip.label <- unname(tree$tip.label)

subset_merged_table <- subset(merged_table, merged_table$taxon_and_author%in%tree$tip.label)
subset_merged_table <- subset(subset_merged_table, !duplicated(subset_merged_table$taxon_and_author))
pruned_tree <- keep.tip(tree, intersect(subset_merged_table$taxon_and_author, tree$tip.label))
pruned_tree$tip.label <- paste0(unlist(lapply(strsplit(pruned_tree$tip.label, " "), "[[", 1)), " ", unlist(lapply(strsplit(pruned_tree$tip.label, " "), "[[", 2)))
subset_merged_table$taxon_and_author <- paste0(unlist(lapply(strsplit(subset_merged_table$taxon_and_author, " "), "[[", 1)), " ", unlist(lapply(strsplit(subset_merged_table$taxon_and_author, " "), "[[", 2)))

# phyloanova
aridity <- subset_merged_table$mean_aridity
mechanisms <- subset_merged_table$CM
names(aridity) <- names(mechanisms) <- subset_merged_table$taxon_and_author

mechanisms[which(mechanisms==1)] <- "Twining"
mechanisms[which(mechanisms==2)] <- "Tendrils"
mechanisms[which(mechanisms==3)] <- "Simple_scrambling"
mechanisms[which(mechanisms==4)] <- "Adhesive_roots"
mechanisms[which(mechanisms==5)] <- "Prehensile_branches"
mechanisms[which(mechanisms==6)] <- "Prehensile_petioles"
mechanisms[which(mechanisms==7)] <- "Hooks_or_grapnels"

phyanova <- phylANOVA(pruned_tree, mechanisms, aridity, nsim=10000,p.adj="bonferroni")

sink("phylanova_results.txt")
print(phyanova)
sink()

# Figure
boxplot(exp(subset_merged_table$mean_aridity) ~ subset_merged_table$CM)

#pdf(paste0(wd, "/figures/AICcWts/aicw.pdf"), width=8, height=5)

table_for_plot <- subset_merged_table
table_for_plot$CM <- as.character(table_for_plot$CM)
table_for_plot$mean_aridity <- exp(table_for_plot$mean_aridity)

table_for_plot$CM[which(table_for_plot$CM==1)] <- "Twining"
table_for_plot$CM[which(table_for_plot$CM==2)] <- "Tendrils"
table_for_plot$CM[which(table_for_plot$CM==3)] <- "Simple scrambling"
table_for_plot$CM[which(table_for_plot$CM==4)] <- "Adhesive roots"
table_for_plot$CM[which(table_for_plot$CM==5)] <- "Prehensile branches"
table_for_plot$CM[which(table_for_plot$CM==6)] <- "Prehensile petioles"
table_for_plot$CM[which(table_for_plot$CM==7)] <- "Hooks or grapnels"


library(ggplot2)
library(ggridges)

# possible plots
#plot_comparisons <- ggplot(table_for_plot, aes(x=mean_aridity, y=CM, group=CM)) + 
#  geom_density_ridges() 

pdf("aridity_climbers.pdf",height=3,width=5.5)
plot_comparisons <- ggplot(table_for_plot, aes(x=CM, y=mean_aridity, fill=CM)) + 
  geom_boxplot(lwd=0.3, outlier.size=0.25, alpha=0.8) + 
  geom_jitter(colour = 2, alpha=0.3, size=0.9) +
  coord_flip()  + 
  scale_fill_brewer(palette="Pastel1")+
  theme_bw(base_size = 8) + 
  theme(legend.position = "none") + 
  annotate(geom="text", x=.9, y=0.9, size=3, hjust=0, label="") + 
  xlab("") +
  ylab("Mean aridity index") +
  theme(axis.text.y = element_text(colour = 'black', size = 10),
        axis.text.x = element_text(colour = 'black', size = 10),
        axis.title.x = element_text(colour = 'black', size = 10)) 
plot_comparisons
dev.off()

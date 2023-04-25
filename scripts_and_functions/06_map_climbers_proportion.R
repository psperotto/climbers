#setwd("~/Desktop/WCVP_special_issue/Patricia_Climbers/climbers")
#setwd("C:/Users/patri/Desktop/temp")
#rm(list=ls())

library(data.table)
library(monographaR)
library(maptools)
library(raster)

##### Reading data #####
all_angios <- fread("neotropics_tracheophyte_filtered_gbif.csv")
all_angios <- all_angios[,2:4]

<<<<<<< HEAD
climbers_final <- readRDS("climbers_final.Rdata")

# Expert-defined Neotropic areas with 7 areas - Antonelli et al. 2018
NTshp<-shapefile("C:/Users/patri/Desktop/temp/Zenodo_scripts_and_data_Antonelli-et-al_PNAS/input/shapes/neotropics_areas.shp")
# Hybrid Neotropics areas (mixed WWF biomes and expert-defined areas) - Antonelli et al 2018
hybridshp<-shapefile("C:/Users/patri/Desktop/temp/Zenodo_scripts_and_data_Antonelli-et-al_PNAS/input/shapes/hibrid_wwf_experts.shp")
=======
wcvp_climbers <- readRDS("wcvp_nt_climbers_final.Rdata")
all_climbers <- paste(wcvp_climbers$taxon_name, wcvp_climbers$taxon_authors, sep=" ")

>>>>>>> 6e9ab718bea0c1fb821d2bb21cd2ff94f0d81e3f

##### Generating distribution points and basic maps #####

<<<<<<< HEAD
pal <- hcl.colors(1000, palette = "YlGnBu", rev=T)
=======
climber_points <- as.data.frame(climber_points)
all_angios <- as.data.frame(all_angios)
>>>>>>> 6e9ab718bea0c1fb821d2bb21cd2ff94f0d81e3f

## All angiosperms
all_angios<-as.data.frame(all_angios)
neotropics_diversity_raster <- mapDiversity(all_angios, resolution = 1.2, plot=F,plot.with.grid = F)

## All climbers
all_climbers <- paste(climbers_final$taxon_name, climbers_final$taxon_authors, sep=" ")
climber_points <- as.data.frame(subset(all_angios, all_angios$species %in% all_climbers))
climber_diversity_raster <- mapDiversity(climber_points, resolution = 1.2, plot=F,plot.with.grid = F)
#plot(climber_diversity_raster, col=pal)
#plot(hybridshp, add=T)

## Climbing mechanism

# Twining
tw<-subset(climbers_final, CM=="1")
tw_names<-paste(tw$taxon_name, tw$taxon_authors, sep = " ")
twining_points<-as.data.frame(subset(all_angios, all_angios$species %in% tw_names))
twining_div_raster<-mapDiversity(twining_points,resolution = 1.2, plot=F,plot.with.grid = F)
#plot(twining_div_raster, col=pal)
#plot(hybridshp, add=T)

# Tendrils
te<-subset(climbers_final, CM=="2")
te_names<-paste(te$taxon_name, te$taxon_authors, sep=" ")
tendrils_points<-as.data.frame(subset(all_angios, all_angios$species %in% te_names))
tendrils_div_raster<-mapDiversity(tendrils_points,resolution = 1.2, plot=F,plot.with.grid = F)
#plot(tendrils_div_raster, col=pal)
#plot(hybridshp, add=T)

# Scrambling
sc<-subset(climbers_final, CM=="3")
sc_names<-paste(sc$taxon_name, sc$taxon_authors, sep=" ")
scrambling_points<-as.data.frame(subset(all_angios, all_angios$species %in% sc_names))
scrambling_div_raster<-mapDiversity(scrambling_points, resolution = 1.2, plot=F,plot.with.grid = F)
#plot(scrabling_div_raster, col=pal)
#plot(hybridshp, add=T)

# Adhesive roots
ar<-subset(climbers_final, CM=="4")
ar_names<-paste(ar$taxon_name, ar$taxon_authors, sep=" ")
roots_points<-as.data.frame(subset(all_angios, all_angios$species %in% ar_names))
roots_div_raster<-mapDiversity(roots_points, resolution = 1.2, plot = F, plot.with.grid = F)
#plot(roots_div_raster, col=pal)
#plot(hybridshp, add=T)

# Prehensile branches
pb<-subset(climbers_final, CM=="5")
pb_names<-paste(pb$taxon_name, pb$taxon_authors, sep=" ")
branches_points<-as.data.frame(subset(all_angios, all_angios$species %in% pb_names))
branches_div_raster<-mapDiversity(branches_points, resolution = 1.2, plot = F, plot.with.grid = F)
#plot(branches_div_raster,col=pal)
#plot(hybridshp, add=T)

# Prehensile leaves
pl<-subset(climbers_final, CM=="6")
pl_names<-paste(pl$taxon_name, pl$taxon_authors, sep=" ")
leaves_points<-as.data.frame(subset(all_angios, all_angios$species %in% pl_names))
leaves_div_raster<-mapDiversity(leaves_points, resolution = 1.2, plot = F, plot.with.grid = F)
#plot(leaves_div_raster,col=pal)
#plot(hybridshp, add=T)

# Hooks/Grapnels
hg<-subset(climbers_final, CM=="7")
hg_names<-paste(hg$taxon_name, hg$taxon_authors, sep=" ")
hooks_points<-as.data.frame(subset(all_angios, all_angios$species %in% hg_names))
hooks_div_raster<-mapDiversity(hooks_points, resolution = 1.2, plot = F, plot.with.grid = F)
#plot(hooks_div_raster,col=pal)
#plot(hybridshp, add=T)

##### Making the maps #####

## Proportion of climbers to all angiosperms
climb_neot_map<-climber_diversity_raster/neotropics_diversity_raster
plot(climb_neot_map,col=pal,zlim=c(0,.3))
plot(hybridshp, add=T)

## Twining/all climbers
#twining_points<-as.data.frame(twining_points)
tw_climb_map<-twining_div_raster/climber_diversity_raster
plot(tw_climb_map,col=pal,zlim=c(0,0.8),main="Proportion of twining in relation to all climbing mechanisms")
plot(hybridshp, add=T)

## Tendrils/all climbers
te_climb_map<-tendrils_div_raster/climber_diversity_raster
plot(te_climb_map,col=pal,zlim=c(0,0.7),main="Proportion of tendrils in relation to all climbing mechanisms")
plot(hybridshp, add=T)

## Scrambling/all climbers
sc_climb_map<-scrambling_div_raster/climber_diversity_raster
#sc_climb_map[which(sc_climb_map[]>=0.5)]
plot(sc_climb_map,col=pal,zlim=c(0,0.5),main="Proportion of scrambling in relation to all climbing mechanisms")
plot(hybridshp, add=T)

## Adhesive roots/all climbers
ar_climb_map<-roots_div_raster/climber_diversity_raster
#ar_climb_map[which(ar_climb_map[]>=0.5)]
plot(ar_climb_map,col=pal,zlim=c(0,0.5),main="Proportion of adhesive roots in relation to all climbing mechanisms")
plot(hybridshp, add=T)

## Prehensile branches/all climbers
pb_climb_map<-branches_div_raster/climber_diversity_raster
plot(pb_climb_map,col=pal,zlim=c(0,0.5),main="Proportion of prehensile branches in relation to all climbing mechanisms")
plot(hybridshp, add=T)

## Prehensile leaves/all climbers
pl_climb_map<-leaves_div_raster/climber_diversity_raster
#pl_climb_map[which(pl_climb_map[]>=.1)]
plot(pl_climb_map,col=pal,zlim=c(0,.08),main="Proportion of prehensile leaves in relation to all climbing mechanisms")
plot(hybridshp, add=T)

## Hooks/all climbers
hg_climb_map<-hooks_div_raster/climber_diversity_raster
hg_climb_map[which(hg_climb_map[]>=.08)]
plot(hg_climb_map,col=pal,zlim=c(0,.08),main="Proportion of hooks/grapnels in relation \nto all climbing mechanisms")
plot(hybridshp, add=T)

## Plotting all maps together
pdf(file="all_maps.pdf")
par(mfrow=c(4,2))
plot(tw_climb_map,col=pal,zlim=c(0,0.8),main="Proportion of twining in relation \nto all climbing mechanisms")
plot(hybridshp, add=T)
plot(te_climb_map,col=pal,zlim=c(0,0.7),main="Proportion of tendrils in relation \nto all climbing mechanisms")
plot(hybridshp, add=T)
plot(sc_climb_map,col=pal,zlim=c(0,0.5),main="Proportion of scrambling in relation \nto all climbing mechanisms")
plot(hybridshp, add=T)
plot(ar_climb_map,col=pal,zlim=c(0,0.5),main="Proportion of adhesive roots in relation \nto all climbing mechanisms")
plot(hybridshp, add=T)
plot(pb_climb_map,col=pal,zlim=c(0,0.5),main="Proportion of prehensile branches in relation \nto all climbing mechanisms")
plot(hybridshp, add=T)
plot(pl_climb_map,col=pal,zlim=c(0,.08),main="Proportion of prehensile leaves in relation \nto all climbing mechanisms")
plot(hybridshp, add=T)
plot(hg_climb_map,col=pal,zlim=c(0,.08),main="Proportion of hooks/grapnels in relation \nto all climbing mechanisms")
plot(hybridshp, add=T)

dev.off()

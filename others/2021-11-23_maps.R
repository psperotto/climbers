
# setwd("~/Desktop/Colabs/Patricia_Climbers/climbers")
# setwd("G:/Meu Drive/Papers/Diversifica??o/climbers")

# rm(list=ls())

library(raster)
library(maptools)
data("wrld_simpl")

############################
# Descriptive results
############################

# PD maps
data_dir <- paste0(getwd(), "/Data")

# (1) Loading main database of Neotropical climbing plants
all_climbers <- read.csv(paste0(data_dir, "/climber_database.csv"))

# (2) Loading maps of all Neotropical plants
all_neotropical <- readRDS(paste0(data_dir, "/full_neotropical_diversity.Rdata"))

# plot(all_neotropical)
# plot(wrld_simpl, add=T)

############################
# Let's make two types of maps.
############################
# First, let's make maps of species richness of climbing plants relative to total
# species richness in an area

############################
# Loading all mechanisms

all_Rdata <- list.files(data_dir)[grep(".Rdata", list.files(data_dir))]
all_mechs <- all_Rdata[grep("Mechanism", all_Rdata)]
mech_maps <- lapply(paste0(data_dir, "/", all_mechs), readRDS)

template.map <- all_neotropical
template.map[!is.na(template.map[])] <- 0
tmp.raster.list <- list()
  for (i in 1:length(mech_maps)) {
    r1 <- mech_maps[[i]]
    r1 <- raster::resample(r1, template.map)
    r1[is.na(r1)] <- 0
    tmp.raster.list[[i]] <- raster::mask(r1, template.map)
    print(i)
  }
  sprichness_map <- raster::calc(raster::stack(tmp.raster.list), sum)
  sprichness_map <- raster::crop(sprichness_map, raster::extent(all_neotropical))
  
  #saveRDS(tmp.raster.list, file=paste0("~/Desktop/MiSSEgradient/MiSSEGradient/regressions_plan/Data/1_rate_rasters/all_species_stack.Rdata"))
  
 raster1 = sprichness_map
 raster2 = all_neotropical
 

    raster1[raster1[]==0] <- NA
    raster2[raster2[]==0] <- NA
    l.model <- stats::lm(raster::getValues(raster1) ~ raster::getValues(raster2), na.action = na.exclude)
    res.raster <- template.map
    res.raster[] <- as.numeric(stats::residuals.lm(l.model))

    plot(res.raster)
    
  plot(sprichness_map)
  plot(all_neotropical)

# ()



  

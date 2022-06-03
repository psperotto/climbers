
##### Load functions ####
# Creates a list of all GBIF records
read.full.gbif <- function(gbif_files, gbif_dir) {
  allfiles <- list() 
  for(file_index in 1:length(gbif_files)) {
    tmp_file <- read.csv(paste0(gbif_dir, "z_filtered_gbif/", gbif_files[file_index]), stringsAsFactors = F)[,2:7]
    if(nrow(tmp_file)>0) {
    tmp_file$order <- strsplit(gbif_files[1], "_")[[1]][1]
    tmp_file$family <- strsplit(gbif_files[1], "_")[[1]][2]
    allfiles[[file_index]] <- tmp_file
    print(file_index)
    } else { next }
  }
  full_list <- data.table::rbindlist(allfiles)
  write.csv(full_list, file=paste0(gbif_dir, "full_tracheophyte_filtered_gbif.csv"))
return(full_list)  
}

save.gbif.neotropics <- function(full_list) {
  centroids <- read.csv(paste0(getwd(),"/centroids.neotropics.csv"), stringsAsFactors = F)[,2:3] # load list of centroids
  table_to_map <- as.data.frame(full_list[,c(1,4,3)])
  cleaned_points <- table_to_map %>% filter(V4 < -20, V4 > -130, V3 < 23, V3 > -23, !V4 %in% centroids$lon, !V3 %in% centroids$lat)   # restricting to the neotropics and removing centroids 
  cleaned_points <- cc_sea(cleaned_points, lon="V4", lat="V3")
  colnames(cleaned_points) <- c("species", "lon", "lat")
  write.csv(cleaned_points, file=paste0(getwd(), "/neotropics_tracheophyte_filtered_gbif.csv"))
}


# Creates a raster of species-richness grids for the Neotropics (also gets rid of annoying centroids)
run.mapDiversity.neotropics <- function(full_list, filename="full_neotropical_diversity", dir) {
 # centroids <- read.csv(paste0(getwd(),"/centroids.neotropics.csv"), stringsAsFactors = F)[,2:3] # load list of centroids
 # table_to_map <- as.data.frame(full_list[,c(1,4,3)])
 # cleaned_points <- table_to_map %>% filter(V4 < -20, V4 > -130, V3 < 23, V3 > -23, !V4 %in% centroids$lon, !V3 %in% centroids$lat)   # restricting to the neotropics and removing centroids 
 # cleaned_points <- cc_sea(cleaned_points, lon="V4", lat="V3")
  full_neotropics <- monographaR::mapDiversity(as.data.frame(full_list), export=T, filename=filename, plot.with.grid = F, plot=T)
  saveRDS(full_neotropics, file=paste0(dir, "/",filename, ".Rdata"))
  return(full_neotropics)
}

# Function to taxize species list of climbing plants according to GBIF taxonomy
gbif.taxize <- function (species) {
  sources <- taxize::gnr_datasources()
  gbif <- sources$id[sources$title == 'GBIF Backbone Taxonomy']
  gnr_resolve_x <- function(x, data_source_ids=gbif, best_match_only = TRUE) {
    new.name <- suppressWarnings(taxize::gnr_resolve(names=x, data_source_ids=data_source_ids, best_match_only=best_match_only)$matched_name)
    if(is.null(new.name)) {
      new.name <- paste0("UNMATCHED_",x)
    }
    return(new.name)
  }
  new.names <- pbapply::pbsapply(species, gnr_resolve_x, cl=48)
  return(new.names)
}

# (it should also adjust for phylogenetic dependency using phyloweights)

mapDiversity.pw <- function (data, phy, resolution = 1) {
  wrld_simpl = NULL
  message("Assuming the columns are ordered as: species, longitude and latitude")
  geo <- data
  colnames(geo) <- c("Species", "x", "y")
  coordinates(geo) = ~x + y
  r0 <- raster(resolution = resolution)
  r0[] <- NA
  r0 <- crop(r0, extent(geo) + +(resolution * 2))
  cells <- data.frame(spp = data[, 1], cells = cellFromXY(r0, 
                                                          data[, 2:3]))
  cells <- unique(cells)
  ### phyloweight ###
  
  phylo.weights.fast <- function(phy) {
    tree <- phy
    Ones <- matrix(1, Ntip(tree), 1)
    C <- vcv.phylo(tree)
    phylo.weights <- t(Ones) %*% solve(C) / sum(t(Ones) %*% solve(C))
    return(phylo.weights)
  }
  tip.weight <- phylo.weights.fast(phy) * 100
  species <- colnames(tip.weight)
  tip.weight <- as.numeric(tip.weight)
  tip.weight <- cbind(species, tip.weight)
  
  cells.pw <- merge(cells, tip.weight, by="species")
  
  unique.cells <- unique(cells.pw$cells)
  sum.pw <- c()
  for(i in 1:length(unique.cells)) {
  tmp <- subset(cells.pw, cells.pw$cells %in% unique.cells[i])
  sum.pw[i] <- sum(as.numeric(tmp$tip.weight))
  }
  names(sum.pw) <- unique.cells
  ########## 

  r0[as.numeric(names(sum.pw))] <- as.numeric(sum.pw)
  dev.off()
  plot(r0)
  return(r0)
}

# Function to plot residuals from linear regression between rasters 
### Regression
plot.res <- function(trait_raster, full_raster, dir=getwd(), pal.name = "RdBu", output = "residuals_sprich_pw", file.name="res") {
  template_background <- full_raster
  template_background[!is.na(template_background)] <- 0
  # set pallete 
  pal <- hcl.colors(30, palette = pal.name, alpha = 1)
  # plot residuals for regressions between sp_rich and trait
    tmp <- crop(full_raster, extent(trait_raster))
    tmp[tmp[]==0] <- NA
    l.model <- lm(getValues(trait_raster) ~ getValues(tmp), na.action = na.exclude)
    template_trait <- trait_raster
    template_trait[!is.na(template_trait)] <- 0
    template_trait[] <- as.numeric(residuals.lm(l.model)) 
    # first you have to make the 0s NAs, otherwise you will have a lot of 0s
    teste_pos <- teste_neg <- template_trait
    max.v <- max(template_trait[!is.na(getValues(template_trait))])
    teste_pos[teste_pos < 0.001] <- NA
    teste_neg[teste_neg > -0.001] <- NA
    template_trait <- do.call(merge, list(teste_pos, teste_neg))
    plot(template_background, col="white", legend=FALSE)
    plot(template_trait, col=rev(pal[c(1:13,17:30)]), add=T, zlim=c(max.v*-1,max.v))
    data("wrld_simpl")
    plot(wrld_simpl, add=T)
    writeRaster(template_trait, file = paste0(getwd(),"/", file.name,".tif"), overwrite=TRUE)
  return(template_trait)
}  

##
phy.list <- function (input.dir=input.dir, names=names, search.for=".taxized.tre") {
  l0 <- list.files(input.dir, search.for)
  l1 <- sub(search.for, "", l0)
  tree_names <- l0[which(l1 %in% names)]
  # Read tree files and .txt
  tmp.files1 <- lapply(paste0(input.dir,"/",tree_names), readLines)
  # Removing single quotes from authority names
  tmp.files2 <- lapply(tmp.files1, function(x) gsub("\\'","",as.character(x)))
  tree_list <- list()
  for(tree_index in sequence(length(tmp.files2))){
    # Read tree files and phylo 
    t0 <- ape::read.tree(text=tmp.files2[[tree_index]])
    t0 <- ape::ladderize(t0)
    # Removing node labels
    t0$node.label <- NULL
    # Removing duplicated tips and tips identified only to genus level
    name <- strsplit(t0$tip.label,"_",fixed = TRUE) 
    t0 <- ape::drop.tip(t0, which(unlist(lapply(1:length(name), function(x) stringr::str_detect(substr(name[[x]][2], 1, 1), "^[:upper:]+$")))))
    t0 <- ape::drop.tip(t0, which(duplicated(t0$tip.label)))
    tree_list[[tree_index]] <- t0
    names(tree_list)[tree_index] <- sub(search.for, "", tree_names)[tree_index]
  }
  return(tree_list)
}  


#' Thinning distribution data to smooth sampling bias
#' @param points A data.frame of three columns containing species and coordinates
#' @param species A character string indicating name of column with species names
#' @param lat A character string indicating name of column with latitudes
#' @param lon character string indicating name of column with longitudes
#' @param n A number indicating how many points to keep in each cell after thinning
Thinning <- function(points, species="species", lat = "decimalLatitude", lon="decimalLongitude", n = 1) {
  tmp_points = points
  colnames(tmp_points)[colnames(tmp_points)==lon] <- "x"
  colnames(tmp_points)[colnames(tmp_points)==lat] <- "y"
  spp <- unique(tmp_points[,species])
  results <- list()
  for(species_index in 1:length(spp)) {
    coords <- tmp_points[tmp_points[,species]==spp[species_index],c("y","x")]
    coords <- coords[!duplicated(coords[,"x"]) & !duplicated(coords[,"y"]),]
    if(nrow(coords) > 1) {
      sp::coordinates(coords) <- ~ y + x
      raster::crs(coords) <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
      r0 <- raster::raster(coords)
      raster::res(r0) <- 1 # cell resolution
      r0 <- raster::extend(r0, raster::extent(r0) + 5) 
      res <- cbind(spp[species_index], as.data.frame(dismo::gridSample(coords, r0, n))) # n = maximum number of points per cell
      colnames(res) <- c("species", "lat","lon")
      results[[species_index]] <- res
    } else {
      res <- cbind(spp[species_index],coords)
      colnames(res) <- c("species", "lat","lon")
      results[[species_index]] <- res
    }
  }
  results <- do.call(rbind, results)
  return(results)
}

#' Get values of selected climatic variables from filtered points
#' @param points A data.frame of three columns containing species and coordinates
#' @param species A character string indicating name of column with species names
#' @param lat A character string indicating name of column with latitudes
#' @param lon character string indicating name of column with longitudes
#' @param layerdir name of directory with climate layers
ClimateFromPoint_custom <- function(points, species="species",lon="lon", lat="lat", layerdir = ""){
  tmp_points = points
  colnames(tmp_points)[which(colnames(tmp_points) == lon)] <- "lon"
  colnames(tmp_points)[which(colnames(tmp_points) == lat)] <- "lat"
  colnames(tmp_points)[which(colnames(tmp_points) == species)] <- "species"
  tmp_points <- tmp_points[,c("species","lat","lon")]
  # Load climatic layers
  temp <- raster::raster(paste0(layerdir, "/current_30sec/bio_1.tif"))
  prec <- raster::raster(paste0(layerdir, "/current_30sec/bio_12.tif"))
  pet <- raster::raster(paste0(layerdir, "/et0_yr/et0_yr.tif"))
  aridity <- raster::raster(paste0(layerdir, "/ai_et0/ai_et0.tif"))
  bio <- list(temp, prec, pet, aridity)
  vars <- c(temp, prec, pet, aridity)
  names(vars) <- c("temp", "prec", "pet", "aridity")
  final_matrix <- matrix(nrow=nrow(tmp_points), ncol=length(vars))
  cat("Extracting climatic information of", nrow(tmp_points), "points",  "\n")
  sp::coordinates(tmp_points) <- ~ lon + lat
  for(var_index in 1:length(vars)) {
    layer <- bio[[var_index]]
    cat("\r",names(vars)[var_index])
    cat("","\n")
    values <- raster::extract(layer, tmp_points)
    final_matrix[,var_index] <- values
  }
  colnames(final_matrix) <- names(vars)
  result <- cbind(tmp_points, final_matrix)
  return(as.data.frame(result))
}

#' Get summary statistics for selected climatic variables
#' @param points A data.frame of three columns containing species and coordinates
#' @param species A character string indicating name of column with species names
#' @param lat A character string indicating name of column with latitudes
#' @param lon character string indicating name of column with longitudes
#' @param n A number indicating how many points to keep in each cell after thinning
GetClimateSummStats_custom <- function (points, type=c("raw","transformed")) {
  tmp_points <- points[,-which(colnames(points) %in% c("lon","lat"))]
  vars <- c("temp","prec", "pet", "aridity")
  allclimatevars <- list()
  spp <- unique(tmp_points$species)
  for(var_index in 1:length(vars)) {
    cat("\r",vars[var_index])
    cat("","\n")
    n_i <- c()
    sigma2_wi <- c()
    summ_stats <- matrix(nrow=length(spp), ncol=5)
    for(species_index in 1:length(spp)){
      sp1 <- tmp_points[tmp_points$species==spp[species_index],]
      cat("\r","Now doing species", species_index)
      cat("","\n")
      values <- sp1[,vars[var_index]]
      values <- values[!is.na(values)]
      if(type=="raw") {
        if(vars[var_index] %in% c("temp")) {
          values <-  (values / 10) 
        }
        n_i[species_index] <- length(values) # sample size
        sigma2_wi[species_index] <- ifelse(is.na(var(values)),0,var(values))  # sample variance
        
      }
      if(type=="transformed") {
        if(vars[var_index] %in% c("temp")) {
          values <-  (values / 10) + 273.15 # transforms to Kelvin
        }
        values <- log(values) # log
        n_i[species_index] <- length(values) # sample size
        sigma2_wi[species_index] <- ifelse(is.na(var(values)),0,var(values))  # sample variance
        
        if(any(values== -Inf)){
          values <- values[-which(values== -Inf)]
        }
      }
      n0 <- length(values)
      mean0 <- round(mean(values), 6)
      sd0 <- round(stats::sd(values), 6)
      se0 <- round(sd0/ sqrt(n0), 6)
      tmp_summ_stats <- c(n0, mean0, sd0, se0)
      summ_stats[species_index,] <- c(spp[species_index], tmp_summ_stats)
      colnames(summ_stats) <- c("species",paste0("n_",vars[var_index]), paste0("mean_",vars[var_index]),
                                paste0("sd_",vars[var_index]), paste0("se_",vars[var_index]))
    }
    sigma2_w <- sum(sigma2_wi*(n_i - 1)) / sum(n_i - 1)
    within_sp_var <-  round(sigma2_w/n_i, 6)
    summ_stats <- cbind(summ_stats, within_sp_var)
    colnames(summ_stats)[6] <- paste0("within_sp_var_",vars[var_index])
    allclimatevars[[var_index]] <- summ_stats
  }
  return(allclimatevars)
}

library(ape)
library(raster)
library(maptools)
library(picante)
library(raster)
data(wrld_simpl)

# Calculating PD using function picante::pd
#' @param list_of_ranges A list of rasters with the approx. distribution of species
#' @param tree A phylogenetic tree including the same species as in list_of_ranges; tip.labels has to be the same as names(list_of_ranges)
#' @param cut A vector of longitudes and latitudes on where to cut the maps to speed up PD calculation
#' @param include.root Whether to include the root of the tree in the calculations of PD
PDranges <- function(list_of_ranges, tree, cut=NULL, include.root=TRUE, template.map=NULL) { 
  ranges = list_of_ranges
  if(is.null(template.map)){
    template.map <- readRDS("data/template.map.Rdata")
  } 
  if(!is.null(cut)) {
    template.map <- crop(template.map, extent(cut[1],cut[2],cut[3],cut[4])) 
  }
  rasters1 <- list() 
  matrixPD1 <- matrix(0, ncol = 1, nrow = ncell(template.map)) 
  matrixPD1 <- data.frame(matrixPD1)
  for (ranges_index in 1:length(ranges)) { # very ineficient loop to get matrix for PD
    Sys.time() -> start_time
    r1 <- ranges[[ranges_index]]
    r1[is.na(r1)] <- 0
    r1 <- raster::resample(r1, template.map) 
    rasters1[[ranges_index]] <- raster::mask(r1, template.map)
    temp <- as.data.frame(values(rasters1[[ranges_index]]))
    colnames(temp) <- names(ranges)[ranges_index]
    matrixPD1 <- cbind(matrixPD1, temp)
    print(c(ranges_index, "in", length(ranges)))
    Sys.time() -> end_time  
    print(end_time-start_time)
  }
  names(rasters1) <- names(ranges)
  matrixPD1[is.na(matrixPD1)] <- 0 #
  matrixPD1 <- matrixPD1[,2:ncol(matrixPD1)]
  pd_final1 <- picante::pd(matrixPD1, tree, include.root=include.root) 
  pdMAP <- template.map
  values(pdMAP) <- pd_final1[,1]
  pdMAP <- raster::mask(pdMAP, template.map)
  return(pdMAP)
}

data <- as.data.frame(full_list[,-4])[1:10000,]

#----------------
# Get areas with low collection numbers
mapCollecting <- function (data, resolution = 1) 
{
  geo <- data
  colnames(geo) <- c("Species", "x", "y")
  coordinates(geo) = ~x + y
  r0 <- raster(resolution = resolution)
  r0[] <- NA
  r0 <- crop(r0, extent(geo) + +(resolution * 2))
  cells <- data.frame(spp = data[, 1], cells = cellFromXY(r0, 
                                                          data[, 2:3]))
  #cells <- unique(cells)
  t.cells <- table(cells$cells)
  r0[as.numeric(names(t.cells))] <- as.numeric(t.cells)
  return(r0)
}



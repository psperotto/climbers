# Utilities:
library(raster)
library(taxize)
library(data.table)

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

get.node.age <- function (phy) {
  root.node <- length(phy$tip.label)+1
  seq.nodes <- phy$edge
  dists <- phy$edge.length
  res <- numeric(max(phy$edge))
  for (i in seq_len(nrow(seq.nodes))) {
    res[seq.nodes[i, 2]] <- res[seq.nodes[i,1]] + dists[i]
  }
  ages <- abs(round(res,3)-round(max(res),3))
  return(ages)
} # fun?ao pra pegar os node ages


organize.bubble.plot <- function(trait_table, reference_table, all_vars, twgd_data) {
  tmp_reference_table <- subset(reference_table, reference_table$gbif_name %in% unique(trait_table$taxized_names))
  wcvp_subset <- subset(all_vars, all_vars$taxon_name %in% tmp_reference_table$wcvp_name)
  wcvp_subset <- subset(wcvp_subset, wcvp_subset$introduced==0)
  wcvp_subset <- subset(wcvp_subset, wcvp_subset$extinct==0)
  wcvp_subset <- subset(wcvp_subset, wcvp_subset$location_doubtful==0)
  
  focal_areas <- unique(wcvp_subset$area_code_l3)
  results <- matrix(nrow=0, ncol=5)
  for(i in 1:length(focal_areas)) {
    one_area <- focal_areas[i]
    one_subset <- subset(wcvp_subset, wcvp_subset$area_code_l3==one_area)
    sp_rich <- length(unique(one_subset$taxon_name))
    family_rich <- length(unique(one_subset$family))
    area_plus_buffer <- twgd_data[which(as.character(twgd_data$LEVEL3_COD) %in% one_area),]
    centroids <- rgeos::gCentroid(area_plus_buffer, byid=TRUE)
    lon <- extent(centroids)[1]
    lat <- extent(centroids)[3]
    results <- rbind(results, cbind(sp_rich, family_rich, one_area, lon, lat))
  }
  results <- as.data.frame(results)
  results$sp_rich <- as.numeric(results$sp_rich)
  results$family_rich <- as.numeric(results$family_rich)
  results$lon <- as.numeric(results$lon)
  results$lat <- as.numeric(results$lat)
  return(results)
}


#' Thinning distribution data to smooth sampling bias
#' @param points A data.frame of three columns containing species and coordinates
#' @param species A character string indicating name of column with species names
#' @param lat A character string indicating name of column with latitudes
#' @param lon character string indicating name of column with longitudes
#' @param n A number indicating how many points to keep in each cell after thinning
Thinning <- function(points, species="species", lat = "decimalLatitude", lon="decimalLongitude", n = 3) {
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



# setwd("~/Desktop/Colabs/Patricia_Climbers/climbers")
# setwd("G:/Meu Drive/Papers/Diversifica??o/climbers")

# rm(list=ls())
## Script and functions to plot Magallon and Sanderson comparison
get.tail.probs <- function (table=table) {
  if(!is.data.frame(table)){
    stop("Input table is not a data.frame.") }
  
  #if(!c("taxa", "diversity", "node", "age_mean") %in% colnames(table)){
  #  stop("Bad formatted input file. Please make sure table has the same column names as the template table. See help page for more info.") }
  # organize table
  # table <- organize.table(table, clade, diversity, node, age_mean, age_upper, age_lower) # find a better way to pass arguments to nested functions
  background = which(table$taxa =="bg_clade")
  
  if (length(background)==0) {
    stop("No clade was defined as background. Make sure one of the clades is named 'bg_clade' in the 'clade' column of the input table.")
  }
  if (length(background)>1) {
    stop("More than one clade was defined as background. Make sure just one of the clades is named 'bg_clade' in the 'clade' column of the input table.")
  }
  if(which.max(table$diversity)!=background & which.max(table$age_mean)!=background) {
    stop("Background clade is not the oldest and/or the most diverse. Make sure to ajust that in the input table.")
  }
  
  background.clade <- table[background,]
  
  results <- list()
  # define vars
  eps.yule = 0
  eps.high.ext = 0.9
  k = background.clade$diversity
  time = background.clade$age_mean
  
  stem.results <- get.tail.prob.stem(eps.yule, eps.high.ext, k, time)
  crown.results <- get.tail.prob.crown(eps.yule, eps.high.ext, k, time)
  
  results[[1]] <- stem.results
  results[[2]] <- crown.results
  names(results) <- c("stem.results", "crown.results")
  
  return(results)
}
##
plottailprobs <- function (table = example, tail.probs = results) {
  
  # organize table
  #table <- organize.table(table, clade, diversity, node, age_mean, age_upper, age_lower) # find a better way to pass arguments to nested functions
  table=table
  tail.probs=tail.probs
  
  background = which(table$taxa =="bg_clade")
  
  if (length(background)==0) {
    stop("No clade was defined as background. Make sure one of the clades is named 'bg_clade' in the 'clade' column of master table.")
  }
  if (length(background)>1) {
    stop("More than one clade was defined as background. Make sure just one of the clades is named 'bg_clade' in the 'clade' column of master table.")
  }
  if(which.max(table$diversity)!=background & which.max(table$age_mean)!=background) {
    stop("Background clade is not the oldest and/or the most diverse. Make sure to ajust that in the master table.")
  }
  
  background.clade <- table[background,]
  
  results <- list()
  # define vars
  k = background.clade$diversity
  time = background.clade$age_mean
  
  par(mfrow=c(1,2))
  ###### stem plot
  stem.results <- tail.probs[[1]]
  k.lower.yule.stem <- stem.results[[1]]
  k.upper.yule.stem <- stem.results[[2]]
  k.lower.high.ext.stem <- stem.results[[3]]
  k.upper.high.ext.stem <- stem.results[[4]]
  
  plot(k.lower.yule.stem$time, log10(k.lower.yule.stem$k), type="l",
       col="red", lwd=1.5, ylim = c(0, round(log10(k), 0)+1), xlim = c(0,time+time*0.1), xlab="", ylab="", yaxt = "n")
  par(new=TRUE)
  plot(k.upper.yule.stem$time, log10(k.upper.yule.stem$k), type="l",
       col="red", lwd=1.5, ylim = c(0, round(log10(k), 0)+1), xlim = c(0,time+time*0.1), xlab="", ylab="", yaxt = "n")
  par(new=TRUE)
  plot(k.lower.high.ext.stem$time, log10(k.lower.high.ext.stem$k), type="l",lty=2,
       col="red",lwd=1.5, ylim = c(0, round(log10(k), 0)+1), xlim = c(0,time+time*0.1), xlab="", ylab="",  yaxt = "n")
  par(new=TRUE)
  plot(k.upper.high.ext.stem$time, log10(k.upper.high.ext.stem$k), type="l",lty=2,
       col="red",lwd=1.5, ylim = c(0, round(log10(k), 0)+1), xlim = c(0,time+time*0.1), xlab='Age of clade (million years)', ylab='K (expected number of species)',yaxt = "n")
  
  points(x=table$age_mean[table$taxa=="bg_clade"], y=log10(table$diversity[table$taxa=="bg_clade"]), pch=1, col="orange", type="p")
  points(x=table$age_mean[table$taxa!="bg_clade" & table$node=="SG"], y=log10(table$diversity[table$taxa!="bg_clade" & table$node=="SG"]), pch=16, type="p")
  
  axis.par <- c(0:round(log10(k), 0)+1)
  axis(2, at=axis.par, labels=axis.par, las=1)
  legend(1, tail(axis.par, 1), legend=c("eps = 0", "eps = 0.9"), col=c("red", "red"), lty=1:2, cex=1)
  
  title(main="stem node rates")
  
  
  ###### crown plot
  crown.results <- tail.probs[[2]]
  k.lower.yule.crown <- crown.results[[1]]
  k.upper.yule.crown <- crown.results[[2]]
  k.lower.high.ext.crown <- crown.results[[3]]
  k.upper.high.ext.crown <- crown.results[[4]]
  
  plot(k.lower.yule.crown$time, log10(k.lower.yule.crown$k), type="l",
       col="blue", lwd=1.5, ylim = c(0, round(log10(k), 0)+1), xlim = c(0,time+time*0.1), xlab="", ylab="", yaxt = "n")
  par(new=TRUE)
  plot(k.upper.yule.crown$time, log10(k.upper.yule.crown$k), type="l",
       col="blue", lwd=1.5, ylim = c(0, round(log10(k), 0)+1), xlim = c(0,time+time*0.1), xlab="", ylab="", yaxt = "n")
  par(new=TRUE)
  plot(k.lower.high.ext.crown$time, log10(k.lower.high.ext.crown$k), type="l",lty=3,
       col="blue",lwd=1.5, ylim = c(0, round(log10(k), 0)+1), xlim = c(0,time+time*0.1), xlab="", ylab="",  yaxt = "n")
  par(new=TRUE)
  plot(k.upper.high.ext.crown$time, log10(k.upper.high.ext.crown$k), type="l",lty=3,
       col="blue",lwd=1.5, ylim = c(0, round(log10(k), 0)+1), xlim = c(0,time+time*0.1), xlab='Age of clade (million years)', ylab='K (expected number of species)',yaxt = "n")
  par(new=TRUE)
  
  #add points for clades
  points(x=table$age_mean[table$taxa=="bg_clade"], y=log10(table$diversity[table$taxa=="bg_clade"]), pch=1, col="orange", type="p")
  points(x=table$age_mean[table$taxa!="bg_clade" & table$node=="CG"], y=log10(table$diversity[table$taxa!="bg_clade" & table$node=="CG"]), pch=16, type="p")
  
  axis.par <- c(0:round(log10(k), 0)+1)
  axis(2, at=axis.par, labels=axis.par, las=1)
  legend(1, tail(axis.par, 1), legend=c("eps = 0", "eps = 0.9"), col=c("blue", "blue"), lty=c(1,3), cex=1)
  title(main="crown node rates")
  
}
##
get.age.nodes <- function (phy) {
  phy <- ape::ladderize(phy)
  root.node <- length(phy$tip.label)+1
  seq.nodes <- phy$edge
  dists <- phy$edge.length
  res <- numeric(max(phy$edge))
  for (i in seq_len(nrow(seq.nodes))) {
    res[seq.nodes[i, 2]] <- res[seq.nodes[i,1]] + dists[i]
  }
  return(res)
}
##
get.bg.rate.tree <- function (phy, sf=NA, eps=NA) {
  if (!inherits(phy, "phylo")) {
    stop("phy should be an object of class \"phylo\".")
  }
  if(is.na(eps)){
    eps=0
    cat("No eps specified, assuming no extinction...\n\n")
  }
  if(is.na(sf)){
    sf=1
    cat("No sf specified, assuming complete sampling...\n\n")
  }
  
  phy=phy
  sf=sf
  eps=eps
  
  diversity <- length(phy$tip.label)/sf
  time <- get.age.nodes(phy)[1]
  bg.rate <- eq6_and_7(eps, diversity, time)
  
  return(bg.rate)
}

####
get.template <- function(output.dir = getwd(), save.file=TRUE) {
  
  template <- data.frame(taxa=c("bg_clade","Calycanthales","Asterales","Myrtales","Arecaceae_c","Arecaceae_s"),
                         diversity=c(262196, 10, 25996, 10782, 2780, 2780),
                         node=c("bg","SG","CG","CG","CG","SG"),
                         age_mean=c(132, 108.8, 28.75, 88.2, 77.4, 84),
                         stringsAsFactors = FALSE)
  template$age_upper <- template$age_mean + 5
  template$age_lower <- template$age_mean - 5
  
  # if(save.file) {
  #    write.csv(template, file=paste0(output.dir,"/template.csv"))
  #  }
  
  return(template)
}

# Internal function
get.tail.prob.crown <- function (eps.yule, eps.high.ext, k, time) {
  eps.yule=eps.yule
  eps.high.ext=eps.high.ext
  k=k
  time=time
  
  results <- list()
  
  bg.rate.crown.yule <- eq6_and_7(eps.yule, k, time, crown=TRUE) # eq 7
  bg.rate.crown.high.ext <- eq6_and_7(eps.high.ext, k, time, crown=TRUE) #
  
  # define base vars
  
  birth.rate.crown.yule = bg.rate.crown.yule / (1 - eps.yule)
  birth.rate.crown.high.ext = bg.rate.crown.high.ext / (1 - eps.high.ext)
  
  death.rate.crown.yule = eps.yule * birth.rate.crown.yule # always 0
  death.rate.crown.high.ext = eps.high.ext * birth.rate.crown.high.ext
  
  # get CI for expected diversity based on background rates from the stem node (eq 10)
  prob <- 0
  time.bins <- seq(1, time, by=1)
  k.bins <- seq(1, k, by=1)
  
  # yule
  k.upper.yule <- matrix (nrow = 0, ncol = 2)
  colnames(k.upper.yule) <- c("k", "time")
  k.lower.yule <- matrix (nrow = 0, ncol = 2)
  colnames(k.lower.yule) <- c("k", "time")
  
  cat("Estimating confidence intervals of expected species diversity according to age of crown group (yule): time bins 1 to", time,  "\n")
  for (time_index in time.bins) {
    cat("\r", time_index)
    for (k_index in k.bins-1) {
      prob <- CrownNtKmore(k_index, birth.rate.crown.yule, death.rate.crown.yule, time_index) # eq 10
      # prob <- CrownNtKmore(k_index, bg.rate.crown.yule, eps.yule, time_index) # eq 10
      if (0.025 == round(prob,3)) {
        k.lower.yule <- rbind(k.lower.yule, c(k_index, time_index)) }
      else if (0.975 == round(prob,3)) {
        k.upper.yule <- rbind(k.upper.yule, c(k_index, time_index)) }
    }
  }
  cat("","\n")
  
  # high extinction
  k.upper.high.ext <- matrix (nrow = 0, ncol = 2)
  colnames(k.upper.high.ext) <- c("k", "time")
  k.lower.high.ext <- matrix (nrow = 0, ncol = 2)
  colnames(k.lower.high.ext) <- c("k", "time")
  
  cat("Estimating confidence intervals of expected species diversity according to age of crown group (high ext): time bins 1 to", time, "\n")
  for (time_index in time.bins) {
    cat("\r", time_index)
    for (k_index in k.bins-1) {
      prob <- CrownNtKmore(k_index, birth.rate.crown.high.ext, death.rate.crown.high.ext, time_index) # eq 10
      # prob <- CrownNtKmore(k_index, bg.rate.crown.high.ext, eps.high.ext, time_index) # eq 10
      if (0.025 == round(prob,3)) {
        k.lower.high.ext <- rbind(k.lower.high.ext, c(k_index, time_index)) }
      else if (0.975 == round(prob,3)) {
        k.upper.high.ext <- rbind(k.upper.high.ext, c(k_index, time_index)) }
    }
  }
  cat("","\n")
  
  results[[1]] <- as.data.frame(k.lower.yule)
  results[[2]] <- as.data.frame(k.upper.yule)
  results[[3]] <- as.data.frame(k.lower.high.ext)
  results[[4]] <- as.data.frame(k.upper.high.ext)
  names(results) <- c("k.lower.yule", "k.upper.yule", "k.lower.high.ext", "k.upper.high.ext")
  
  return(results)
}

# Internal function
get.tail.prob.stem <- function (eps.yule, eps.high.ext, k, time) {
  eps.yule=eps.yule
  eps.high.ext=eps.high.ext
  k=k
  time=time
  results <- list()
  
  bg.rate.stem.yule <- eq6_and_7(eps.yule, k, time, crown=FALSE) # eq 6
  bg.rate.stem.high.ext <- eq6_and_7(eps.high.ext, k, time, crown=FALSE) # eq 6
  
  # define base vars
  death.rate.stem.yule = eps.yule * bg.rate.stem.yule/ (1 - eps.yule) # always 0
  death.rate.stem.high.ext = eps.high.ext * bg.rate.stem.high.ext / (1 - eps.high.ext)
  
  birth.rate.stem.yule = bg.rate.stem.yule / (1 - eps.yule)
  birth.rate.stem.high.ext  = bg.rate.stem.high.ext / (1 - eps.high.ext)
  
  # get CI for expected diversity based on background rates from the stem node (eq 10)
  prob <- 0
  time.bins <- seq(1, time, by=1)
  k.bins <- seq(1, k, by=1)
  
  # yule
  k.upper.yule <- matrix (nrow = 0, ncol = 2)
  colnames(k.upper.yule) <- c("k", "time")
  k.lower.yule <- matrix (nrow = 0, ncol = 2)
  colnames(k.lower.yule) <- c("k", "time")
  
  cat("Estimating confidence intervals of expected species diversity according to the age of stem group (yule): time bins 1 to", time, "\n")
  for (time_index in time.bins) {
    cat("\r", time_index)
    for (k_index in k.bins-1) {
      prob <- StemNtKmore(k_index, birth.rate.stem.yule, death.rate.stem.yule, time_index) # eq 10
      # prob <- StemNtKmore(k_index, bg.rate.stem.yule, eps.yule, time_index) # eq 10
      if (0.025 == round(prob,3)) {
        k.lower.yule <- rbind(k.lower.yule, c(k_index, time_index)) }
      else if (0.975 == round(prob,3)) {
        k.upper.yule <- rbind(k.upper.yule, c(k_index, time_index)) }
    }
  }
  cat("","\n")
  
  # high extinction
  k.upper.high.ext <- matrix (nrow = 0, ncol = 2)
  colnames(k.upper.high.ext) <- c("k", "time")
  k.lower.high.ext <- matrix (nrow = 0, ncol = 2)
  colnames(k.lower.high.ext) <- c("k", "time")
  
  cat("Estimating confidence intervals of expected species diversity according to age of the stem group (high ext): time bins 1 to", time, "\n" )
  for (time_index in time.bins) {
    cat("\r", time_index)
    for (k_index in k.bins-1) {
      prob <- StemNtKmore(k_index, birth.rate.stem.high.ext, death.rate.stem.high.ext, time_index) # eq 10
      # prob <- StemNtKmore(k_index, bg.rate.stem.high.ext, eps.high.ext, time_index) # eq 10
      if (0.025 == round(prob,3)) {
        k.lower.high.ext <- rbind(k.lower.high.ext, c(k_index, time_index)) }
      else if (0.975 == round(prob,3)) {
        k.upper.high.ext <- rbind(k.upper.high.ext, c(k_index, time_index)) }
    }
  }
  cat("","\n")
  
  results[[1]] <- as.data.frame(k.lower.yule)
  results[[2]] <- as.data.frame(k.upper.yule)
  results[[3]] <- as.data.frame(k.lower.high.ext)
  results[[4]] <- as.data.frame(k.upper.high.ext)
  names(results) <- c("k.lower.yule", "k.upper.yule", "k.lower.high.ext", "k.upper.high.ext")
  
  return(results)
}

#-------------------------
# Equations from paper:
# equations 6 and 7
eq6_and_7 <- function(eps, diversity, time, crown=TRUE) {
  if(crown == TRUE) {
    if(eps == 0) {
      net.div.mle <- (log(diversity) - log(2))/time # eq 4
    } else {
      net.div.mle <- 1/time*(log(diversity/2*(1-eps^2)+ 2*eps+1/2*(1-eps)*sqrt(diversity*(diversity*eps^2-8*eps+2*diversity*eps+diversity)))-log(2)) # eq 7
    }
  }else{
    if(eps == 0) {
      net.div.mle <- log(diversity)/time # eq 3
    } else {
      net.div.mle  <- 1/time*log(diversity*(1-eps)+eps) # eq 6
    }
  }
  return(net.div.mle)
}

#Eq. 1a.
ProbNt0 <- function(birth.rate, death.rate, time, crown=TRUE){
  if(crown == TRUE){
    a <- 2
  }else{
    a <- 1
  }
  p0 <- AlphaT(birth.rate=birth.rate, death.rate=death.rate, time=time) ^ a
  return(p0)
}

#Eq. 2a.
AlphaT <- function(birth.rate, death.rate, time){
  beta.t <- BetaT(birth.rate=birth.rate, death.rate=death.rate, time=time)
  alpha.t <- (death.rate/birth.rate) * beta.t
  return(alpha.t)
}

#Eq. 2b.
BetaT <- function(birth.rate, death.rate, time){
  net.diver.rate <- birth.rate - death.rate
  exprt <- exp(net.diver.rate * time)
  beta.t <- (exprt - 1) / (exprt - (death.rate/birth.rate))
  return(beta.t)
}

#Eq 10a
StemNtKmore <- function(k, birth.rate, death.rate, time){
  beta.t <- BetaT(birth.rate=birth.rate, death.rate=death.rate, time=time)
  probNgeNtax <- beta.t ^ (k-1)
  return(probNgeNtax)
}

#Eq 11a
CrownNtKmore <- function(k, birth.rate, death.rate, time) {
  alpha.t <- AlphaT(birth.rate=birth.rate, death.rate=death.rate, time=time)
  beta.t <- BetaT(birth.rate=birth.rate, death.rate=death.rate, time=time)
  #probNgeNtax <- ((beta.t^(k-2))/ (1+alpha.t))* ((k * (1 - alpha.t - beta.t + (alpha.t*beta.t)) + alpha.t + (2*beta.t) - 1))
  probNgeNtax <- (beta.t^(k-2))*(k*(1 - alpha.t - beta.t + alpha.t*beta.t) + alpha.t + 2*beta.t -1)/(1+alpha.t)
  #print(probNgeNtax)
  #probNgeNtax <- (beta.t^(k-2))*(k*(1 - alpha.t - beta.t + alpha.t*beta.t) + alpha.t + 2*beta.t -1)/(1 - alpha.t + 2*alpha.t*beta.t)
  #print(probNgeNtax)
  return(probNgeNtax)
}

# WWFload is taken from speciesgeocodeR; all credit goes to the original authors
WWFload <- function(x = NULL) {
  if (missing(x)) {
    x <- getwd()
  }
  download.file("http://assets.worldwildlife.org/publications/15/files/original/official_teow.zip",
                destfile = file.path(x, "wwf_ecoregions.zip"), quiet=TRUE)
  unzip(file.path(x, "wwf_ecoregions.zip"), exdir = file.path(x, "WWF_ecoregions"))
  file.remove(file.path(x, "wwf_ecoregions.zip"))
  wwf <- maptools::readShapeSpatial(file.path(x, "WWF_ecoregions", "official",
                                              "wwf_terr_ecos.shp"))
  return(wwf)
}


localityToBiome <- function (points,lat="lat",lon="lon") {
  #colnames(points) <- c("acceptedScientificName","key","decimalLatitude","decimalLongitude","basisOfRecord","issues")
  cat("Getting biome from locality data...")
  points[,lat] <-  as.numeric(points[,lat])
  points[,lon] <-  as.numeric(points[,lon])
  locations.spatial <- sp::SpatialPointsDataFrame(coords=points[,c(lon, lat)], data=points)
  wwf <- WWFload(tempdir())
  mappedregions <- sp::over(locations.spatial, wwf)
  biomes <- c("Tropical & Subtropical Moist Broadleaf Forests", "Tropical & Subtropical Dry Broadleaf Forests", "Tropical & Subtropical Coniferous Forests", "Temperate Broadleaf & Mixed Forests", "Temperate Conifer Forests", "Boreal Forests/Taiga", "Tropical & Subtropical Grasslands, Savannas & Shrubland", "Temperate Grasslands, Savannas & Shrublands", "Flooded Grasslands & Savannas", "Montane Grasslands & Shrublands", "Tundra", "Mediterranean Forests, Woodlands & Scrub", "Deserts & Xeric Shrublands", "Mangroves")
  points$eco_name <- mappedregions$ECO_NAME
  points$biome <- biomes[mappedregions$BIOME]
  return(points)
}


# getting biomes for each species
getBiomes <- function (points, species="species") {
  cat("Summarizing biome from locality data...")
  points <- as.data.frame(points) # not sure how to do it without transforming back to data.frame
  points <- subset(points, !is.na(points[,"biome"]))
  categories <- unique(points[,"biome"])
  taxa <- as.character(unique(points[,species]))
  result <- matrix(0, nrow=length(taxa), ncol=length(categories))
  rownames(result) <- taxa
  colnames(result) <- categories
  for (taxon_index in seq_along(taxa)) {
    for (category_index in seq_along(categories)) {
      x0 <- points[,species]==taxa[taxon_index]
      x1 <- points[,"biome"]==categories[category_index]
      result[taxon_index, category_index] <- length(which(x0 & x1))
    }
  }
  return(result)
}

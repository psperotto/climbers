#' Estimates range and range size based on species distribution modeling.
#' @param points A data.frame of distribution points with at least three columns where one column represents species names and the other two decimal coordinates
#' @param species The name of the column in the data.frame with the names of species
#' @param lat The name of the column in the data.frame with latitudes
#' @param lon The name of the column in the data.frame with longitudes
#' @param threshold A number between 0 and 1 for the threshold to make models binary
#' @param buffer The radius in kilometers around points to estimate range size when there are three or less valid points
#' @param res The resolution (2.5, 5 or 10) of climatic variables used for modeling
#' @return A list containing information about the species distribution modeling performance and a raster with the possible range of each species
#' @export
#'
#'
#'
GetRanges <- function(points, species="species", lat="decimalLatitude", lon="decimalLongitude", threshold=0.75, buffer=25, res=10) {
  tmp_points = as.data.frame(points)
  tmp_points = tmp_points[,c(which(colnames(tmp_points)==species),which(colnames(tmp_points)==lon),which(colnames(tmp_points)==lat))]
  colnames(tmp_points) <- c("species","lon","lat")
  spp <- unique(tmp_points[,1])
  list_of_ranges <- list()
  for(species_index in 1:length(spp)) {
    points_for_range <- tmp_points[tmp_points$species==spp[species_index],]
    tmp_list <- GetOneRange(points_for_range, threshold, buffer, res)
    list_of_ranges[[species_index]] <- tmp_list
    names(list_of_ranges)[species_index] <- spp[species_index]
    #cat("","\n")
    #cat(species[species_index], "done.")
    #cat("","\n")
    cat("\r", species_index, " out of ",length(spp))
  }
  try(unlink("wc2-5", recursive = TRUE))
  try(unlink("wc5", recursive = TRUE))
  try(unlink("wc10", recursive = TRUE))
  return(list_of_ranges)
}


#make.new.dir <- function(dir, namedir, overwrite.dir) {
#  new.dir <- paste0(getwd(),"/", namedir)
#  if (overwrite.dir) {
#    unlink(new.dir, recursive=TRUE) # overwrites output directory
#    dir.create(file.path(new.dir))
#  } else { suppressWarnings(dir.create(new.dir)) }
#  return(new.dir)
#}

#' @param points_for_range points_for_range
#' @param threshold threshold
#' @param buffer buffer
#' @param res res
GetOneRange <- function(points_for_range, threshold, buffer, res) {
  #cat("Thinning points...")
  if(nrow(points_for_range) > 2) {
    thinned_points <- Thinning_internal(points_for_range)
  } else { thinned_points = points_for_range}
  #cat("Loading environmental predictors...")
  predictors <- LoadWcLayers(res.layers=res)
  #cat("Creating background polygon...")
  bg <- BgPolygon(thinned_points)
  if(nrow(thinned_points) < 4) { # If three or fewer valid points, the range will be retrieved from a circle around these points
    list_of_model_results <- RangeFromFewPoints(thinned_points, predictors, buffer)
  } else {
    #cat("Removing predictors with colinearity problems ...")
    predictors_final <- ColinearityTest(bg, predictors)
    #cat("Performing species distribution modeling...")
    list_of_model_results <- RangeFromSDM_dismo(thinned_points, predictors_final, bg)
  }
  #cat("Making models binary based on threshold...")
  fullresults <- RangeSize(list_of_model_results, threshold)
  #cat("Adding alerts...")
  fullresults_w_alerts <- AddAlerts(fullresults, bg)
  return(fullresults_w_alerts)
}


#' Performs sdm using dismo's maxent
#' For great maxent and sdm tutorials see:
#' jcoliver.github.io/learn-r/011-species-distribution-models.html
#' cran.r-project.org/web/packages/dismo/vignettes/sdm.pdf
#' etc
#' @importFrom dismo randomPoints kfold maxent predict evaluate
#' @importFrom raster raster
#' @param thinned_points thinned_points
#' @param predictors_final predictors_final
#' @param bg bg
RangeFromSDM_dismo <- function (thinned_points, predictors_final, bg) {
  points_for_sdm <- thinned_points[,c("lon","lat")]
  species_name <- as.character(thinned_points[1,1])
  mask <- raster::raster(predictors_final[[1]])
  background <- dismo::randomPoints(mask = mask, n = abs((bg[1] - bg[2]) * (bg[3] - bg[4])), # Number of random background points
                                    ext = bg, extf = 1.25)
  background <- as.data.frame(background)
  testing_group <- 1 # arbitrarily assign group 1 as the testing data group
  group_presence <- dismo::kfold(points_for_sdm, k = 2) # randomly divide dataset into two groups
  presence_train <- points_for_sdm[group_presence != testing_group, ]
  presence_test <- points_for_sdm[group_presence == testing_group, ]

  group_background <- dismo::kfold(background, k = 2) # background points
  background_train <- background[group_background != testing_group, ]
  background_test <- background[group_background == testing_group, ]

  sdm_raw <- dismo::maxent(x = predictors_final, p = presence_train, nbg=1000)
  predict <- dismo::predict(object = sdm_raw, x = predictors_final, ext = bg)

  eval <- dismo::evaluate(p = presence_test, a = background_test, model = sdm_raw, x = predictors_final)
  auc0 <- methods::slot(eval, "auc")

  results <- list()
  results[[1]] <- paste0(species_name)
  results[[2]] <- predict
  results[[3]] <- paste("maxent (dismo)")
  results[[4]] <- names(predictors_final)
  results[[5]] <- round(auc0, 2)
  names(results) <- c("species_name","original_model", "sdm_method", "predictors","auc")

  return(results)
}

#' Detecting collinearity in predictors
#' @importFrom raster crop stack
#' @importFrom usdm exclude vifcor
#' @param bg bg
#' @param predictors predictors
ColinearityTest <- function(bg, predictors) {
  layers <- raster::crop(predictors, raster::extent(bg))
  v0 <- suppressWarnings(usdm::vifcor(layers, th=0.8)) # stablished threshold
  predictors_final <- usdm::exclude(layers, v0) # excludes variables that have collinearity problems
  predictors_final <- raster::stack(predictors_final)
  return(predictors_final)
}

#' Internal function to add alerts to the results
#' @param fullresults fullresults
#' @param bg bg
AddAlerts <- function(fullresults, bg) {
  if(length(fullresults)==5){
    # is the number of good points three or less? (data deficiency alert)
    fullresults[[6]] <- "Data deficiency alert: species with three or less valid points."
    names(fullresults)[6] <- "alerts"
    return(fullresults)
  } else {
    alerts <- c()
    #cross check with list of crops and invasive species
    #crops <- readRDS("R/crops_taxized_gbif.Rdata")
    #invasive <- readRDS("R/invasive_taxized_gbif.Rdata")
    #if(fullresults$species_name %in% crops) {
      #crop_alert <- "Crop alert: this species is listed as a crop at either fao.org, hort.purdue.edu/newcrop/ or Meyer et al. (2012)."
      # Meyer, R. S., DuVal, A. E. & Jensen, H. R. Patterns and processes in crop domestication: an historical review and quantitative analysis of 203 global food crops. New Phytol. 196, 29–48 (2012).
     # alerts <- c(alerts, crop_alert)
    #}
    #if(fullresults$species_name %in% invasive) {
    #  invasive_alert <- "Invasive alert: this species is listed as invasive at GISD. Check www.iucngisd.org for more information."
    #  alerts <- c(alerts, invasive_alert)
    #}
    # is AUC very low?
    if(fullresults$auc < 0.75) {
      low_auc_alert <- "Alert of low AUC: model with AUC lower than 0.75."
      alerts <- c(alerts, low_auc_alert)
    }
    # is range too wide? (proxy for non-natural distribution)
    lims <- bg[]
    if(lims[2] > lims[1] + 100 & lims[4] > lims[3] + 50) {
      wide_range_alert <- "Alert of wide range: range seems too wide for a natural distribution. Check quality of distribution data."
      alerts <- c(alerts, wide_range_alert)
    } else if(fullresults$range_size > 10000000) {
      wide_range_alert <- "Alert of wide range: range seems too wide for a natural distribution. Check quality of distribution data."
      alerts <- c(alerts, wide_range_alert)
    }
    if(is.null(alerts)) {
      fullresults[[9]] <- "No alerts returned."
      names(fullresults)[9] <- "alerts"
    } else {
      fullresults[[9]] <- alerts
      names(fullresults)[9] <- "alerts"
    }
    return(fullresults)
  }
}

#' Performs spatial thinning before modeling by selecting one occurence per grid cell
#' @importFrom dismo gridSample
#' @importFrom raster crs raster res extend extent
#' @importFrom sp coordinates
#' @param points_for_range points_for_range
Thinning_internal <- function(points_for_range) {
  coords <- points_for_range[,c(3:2)]
  sp::coordinates(coords) <- ~ lat + lon
  raster::crs(coords) <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
  coords@bbox[1:4] <- c(coords@bbox[1] - 1, coords@bbox[2] - 1, coords@bbox[3] + 1, coords@bbox[4] + 1) # that's to avoid errors when min and max lon and lat are the same
  r0 <- raster::raster(coords)
  raster::res(r0) <- 1 # cell resolution
  r0 <- raster::extend(r0, raster::extent(r0) + 5) # expand the extent of the RasterLayer a little
  thinned_points <- as.data.frame(dismo::gridSample(coords, r0, n = 1)) # n = maximum number of points per cell
  species <- points_for_range$species[sequence(nrow(thinned_points))]
  thinned_points <- cbind(species, thinned_points)
  return(thinned_points)
}

#' Getting a proxy for range size when you have three points or less.
#' When there are three or less valid points, the distribution will be estimated from a circle around each point of radius equal to the argument buffer.
#' @importFrom sp coordinates polygons
#' @importFrom raster raster crop mask
#' @importFrom dismo circles
#' @param thinned_points thinned_points
#' @param predictors predictors
#' @param buffer buffer
RangeFromFewPoints <- function(thinned_points, predictors, buffer) {
  point <- thinned_points
  species_name <- as.character(point[1,1])
  point[,1] <- 1
  sp::coordinates(point) <- ~ lon + lat
  area_mask <- raster::raster(res=res(predictors))
  area_mask[] <- 1
  buffer_mask = buffer # Radius of the circle in kilometers
  circle_around_point <- dismo::circles(point, d=buffer_mask*1000, lonlat=TRUE)
  circle_around_point <- sp::polygons(circle_around_point)
  area_mask <- crop(area_mask, circle_around_point)
  area_mask <- mask(area_mask, circle_around_point)
  results <- list()
  results[[1]] <- paste0(species_name)
  results[[2]] <- area_mask
  results[[3]] <- "no AUC"
  names(results) <- c("species_name","range","auc")
  return(results)
}

#' Estimates range size from binary model
#' @importFrom raster area
#' @param list_of_model_results list_of_model_results
#' @param threshold threshold
RangeSize <- function (list_of_model_results, threshold) {
  if(length(list_of_model_results) == 3) {
    range_size <- round(sum(raster::area(list_of_model_results$range, na.rm=T)[], na.rm=T), 2)
    list_of_model_results[[4]] <- range_size
    list_of_model_results[[5]] <- paste("Range estimated from three points or less.")
    names(list_of_model_results)[c(4,5)] <- c("range_size", "note")
    return(list_of_model_results)
  } else {
    model <- list_of_model_results$original_model
    model[model[] < threshold] <- NA
    model[!is.na(model)] <- 1
    range_size <- round(sum(raster::area(model, na.rm=T)[], na.rm=T), 2)
    list_of_model_results[[6]] <- model
    list_of_model_results[[7]] <- paste(threshold)
    list_of_model_results[[8]] <- range_size
    names(list_of_model_results)[c(6,7,8)] <- c("range", "threshold","range_size")
    return(list_of_model_results)
  }
}

#' Internal function -- gets polygon around distribution to crop predictors before modeling
#' @importFrom raster extent
#' @param thinned_points thinned_points
#' @param buffer.polygon buffer.polygon
BgPolygon <- function (thinned_points, buffer.polygon=c(5, 10)) {
  max.lat <- ceiling(max(thinned_points[,2])) + buffer.polygon[1]
  min.lat <- floor(min( thinned_points[,2])) - buffer.polygon[1]
  max.lon <- ceiling(max(thinned_points[,3])) + buffer.polygon[2]
  min.lon <- floor(min(thinned_points[,3])) - buffer.polygon[2]
  #if(max.lat > min.lat + 50 & max.lon > min.lon + 100){
  #  warning("Very wide background detected - this may not correspond to a natural distribution.")
  #}
  bg.area <- raster::extent(x = c(min.lon, max.lon, min.lat, max.lat))
  if(bg.area[1] > 180 | bg.area[2] > 180 |  bg.area[2] < -180 | bg.area[1] < -180){
    bg.area[which(bg.area[] < -180)] <- -179.9
    bg.area[which(bg.area[] > 180)] <- 179.9
  }
  #if(plot.polygons) {
  #  polygon.dir <- paste0(getwd(),"/polygon")
  #  if (overwrite.dir) {
  #    unlink(polygon.dir, recursive=TRUE) # overwrites output directory for models
  #    dir.create(file.path(polygon.dir))
  #  } else { suppressWarnings(dir.create(polygon.dir)) }
  #  pdf(file=paste0(polygon.dir,"/", thinned_points[1,1],"__background_polygon.pdf"), height=6, width=10)
  #  plot(wrld_simpl)
  #  plot(bg.area, col="orange", add=T)
  #  graphics::title(main=paste0(thinned_points[1,1]))
  #  dev.off()
  #}
  return(bg.area)
}

#' Internal function for now -- loads climatic layers from Worldclim
#' @importFrom raster getData
#' @param res.layers res.layers
LoadWcLayers <- function (res.layers) {
  bio <- raster::getData("worldclim", var="bio", res=res.layers, path=getwd()) # climatic layers
  return(bio)
}


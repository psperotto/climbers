
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
run.mapDiversity.neotropics <- function(full_list, filename="full_neotropical_diversity") {
 # centroids <- read.csv(paste0(getwd(),"/centroids.neotropics.csv"), stringsAsFactors = F)[,2:3] # load list of centroids
 # table_to_map <- as.data.frame(full_list[,c(1,4,3)])
 # cleaned_points <- table_to_map %>% filter(V4 < -20, V4 > -130, V3 < 23, V3 > -23, !V4 %in% centroids$lon, !V3 %in% centroids$lat)   # restricting to the neotropics and removing centroids 
 # cleaned_points <- cc_sea(cleaned_points, lon="V4", lat="V3")
  full_neotropics <- monographaR::mapDiversity(cleaned_points, export=T, filename=filename, plot.with.grid = F, plot=T)
  saveRDS(full_neotropics, file=paste0(filename, ".Rdata"))
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
  new.names <- pbapply::pbsapply(species, gnr_resolve_x)
  
  return(new.names)
}

# Function to plot residuals from linear regression between rasters 
# (it should also adjust for phylogenetic dependency using phyloweights)
run.mapDiversity.neotropics.pw <- function(full_list, filename="full_neotropical_diversity") {
  centroids <- read.csv(paste0(getwd(),"/centroids.neotropics.csv"), stringsAsFactors = F)[,2:3] # load list of centroids
  table_to_map <- as.data.frame(full_list[,c(1,4,3)])
  cleaned_points <- table_to_map %>% filter(V4 < -20, V4 > -130, V3 < 23, V3 > -23, !V4 %in% centroids$lon, !V3 %in% centroids$lat)   # restricting to the neotropics and removing centroids 
  cleaned_points <- cc_sea(cleaned_points, lon="V4", lat="V3")

  taxonomy_tree <- phy.list(input.dir=getwd(), names="GBMB", search.for=".taxized.tre")[[1]]
  taxonomy_tree$tip.label <- gsub("_", " ", taxonomy_tree$tip.label)
  
  tips_trait <- intersect(taxonomy_tree$tip.label, cleaned_points$V1)
  trait_tree <- keep.tip(taxonomy_tree, tips_trait)
  
  full_neotropics <- monographaR::mapDiversity(cleaned_points, export=T, filename=filename, plot.with.grid = F, plot=T)
  saveRDS(full_neotropics, file=paste0(filename, ".Rdata"))
  
  ### Phyloweight ###

  phylo.weight.rasters <- tmp.raster.list[species.in.common]
  
  phylo.weights.fast <- function(phy){
    tree <- phy
    Ones <- matrix(1, Ntip(tree), 1)
    C <- vcv.phylo(tree)
    phylo.weights <- t(Ones) %*% solve(C) / sum(t(Ones) %*% solve(C))
    return(phylo.weights)
  }
  tip.weight <- phylo.weights.fast(tree) * 1000
  n0 <- colnames(tip.weight)
  tip.weight <- as.numeric(tip.weight)
  names(tip.weight) <- n0
  
  mean.tip.rates$taxon <- sub("_"," ", mean.tip.rates$taxon)
  tmp.raster.rates <- tmp.raster.list[which(names(tmp.raster.list) %in% mean.tip.rates$taxon)]
  tip.weight <- tip.weight[names(tmp.raster.rates)]
  
  # sp richness pw
  sp.rich <- tmp.raster.rates
  for(species_index in 1:length(sp.rich)){
    tmp.tip.weight = as.numeric(tip.weight[names(sp.rich)[species_index]])
    values(sp.rich[[species_index]]) <- values(sp.rich[[species_index]]) * tmp.tip.weight
  }
  sp.rich <- as.list(calc(stack(sp.rich), sum))
  names(sp.rich) <- "species_richness_pw"
}  

  ### Regression
plot.res <- function(trait_raster, full_raster, dir=getwd(), pal.name = "RdBu", output = "residuals_sprich_pw") {
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
    plot(template_background, col=hcl.colors(5, palette = "Light Grays", alpha = 1)[2], legend=FALSE)
    plot(template_trait, col=rev(pal[c(1:13,17:30)]), add=T, zlim=c(max.v*-1,max.v))
    data("wrld_simpl")
    plot(wrld_simpl, add=T)
    writeRaster(template.raster, file = paste0(getwd(),"/", file.name,".tif"), overwrite=TRUE)
  return(residuals)
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



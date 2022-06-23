
library(ape)

load.trees <- function(tree.dir) {
  tree_files <- list.files(tree.dir, full.names = T)
  all_trees <- list()
  for(i in 1:length(tree_files)) {
    load(tree_files[i])
    if(exists("one_tree")) {
      all_trees[[i]] <- one_tree
      names(all_trees)[i] <- gsub(paste0(c(paste0(tree.dir,"/"), ".Rsave"), collapse="|"),"", tree_files[i])
      rm("one_tree")
    }
  }
  return(all_trees)
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

# most recent common ancestor of a bunch of tips
big_tree <- readRDS("trees/taxized_GBMB.Rdata")
big_tree$tip.label <- unname(big_tree$tip.label)

climber_clades <- readRDS("bg.clades.list.Rdata")
all_climbers_genera <- climber_clades[[1]]
all_climbers_genera <- all_climbers_genera$taxa
all_climbers_genera <- all_climbers_genera[-1]

tips_to_keep <- c()
for(i in 1:length(all_climbers_genera)){
  one_label <- all_climbers_genera[i]
  tips_to_keep_tmp <- grep(one_label, big_tree$tip.label)
  big_tree$tip.label[tips_to_keep_tmp] <- one_label
  tips_to_keep <- c(tips_to_keep, tips_to_keep_tmp)
}
pruned_big_tree <- keep.tip(big_tree, tips_to_keep)
pruned_big_tree$node.label <- NULL
pruned_big_tree <- keep.tip(pruned_big_tree, which(!duplicated(pruned_big_tree$tip.label)))
write.tree(pruned_big_tree, file="climber_genera_tree.tre")
#plot(ladderize(pruned_big_tree), cex=0.3)


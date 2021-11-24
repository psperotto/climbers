# setwd("~/Desktop/Colabs/Patricia_Climbers/climbers")
# setwd("~/Desktop/climbers")
setwd("C:/Users/patri/Google Drive/Papers/Diversifica??o/climbers")
getwd()

library(phytools)
library(geiger)
library(ape)
library(magrittr)
library(diversitree)
library(dismo)
library(hisse)
library(CoordinateCleaner)
library(monographaR)
library(dplyr)
library(devtools)
library(lcvplants)
library(beepr)
library(tidyverse)
library(taxize)


  ## creating a list of only the genera of neotropical climbers
climbers <- read.csv("climber_database.csv", stringsAsFactors = F)
x<-as.data.frame(climbers[,2])
genera2<-as.character(x[!duplicated(x$`climbers[, 2]`), ]) #this thing between '' is the name of the column

  ## creating a list with every species under every genera
spppergenera <- list()
for (i in 1:length(genera2))  {
  spppergenera[[i]]<-LCVP(genera2[[i]], genus_tab = T)
}

  ## glueing together all species
speciespergenera<-do.call(rbind,spppergenera)
      # do.call(): do.call constructs and executes a function 
      # call from a name or a function and 
      # a list of arguments to be passed to it.

  ## leaving only accepted SPECIES names
bb<-subset(speciespergenera, Status=="accepted")
acceptedspecies<-subset(bb, Infrasp=="species")

  ## counting spp per genera
sppclimbers<-as.data.frame(table(unlist(climbers$Genus)))
spptotal<-as.data.frame(table(unlist(acceptedspecies$Genus)))
#write.csv(sppclimbers, file="sppclimbers.csv")
#write.csv(spptotal, file="spptotal.csv")
  
  ## checking the acceptance of names in my climbing species list
a<-as.matrix(climbers$Species)
tipLabel<-strsplit(a,"_")
new.names<-matrix() 
for(i in 1:length(tipLabel)){
  new.names[i]<-paste(tipLabel[[i]][1:2], collapse=" ")
}
names<-as.vector(new.names)
    # these steps here were to create the vector the function LCVP needs

checknames<-LCVP(names)
beep("mario")

b<-as.data.frame(checknames$LCVP_Accepted_Taxon)
updatedclimbernames<-as.data.frame(b[!duplicated(b$`checknames$LCVP_Accepted_Taxon`), ])
names(updatedclimbernames)[1]<-"Species" # renaming the only column in the 
                                         # 'updatedclimbernames' dataframe
updatedclimbernames<-subset(updatedclimbernames, Species!="unresolved")
updatedclimbernames<-subset(updatedclimbernames, Species!="")
write.csv(updatedclimbernames, file = "updatedclimbernames.csv")

  ## recounting spp per genus in 'updatedclimbernames' 
updatedcn<-read.csv(file = "updatedclimbernames.csv", sep=";") # 'cn' stands for 'climber names'
spptotal_recount<-as.data.frame(table(unlist(updatedcn$Genus)))
saveRDS(tree, file="treegenera.RDS")

      ######## taxize ######
sources <- taxize::gnr_datasources()
sources$title

h<-as.data.frame(taxize::gnr_resolve(names="Justicia schomburgkiana", 
                    data_source_ids=sources$id[sources$title == "The International Plant Names Index"], 
                    best_match_only=FALSE))
            ## esse names=x ? um nome de esp?cie por vez, d? pra fazer um loop iterando pra uma lista de nomes
    # $matched_name

# ?gnr_resolve() #but this function doesn't seem to give me the accepted names back... it just seems if they exist at all, from what i understood.
h <- gnr_resolve(names = c("Justicia schomburgkiana"))
#=======

h<-gnr_resolve(names = c("Asteraceae", "Plantae"))


#-----------------------
# getting age of stem and crown node

# load tree
s_b_tree <- readRDS("treegenera.RDS")
example <- "Justicia"

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
}


library(phangorn)

tip_numbers <- grep(example, s_b_tree$tip.label)
crown_node <- mrca.phylo(s_b_tree, tip_numbers)
stem_node <- Ancestors(s_b_tree, crown_node, type = "parent")
crown_age <- get.node.age(s_b_tree)[crown_node]
stem_age <- get.node.age(s_b_tree)[stem_node]





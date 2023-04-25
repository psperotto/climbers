# Script for preparing the objects for the diversification and spatial analyses

## Reading data
#library(xlsx)
#APGIV <- read.xlsx("C:/Users/patri/Desktop/temp/APGIV.xlsx", sheetIndex = 1) # list of angiosperm families, orders and clades following APG IV
APGIV <- read.csv("datasets/APGIV.csv") # list of angiosperm families, orders and clades following APG IV

#wcvp_names <- read.table("C:/Users/patri/Desktop/temp/wcvp_names_and_distribution_special_edition_2022/wcvp_names.txt", sep="|", header=T, quote = "", fill=TRUE, encoding = "UTF-8")
#climbers_final <- read.csv(file = "climbers_final.csv")[,2:10]
climbers_final <- read.csv(file = "datasets/database_corrected.csv")
# alternatively: climbers_final <- readRDS(file="climbers_final.Rdata")
library(ape)
tree <- read.tree(file="trees/GBMB.tre")

## Counting species of climbers per family
spp.family<-as.data.frame(table(unlist(climbers_final$family)))
colnames(spp.family)<-c("Family", "Nr")
spp.family<-spp.family[order(spp.family$Nr, decreasing=T),]
rownames(spp.family)<-c(1:98)

## Counting species of climbers per genera
spp.genera<-as.data.frame(table(unlist(climbers_final$genus)))
colnames(spp.genera)<-c("Genus", "Nr")
spp.genera<-spp.genera[order(spp.genera$Genus),]

saveRDS(spp.family, file="spp.family.Rdata")
saveRDS(spp.genera, file="spp.genera.Rdata")

## Counting all species per genera of climbers 
spp.all <- subset(wcvp_names, taxon_status=="Accepted" & taxon_rank=="Species") %>% filter(genus %in% spp.genera$Genus)
spp.all <- as.data.frame(table(unlist(spp.all$genus)))
colnames(spp.all)<-c("Genus", "Nr") 
spp.all <- spp.all[order(spp.all$Genus),]
length(spp.all$Genus) # 786
rownames(spp.all) <-c (1:786)

saveRDS(spp.all, file="spp.all.Rdata")

spp.all <- readRDS("spp.all.Rdata")
## checking if the genera in spp.genera and spp.all match
spp.genera[,1] == spp.all[,1] # checking if the genera in both table are the same
# all are TRUE, ok!

## Selecting genera to be included in the analyses (for a genus to be considered in the analyses, we chose a threshold of climbing at least 75% of its diversity to be of climbing neotropical species)
percent<-data.frame(spp.all$Genus, c(spp.genera$Nr/spp.all$Nr))
colnames(percent)<-c("Genus","percent_nt_climbers")
percent<-percent[order(-percent$percent_nt_climbers),]
genera75<-subset(percent, percent_nt_climbers>=0.75) # 244 genera

saveRDS(genera75, file="genera75.Rdata")

## Renaming the tips of Smith & Brown's (2018) supertree to contain only genera
n<-tree$tip.label 
n<-strsplit(n,"_") 
list<-list()
for (i in 1:length(n)){
  if (length(n[[i]]) > 2) {
    list[i]<-paste0(n[[i]][1:2])
  } else {
    list[i]<-paste0(n[[i]])
  }
}
list<-do.call(rbind,list) 
list<-as.character(list)
tree$tip.label<-list

saveRDS(tree, file="treegenera.Rdata")

# function to get node ages -> THAIS: n?o sei se deixo isso aqui ja q tem um script s? com as fun??es
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

## Getting the ages of genera
treegenera <- readRDS("treegenera.Rdata")
library(phangorn)
genera <- as.character(genera75$Genus)
ages <- data.frame(nrow=length(genera),ncol=3)
for(i in 1:length(genera)){
  
  tip_numbers <- grep(genera[i], treegenera$tip.label) 
  if (length(tip_numbers)==0) {
    ages[i,1]<-genera[i]
    ages[i,2:3]<-"/"
    next
  } 
  
  if (length(tip_numbers)==1){
    stem_node <- Ancestors(treegenera, tip_numbers, type = "parent")
    stem_age <- get.node.age(treegenera)[stem_node]
    ages[i,1]<-genera[i]
    ages[i,2]<-stem_age
    ages[i,3]<-"/"
    next 
  }
  if (length(tip_numbers)>1){ 
    crown_node <- mrca.phylo(treegenera, tip_numbers)
    stem_node <- Ancestors(treegenera, crown_node, type = "parent")
    stem_age <- get.node.age(treegenera)[stem_node]
    crown_age <- get.node.age(treegenera)[crown_node]
    
    ages[i,1]<-genera[i]
    ages[i,2]<-stem_age
    ages[i,3]<-crown_age
  }
}
colnames(ages)<-c("Genus","Stem_Age", "Crown_Age")
ages<-ages[order(ages$Genus),]
ages<-subset(ages, Stem_Age!="/")
length(ages$Genus) # 176
rownames(ages)<-c(1:176)
#saveRDS(ages, file = "ages.Rdata")

##############################################################################
##### Table for the Magall?n & Sanderson (2001) analyses with crown ages #####

## creating pre-table with the genera, their crown ages and number of species
library(dplyr)
f <- spp.all %>% filter(Genus %in% ages$Genus)
f$Genus == ages$Genus # checking if the genera were correctly selected
                      # if all say TRUE, ok!
t<-cbind(ages[,-2],f$Nr)
t<-subset(t, Crown_Age!="/")
colnames(t)<-c("Genus", "Crown_Age", "Nr")

## Excluding genera that present more than one climbing mechanism
genera_list <- t$Genus 
include <- data.frame(matrix(nrow=0,ncol=2))
for(i in 1:length(genera_list)){
  subset <- climbers_final[climbers_final$genus==genera_list[i],]
  if(length(table(subset$CM))==1) {
    include <- rbind(include, c(genera_list[i], subset$CM[1]))
  } 
}
colnames(include) <- c("Genus","CM")
table <- merge(t, include, by="Genus") 

## Assiging genera to major angiosperm clades
for (k in 1:length(table$Genus)) {
  i <- which(climbers_final$genus == table$Genus[k])
  for (p in 1:length(APGIV$family)) {
    m<-unique(climbers_final$family[i])
    if (m == APGIV$family[p]) {
      table[k,"Family"] <- APGIV$family[p]
      table[k,"Order"] <- APGIV$order[p]
      table[k,"Clade"] <- APGIV$clade[p]
    }
  }
}

## Adjusting major clades of the table to be just "Magnoliids", "Monocots", "Superrosids" and "Superasterids" (and 'Ranunculales') sensu Magall?n et al. (2015)
table[,"Clade_M2015"]<-NA # 'M2015' is in reference to Magall?n et al. (2015)
for (i in 1:length(table$Clade)) {
  if (table$Clade[i] == "Asterids") {
    table$Clade_M2015[i] <- "Superasterids"
    next  
  }
  if (table$Clade[i] == "Rosids") {
    table$Clade_M2015[i] <- "Superrosids"
    next
  }
  if (table$Clade[i] == "Core_eudicots") {
    table$Clade_M2015[i] <- "Superasterids"
    next
  }
  if (table$Clade[i] == "Eudicots") {
    table$Clade_M2015[i] <- "Ranunculales"
  } else {
    table$Clade_M2015[i] <- table$Clade[i]
  }
}

## Creating final tables for Magall?n & Sanderson's (2001) analyses
a<-data.frame(taxa=table$Genus,
                diversity=table$Nr,
                node=c("CG"),
                age_mean=as.numeric(table$Crown_Age),
                clade=table$Clade_M2015,
                CM=table$CM,
                stringsAsFactors = FALSE)

d<-data.frame(taxa=c("bg_clade","Magnoliids", "Monocots", "Superrosids", "Superasterids", "Ranunculales"),
              diversity=c(352000,10293,69335,87302,111856, 4510), 
              node=c("bg","CG","CG","CG","CG","CG"),
              age_mean=as.numeric(c(139.4,132.4,133.2,122.4,122.6,114.8)),
              clade=NA,
              CM=NA,
              stringsAsFactors = FALSE)

bg.ang<-rbind(d[1,],a) # 'ang' refers to Angiosperms

bg.mag<-rbind(d[2,],subset(a, clade=="Magnoliids"))
bg.mag[1,1]<-"bg_clade"
bg.mag[1,3]<-"bg"

bg.mon<-rbind(d[3,],subset(a, clade=="Monocots"))
bg.mon[1,1]<-"bg_clade"
bg.mon[1,3]<-"bg"

bg.sros<-rbind(d[4,],subset(a, clade=="Superrosids"))
bg.sros[1,1]<-"bg_clade"
bg.sros[1,3]<-"bg"

bg.sast<-rbind(d[5,],subset(a, clade=="Superasterids"))
bg.sast[1,1]<-"bg_clade"
bg.sast[1,3]<-"bg"

bg.ranunc<-rbind(d[6,],subset(a, clade=="Ranunculales"))
bg.ranunc[1,1]<-"bg_clade"
bg.ranunc[1,3]<-"bg"

bg.all<-rbind(d,a)

## list with all final tables
bg.clades<-vector(mode = "list", length = 7)
bg.clades[[1]]<-bg.ang
bg.clades[[2]]<-bg.mag
bg.clades[[3]]<-bg.mon
bg.clades[[4]]<-bg.sros
bg.clades[[5]]<-bg.sast
bg.clades[[6]]<-bg.ranunc
bg.clades[[7]]<-bg.all
saveRDS(bg.clades, file = "bg.clades.crown.Rdata")

##############################################################################
##### Table for the Magall?n & Sanderson (2001) analyses with stem ages #####

## creating pre-table with the genera, their stem ages and number of species
library(dplyr)
f <- spp.all %>% filter(Genus %in% ages$Genus)
f$Genus == ages$Genus # checking if the genera were correctly selected
                      # if all say TRUE, ok!
t<-cbind(ages[,1:2],f$Nr)
colnames(t)<-c("Genus", "Stem_Age", "Nr")

## Excluding genera that present more than one climbing mechanism
genera_list <- t$Genus 
include <- data.frame(matrix(nrow=0,ncol=2))
for(i in 1:length(genera_list)){
  subset <- climbers_final[climbers_final$Genus==genera_list[i],]
  if(length(table(subset$CM))==1) {
    include <- rbind(include, c(genera_list[i], subset$CM[1]))
  } 
}
colnames(include) <- c("Genus","CM")
table <- merge(t, include, by="Genus") 

## Assiging genera to major angiosperm clades
for (k in 1:length(table$Genus)) {
  i <- which(climbers_final$Genus == table$Genus[k])
  for (p in 1:length(APGIV$family)) {
    m<-unique(climbers_final$Family[i])
    if (m[1] == APGIV$family[p]) {
      table[k,"Family"] <- APGIV$family[p]
      table[k,"Order"] <- APGIV$order[p]
      table[k,"Clade"] <- APGIV$clade[p]
    }
  }
}

## Adjusting major clades of the table to be just "Magnoliids", "Monocots", "Superrosids" and "Superasterids" (and 'Ranunculales') sensu Magall?n et al. (2015)
table[,"Clade_M2015"]<-NA # 'M2015' is in reference to Magall?n et al. (2015)
for (i in 1:length(table$Clade)) {
  if (table$Clade[i] == "Asterids") {
    table$Clade_M2015[i] <- "Superasterids"
    next  
  }
  if (table$Clade[i] == "Rosids") {
    table$Clade_M2015[i] <- "Superrosids"
    next
  }
  if (table$Clade[i] == "Core_eudicots") {
    table$Clade_M2015[i] <- "Superasterids"
    next
  }
  if (table$Clade[i] == "Eudicots") {
    table$Clade_M2015[i] <- "Ranunculales"
  } else {
    table$Clade_M2015[i] <- table$Clade[i]
  }
}

## Creating final tables for Magall?n & Sanderson's (2001) analyses
a<-data.frame(taxa=table$Genus,
              diversity=table$Nr,
              node=rep("SG",length(table$Nr)),
              age_mean=as.numeric(table$Stem_Age),
              clade=table$Family,
              CM=table$CM,
              stringsAsFactors = FALSE)

d<-data.frame(taxa=c("bg_clade","Magnoliids", "Monocots", "Superrosids", "Superasterids", "Ranunculales"),
              diversity=c(352000,10293,69335,87302,111856, 4510), 
              node=c("bg","SG","SG","SG","SG","SG"),
              age_mean=as.numeric(c(325.1, 134.6, 135.8,123.7,123.7,137.7)),
              clade=NA,
              CM=NA,
              stringsAsFactors = FALSE)

bg.ang<-rbind(d[1,],a) # 'ang' refers to Angiosperms

bg.mag<-rbind(d[2,],subset(a, clade=="Magnoliids"))
bg.mag[1,1]<-"bg_clade"
bg.mag[1,3]<-"bg"

bg.mon<-rbind(d[3,],subset(a, clade=="Monocots"))
bg.mon[1,1]<-"bg_clade"
bg.mon[1,3]<-"bg"

bg.sros<-rbind(d[4,],subset(a, clade=="Superrosids"))
bg.sros[1,1]<-"bg_clade"
bg.sros[1,3]<-"bg"

bg.sast<-rbind(d[5,],subset(a, clade=="Superasterids"))
bg.sast[1,1]<-"bg_clade"
bg.sast[1,3]<-"bg"

bg.ranunc<-rbind(d[6,],subset(a, clade=="Ranunculales"))
bg.ranunc[1,1]<-"bg_clade"
bg.ranunc[1,3]<-"bg"

bg.all<-rbind(d,a)

## List with all final tables
bg.clades<-vector(mode = "list", length = 7)
bg.clades[[1]]<-bg.ang
bg.clades[[2]]<-bg.mag
bg.clades[[3]]<-bg.mon
bg.clades[[4]]<-bg.sros
bg.clades[[5]]<-bg.sast
bg.clades[[6]]<-bg.ranunc
bg.clades[[7]]<-bg.all

saveRDS(bg.clades, file = "REVIEW_bg.clades.stem.Rdata")





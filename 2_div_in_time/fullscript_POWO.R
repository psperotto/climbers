#setwd("G:/Meu Drive/Papers/Diversificação/climbers")

# readig packages
library(beepr)
library(phangorn)
library(ape)
library(dplyr)
library(arsenal)
library(tidyverse)
library(geiger)

# reading the data
# setwd("G:/Meu Drive/Papers/Diversificação/climbers/Data")
climbers<-read.csv("G:/Meu Drive/Papers/Diversificação/climbers/Data/climber_database.csv", stringsAsFactors = F)
p<- read.delim("G:/Meu Drive/Papers/Diversificação/wcvp_v7_dec_2021.txt", sep="|", header=T) 
powo<-subset(p, rank=="SPECIES") # coluna 'taxon_name' pra nomes
#rm(p, powo)

# filtering powo to climbing accepted species
c <- powo %>% filter(taxon_name %in% climbers$Species) # works fine
climbers_acc <- subset(c, taxonomic_status=="Accepted") # accepted names
climbers_exc <- subset(c, taxonomic_status!="Accepted") # non accepted names

# counting species of climbers per family and per genera
spp.family<-as.data.frame(table(unlist(climbers_acc$family)))
colnames(spp.family)<-c("Family", "Nr")
spp.family<-spp.family[order(-spp.family$Nr),]

spp.genera<-as.data.frame(table(unlist(climbers_acc$genus)))
colnames(spp.genera)<-c("Genus", "Nr")
#spp.genera<-spp.genera[order(-spp.genera$Nr),]

#saveRDS(spp.family, file="G:/Meu Drive/Papers/Diversificação/climbers/Data/spp.family.POWO.Rdata")
#saveRDS(spp.genera, file="G:/Meu Drive/Papers/Diversificação/climbers/Data/spp.genera.POWO.Rdata")

# finding all species per genera of climbers 
sppall.POWO <- powo %>% filter(genus %in% spp.genera$Genus)
sppall.POWO <- subset(sppall.POWO,taxonomic_status=="Accepted")
sppall.POWO <- as.data.frame(table(unlist(sppall.POWO$genus)))
colnames(sppall.POWO)<-c("Genus", "Nr") 
#sppall.POWO <- sppall.POWO[order(-sppall.POWO$Nr),] # agora foi!!

# calculating the proportions
percent.POWO<-data.frame(sppall.POWO$Genus, c(spp.genera$Nr/sppall.POWO$Nr))
# isso ta mto estranho, olhar melhor
colnames(percent.POWO)<-c("Genus","perc_NT_spp")
percent.POWO<-percent.POWO[order(-percent.POWO$perc_NT_spp),]

# selecting only genera with >75%
genera75.POWO<-subset(percent.POWO, perc_NT_spp>=0.75) # 184 genera
# saveRDS(genera75.POWO, file="G:/Meu Drive/Papers/Diversificação/climbers/Data/genera75.POWO.Rdata")

# read tree
tree <- read.tree(file="GBMB.tre")
x<-as.data.frame(tree$node.label)

# renaming the tips of the tree to contain only the genus
n<-tree$tip.label
n<-strsplit(n,"_") #n is a list

# selecting only the first object in the vectors within the list 
list<-list()
for (i in 1:length(n)){
  if (length(n[[i]]) > 2) {
    list[i]<-paste0(n[[i]][1:2])
  } else {
    list[i]<-paste0(n[[i]])
  }
}
list<-do.call(rbind,list) # returns a dataframe with
list<-as.character(list)
tree$tip.label<-list

# saving the tree with only genera as tips
saveRDS(tree, file="treegenera.RDS")

# getting age of stem and crown node
treegenera <- readRDS("G:/Meu Drive/Papers/Diversificação/climbers/Data/treegenera.Rdata")

library(phangorn)
library(ape)

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

genera<-as.character(genera75.POWO$Genus)
results <- data.frame(nrow=length(genera),ncol=3)


for(i in 1:length(genera)){
  
  tip_numbers <- grep(genera[i], treegenera$tip.label) 
  if (length(tip_numbers)==0) {
    results[i,1]<-genera[i]
    results[i,2:3]<-"/"
    # aqui eu disse q se tip_numbers for um integer de tamanho 0, ? pra escrever 
    # o nome do genero e / nas duas colunas -> FUNCIONA!
    # nos 'if' statements, se a condi??o n?o ? fulfilled, o loop ignora e segue pro proximo passo
    next
  } 
  
  if (length(tip_numbers)==1){
    stem_node <- Ancestors(treegenera, tip_numbers, type = "parent")
    stem_age <- get.node.age(treegenera)[stem_node]
    results[i,1]<-genera[i]
    results[i,2]<-"/"
    results[i,3]<-stem_age
    next # precisa desse 'next' aqui pra mandar o loop NAO seguir adiante 
    # e ir pro proximo 'i' de 'genera[i]' em vez de seguir adiante no codigo
  }
  if (length(tip_numbers)>1){ # padronizei o codigo aqui pra 'girar em torno' dos lengths
    crown_node <- mrca.phylo(treegenera, tip_numbers)
    stem_node <- Ancestors(treegenera, crown_node, type = "parent")
    crown_age <- get.node.age(treegenera)[crown_node]
    stem_age <- get.node.age(treegenera)[stem_node]
    
    results[i,1]<-genera[i]
    results[i,2]<-crown_age
    results[i,3]<-stem_age
  }
}

########## preparando a tabela q vai servir de base pros calculos de net.div ##########

colnames(results)<-c("Genus", "Crown_Age","Stem_Age")
#results<-results[order(results$Genus),]
results<-subset(results, Stem_Age!="/") # deixando s? generos que tem Stem Age
length(results$Genus) # 143
#rownames(results)<-c(1:119)
#saveRDS(results, file = "G:/Meu Drive/Papers/Diversificação/climbers/Data/crown_stem_ages.POWO.Rdata")

f<-sppall.POWO$Genus %in% results$Genus
f<-as.numeric(which(f==T))
d<-sppall.POWO[f,]
d<-d[order(d$Genus),]
d
results2<-cbind(results,d$Nr)
colnames(results2)[4]<-"Nr"
#rm(d,f,results2)

setwd("C:/Users/patri/Google Drive/Papers/Diversificação/climbers")

# loading packages
library(lcvplants)
library(beepr)
library(phangorn)
library(ape)
library(dplyr)
library(arsenal)
library(tidyverse)

########## preparing the data ##########

# reading the data
climbers<-read.csv("climber_database.csv", stringsAsFactors = F)

# checking status of species with LCVP
h<-LCVP(climbers$Species)
beep("mario")
saveRDS(h,file="LCVP_climbers.Rdata")
h<-readRDS(file="LCVP_climbers.Rdata")

accepted<-subset(h, Status=="accepted")
problems<-subset(h, Status!="accepted")

# counting genera of climbers
genera<-unique(accepted$Genus) # 689 generos, 100 a menos que o database original...

# counting species of NT climbers in every genus
sppclimbers<-as.data.frame(table(unlist(accepted$Genus)))
colnames(sppclimbers)<-c("Genus", "Nr")

# counting all species for every genus of NT climbers
list <- list()
for (i in 1:length(genera))  {
  list[[i]]<-LCVP(genera[[i]], genus_tab = T)
}
beep("mario")

# glueing together all species for every genus of climbers and counting them
x<-do.call(rbind,list)
x<-subset(x, Status=="accepted")
sppall<-as.data.frame(table(unlist(x$Genus)))
      # sppall has 705 genera, how?? try to solve it tomorrow

colnames(sppall)<-c("Genus", "Nr")
saveRDS(sppall, file="sppall.Rdata")
#sppall<-readRDS(file="sppall.Rdata")

# finding the genera in sppall but not in sppclimbers
z<-sppall$Genus %in% sppclimbers$Genus
z<-as.numeric(which(z==F))

# excluding these genera from sppall
sppall2<-sppall[-z,]
sppall2<-sppall2[order(sppall2$Genus),]
rownames(sppall2)<-c(1:689)
sppclimbers<-sppclimbers[order(sppclimbers$Genus),]

########## choosing which genera will be included in the analysis ##########

# calculating the proportions
percent<-data.frame(sppall2$Genus, c(sppclimbers$Nr/sppall2$Nr))
colnames(percent)<-c("Genus","perc_NT_spp")

# selecting only genera with >75%
genera75<-subset(percent, perc_NT_spp>=0.75)
genera75<-genera75[order(genera75$perc_NT_spp, decreasing=T),]
rownames(genera75)<-c(1:155)
saveRDS(genera75, file="genera75.Rdata")

########## extracting clade/genera ages from S&B 2018 tree #########

# read tree
tree <- read.tree(file="GBMB.tre")

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
treegenera <- readRDS("treegenera.RDS")

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
} # funçao pra pegar os node ages

genera<-as.character(genera75$Genus)
results <- matrix(nrow=length(genera),ncol=3)

for(i in 1:length(genera)){
  
  tip_numbers <- grep(genera[i], treegenera$tip.label) 
  if (length(tip_numbers)==0) {
    results[i,1]<-genera[i]
    results[i,2:3]<-"/"
    # aqui eu disse q se tip_numbers for um integer de tamanho 0, é pra escrever 
    # o nome do genero e / nas duas colunas -> FUNCIONA!
    # nos 'if' statements, se a condição não é fulfilled, o loop ignora e segue pro proximo passo
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

#library(beepr)
beep("mario")

########## preparando a tabela q vai servir de base pros calculos de net.div ##########

# ajeitando o objeto 'results'
results<-as.data.frame(results)
colnames(results)<-c("Genus", "Crown_Age","Stem_Age")
results<-results[order(results$Genus),]
results<-subset(results, Stem_Age!="/") # deixando só generos que tem Stem Age
length(results$Genus) # 119
rownames(results)<-c(1:119)
#saveRDS(results, file = "crown_stem_ages.Rdata")

# adicionando colunas

# selecionar do sppall2 as linhas de Genus q sao iguais aos genera in 'results'
f<-sppall2$Genus %in% results$Genus
f<-as.numeric(which(f==T))

sppall3<-sppall2[f,]
sppall3<-sppall3[order(sppall3$Genus),]
rownames(sppall3)<-c(1:119)

# colar a coluna 'Nr' de sppall3 com results
results2<-cbind(results,sppall3$Nr)
colnames(results2)[4]<-"Nr"

saveRDS(results2, file = "partialresults.Rdata")

# criar mais uma coluna com o numero total de spp de cada genero -> OK!!
# criar outra coluna dizendo qual grande grupo (monocots, asteridae, rosidae, etc)
# criar outra coluna com o mecanismo de escalada de cada genero


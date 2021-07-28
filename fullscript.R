# setwd("~/Desktop/Colabs/Patricia_Climbers/climbers")
setwd("C:/Users/patri/Google Drive/Papers/Diversificação/climbers")

# loading packages
library(lcvplants)
library(beepr)
library(phangorn)
library(ape)
library(dplyr)
library(arsenal)
library(tidyverse)
library(geiger)

########## preparing the data ##########

# reading the data
climbers<-read.csv("climber_database.csv", stringsAsFactors = F)

# checking status of species with LCVP
h<-LCVP(climbers$Species)
beep("mario")
#saveRDS(h,file="LCVP_climbers.Rdata")
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
#saveRDS(sppall, file="sppall.Rdata")
sppall<-readRDS(file="sppall.Rdata")

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
genera75 <- readRDS("genera75.Rdata")

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
} # fun?ao pra pegar os node ages

genera<-as.character(genera75$Genus)
results <- matrix(nrow=length(genera),ncol=3)

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

#library(beepr)
beep("mario")

########## preparando a tabela q vai servir de base pros calculos de net.div ##########

# ajeitando o objeto 'results'
results<-as.data.frame(results)
colnames(results)<-c("Genus", "Crown_Age","Stem_Age")
results<-results[order(results$Genus),]
results<-subset(results, Stem_Age!="/") # deixando s? generos que tem Stem Age
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
#saveRDS(results2, file = "partialresults.Rdata")

# criar outra coluna com o mecanismo de escalada de cada genero #
climbers <- read.csv("climber_database.csv")

# excluding genera with more than one mechanism
genera_list <- unique(climbers$Genus)
to_exclude <- c()
to_include <- data.frame(matrix(nrow=0,ncol=2))
for(i in 1:length(genera_list)){
  subset <- climbers[climbers$Genus==genera_list[i],]
  if(length(table(subset$CM))==1) {
    to_include <- rbind(to_include, c(genera_list[i], subset$CM[1]))
  } else {
    to_exclude <- c(to_exclude, genera_list[i]) 
  }
}
print(to_exclude) #check those to exclude
colnames(to_include) <- c("Genus","CM")
table_final <- merge(results2, to_include, by="Genus")

# criar coluna dizendo qual grande grupo (monocots, asteridae, rosidae, etc) #
orders<-accepted$Genus %in% table_final$Genus
orders<-as.numeric(which(orders==T))
orders<-accepted[orders,]
orders<-distinct(orders,Genus,.keep_all = T)
orders<-orders[,c("Order","Genus")]
orders<-orders[order(orders$Genus),]
rownames(orders)<-c(1:117)
orders<-orders %>% add_column(Group =NA)

for (i in 1:length(orders$Order)) {
  if (orders[i,1]=="Piperales"){
    orders[i,3]<-"Magnoliids"
    next
  }
  if (orders[i,1] %in% c("Alismatales", "Asparagales", "Pandanales")==T){
    orders[i,3]<-"Monocots"
    next
  }
  if (orders[i,1] %in% c("Celastrales","Malpighiales","Cucurbitales","Fabales","Rosales","Sapindales")==T){
    orders[i,3]<-"Rosids"
    next
  }
  if (orders[i,1] %in% c("Gentianales","Asterales","Lamiales","Ericales","Cornales","Solanales","Icacinales")==T){
    orders[i,3]<-"Asterids"
    next
  }  
  if (orders[i,1] %in%  c("Dilleniales","Caryophyllales","Ranunculales")==T){
    orders[i,3]<-"Other"
  }
}
# it worked!
 #use rbind

########## calculando net.div ##########
# library(geiger)
net_div<-data.frame(matrix(nrow=length(table_final$Genus),ncol=4))
for (i in 1:length(table_final$Genus)) {
  net_div[i,1]<-table_final$Genus[i]
  net_div[i,2]<-bd.ms(phy=NULL, time = as.numeric(table_final[i,3]), 
                      n = as.numeric(table_final[i,4]),crown=F, epsilon = 0)
  net_div[i,3]<-bd.ms(phy=NULL, time = as.numeric(table_final[i,3]), 
                      n = as.numeric(table_final[i,4]),crown=F, epsilon = 0.5)
  net_div[i,4]<-bd.ms(phy=NULL, time = as.numeric(table_final[i,3]), 
                      n = as.numeric(table_final[i,4]),crown=F, epsilon = 0.9)
  }
colnames(net_div)<-c("Genus","e=0","e=0.5","e=0.9")
#saveRDS(net_div,file="net_div.Rdata")

# calculando limites de stem age com CI=95% #
limits<-data.frame(matrix(nrow=length(table_final$Genus),ncol=7))
for (i in 1:length(table_final$Genus)){
  limits[i,1]<-table_final$Genus[i]
  limits[i,2]<-stem.limits(time = as.numeric(table_final[i,3]), 
                           r=net_div[i,2], epsilon = 0, CI=.95)[1]
  limits[i,3]<-stem.limits(time = as.numeric(table_final[i,3]), 
                           r=net_div[i,2], epsilon = 0, CI=.95)[2]
  limits[i,4]<-stem.limits(time = as.numeric(table_final[i,3]), 
                           r=net_div[i,3], epsilon = 0.5, CI=.95)[1]
  limits[i,5]<-stem.limits(time = as.numeric(table_final[i,3]), 
                           r=net_div[i,3], epsilon = 0.5, CI=.95)[2]
  limits[i,6]<-stem.limits(time = as.numeric(table_final[i,3]), 
                           r=net_div[i,4], epsilon = 0.9, CI=.95)[1]
  limits[i,7]<-stem.limits(time = as.numeric(table_final[i,3]), 
                           r=net_div[i,4], epsilon = 0.9, CI=.95)[2]
}
colnames(limits)<-c("Genus","lb_0","ub_0","lb_0.5","ub_0.5","lb_0.9","ub_0.9")
#saveRDS(limits,file="stem_ages_CI95.Rdata")

###### graficos Magallon & Sanderson ####

table <- data.frame(taxa=table_final$Genus,
                       diversity=table_final$Nr,
                       node=c("SG"),
                       age_mean=table_final$Stem_Age,
                       stringsAsFactors = FALSE)

#teste com ages of angiosperms, monocots, etc#
b<-data.frame(taxa=c("bg_clade","Monocots", "Rosids", "Asterids"),
              diversity=c(295383,74273,10,10),
              node=c("CG","SG","SG","SG"),
              age_mean=c(209, 154,130,125),
              stringsAsFactors = FALSE)
table<-rbind(b,table)
## quando for fazer os graficos tipo da magallon e sanderson fazer pra angiospermas no geral e tb pros grandes grupos asteridae, rosidae, etc
# possible mean ages of angiosperm crown nodes: Magallon et al. 2015: 139.4 ; Li et al. 2019: 209 myr

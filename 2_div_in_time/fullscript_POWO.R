#setwd("G:/Meu Drive/Papers/Diversificação/climbers")
# "G:/Meu Drive/Papers/Diversificação/climbers/Data/

# reading packages
library(beepr)
library(phangorn)
library(ape)
library(dplyr)
library(arsenal)
library(tidyverse)
library(geiger)

##### reading the data #####
# setwd("G:/Meu Drive/Papers/Diversificação/climbers/Data")
APGIV<-read.csv("G:/Meu Drive/Papers/Diversificação/climbers/Data/APGIV.csv", sep=",") # dataset with APG IV families and orders
climbers<-read.csv("G:/Meu Drive/Papers/Diversificação/climbers/Data/climber_database.csv", stringsAsFactors = F)
p<- read.delim("G:/Meu Drive/Papers/Diversificação/wcvp_v7_dec_2021.txt", sep="|", header=T) 
p2<-subset(p, rank=="SPECIES") # coluna 'taxon_name' pra nomes
powo<-subset(p2, taxonomic_status=="Accepted")
# rm(p, p2, powo) # -> apagar esses 2 quando for fechar o R pra nao ficar com um R space gigante

##### generating objects #####
# filtering powo to climbing accepted species
c <- p2 %>% filter(taxon_name %in% climbers$Species) # works fine
climbers_acc <- subset(c, taxonomic_status=="Accepted") # accepted names
#climbers_exc <- subset(c, taxonomic_status!="Accepted") # non accepted names

# counting species of climbers per family
spp.family<-as.data.frame(table(unlist(climbers_acc$family)))
colnames(spp.family)<-c("Family", "Nr")
#write.csv(spp.family, file="climberfamilies.csv")

# counting species of climbers per genera
spp.genera<-as.data.frame(table(unlist(climbers_acc$genus)))
colnames(spp.genera)<-c("Genus", "Nr")
#spp.genera<-spp.genera[order(spp.genera$Nr),]

#saveRDS(spp.family, file="G:/Meu Drive/Papers/Diversificação/climbers/Data/spp.family.POWO.Rdata")
#saveRDS(spp.genera, file="G:/Meu Drive/Papers/Diversificação/climbers/Data/spp.genera.POWO.Rdata")

# finding all species per genera of climbers 
spp.all <- powo %>% filter(genus %in% spp.genera$Genus)
spp.all <- as.data.frame(table(unlist(spp.all$genus)))
colnames(spp.all)<-c("Genus", "Nr") 
#sppall.POWO <- sppall.POWO[order(-sppall.POWO$Nr),] # agora foi!!

##### selecting genera for the analyses #####
# calculating the proportions
spp.genera[,1] == spp.all[,1] # checking if the genera in both table are the same
percent<-data.frame(spp.all$Genus, c(spp.genera$Nr/spp.all$Nr))
colnames(percent)<-c("Genus","perc_NT_spp")
percent<-percent[order(-percent$perc_NT_spp),]

# selecting only genera with >75%
genera75<-subset(percent.POWO, perc_NT_spp>=0.75) # 184 genera
# saveRDS(genera, file="G:/Meu Drive/Papers/Diversificação/climbers/Data/genera75.POWO.Rdata")

##### extracting ages of genera from S&B tree for the analyses #####
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

# getting ages of stem nodes
treegenera <- readRDS("G:/Meu Drive/Papers/Diversificação/climbers/Data/treegenera.Rdata")

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

# getting the ages
genera<-as.character(genera75$Genus)
results <- data.frame(nrow=length(genera),ncol=2)
for(i in 1:length(genera)){
  
  tip_numbers <- grep(genera[i], treegenera$tip.label) 
  if (length(tip_numbers)==0) {
    results[i,1]<-genera[i]
    results[i,2]<-"/"
    # aqui eu disse q se tip_numbers for um integer de tamanho 0, ? pra escrever 
    # o nome do genero e / nas duas colunas -> FUNCIONA!
    # nos 'if' statements, se a condi??o n?o ? fulfilled, o loop ignora e segue pro proximo passo
    next
  } 
  
  if (length(tip_numbers)==1){
    stem_node <- Ancestors(treegenera, tip_numbers, type = "parent")
    stem_age <- get.node.age(treegenera)[stem_node]
    results[i,1]<-genera[i]
    results[i,2]<-stem_age
    next # precisa desse 'next' aqui pra mandar o loop NAO seguir adiante 
    # e ir pro proximo 'i' de 'genera[i]' em vez de seguir adiante no codigo
  }
  if (length(tip_numbers)>1){ # padronizei o codigo aqui pra 'girar em torno' dos lengths
    crown_node <- mrca.phylo(treegenera, tip_numbers)
    stem_node <- Ancestors(treegenera, crown_node, type = "parent")
    stem_age <- get.node.age(treegenera)[stem_node]
    
    results[i,1]<-genera[i]
    results[i,2]<-stem_age
  }
}

##### preparing the base table for the any analyses #####
colnames(results)<-c("Genus","Stem_Age")
results<-results[order(results$Genus),]
results<-subset(results, Stem_Age!="/") # deixando s? generos que tem Stem Age
length(results$Genus) # 143
#saveRDS(results, file = "G:/Meu Drive/Papers/Diversificação/climbers/Data/crown_stem_ages.POWO.Rdata")

# selecting climbing genera and their total diversity
f<-spp.all$Genus %in% results$Genus
f<-as.numeric(which(f==T))
d<-spp.all[f,]
d<-d[order(d$Genus),]
d$Genus == results$Genus # checking if the genera were correctly selected
t<-cbind(results,d$Nr)
colnames(t)<-c("Genus", "Stem_Age", "Nr")

# creating a column with the climbing mechanisms of genera #
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
tt <- merge(t, to_include, by="Genus") # 139 genera with >=75% of species being neotropical climbers, all with the same mechanism

##### table for the M&S analyses #####
a <- data.frame(taxa=tt$Genus,
                diversity=tt$Nr,
                node=c("SG"),
                age_mean=as.numeric(tt$Stem_Age),
                stringsAsFactors = FALSE)

p3 <- subset(p,rank=="GENUS")
p3 <- subset(p3,taxonomic_status=="Accepted")
orders<-p3$genus %in% a$taxa
orders<-as.numeric(which(orders==T))
orders<-p3[orders,]
orders<-distinct(orders,genus,.keep_all = T)
orders<-orders[,2:3]
orders<-orders[order(orders$genus),]
rownames(orders)<-c(1:139)
#orders<-orders %>% add_column(order = NA)
orders<-orders %>% add_column(clade = NA)

for (i in 1:length(orders$genus) {
  #if ( ) ### continue from here tomorrow
}


for (i in 1:length(orders$order)) {
  if (orders[i,1]=="Piperales"){
    orders[i,3]<-"Magnoliids"
    next
  }
  if (orders[i,1] %in% c("Alismatales", "Asparagales", "Pandanales")==T){
    orders[i,3]<-"Monocots"
    next
  }
  if (orders[i,1] %in% c("Celastrales","Malpighiales","Cucurbitales","Fabales","Rosales","Sapindales")==T){
    orders[i,3]<-"Superrosids"
    next
  }
  if (orders[i,1] %in% c("Gentianales","Asterales","Dilleniales","Lamiales","Ericales","Cornales","Caryophyllales","Solanales","Icacinales")==T){
    orders[i,3]<-"Superasterids"
    next
  }  
  if (orders[i,1]=="Ranunculales"){
    orders[i,3]<-"Ranunculales"
  }
}
# it worked!
#use rbind

# teste com ages of angiosperms, monocots, etc #
b<-data.frame(taxa=c("bg_clade","Magnoliids", "Monocots", "Superrosids", "Superasterids"),
              diversity=c(352000,10293,69335,87302, 111856), 
              node=c("bg","SG","SG","SG", "SG"),
              age_mean=as.numeric(c(325.1, 134.6, 135.8,123.7,123.7)), 
              stringsAsFactors = FALSE)

#table<-rbind(b,a) 
write.csv(table, file = "table.csv")

# tables for each background clade

bg.ang<-rbind(b[1,],a)
saveRDS(bg.ang, file="bg.ang.Rdata")

bg.mag<-rbind(b[2,],a)
bg.mag[1,1]<-"bg_clade"
bg.mag[1,3]<-"bg"
saveRDS(bg.mag, file="bg.mag.Rdata")

bg.mon<-rbind(b[3,],a)
bg.mon[1,1]<-"bg_clade"
bg.mon[1,3]<-"bg"
saveRDS(bg.mon, file="bg.mon.Rdata")

bg.sros<-rbind(b[4,],a)
bg.sros[1,1]<-"bg_clade"
bg.sros[1,3]<-"bg"
saveRDS(bg.sros, file="bg.sros.Rdata")

bg.sast<-rbind(b[5,],a)
bg.sast[1,1]<-"bg_clade"
bg.sast[1,3]<-"bg"
saveRDS(bg.sast, file="bg.sast.Rdata")

# lista pra rodar os graficos mais facilmente depois
bg.clades<-vector(mode = "list", length = 5)
bg.clades[[1]]<-bg.ang
bg.clades[[2]]<-bg.mag
bg.clades[[3]]<-bg.mon
bg.clades[[4]]<-bg.sros
bg.clades[[5]]<-bg.sast


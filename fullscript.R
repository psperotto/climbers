# setwd("~/Desktop/WCVP_special_issue/Patricia_Climbers/climbers")
setwd("G:/Meu Drive/Papers/Diversifica??o/climbers")

# FULL SCRIPT FOR COUNTING CLIMBER SPECIES
# AND FOR NETDIV ANALYSES

# loading packages
library(lcvplants)
library(beepr)
library(phangorn)
library(ape)
library(dplyr)
library(arsenal)
library(tidyverse)
library(geiger)

# Source useful functions:
source("utility_functions.R")

################################
# Preparing the data ##########
################################

# Reading the data
# -------------------- Make sure this table is always up to date ------------------ #
climbers <- read.csv("Data/climber_database.csv", stringsAsFactors = F)
climbers$taxized_names <- fix.names.taxize(climbers$taxized_names)

# checking status of species with LCVP
h <- LCVP(climbers$Species)
beep("mario")
# saveRDS(h,file="LCVP_climbers.Rdata")

#------------------------------- checkpoint
h <- readRDS(file="LCVP_climbers.Rdata")

accepted <- subset(h, Status=="accepted")
problems <- subset(h, Status!="accepted")

# counting genera of climbers
genera <- unique(accepted$Genus) # 689 generos, 100 a menos que o database original...

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
x <- do.call(rbind,list)
x <- subset(x, Status=="accepted")
sppall<-as.data.frame(table(unlist(x$Genus)))
      # sppall has 705 genera, how?? try to solve it tomorrow

colnames(sppall)<-c("Genus", "Nr")
#saveRDS(sppall, file="sppall.Rdata")

#------------------------------- checkpoint
sppall <- readRDS(file="sppall.Rdata")

# finding the genera in sppall but not in sppclimbers
z <- sppall$Genus %in% sppclimbers$Genus
z <- as.numeric(which(z==F))

# excluding these genera from sppall
sppall2 <- sppall[-z,]
sppall2 <- sppall2[order(sppall2$Genus),]
rownames(sppall2) <- 1:nrow(sppall2)
sppclimbers <- sppclimbers[order(sppclimbers$Genus),]

########## choosing which genera will be included in the analysis ##########

# calculating the proportions
percent<-data.frame(sppall2$Genus, c(sppclimbers$Nr/sppall2$Nr))
colnames(percent)<-c("Genus","perc_NT_spp")

# selecting only genera with >75%
genera75<-subset(percent, perc_NT_spp>=0.75)
genera75<-genera75[order(genera75$perc_NT_spp, decreasing=T),]
rownames(genera75)<- 1:nrow(genera75)
#saveRDS(genera75, file="genera75.Rdata")

#------------------------------- checkpoint
genera75 <- readRDS("genera75.Rdata")

########## extracting clade/genera ages from S&B 2018 tree #########

# read tree
tree <- read.tree(file="Data/GBMB.tre")

# renaming the tips of the tree to contain only the genus
n <- tree$tip.label
n <- strsplit(n,"_") #n is a list

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

#------------------------------- checkpoint
# getting age of stem and crown node
treegenera <- readRDS("treegenera.RDS")

library(phangorn)
library(ape)

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
rownames(results)<-1:nrow(results)
#saveRDS(results, file = "crown_stem_ages.Rdata")

# adicionando colunas

# selecionar do sppall2 as linhas de Genus q sao iguais aos genera in 'results'
f<-sppall2$Genus %in% results$Genus
f<-as.numeric(which(f==T))

sppall3<-sppall2[f,]
sppall3<-sppall3[order(sppall3$Genus),]
rownames(sppall3)<-1:nrow(sppall3)

# colar a coluna 'Nr' de sppall3 com results
results2<-cbind(results,sppall3$Nr)
colnames(results2)[4]<-"Nr"
#saveRDS(results2, file = "partialresults.Rdata")

# criar outra coluna com o mecanismo de escalada de cada genero #
climbers <- read.csv("Data/climber_database.csv")

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
write.csv(table_final, file="Data/table_final_diver.csv", row.names=F)


#
table_final <- read.csv("Data/table_final_diver.csv")
table_final$node <- "SG"
table_final <- table_final[c("Genus","Nr","node","Stem_Age","CM")]
colnames(table_final) <- c("taxa","diversity","node","age_mean","CM")

results_ms <- table_final
results_ms$bd_0 <- NA
results_ms$bd_09 <- NA
for(i in 1:nrow(results_ms)){
  time <- results_ms[i,"age_mean"]
  n <- results_ms[i,"diversity"]
  results_ms$bd_0[i] <- geiger::bd.ms(phy=NULL, time, n, crown=F, epsilon = 0)
  results_ms$bd_09[i] <- geiger::bd.ms(phy=NULL, time, n, crown=F, epsilon = 0.9)
}

#results_ms <- subset(results_ms, results_ms$bd_0!=0)
#results_ms <- subset(results_ms, results_ms$bd_09!=0)


#

boxplot(results_ms$bd_0~results_ms$CM, ylim=c(0,3))
boxplot(results_ms$bd_09~results_ms$CM, ylim=c(0,1))

library(PMCMR)

results_ms <- subset(results_ms, results_ms$CM %in% c(1:4))

table(results_ms$CM)

result <- posthoc.kruskal.conover.test(results_ms$bd_09,results_ms$CM)
result <- posthoc.kruskal.conover.test(results_ms$bd_0,results_ms$CM)




#

example <- get.template(save.file=F)

mech1 <- table_final[which(table_final$CM==1),]
mech1 <- rbind(example[1,1:4], mech1[,-5])
results_mech1 <- get.tail.probs(table=mech1)
plottailprobs(mech1, results_mech1)
View(mech2)

mech2 <- table_final[which(table_final$CM==2),]
mech2 <- rbind(example[1,1:4], mech2[,-5])
results_mech2 <- get.tail.probs(table=mech2)
plottailprobs(mech2, results_mech2)
#dev.off()

mech3 <- table_final[which(table_final$CM==3),]
mech3 <- rbind(example[1,1:4], mech3[,-5])
results_mech3 <- get.tail.probs(table=mech3)
plottailprobs(mech3, results_mech3)

mech4 <- table_final[which(table_final$CM==4),]
mech4 <- rbind(example[1,1:4], mech4[,-5])
results_mech4 <- get.tail.probs(table=mech4)
plottailprobs(mech4, results_mech4)

mech5 <- table_final[which(table_final$CM==5),]
mech5 <- rbind(example[1,1:4], mech5[,-5])
results_mech5 <- get.tail.probs(table=mech5)
plottailprobs(mech5, results_mech5)

mech6 <- table_final[which(table_final$CM==6),]
mech6 <- rbind(example[1,1:4], mech6[,-5])
results_mech6 <- get.tail.probs(table=mech6)
plottailprobs(mech6, results_mech6)
beepr::beep("complete")


#-----------------

######## PLOT FIG 3 & 4 ######

results <- get.tail.probs(table=example)


#taxa diversity node age_mean age_upper age_lower
#1      bg_clade    262196   bg   132.00    137.00    127.00



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
              diversity=c(295383,74273,10,10), # procurar o numero real de spp de rosids e asterids
              node=c("CG","SG","SG","SG"),
              age_mean=c(209, 154, 130,125), # mudar esses valores depois de determinar as idades correta a serem usadas
              stringsAsFactors = FALSE)
table<-rbind(b,table)
## quando for fazer os graficos tipo da magallon e sanderson fazer pra angiospermas no geral e tb pros grandes grupos asteridae, rosidae, etc
# possible mean ages of angiosperm crown nodes: Magallon et al. 2015: 139.4 ; Li et al. 2019: 209 myr

write.csv(table, file = "table.csv")

table$age_mean<-as.numeric(table$age_mean) # yaaaay :D

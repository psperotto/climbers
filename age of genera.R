setwd("C:/Users/patri/Google Drive/Papers/Diversificação/climbers")

######### extracting clade/genera ages from S&B 2018 tree

## read tree
tree <- read.tree(file="GBMB.tre")

## renaming the tips of the tree to contain only the genus
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


#-----------------------
# getting age of stem and crown node

# load tree
# tree <- readRDS("treegenera.RDS")

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

genera <- as.character(sppclimbers_recount[,1])
#saveRDS(sppclimbers_recount, file = "sppclimbers.Rdata")

results <- matrix(nrow=length(genera),ncol=3)
# tem generos no genera q nao estao amostrados na tree, 
# fazer um if else pra pular pro proximo qdo encontrar 
# um desses, e tb um if else pra quando tem um tip só

for(i in 1:length(genera)){
  
tip_numbers <- grep(genera[i], tree$tip.label) 
 if (length(tip_numbers)==0) {
    results[i,1]<-genera[i]
    results[i,2:3]<-"/"
    # aqui eu disse q se tip_numbers for um integer de tamanho 0, é pra escrever 
    # o nome do genero e / nas duas colunas -> FUNCIONA!
    # nos 'if' statements, se a condição não é fulfilled, o loop ignora e segue pro proximo passo
    next
    } 
  
  if (length(tip_numbers)==1){
    stem_node <- Ancestors(tree, tip_numbers, type = "parent")
    stem_age <- get.node.age(tree)[stem_node]
    results[i,1]<-genera[i]
    results[i,2]<-"/"
    results[i,3]<-stem_age
    next # precisa desse 'next' aqui pra mandar o loop NAO seguir adiante 
         # e ir pro proximo 'i' de 'genera[i]' em vez de seguir adiante no codigo
    }
  if (length(tip_numbers)>1){ # padronizei o codigo aqui pra 'girar em torno' dos lengths
    crown_node <- mrca.phylo(tree, tip_numbers)
    stem_node <- Ancestors(tree, crown_node, type = "parent")
    crown_age <- get.node.age(tree)[crown_node]
    stem_age <- get.node.age(tree)[stem_node]

    results[i,1]<-genera[i]
    results[i,2]<-crown_age
    results[i,3]<-stem_age
  }
}

#library(beepr)
beep("mario")
# WE DIDI IT BOYSSSSSSSSSSSSSSS!!!!!! IT FUCKING WORKED!!!!!!!

colnames(results)<-c("Genus", "Crown Age","Stem Age")
#as.data.frame(results)
saveRDS(results, file = "crown_stem_ages.Rdata")

  # criar mais uma coluna com o numero total de spp de cada genero
  # criar outra coluna dizendo qual grande grupo (monocots, asteridae, etc)
  # criar outra coluna com o mecanismo de escalada de cada genero

######### funçoes do phangorn
tree2 <- rtree(10)
plot(tree2, show.tip.label = FALSE)
nodelabels()
tiplabels()
Ancestors(tree2, 1:3, "parent") # te dá todos os ancestrais dos tips indicados 
                             # (no caso, o '1:3' ali)

Children(tree2, 11) # te dá os nós descendentes do nó indicado 
                    # (nó 11 tem os nós 12 e 15 como descendentes)

Descendants(tree2, 11, "tips") # te dá os tips descendetes do nó indicado
                               # (nó 11 tem os tips 1:10 como descendentes)
Siblings(tree2, 3) # te dá os tips irmãos do tip indicado
                   # tip 3 tem otip 2 como irmão)

# Siblings of all nodes
Siblings(tree2) # te dá os nós irmãos dos nós indicados
mrca.phylo(tree2, 1:3) # te dá o mrca dos tips indicados (tips 1 a 3 tem o nó 13 como mrca)
mrca.phylo(tree2, match(c("t1", "t2", "t3"), tree2$tip)) # te dá o mrca dos tips (unlabeled) indicados
                                                         # (mrca dos tips t1, t2 e t3 é o nó 15)
mrca.phylo(tree2)
# same as mrca(tree), but faster for large trees
# te retorna uma matriz com todos mrca entre 2 tips (unlabeled tips)
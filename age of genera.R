setwd("C:/Users/patri/Google Drive/Papers/Diversificação/climbers")

######### extracting clade/genera ages from S&B 2018 tree

## read tree
tree <- read.tree(file="GBMB.tre")

## renaming the tips of the tree to contain only the genus
n<-tree$tip.label
n<-strsplit(n,"_") #n is a list

  # selecting only the first object un the vectors within the list 
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






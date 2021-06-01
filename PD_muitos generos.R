### Phylogenetic diversity maps ####

##Vamos criar o PD para todo o campo rupestre! Para isso, faremos em 3 etapas: 3 vezes
#calcularemos o PD


#rm(list=ls()) # uncomment to clean your environment

library(raster)
library(dplyr)
library(maptools)
library(sf)
library(rgdal)
library(picante)
data(wrld_simpl)


#wd <- "C:/Users/raque/OneDrive/Documentos/Lista de espécies - CR/"
#setwd(wd)

# Carregando a arvore de genero
tree <- read.tree("tree_cr_genus.txt")

# vamos carregar os poligonos que serão usados para criar a matriz e os modelos binários
wd_polygon <- "polygons/shape_files/" 
wd_current <- "rasters/CHELSA_cur_V1_2B_r2_5m/2_5min/"

setwd(paste(wd, wd_polygon, sep=""))
list.files(pattern="_current.shp") -> files_cur #pegando todos os polígonos de distribuição que criamos para cada genero
sapply(files_cur, readShapePoly, simplify = F) -> modelos
sub("_current.shp", "", files_cur) -> labs
labs <- sub(" ","_", labs)
names(modelos) <- labs

###############################################
### PARTE 1: primeiro PD - até latitude -14 ###
###############################################

# criando o raster vazio para os modelos binários
setwd(paste(wd, wd_current, sep=""))
mask <- raster("bio_1.tif")
mask <- crop(mask, extent(-60, -30, -14, -8.5)) # cortando para a faixa desejada (no meu caso, norte do campo rupestre)
mask[!is.na(mask)] <- 0 # mudando todos os valores para 0
plot(mask, axes = T, box = T, legend = T, main = "Neotropics")

rasters1<-list() # criando lista vazia para empilhar os rasters
matrixPD1 <- matrix(0, ncol = 1, nrow = 95040) # criando matriz vazia para fazer matrixPD (nessa hora vc provavelmente terá que mudar o numero de linhas)
matrixPD1 <- data.frame(matrixPD1)

for (i in 1:length(modelos)) { # criando matrixPD e rasterizando poligonos dos modelos para empilhar
  Sys.time() -> start_time
  r1 <- rasterize(modelos[[i]], mask, field=1)
  r1[is.na(r1)] <- 0
  resample(r1, mask) -> r1
  mask(r1, mask) -> rasters1[[i]]
  meh <- as.data.frame(values(rasters1[[i]]))
  tip <- sub(" ", "_", labs[i])
  colnames(meh) <- tip
  matrixPD1 <- cbind(matrixPD1, meh)
  print(c(i, "in", length(modelos)))
  Sys.time() -> end_time  
  print(end_time-start_time)
}

names(rasters1) <- names(modelos)
matrixPD1[is.na(matrixPD1)] <- 0 # substituindo os NAs por 0s
matrixPD1 <- matrixPD1[,2:614]


# Fazendo matrixPD e arvore conversarem #
tree1 <- drop.tip(tree, setdiff(tree$tip.label, colnames(matrixPD1))) # cortando a arvore para as colunas da matrixPD (caso precise)
matrixPD1 <- matrixPD1[,c(which(colnames(matrixPD1) %in% tree1$tip.label))] # cortando a matrixPD de acordo com os terminais da arvore (caso precise)

setwd("Results")
write.csv(matrixPD1, file="matrixPD1.csv")

# Rodando analise de PD #
pd_final1 <- pd(matrixPD1, tree1, include.root=TRUE) # dá um valor de PD por célula (SR é o species richness?)


# Empilhando os rasters e criando mapa de SR#
rasters1 <- rasters1[which(names(rasters1) %in% tree$tip.label)]
sr1 <- calc(stack(rasters1), sum)
plot(sr1, main = "Current SR")

# criando mapas de PD #
pdMAP <- mask
values(pdMAP) <- pd_final1[,1]

plot(pdMAP, main = "PD1")
plot(wrld_simpl, add = T)

writeRaster(pdMAP, "PD1.tif")
writeRaster(sr1, "SR1.tif")


# which cells have higher or lower PD values as expected by their SR?
plot(pd_final1[,2], pd_final1[,1], xlab = "species richness", ylab = "pd values")

# using linear regression (not sure if the best method...)
linearRegression <- lm(pd_final1[,1]~pd_final1[,2])
abline(linearRegression, col="red") # adding trend line
summary(linearRegression)
res <- as.numeric(linearRegression$residuals)
resMap1 <- mask
values(resMap1) <- res
plot(resMap1, main="residual") # positive values = more pd than expected by species richness/ negative values = less pd than expected by species richness
plot(wrld_simpl, add = T)
writeRaster(resMap1, "Res1.tif")


###############################################
### PARTE 2: segundo PD - até latitude -19.5 ###
###############################################

#criando o raster vazio para os modelos binários
setwd(paste(wd, wd_current, sep=""))
mask <- raster("bio_1.tif")
mask <- crop(mask, extent(-60, -30, -19.5, -14)) # cortando para a faixa desejada
mask[!is.na(mask)] <- 0 # mudando todos os valores para 0
plot(mask, axes = T, box = T, legend = T, main = "Neotropics")

rasters2<-list() # criando lista vazia para empilhar os rasters
matrixPD2 <- matrix(0, ncol = 1, nrow = 95040) # criando matriz vazia para fazer matrixPD
matrixPD2 <- data.frame(matrixPD2)

for (i in 1:length(modelos)) { # criando matrixPD e rasterizando poligonos dos modelos para empilhar
  Sys.time() -> start_time
  r1 <- rasterize(modelos[[i]], mask, field=1)
  r1[is.na(r1)] <- 0
  resample(r1, mask) -> r1
  mask(r1, mask) -> rasters2[[i]]
  meh <- as.data.frame(values(rasters2[[i]]))
  tip <- sub(" ", "_", labs[i])
  colnames(meh) <- tip
  matrixPD2 <- cbind(matrixPD2, meh)
  print(c(i, "in", length(modelos)))
  Sys.time() -> end_time  
  print(end_time-start_time)
}

names(rasters2) <- names(modelos)
matrixPD2[is.na(matrixPD2)] <- 0 # substituindo os NAs por 0s
matrixPD2 <- matrixPD2[,2:614]


# Fazendo matrixPD e arvore conversarem #
tree2 <- drop.tip(tree, setdiff(tree$tip.label, colnames(matrixPD2))) # cortando a arvore para as colunas da matrixPD (caso precise)
matrixPD2 <- matrixPD2[,c(which(colnames(matrixPD2) %in% tree$tip.label))] # cortando a matrixPD de acordo com os terminais da arvore (caso precise)

setwd("Results")
write.csv(matrixPD2, file="matrixPD2.csv")

# Rodando analise de PD #
pd_final2 <- pd(matrixPD2, tree2, include.root=TRUE) # dá um valor de PD por célula (SR é o species richness?)


# Empilhando os rasters e criando mapa de SR#
rasters2 <- rasters2[which(names(rasters2) %in% tree2$tip.label)]
sr2 <- calc(stack(rasters2), sum)
plot(sr2, main = "Current SR")

# criando mapas de PD #
pdMAP <- mask
values(pdMAP) <- pd_final2[,1]

plot(pdMAP, main = "PD2")
plot(wrld_simpl, add = T)

writeRaster(pdMAP, "PD2.tif")
writeRaster(sr2, "SR2.tif")

# which cells have higher or lower PD values as expected by their SR?
plot(pd_final2[,2], pd_final2[,1], xlab = "species richness", ylab = "pd values")

# using linear regression (not sure if the best method...)
linearRegression <- lm(pd_final2[,1]~pd_final2[,2])
abline(linearRegression, col="red") # adding trend line
summary(linearRegression)
res <- as.numeric(linearRegression$residuals)
resMap2 <- mask
values(resMap2) <- res
plot(resMap2, main="residual") # positive values = more pd than expected by species richness/ negative values = less pd than expected by species richness
plot(wrld_simpl, add = T)
writeRaster(resMap2, "Res2.tif")

###############################################
### PARTE 3: segundo PD - até latitude -25 ###
###############################################

#criando o raster vazio para os modelos binários
setwd(paste(wd, wd_current, sep=""))
mask <- raster("bio_1.tif")
mask <- crop(mask, extent(-60, -30, -25, -19.5)) # cortando para a faixa desejada
mask[!is.na(mask)] <- 0 # mudando todos os valores para 0
plot(mask, axes = T, box = T, legend = T, main = "Neotropics")

rasters3<-list() # criando lista vazia para empilhar os rasters
matrixPD3 <- matrix(0, ncol = 1, nrow = 95040) # criando matriz vazia para fazer matrixPD
matrixPD3 <- data.frame(matrixPD3)

for (i in 1:length(modelos)) { # criando matrixPD e rasterizando poligonos dos modelos para empilhar
  Sys.time() -> start_time
  r1 <- rasterize(modelos[[i]], mask, field=1)
  r1[is.na(r1)] <- 0
  resample(r1, mask) -> r1
  mask(r1, mask) -> rasters3[[i]]
  meh <- as.data.frame(values(rasters3[[i]]))
  tip <- sub(" ", "_", labs[i])
  colnames(meh) <- tip
  matrixPD3 <- cbind(matrixPD3, meh)
  print(c(i, "in", length(modelos)))
  Sys.time() -> end_time  
  print(end_time-start_time)
}

names(rasters3) <- names(modelos)
matrixPD3[is.na(matrixPD3)] <- 0 # substituindo os NAs por 0s
matrixPD3 <- matrixPD3[,2:614]


# Fazendo matrixPD e arvore conversarem #
tree3 <- drop.tip(tree, setdiff(tree$tip.label, colnames(matrixPD3))) # cortando a arvore para as colunas da matrixPD (caso precise)
matrixPD3 <- matrixPD3[,c(which(colnames(matrixPD3) %in% tree$tip.label))] # cortando a matrixPD de acordo com os terminais da arvore (caso precise)

setwd(wd)
setwd("Results")
write.csv(matrixPD3, file="matrixPD3.csv")

# Rodando analise de PD #
pd_final3 <- pd(matrixPD3, tree3, include.root=TRUE) # dá um valor de PD por célula (SR é o species richness?)


# Empilhando os rasters e criando mapa de SR#
rasters3 <- rasters3[which(names(rasters3) %in% tree$tip.label)]
sr3 <- calc(stack(rasters3), sum)
plot(sr3, main = "Current SR")

# criando mapas de PD #
pdMAP <- mask
values(pdMAP) <- pd_final3[,1]

plot(pdMAP, main = "PD2")
plot(wrld_simpl, add = T)

writeRaster(pdMAP, "PD3.tif")
writeRaster(sr3, "SR3.tif")

# which cells have higher or lower PD values as expected by their SR?
plot(pd_final3[,2], pd_final3[,1], xlab = "species richness", ylab = "pd values")

# using linear regression (not sure if the best method...)
linearRegression <- lm(pd_final3[,1]~pd_final3[,2])
abline(linearRegression, col="red") # adding trend line
summary(linearRegression)
res <- as.numeric(linearRegression$residuals)
resMap3 <- mask
values(resMap3) <- res
plot(resMap3, main="residual") # positive values = more pd than expected by species richness/ negative values = less pd than expected by species richness
plot(wrld_simpl, add = T)
writeRaster(resMap3, "Res3.tif")




#### Phylogenetic diversity maps ####

#rm(list=ls()) # uncomment to clean your environment

library(raster)
library(dplyr)
library(maptools)
library(sf)
library(rgdal)
library(picante)
data(wrld_simpl) # puxa dados do shapefile do world map

# setting wd #
#wd <- "C:/Users/raque/OneDrive/Documentos/EDGE_MANUAL 3/pd_maps/"
#setwd(wd)

# Carregando a arvore #
tree <- read.tree("final_tree.txt") # ja com a C. multinervia inserida
  #no meu caso (tita) essa arvore seria a arvore de generos?

# Mudando diretorio para puxar os shapes #
wd_polygon <- "shape_files/"

# Loading polygons #
# Vamos carregar os poligonos dos modelos #
setwd(paste(wd, wd_polygon, sep=""))
list.files(pattern="_current.shp") -> files_cur
sapply(files_cur, readShapePoly, simplify = F) -> modelos
sub("_current.shp", "", files_cur) -> labs
labs <- sub(" ","_", labs)
names(modelos) <- labs

# Agora vamos carregar um raster "vazio" do leste do Brasil pra projetar os valores
setwd(wd)
mask <- raster("bio_1.tif")
mask <- crop(mask, extent(-60, -20, -25, 0)) # cortando para o leste do Brasil
mask[!is.na(mask)] <- 0 # mudando todos os valores para 0
plot(mask, axes = T, box = T, legend = T, main = "Neotropics")

##
rasters1<-list() # criando lista vazia para empilhar os rasters
matrixPD <- matrix(0, ncol = 1, nrow = 36000) # criando matriz vazia para fazer matrixPD
matrixPD <- data.frame(matrixPD)

for (i in 1:length(modelos)) { # criando matruxPD e rasterizando poligonos dos modelos para empilhar
  Sys.time() -> start_time
  r1 <- rasterize(modelos[[i]], mask, field=1)
  r1[is.na(r1)] <- 0
  resample(r1, mask) -> r1
  mask(r1, mask) -> rasters1[[i]]
  meh <- as.data.frame(values(rasters1[[i]]))
  tip <- sub(" ", "_", labs[i])
  colnames(meh) <- tip
  matrixPD <- cbind(matrixPD, meh)
  print(c(i, "in", length(modelos)))
  Sys.time() -> end_time  
  print(end_time-start_time)
}
matrixPD[is.na(matrixPD)] <- 0 # substituindo os NAs por 0s

# Empilhando os rasters #
names(rasters1) <- names(modelos)
rasters1 <- rasters1[which(names(rasters1) %in% tree$tip.label)]
res1 <- Reduce("+", rasters1, accumulate = TRUE)
l1<-length(res1)
speciesRichness <- res1[[l1]] # criando um raster de species richness a partir dos modelos binarizados empilhados
plot(speciesRichness, main = "Current SR")

# Fazendo matrixPD e arvore conversarem #
tree <- drop.tip(tree, setdiff(tree$tip.label, colnames(matrixPD))) # cortando a arvore para as colunas da matrixPD (caso precise)
matrixPD <- matrixPD[,c(which(colnames(matrixPD) %in% tree$tip.label))] # cortando a matrixPD de acordo com os terminais da arvore (caso precise)

# Rodando analise de PD #
pd_final <- pd(matrixPD, tree, include.root=TRUE) # daí dá um valor de PD por célula (SR é o species richness?)

# Comparando resultados
pdMAP <- mask
values(pdMAP) <- pd_final[,1]
srMAPfromPD <- mask
values(srMAPfromPD) <- pd_final[,2]

par(mfrow=c(1,3))
plot(speciesRichness, main = "SR")
plot(wrld_simpl, add = T)
plot(pdMAP, main = "PD")
plot(wrld_simpl, add = T)
plot(srMAPfromPD, main = "SR from PD")
plot(wrld_simpl, add = T) # checar se primeiro e terceiro mapa sao iguais (devem ser)
dev.off()
writeRaster(pdMAP, "PD.tif")
writeRaster(speciesRichness, "SR.tif")
writeRaster(srMAPfromPD, "SR_PD.tif")

# which cells have higher or lower PD values as expected by their SR?
plot(pd_final[,2], pd_final[,1], xlab = "species richness", ylab = "pd values")

# using linear regression (not sure if the best method...)
linearRegression <- lm(pd_final[,1]~pd_final[,2])
abline(linearRegression, col="red") # adding trend line
summary(linearRegression)
res <- as.numeric(linearRegression$residuals)
resMap <- mask
values(resMap) <- res
plot(resMap, main="residual") # positive values = more pd than expected by species richness/ negative values = less pd than expected by species richness
plot(wrld_simpl, add = T)
writeRaster(resMap, "Res.tif")

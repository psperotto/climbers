
  ## distribution of neotropical climbers ##
setwd("C:/Users/patri/Google Drive/MESTRADO/MACROEVO/SSE") # novo laptop
setwd("C:/Users/Patricia/Google Drive/MESTRADO/MACROEVO/SSE")
getwd()

library(phytools)
library(geiger)
library(ape)
library(magrittr)
library(diversitree)
library(dismo)
library(hisse)
library(CoordinateCleaner)
library(monographaR)
library(dplyr)
library(viridis)
library(schoolmath)
library(RColorBrewer)
library(taxize)
data(wrld_simpl)

?saveRDS() #-> funcao pra salvar objetos grandes no directory
?readRDS() #-> funcao pra chamar o objeto salvo

# script pra pegar spp de cada CM quando tiver a lista pronta:
 g <- read.csv("climbers dez 2019.csv")
 # voluvel (1)
 f <- as.character(subset(g, CM5==1)[,3])
 climbers1 <- sub("_", " ", f)
 # gavinha (2)
 a <- as.character(subset(g, CM5==2)[,3])
 climbers2 <- sub("_", " ", a)
 # apoiante (3)
 b <- as.character(subset(g, CM5==3)[,3])
 climbers3 <- sub("_", " ", b)
 # raizes (4)
 c <- as.character(subset(g, CM5==4)[,3])
 climbers4 <- sub("_", " ", c)
 # ramos (5)
 d <- as.character(subset(g, CM5==5)[,3])
 climbers5 <- sub("_", " ", d)
 # petioles (6)
 e <- as.character(subset(g, CM5==6)[,3])
 climbers6 <- sub("_", " ", e)
 # hooks (6)
 h <- as.character(subset(g, CM5==7)[,3])
 climbers7 <- sub("_", " ", h)
 
# script pra pegar todas spp quando tiver a lista pronta:
# g <- read.csv("trepadeiras outubro 2019.csv")
# f <- as.character(g$Species)
# allclimbers <- sub("_", " ", f)
# write.csv(subset(g, CM==1)[,1], file="twiners.csv")

dd<-c("Pacouria boliviensis","Pacouria guianensis","Pacouria paraensis") 
 
########## Load distribution points into dataframe ##########
 
# "for(i in xx:length(lista){..." ("xx" numero daonde começar)
w<-as.data.frame(climbers1)
q<- w$climbers1==  "Fallopia convolvulus" 
r<- as.data.frame(q[T])
    # essas 3 linhas sao pra quando der pau na busca no GBIF e eu 
    # descorbir o numero da spp em q parou a busca pra ajustar o loop

#library(dismo) # pacote para acessar gbif => tem um limite de registros 
                                  # q ele consegue pegar por aqui (200.000)
datalist = list()
for (z in 4407:length(climbers1)){
  tt <- strsplit(climbers1," ",fixed = TRUE)[[z]]
  tt[1] -> x
  tt[2] -> y
  print(c(x,y)) # to check names
  dat <- gbif(x, y) # download points from gbif
  datalist[[z]] <- dat 
}
big_data <- data.table::rbindlist(datalist, fill=T)
big_data <- as.data.frame(big_data[,c("species","lon","lat")])
sub(" ", "_", big_data$species) -> names
big_data[,1] <- names 
x6<-big_data

  # salvando blocos de ocorrencias como objetos
saveRDS(big_data,file = "bigdata_roots.Rdata") 
  
  # lendo os blocos
roots <- readRDS("bigdata_roots.Rdata")
saveRDS(roots, file="roots.Rdata")
hooks <- readRDS("bigdata_hooks.Rdata")

twiners<-readRDS("twiners.Rdata") 
tendrils<-readRDS("tendrils.Rdata") 
roots<-readRDS("roots.Rdata")
scramblers<-readRDS("scramblers.Rdata")
petioles<-readRDS("petioles.Rdata")
hooks<-readRDS("hooks.Rdata")
branches<-readRDS("branches.Rdata")

allclimbers<-rbind(twiners,tendrils,scramblers,petioles,branches,roots,hooks)

  # juntando os blocos
twiners <- rbind(x1, x2, x3, x4, x5, x6)

  # renaming columns in a dataframe
hooks <- hooks %>% 
  rename(
    Species = species,
    Longitude = lon,
    Latitude = lat
  )

  # salvando o 'raw file' de ocorrencia pra cada mecanismo 
#saveRDS(twiners,file="twiners2.Rdata")
#twiners<-readRDS("twiners2.Rdata")

########## Data cleaning ##########

#library(CoordinateCleaner) # pacote que limpa coordenadas
points <- subset(allclimbers, !is.na(Longitude) & !is.na(Latitude)) 
                              # talvez criar aqui o subset com numero de 
                              # coletor tb p/ evitar exclusão de coletas 
                              # diversas no mesmo lugar usando o cc_dup
cleaned_points <- cc_cap(points, lon="Longitude", lat="Latitude", species ="Species")
cleaned_points <- cc_cen(cleaned_points, lon = "Longitude", lat = "Latitude", species ="Species") 
cleaned_points <- cc_dupl(cleaned_points, lon="Longitude", lat="Latitude", species ="Species")
cleaned_points <- cc_equ(cleaned_points, lon="Longitude", lat="Latitude")
cleaned_points <- cc_inst(cleaned_points, lon="Longitude", lat="Latitude", species ="Species")
cleaned_points <- cc_val(cleaned_points, lon="Longitude", lat="Latitude")
cleaned_points2 <- cc_sea(cleaned_points, lon="Longitude", lat="Latitude")


# Filtering centroids #
# centroids for: 
# Brazil, Bolivia, Peru, Colombia, Venezuela, Chile... 
# (should probably add others..)
#loncent <- c(-52.87310,-64.435,-74.14,-72.8667,-65.9119,-70.8906,
#             -64.9208, -79.016069, -78.7520357, -58.1689,	-56.0180680,
 #            -55.625, -52.9708, -58.7017,	-83.9519, -90.1792, -89.35083,
  #           -89.93333)
#latcent <- c(-10.833900,-16.7261,-9.1839,3.8811,7.0758,-35.8156,
#             -35.3869, 21.6229002, -1.4238195, -23.2025, -32.7995198,
 #            4.1006, 3.8558, 4.7364, 10.0114, 15.7422, 16.08111,
  #           13.81667)

# flaging and removing centroids
remove_centroid <- c()
for(i in 1:length(cleaned_points2$Species)) {
  for(u in 1:length(latcent))
    if(cleaned_points2[i,2] == loncent[u] && cleaned_points2[i,3] == latcent[u]){
      print(paste(i, "is centroid!", sep=" "))
      remove_centroid <- c(remove_centroid, i)
    }
}
hooks <-cleaned_points2[-remove_centroid, ]

#write.csv(tendrils, file="tendrils.csv")
#read.csv("hooks.csv")->hooks

# 3.4. Removing inaccurate points #
# removing points without decimal cases #
keep_accurate <- c()
for(i in 1:length(hooks$Latitude)) {
  if(is.decimal(hooks[i,2]) && is.decimal(hooks[i,3])){
    print(paste(i, "is accurate!", sep=" "))
    keep_accurate <- c(keep_accurate, i)
  }
}
hooks <- hooks[keep_accurate, ]

# 3.5. Removing outliers in distribution #
names(table(hooks$Species))->species
cleaned_point <- data.frame()
for(u in 1:length(species)){
  sp0 <- hooks[hooks$Species==species[u],]
  out_lat <- boxplot.stats(sp0$Latitude)$out
  out_lon <- boxplot.stats(sp0$Longitude)$out
  sp <- sp0[ ! sp0$Latitude %in% out_lat, ]
  sp <- sp[ ! sp$Longitude %in% out_lon, ]
  cleaned_point <- rbind(cleaned_point, sp)
}


########## restricting to the Neotropics ########## 
#library(dplyr)
allclimbers %>% filter(Latitude> -14)->allclimbers2
allclimbers2 %>% filter(Latitude< -9)->allclimbers3
allclimbers3 %>% filter(Longitude> -44)->allclimbers4
allclimbers4 %>% filter(Longitude< -41)->allclimbers_chapada
saveRDS(allclimbers_chapada, file="allclimbers_chapada.Rdata")
write.csv(allclimbers, file = "allclimbers.csv")


twiners %>% filter(Latitude> -14)->twiners2
twiners2 %>% filter(Latitude< -9)->twiners3
twiners3 %>% filter(Longitude> -43)->twiners4
twiners4 %>% filter(Longitude< -40)->twiners_chapada
saveRDS(twiners_chapada, file="twiners_chapada2.Rdata")
write.csv(twiners_chapada, file = "twiners_chapada2.csv")

    ## tabela com pontos - pode ir como material suplementar, etc ##
#write.csv(cleaned_points2, file="cleaned_points2.csv") 

########## checar a distribuicao ##########
#library(rgeos)
#library(monographaR)
#plot(cleaned_points2$lon, cleaned_points2$lat, asp=1, add=T)
nameCol <- c("sp", "Longitude","Latitude")
colnames(allclimbers_chapada) <- nameCol
mapallclimbers_chapada<-
  mapDiversity(allclimbers, resolution = 1, plot = TRUE, 
                 plot.with.grid = F, legend = T, 
                 export = F, alpha=1,
                 col=hcl.colors(50, palette = "inferno",rev=T))
#saveRDS(maphooks, file="maphooks.Rdata")

#tiff("teste.tiff", width = 10, height = 16, units = 'in', res = 300)
pdf(file="mapdiversity.pdf", width=10, height=16)
#jpeg(file="mapdiversity.jpeg", width = 10, height = 16, units = 'in', res = 300)
par(mfrow=c(4,2))

plot(crop(mapallclimbers, extent(-130, -20, -30, 30)),
     main = "All climbers (n=9822)", zlim=c(1,850), cex.main=2,cex.axis=2,
     col=hcl.colors(50, palette = "inferno",rev=T))
plot(wrld_simpl, add=T,  border="black", lwd=0.3)

plot(crop(maptwiners, extent(-130, -20, -30, 30)),
     main = "Twining (n=4680)", zlim=c(1,300), cex.main=2,cex.axis=2,
     col=hcl.colors(50, palette = "inferno",rev=T))
plot(wrld_simpl, add=T,  border="black", lwd=0.3)

plot(crop(maptendrils, extent(-130, -20, -30, 30)),
     main = "Tendrils (n=2126)", zlim=c(1,300), cex.main=2,cex.axis=2,
     col=hcl.colors(50, palette = "inferno",rev=T))
plot(wrld_simpl, add=T,  border="black", lwd=0.3)

plot(crop(mapscramblers, extent(-130, -20, -30, 30)),
     main = "Simple scrambling (n=1404)", zlim=c(1,300), cex.main=2,cex.axis=2,
     col=hcl.colors(50, palette = "inferno",rev=T))
plot(wrld_simpl, add=T,  border="black", lwd=0.3)

plot(crop(maproots, extent(-130, -20, -30, 30)),
     main = "Adhesive roots (n=1124)", zlim=c(1,300), cex.main=2,cex.axis=2,
     col=hcl.colors(50, palette = "inferno",rev=T))
plot(wrld_simpl, add=T,  border="black", lwd=0.3)

plot(crop(mapbranches, extent(-130, -20, -30, 30)),
     main = "Prehensile branches (n=349)", zlim=c(1,300), cex.main=2,cex.axis=2,
     col=hcl.colors(50, palette = "inferno",rev=T))
plot(wrld_simpl, add=T,  border="black", lwd=0.3)

plot(crop(mappetioles, extent(-130, -20, -30, 30)),
     main = "Prehensile petioles (n=105)", zlim=c(1,300), cex.main=2,cex.axis=2,
     col=hcl.colors(50, palette = "inferno",rev=T))
plot(wrld_simpl, add=T,  border="black", lwd=0.3)

plot(crop(maphooks, extent(-130, -20, -30, 30)),
     main = "Hooks or grapnels (n=31)", zlim=c(1,300), cex.main=2,cex.axis=2,
     col=hcl.colors(50, palette = "inferno",rev=T))
plot(wrld_simpl, add=T,  border="black", lwd=0.3)

dev.off()
                 
#title(main="Hooks or grapnels")
# salvar pdf tamanho 10x6

??plot()
  
?mapDiversity()




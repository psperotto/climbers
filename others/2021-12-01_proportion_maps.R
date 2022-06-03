setwd("~/Desktop/Colabs/Patricia_Climbers/climbers")

library(data.table)
library(monographaR)
library(tidyverse)
library(raster)
library(CoordinateCleaner)
library(ape)
library(maptools)
data("wrld_simpl")

source("neotropical_climbers_functions.R") 
# source das funcoes q vamos usar

#rm(list=ls())
# First get grid cell values from filtered GBIF data 
#### from TV local
#gbif_dir <- paste0(getwd(), "/full_gbif_quest")
#gbif_files <- list.files(paste0(gbif_dir, "/z_filtered_gbif"), ".csv")
#full_list <- read.full.gbif(gbif_files, gbif_dir)
#write.csv(full_list, file=paste0(gbif_dir, "full_tracheophyte_filtered_gbif.csv"))

# load table back
#full_list <- fread(paste0(gbif_dir, "/full_tracheophyte_filtered_gbif.csv"))
#full_list <- full_list[,-1]
#save.gbif.neotropics(full_list)

# Load master table
full_list <- fread("full_gbif_quest/master_table.csv")

##############
# Directory for descriptive results
descriptive_dir <- paste0(getwd(), "/1_descriptive")

#full_list <- fread(paste0(gbif_dir, "/neotropics_tracheophyte_filtered_gbif.csv"))
#full_list <- full_list[,-1]
full_map <- run.mapDiversity.neotropics(full_list[,-4], filename="full_neotropical_diversity", dir=descriptive_dir) # function modified from monographaR


#-------
all_climbers <- subset(full_list, full_list$mechanism != "not_a_climber")
#write.csv(as.data.frame(all_climbers[,1:3]), file="climbers_points.csv", row.names=F)

all_climbers_map <- run.mapDiversity.neotropics(all_climbers[,-4], filename="all_climbers_neotropical_diversity", dir=descriptive_dir)


proportion <- (all_climbers_map /full_map)
proportion[which(proportion[]>0.25)] <- NA 
proportions_coordinates <- rasterToPoints(proportion)
proportions_coordinates <- as.data.frame(proportions_coordinates)

bio12 <- raster("climate_layers/current_30sec/bio_12.tif")
bio12 <- crop(bio12, extent(proportion))
bio12 <- rasterToPoints(bio12)
bio12 <- as.data.frame(bio12)
bio12[,1:2] <- round(bio12[,1:2],1)

proportions_coordinates$merge <- paste0(proportions_coordinates[,1], "_",proportions_coordinates[,2])
bio12$merge <- paste0(bio12[,1], "_",bio12[,2])

merged_table <- merge(bio12, proportions_coordinates, by="merge")
merged_table <- merged_table[,c("x.x","y.x","bio_12","layer")]
colnames(merged_table) <- c("x","y","bio_12","proportion")

plot(merged_table$bio_12, merged_table$proportion)


points <- SpatialPoints(merged_table[,1:2])
points <- SpatialPointsDataFrame(points, merged_table[,3:4])

plot(points$proportion)

ll = "+proj=longlat +datum=WGS84" 
proj4string(points) = CRS(ll)

require(spdep)
hexnb = poly2nb(points)



slot(points, "proj4string") 

lcc = "+proj=lcc +lat_1=24 +lat_2=-24 +lon_0=-70" 
require(rgdal)
W.sdf = spTransform(points, CRS(ll))
bbox(W.sdf)

hpt = spsample(W.sdf, type="hexagonal", n=250, bb=bbox(W.sdf) * 1.2, offset=c(1, -1))
hpg = HexPoints2SpatialPolygons(hpt)
hexid = over(x=W.sdf, y=hpg)
hexid <- subset(hexid, !is.na(hexid))
hpg = hpg[unique(hexid)]
int = over(x=hpg, y=W.sdf, fn=max) 
colnames(int) = c("WmaxS")
head(int)
hspdf = SpatialPolygonsDataFrame(hpg, int, match.ID = TRUE)
plot(hspdf)
names(hspdf) <- c("climate","proportion")
plot(W.sdf, pch=20, cex=.3, add=TRUE)

spplot(hspdf, "climate")


library(spdep)
queen.nb = poly2nb(hspdf)
listw1 = nb2listw(queen.nb)
reglm <- lm(climate~proportion, data=hspdf)
summary(reglm)
lm.morantest(reglm, listw1)



W.sdf@data = data.frame(num=rep(1, ch))

head(slot(hspdf, "data"))
require(maps)
require(maptools)
cl = map("world", xlim=c(-120, 20), ylim=c(-10, 70), plot=FALSE)
clp = map2SpatialLines(cl, proj4string=CRS(ll))
clp = spTransform(clp, CRS(ll))
l2 = list("sp.lines", clp, col="gray")
require(colorRamps)
cr = blue2yellow(20)

spplot(hspdf, "proportion", col="white", col.regions=blue2yellow(20),sp.layout=list(ll, l2),
       colorkey=list(space="bottom"), sub="climate")

require(spdep)
hexnb = poly2nb(hspdf)
wts = nb2listw(hexnb, style="W")
summary(wts)
m = length(hspdf$proportion)
s = Szero(wts)
moran(hspdf$proportion, wts, n=m, S0=s)

install.packages("spgwr")
require(spgwr)

bw = gwr.sel(proportion ~ climate, data=hspdf)
bw * .001
model = gwr(proportion ~ climate, data=hspdf, bandwidth=bw)


#
#CRS arguments:
#  +proj=longlat +datum=WGS84 +ellps=WGS84
#+towgs84=0,0,0


(Sr1 = Polygon(points[,1:2]))
(Srs1 = Polygons(list(Sr1), "s1"))
?SpatialPolygons
(SpP = SpatialPolygons(list(Srs1), 1:1, proj4string= crs("+proj=utm +zone=23 +south +datum=WGS84 +units=m +no_defs"))) 
plot(SpP, col = 3:3, pbg="white") 
SpP ### can not write as shapefile

### Convert the SpatialPolygons to SpatialPolygonsDataFrame
shape_pol <- SpatialPolygonsDataFrame(SpP, match.ID=F, data= data.frame(x=spdf[1:1,1], y=spdf[1:1,2]))
shape_pol ### can be write as shapefile
plot(shape_pol, col = 4)

#----------------------
library(sp)
library(raster)

### Example data: creating a SpatialPointsDataFrame object
x = c(1,2,5,4,3)
y = c(3,2,3,6,6)
df_points <- as.data.frame(cbind(x,y))
S <- SpatialPoints(cbind(x,y))
# S <- SpatialPoints(list(x,y))
# S <- SpatialPoints(data.frame(x,y))
S
plot(S)
spdf <- SpatialPointsDataFrame(S, df_points)
spdf
plot(spdf)
# crs(spdf) <- ("+proj=utm +zone=23 +south +datum=WGS84 +units=m +no_defs") ### add a crs

### Convert the SpatialPointsDataFrame to SpatialPolygons
(Sr1 = Polygon(spdf[,1:2]))
(Srs1 = Polygons(list(Sr1), "s1"))
(SpP = SpatialPolygons(list(Srs1), 1:1, proj4string= crs("+proj=utm +zone=23 +south +datum=WGS84 +units=m +no_defs"))) 
plot(SpP, col = 3:3, pbg="white", add=T) 
SpP ### can not write as shapefile

### Convert the SpatialPolygons to SpatialPolygonsDataFrame
shape_pol <- SpatialPolygonsDataFrame(SpP, match.ID=F, data= data.frame(x=spdf[1:1,1], y=spdf[1:1,2]))
shape_pol ### can be write as shapefile
plot(shape_pol, col = 4, add=T)

### write shapefile
library(rgdal)
writeOGR(shape_pol, paste0(getwd(), "/Output_shapes"), "p_to_shape_pol", driver="ESRI Shapefile")


# this code just creates a sample SpatialPointsDataFrame 
x <- c(-10,-10,10,10,-10)
y <- c(-10,10,10,-10,-10)
df.1 <- data.frame(x,y,id=47, order=1:5,hole=F,piece=1,group=47.1,box_id=1)
coordinates(df.1)=c("x","y")
x <- c(+15+3*cos(2*pi/5*(0:5)))
y <- c(-15+3*sin(2*pi/5*(0:5)))
df.2 <- data.frame(x,y,id=48, order=1:6,hole=F,piece=1,group=48.1,box_id=2)
coordinates(df.2)=c("x","y")
x <- c(-15+2*cos(2*pi/8*(0:8)))
y <- c(+15+2*sin(2*pi/8*(0:8)))
df.3.1 <- data.frame(x,y,id=20, order=1:9,hole=F,piece=1,group=20.1,box_id=3)
coordinates(df.3.1)=c("x","y")
x <- c(0+2*cos(2*pi/8*(0:8)))
y <- c(+15+2*sin(2*pi/8*(0:8)))
df.3.2 <- data.frame(x,y,id=20, order=1:9,hole=F,piece=1,group=20.2,box_id=3)
coordinates(df.3.2)=c("x","y")
x <- c(+15+2*cos(2*pi/8*(0:8)))
y <- c(+15+2*sin(2*pi/8*(0:8)))
df.3.3 <- data.frame(x,y,id=20, order=1:9,hole=F,piece=1,group=20.3,box_id=3)
coordinates(df.3.3)=c("x","y")
df <- rbind(df.1,df.2,df.3.1,df.3.2,df.3.3)
df
#              coordinates id order  hole piece group box_id
# 1             (-10, -10) 47     1 FALSE     1  47.1      1
# 2              (-10, 10) 47     2 FALSE     1  47.1      1
# 3               (10, 10) 47     3 FALSE     1  47.1      1
# 4              (10, -10) 47     4 FALSE     1  47.1      1
# 5             (-10, -10) 47     5 FALSE     1  47.1      1
# 6              (18, -15) 48     1 FALSE     1  48.1      2
# 7  (15.92705, -12.14683) 48     2 FALSE     1  48.1      2
# 8  (12.57295, -13.23664) 48     3 FALSE     1  48.1      2
# ...
data <- data.frame(box_id=unique(df$box_id),row.names=unique(df$id))

points2polygons <- function(df,data) {
  get.grpPoly <- function(group,ID,df) {
    Polygon(coordinates(df[df$id==ID & df$group==group,]))
  }
  get.spPoly  <- function(ID,df) {
    Polygons(lapply(unique(df[df$id==ID,]$group),get.grpPoly,ID,df),ID)
  }
  spPolygons  <- SpatialPolygons(lapply(unique(df$id),get.spPoly,df))
  SpatialPolygonsDataFrame(spPolygons,match.ID=T,data=data)
}

spDF <- points2polygons(df,data)
plot(spDF,col=spDF$box_id+1)


library(rgdal)
writeOGR(spDF,dsn=".",layer="myShapefile", driver="ESRI Shapefile")


class(points)

poly2nb(points)

plot(merged_table$bio_12, merged_table$proportion, pch=19, cex=0.1)

model <- lm(merged_table$bio_12~merged_table$proportion)
lm_col <- lm.morantest(model, list)

?nb2listw

?lm.morantest

install.packages("spatialreg")
library(spatialreg)
?lagsarlm


# Diversification analyses
setwd("~/Desktop/WCVP_special_issue/Patricia_Climbers/climbers")
#rm(list=ls())

source("00_utility_functions.R")

#--------------------------------------
# add information on main biome
climber_points <- as.data.frame(fread("gbif_climbers/climbers_cleaned_points.csv"))
biomes_for_points <- localityToBiome(climber_points, lat="decimalLatitude",lon="decimalLongitude")
biomes_for_genera <- getBiomes(biomes_for_points, species="genus") #we want summaries for each genus
biomes_for_genera <- as.data.frame(biomes_for_genera)

closed_canopy <- c("Tropical & Subtropical Moist Broadleaf Forests","Tropical & Subtropical Dry Broadleaf Forests","Tropical & Subtropical Coniferous Forests","Mangroves","Mediterranean Forests, Woodlands & Scrub")
open_canopy <- c("Tropical & Subtropical Grasslands, Savannas & Shrubland","Deserts & Xeric Shrublands",
                 "Montane Grasslands & Shrublands","Deserts & Xeric Shrublands","Flooded Grasslands & Savannas")

summary_biome <- data.frame(genera=row.names(biomes_for_genera),habitat=NA)
for(genus_index in 1:nrow(biomes_for_genera)) {
  one_genus <- rownames(biomes_for_genera)[genus_index]
  subset_biomes <- biomes_for_genera[genus_index,]
  closed_canopy_points <- sum(subset_biomes[which(names(subset_biomes)%in%closed_canopy)])
  open_canopy_points <- sum(subset_biomes[which(names(subset_biomes)%in%open_canopy)])
  summary_biome[genus_index,2] <- ifelse(closed_canopy_points>=open_canopy_points, "closed_canopy","open_canopy")
}

#--------------------------------------
climber_clades <- readRDS("bg.clades.list.Rdata")

for (i in 1:length(climber_clades)) {
  test <- climber_clades[[i]]
  
  # tendrils against all others
  CMs <- as.numeric(test$CM[-1])
  CMs[which(CMs==1)] <- "Other"
  CMs[which(CMs==2)] <- "Tendrils"
  CMs[which(CMs==3)] <- "Other"
  CMs[which(CMs==4)] <- "Other"
  CMs[which(CMs==5)] <- "Other"
  CMs[which(CMs==6)] <- "Other"
  CMs[which(CMs==7)] <- "Other"
  CMs[which(CMs==8)] <- "Other"
  
  test$CM <- c(NA, CMs)
  
  #example <- get.template(save.file=F)
  
  ######## PLOT FIG 3 & 4 ######
  test$node <- c("bg",rep("CG",length(test$node)-1))
  test <- subset(test, test$diversity>1)
  
  if(nrow(test)>1) {
    results <- get.tail.probs(table=test)
    
    k = test$diversity[which(test$taxa=="bg_clade")]
    time = test$age_mean[which(test$taxa=="bg_clade")]
    
    tail.probs=results
    crown.results <- tail.probs[[2]]
    k.lower.yule.crown <- crown.results[[1]]
    k.upper.yule.crown <- crown.results[[2]]
    k.lower.high.ext.crown <- crown.results[[3]]
    k.upper.high.ext.crown <- crown.results[[4]]
    
    max_x <- head(sort(test$age_mean, decreasing=T))[2]
    max_y <- head(sort(test$diversity, decreasing=T))[2]
      
    pdf(file=paste0(test$clade[2],".pdf"))
    #dev.off()
    plot(k.lower.yule.crown$time, log10(k.lower.yule.crown$k), type="l",
         col="gray", lwd=1.5, ylim = c(0, round(log10(max_y), 0)+1), xlim = c(0,max_x+max_x*0.1), xlab="", ylab="", yaxt = "n")
    par(new=TRUE)
    
    #plot(k.lower.yule.crown$time, log10(k.lower.yule.crown$k), type="l",
    #     col="gray", lwd=1.5, ylim = c(0, round(log10(k), 0)+1), xlim = c(0,time+time*0.1), xlab="", ylab="", yaxt = "n")
    #par(new=TRUE)
    
    plot(k.upper.yule.crown$time, log10(k.upper.yule.crown$k), type="l",
         col="gray", lwd=1.5, ylim = c(0, round(log10(max_y), 0)+1), xlim = c(0,max_x+max_x*0.1), xlab="", ylab="", yaxt = "n")
    par(new=TRUE)
    plot(k.lower.high.ext.crown$time, log10(k.lower.high.ext.crown$k), type="l",lty=3,
         col="gray",lwd=1.5, ylim = c(0, round(log10(max_y), 0)+1), xlim = c(0,max_x+max_x*0.1), xlab="", ylab="",  yaxt = "n")
    par(new=TRUE)
    plot(k.upper.high.ext.crown$time, log10(k.upper.high.ext.crown$k), type="l",lty=3,
         col="gray",lwd=1.5, ylim = c(0, round(log10(max_y), 0)+1), xlim = c(0,max_x+max_x*0.1), xlab='Age of clade (million years)', ylab='K (expected number of species)',yaxt = "n")
    par(new=TRUE)
    
    # defining color vector
    mode <- test$CM[-1]
    colors_states <- c("midnightblue", "goldenrod")
    mode_cols <- mode
    mode_cols[mode_cols=="Tendrils"] <- "goldenrod"
    mode_cols[mode_cols=="Other"] <- "midnightblue"
    
    #add points for clades
    points(x=test$age_mean[test$taxa=="bg_clade"], y=log10(test$diversity[test$taxa=="bg_clade"]), pch=3, col="black", type="p")
    points(x=test$age_mean[test$taxa!="bg_clade" & test$node=="CG"], y=log10(test$diversity[test$taxa!="bg_clade" & test$node=="CG"]), pch=16, type="p", col=mode_cols)
    
    #axis.par <- c(0:round(log10(k), 0)+1)
    #axis(2, at=axis.par, labels=axis.par, las=1)
    #legend(1, tail(axis.par, 1), legend=c("eps = 0", "eps = 0.9"), col=c("gray", "gray"), lty=c(1,3), cex=1)
    title(main=table$clade[2])
    dev.off()
  }
}


# save(results, file="test.Rsave")

pdf(file="plot.tail.probs.pdf", height=6, width=8)
plottailprobs(test, results)
dev.off()

# rm(list=ls())
#####################################################################
#####################################################################
#####################################################################
#####################################################################
#####################################################################
#####################################################################
# Now with some families as background rate

setwd("~/Desktop/climbers")

#climber_clades_stem <- readRDS("REVIEW_bg.clades.stem.Rdata")
#all_climbers_stem <- climber_clades_stem[[7]]
#diverse_families <- sort(table (all_climbers_stem$clade), decreasing=T)[1:5]
#climber_families_stem <- subset(all_climbers_stem, all_climbers_stem$clade%in%names(diverse_families))
#climber_families_stem <- rbind(all_climbers_stem[1:5,],climber_families_stem)

#write.csv(climber_families_stem, file="climber_families_stem_review.csv", row.names=F)

all_climbers_stem <- read.csv("climber_families_stem_review_complete.csv")

# now let's exclude everything based on specific cutoffs:
climber_clades_stem <- readRDS("stem.bg.clades.list.Rdata")
root.age <- as.numeric(climber_clades_stem[[7]][1,"age_mean"])

lsca.df.sum1 <- all_climbers_stem[1:5,]
lsca.df.sum1 <- rbind(climber_clades_stem[[7]][1,],lsca.df.sum1)

cutoff <- c(0.75)
results_folder <- "/revision_results/revision_results_background_families/"
families <- c("Malpighiaceae","Bignoniaceae","Fabaceae","Cucurbitaceae","Apocynaceae")
pal <- hcl.colors(5, palette = "Viridis", alpha = 0.7)
names(pal) <- families

for(family_index in 1:length(families)) {
  
  #####################################################################
  # Stem node and high extinction
  
  # get diversification rates
  epsilon_high <- 0.9 # 
  
  r_all_angios <- bd.ms(time = as.numeric(lsca.df.sum1[1,'age_mean']), n = as.numeric(lsca.df.sum1[1,'diversity']), crown = FALSE, epsilon = epsilon_high)
  r_malpighiaceae <- bd.ms(time = as.numeric(lsca.df.sum1[2,'age_mean']), n = as.numeric(lsca.df.sum1[2,'diversity']), crown = FALSE, epsilon = epsilon_high)
  r_bignoniaceae <- bd.ms(time = as.numeric(lsca.df.sum1[3,'age_mean']), n = as.numeric(lsca.df.sum1[3,'diversity']), crown = FALSE, epsilon = epsilon_high)
  r_fabaceae <- bd.ms(time = as.numeric(lsca.df.sum1[4,'age_mean']), n = as.numeric(lsca.df.sum1[4,'diversity']), crown = FALSE, epsilon = epsilon_high)
  r_cucurbitaceae <- bd.ms(time = as.numeric(lsca.df.sum1[5,'age_mean']), n = as.numeric(lsca.df.sum1[5,'diversity']), crown = FALSE, epsilon = epsilon_high)
  r_apocynaceae <- bd.ms(time = as.numeric(lsca.df.sum1[6,'age_mean']), n = as.numeric(lsca.df.sum1[6,'diversity']), crown = FALSE, epsilon = epsilon_high)
  
  tot.div1 <- c(as.numeric(lsca.df.sum1[1,'diversity']),
                as.numeric(lsca.df.sum1[2,'diversity']),
                as.numeric(lsca.df.sum1[3,'diversity']),
                as.numeric(lsca.df.sum1[4,'diversity']),
                as.numeric(lsca.df.sum1[5,'diversity']),
                as.numeric(lsca.df.sum1[6,'diversity']))
  all_r <- c(r_all_angios, r_malpighiaceae, r_bignoniaceae, r_fabaceae, r_cucurbitaceae, r_apocynaceae)
  
  getStemLimits <- function(time, e=epsilon, r) {
    t(sapply(time, stem.limits, epsilon=e, r=r))
  }
  
  limits <- mclapply(all_r, getStemLimits, time=2:root.age, e=epsilon_high, mc.cores = detectCores()-1)
  limits.df_high <- data.frame(rbind(limits[[1]], limits[[2]],
                                     limits[[3]], limits[[4]],
                                     limits[[5]], limits[[6]]),
                               Diversity=factor(rep(tot.div1, each=nrow(limits[[1]]))),
                               r=rep(all_r, each=nrow(limits[[1]])))
  colnames(limits.df_high)[1:2] <- c("lb", "ub")
  
  #####
  
  epsilon_low <- 0 # 
  root.age <- as.numeric(lsca.df.sum1[1,'age_mean'])
  
  r_all_angios <- bd.ms(time = as.numeric(lsca.df.sum1[1,'age_mean']), n = as.numeric(lsca.df.sum1[1,'diversity']), crown = FALSE, epsilon = epsilon_low)
  r_malpighiaceae <- bd.ms(time = as.numeric(lsca.df.sum1[2,'age_mean']), n = as.numeric(lsca.df.sum1[2,'diversity']), crown = FALSE, epsilon = epsilon_low)
  r_bignoniaceae <- bd.ms(time = as.numeric(lsca.df.sum1[3,'age_mean']), n = as.numeric(lsca.df.sum1[3,'diversity']), crown = FALSE, epsilon = epsilon_low)
  r_fabaceae <- bd.ms(time = as.numeric(lsca.df.sum1[4,'age_mean']), n = as.numeric(lsca.df.sum1[4,'diversity']), crown = FALSE, epsilon = epsilon_low)
  r_cucurbitaceae <- bd.ms(time = as.numeric(lsca.df.sum1[5,'age_mean']), n = as.numeric(lsca.df.sum1[5,'diversity']), crown = FALSE, epsilon = epsilon_low)
  r_apocynaceae <- bd.ms(time = as.numeric(lsca.df.sum1[6,'age_mean']), n = as.numeric(lsca.df.sum1[6,'diversity']), crown = FALSE, epsilon = epsilon_low)
  
  tot.div1 <- c(as.numeric(lsca.df.sum1[1,'diversity']),
                as.numeric(lsca.df.sum1[2,'diversity']),
                as.numeric(lsca.df.sum1[3,'diversity']),
                as.numeric(lsca.df.sum1[4,'diversity']),
                as.numeric(lsca.df.sum1[5,'diversity']),
                as.numeric(lsca.df.sum1[6,'diversity']))
  all_r <- c(r_all_angios, r_malpighiaceae, r_bignoniaceae, r_fabaceae, r_cucurbitaceae, r_apocynaceae)
  
  getStemLimits <- function(time, e=epsilon, r) {
    t(sapply(time, stem.limits, epsilon=e, r=r))
  }
  
  limits <- mclapply(all_r, getStemLimits, time=2:root.age, e=epsilon_low, mc.cores = detectCores()-1)
  limits.df_low <- data.frame(rbind(limits[[1]], limits[[2]],
                                    limits[[3]], limits[[4]],
                                    limits[[5]], limits[[6]]),
                              Diversity=factor(rep(tot.div1, each=nrow(limits[[1]]))),
                              r=rep(all_r, each=nrow(limits[[1]])))
  colnames(limits.df_low)[1:2] <- c("lb", "ub")
  
  
  
  ################################################
  subset_stem <- subset(all_climbers_stem, !is.na(all_climbers_stem$clade))
  subset_stem_t <- subset(subset_stem, subset_stem$CM==2)
  subset_stem_other <- subset(subset_stem, subset_stem$CM!=2)
  
  #nrow(subset_stem_t)
  #nrow(subset_stem_other[subset_stem_other$clade=="Malpighiaceae",]) + nrow(subset_stem_other[subset_stem_other$clade=="Malpighiaceae",])

  ####
  pdf(file=paste0("revision_results/revision_results_background_families/", families[family_index],".pdf"), height=5, width=3)
  ########################################################################################
  par(mar=c(3,2,1,0.5))
  par(lwd=.3)
  #par(mfrow=c(2,1))
  

  plot(x=subset_stem_t$age_mean[subset_stem_t$clade==families[family_index]], y=subset_stem_t$diversity[subset_stem_t$clade==families[family_index]], ylim=c(1,1000), xlim=c(0,30), log="y", xaxt="n", yaxt="n", bty="n", xlab="", ylab="", pch=24, bg=pal[families[family_index]], cex=1, lwd=0.5)
  points(x=subset_stem_other$age_mean[subset_stem_other$clade==families[family_index]], y=subset_stem_other$diversity[subset_stem_other$clade==families[family_index]], pch=21, bg=pal[families[family_index]], cex=1, lwd=0.5)
  
  # adding axes
  par(lwd=.6)
  axis(side=1, at = seq(0,260,10), labels = seq(0,260,10), tick = TRUE, line = 0, lwd = .4, cex.axis=.6, tcl=NA, mgp=c(0,0.1,0),las=1, cex.lab=.4)
  axis(side=2, at = 10^{0:5}, labels = 10^{0:5}, tick = TRUE, line = 0, lwd = .4, cex.axis=.6, tcl=NA, mgp=c(2.5,0.2,0),las=1)
  
  mtext(text="Stem age (Mya)", side=1, line=1, cex=.5)
  par(las=0)
  mtext(text="Species richness", side=2, line=1.5, cex=.5)
  par(lwd=.5)
  
  
  # adding CI lines
  length_x <- length(limits.df_low$lb[limits.df_low$Diversity==  lsca.df.sum1$diversity[lsca.df.sum1$taxa==families[family_index]]])
  lines(y=limits.df_low$lb[limits.df_low$Diversity== lsca.df.sum1$diversity[lsca.df.sum1$taxa==families[family_index]]], x=1:length_x)
  lines(y=limits.df_low$ub[limits.df_low$Diversity== lsca.df.sum1$diversity[lsca.df.sum1$taxa==families[family_index]]], x=1:length_x)
  lines(y=limits.df_high$lb[limits.df_high$Diversity== lsca.df.sum1$diversity[lsca.df.sum1$taxa==families[family_index]]], x=1:length_x,lty=2)
  lines(y=limits.df_high$ub[limits.df_high$Diversity== lsca.df.sum1$diversity[lsca.df.sum1$taxa==families[family_index]]], x=1:length_x,lty=2)
  
  #legend(lty=1:5, c("All angiosperms","Magnoliids","Monocots","Superasterids","Superrosids"), x=150, y=5, bty="n", cex = .7, lwd=.8)
  legend(pch=c(21,24), legend = c("other", "tendrils"), x=0, y=1000, bty="n", cex = .7, pt.cex=1.25)
  
  ################
  test.line_low <- approxfun(x=1:length_x, y=limits.df_low$lb[limits.df_low$Diversity== lsca.df.sum1$diversity[lsca.df.sum1$taxa==families[family_index]]])
  test.line_high <- approxfun(x=1:length_x, y=limits.df_high$ub[limits.df_high$Diversity== lsca.df.sum1$diversity[lsca.df.sum1$taxa==families[family_index]]])
  
  keep <- c()
  for(i in 1:nrow(subset_stem)) {
    if(!is.na(test.line_low(subset_stem$age_mean[i]))){
      if(test.line_low(subset_stem$age_mean[i]) > subset_stem$diversity[i])
        keep <- c(keep,i)    
    }
  }
  
  # clades below lower bound:
  lower_div <- as.data.frame(subset_stem[keep,])
  lower_div <- subset(lower_div, lower_div$clade==families[family_index])
  write.csv(lower_div, file=paste0("revision_results/revision_results_background_families/clades_below_lower_bound_stem_",families[family_index],".csv"))
  
  
  keep <- c()
  for(i in 1:nrow(subset_stem)) {
    if(!is.na(test.line_high(subset_stem$age_mean[i]))){
      if(test.line_high(subset_stem$age_mean[i]) < subset_stem$diversity[i])
        keep <- c(keep,i)    
    }
  }
  
  # clades above upper bound:
  upper_div <- as.data.frame(subset_stem[keep,])
  upper_div <- subset(upper_div, upper_div$clade==families[family_index])
  write.csv(upper_div, file=paste0("revision_results/revision_results_background_families/clades_above_upper_bound_stem_",families[family_index],".csv"))
  
  
  nrow(lower_div) # number below lower bound
  nrow(upper_div) # number above upper bound
  #nrow(subset_stem[subset_stem$clade=="Malpighiaceae",]) - nrow(lower_div) - nrow(upper_div) # number below lower bound
  
  dev.off()
}




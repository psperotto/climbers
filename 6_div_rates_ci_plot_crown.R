#####################################################################
#####################################################################
# Code for confidence interval plots (crown node)
#
# rm(list=ls())

library(ape)
library(phangorn)
library(parallel)
library(geiger)

setwd("~/Desktop/WCVP_special_issue/Patricia_Climbers/climbers")
climber_clades_crown <- readRDS("crown.bg.clades.list.Rdata")
all_climbers_crown <- climber_clades_crown[[7]]
all_climbers_crown <- subset(all_climbers_crown, all_climbers_crown$diversity > 1) #need to remove those with only one species

#####################################################################
# crown node and high extinction

# get diversification rates
lsca.df.sum1 <- all_climbers_crown[1:6,]

epsilon <- 0.9 # 
root.age <- as.numeric(lsca.df.sum1[1,'age_mean'])

r_all_angios <- bd.ms(time = as.numeric(lsca.df.sum1[1,'age_mean']), n = as.numeric(lsca.df.sum1[1,'diversity']), crown = FALSE, epsilon = epsilon)
r_magnoliids <- bd.ms(time = as.numeric(lsca.df.sum1[2,'age_mean']), n = as.numeric(lsca.df.sum1[2,'diversity']), crown = FALSE, epsilon = epsilon)
r_monocots <- bd.ms(time = as.numeric(lsca.df.sum1[3,'age_mean']), n = as.numeric(lsca.df.sum1[3,'diversity']), crown = FALSE, epsilon = epsilon)
r_rosids <- bd.ms(time = as.numeric(lsca.df.sum1[4,'age_mean']), n = as.numeric(lsca.df.sum1[4,'diversity']), crown = FALSE, epsilon = epsilon)
r_asterids <- bd.ms(time = as.numeric(lsca.df.sum1[5,'age_mean']), n = as.numeric(lsca.df.sum1[5,'diversity']), crown = FALSE, epsilon = epsilon)
r_ranunculales <- bd.ms(time = as.numeric(lsca.df.sum1[6,'age_mean']), n = as.numeric(lsca.df.sum1[6,'diversity']), crown = FALSE, epsilon = epsilon)

tot.div1 <- c(as.numeric(lsca.df.sum1[1,'diversity']),
              as.numeric(lsca.df.sum1[2,'diversity']),
              as.numeric(lsca.df.sum1[3,'diversity']),
              as.numeric(lsca.df.sum1[4,'diversity']),
              as.numeric(lsca.df.sum1[5,'diversity']),
              as.numeric(lsca.df.sum1[6,'diversity']))
all_r <- c(r_all_angios, r_magnoliids, r_monocots, r_rosids, r_asterids, r_ranunculales)

getcrownLimits <- function(time, e=epsilon, r) {
  t(sapply(time, crown.limits, epsilon=e, r=r))
}

limits <- mclapply(all_r, getcrownLimits, time=2:root.age, e=epsilon, mc.cores = detectCores()-1)
limits.df <- data.frame(rbind(limits[[1]], limits[[2]],
                              limits[[3]], limits[[4]],
                              limits[[5]], limits[[6]]),
                        Diversity=factor(rep(tot.div1, each=nrow(limits[[1]]))),
                        r=rep(all_r, each=nrow(limits[[1]])))
colnames(limits.df)[1:2] <- c("lb", "ub")


################################################
subset_crown <- subset(all_climbers_crown, !is.na(all_climbers_crown$clade))
subset_crown_t <- subset(subset_crown, subset_crown$CM==2)
subset_crown_other <- subset(subset_crown, subset_crown$CM!=2)

pdf("ci_plot_crown_nodes.pdf",height=5,width=3)
# then mag sand
par(mar=c(3,2,1,0.5))
par(lwd=.3)
par(mfrow=c(2,1))

pal <- hcl.colors(5, palette = "Viridis", alpha = 0.7)
names(pal) <- c("Superasterids","Superrosids","Ranunculales","Monocots","Magnoliids")

plot(x=subset_crown_t$age_mean[subset_crown_t$clade=="Superrosids"], y=subset_crown_t$diversity[subset_crown_t$clade=="Superrosids"], ylim=c(1,1000), xlim=c(0,30), log="y", xaxt="n", yaxt="n", bty="n", xlab="", ylab="", pch=24, bg=pal['Superrosids'], cex=1, lwd=0.5)
points(x=subset_crown_other$age_mean[subset_crown_other$clade=="Superrosids"], y=subset_crown_other$diversity[subset_crown_other$clade=="Superrosids"], pch=21, bg=pal['Superrosids'], cex=1, lwd=0.5)

#points(x=subset_crown_t$age_mean[subset_crown_t$clade=="Superasterids"], y=subset_crown_t$diversity[subset_crown_t$clade=="Superasterids"], pch=24, bg=pal['Superasterids'], cex=1, lwd=0.5)
#points(x=subset_crown_other$age_mean[subset_crown_other$clade=="Superasterids"], y=subset_crown_other$diversity[subset_crown_other$clade=="Superasterids"], pch=21, bg=pal['Superasterids'], cex=1, lwd=0.5)

#points(x=subset_crown_t$age_mean[subset_crown_t$clade=="Magnoliids"], y=subset_crown_t$diversity[subset_crown_t$clade=="Magnoliids"], pch=24, bg=pal['Magnoliids'], cex=1, lwd=0.5)
#points(x=subset_crown_other$age_mean[subset_crown_other$clade=="Magnoliids"], y=subset_crown_other$diversity[subset_crown_other$clade=="Magnoliids"], pch=21, bg=pal['Magnoliids'], cex=1, lwd=0.5)

#points(x=subset_crown_t$age_mean[subset_crown_t$clade=="Monocots"], y=subset_crown_t$diversity[subset_crown_t$clade=="Monocots"], pch=24, bg=pal['Monocots'], cex=1, lwd=0.5)
#points(x=subset_crown_other$age_mean[subset_crown_other$clade=="Monocots"], y=subset_crown_other$diversity[subset_crown_other$clade=="Monocots"], pch=21, bg=pal['Monocots'], cex=1, lwd=0.5)

#points(x=lsca.df.sum$Age_mean, y=lsca.df.sum$Clade_richness, pch=22, bg='white', cex=1.2)
#text(x=lsca.df.sum$Age_mean+c(-10,10,-10), y=lsca.df.sum$Clade_richness-c(1000,2000,1000), c("D", "P", "R"), cex=.8)

# adding axes
par(lwd=.6)
axis(side=1, at = seq(0,260,10), labels = seq(0,260,10), tick = TRUE, line = 0, lwd = .4, cex.axis=.6, tcl=NA, mgp=c(0,0.1,0),las=1, cex.lab=.4)
axis(side=2, at = 10^{0:5}, labels = 10^{0:5}, tick = TRUE, line = 0, lwd = .4, cex.axis=.6, tcl=NA, mgp=c(2.5,0.2,0),las=1)

mtext(text="Crown age (Mya)", side=1, line=1, cex=.5)
par(las=0)
mtext(text="Species richness", side=2, line=1.5, cex=.5)
par(lwd=.5)

# adding CI lines
length_x <- length(limits.df$lb[limits.df$Diversity==352000])
lines(y=limits.df$lb[limits.df$Diversity==352000], x=1:length_x)
lines(y=limits.df$ub[limits.df$Diversity==352000], x=1:length_x)
#lines(y=limits.df$lb[limits.df$Diversity==10293], x=1:length_x,lty=2)
#lines(y=limits.df$ub[limits.df$Diversity==10293], x=1:length_x,lty=2)
#lines(y=limits.df$lb[limits.df$Diversity==69335], x=1:length_x,lty=3)
#lines(y=limits.df$ub[limits.df$Diversity==69335], x=1:length_x,lty=3)
#lines(y=limits.df$lb[limits.df$Diversity==87302], x=1:length_x,lty=4)
#lines(y=limits.df$ub[limits.df$Diversity==87302], x=1:length_x,lty=4)
lines(y=limits.df$lb[limits.df$Diversity==111856], x=1:length_x,lty=5)
lines(y=limits.df$ub[limits.df$Diversity==111856], x=1:length_x,lty=5)

#legend(lty=1:5, c("All angiosperms","Magnoliids","Monocots","Superasterids","Superrosids"), x=150, y=5, bty="n", cex = .7, lwd=.8)
legend(pch=c(21,24), legend = c("other", "tendrils"), x=0, y=1000, bty="n", cex = .7, pt.cex=1.25)

# ASTERIDS
pal <- hcl.colors(5, palette = "Viridis", alpha = 0.7)
names(pal) <- c("Superasterids","Superrosids","Ranunculales","Monocots","Magnoliids")

plot(x=subset_crown_t$age_mean[subset_crown_t$clade=="Superasterids"], y=subset_crown_t$diversity[subset_crown_t$clade=="Superasterids"], ylim=c(1,1000), xlim=c(0,30), log="y", xaxt="n", yaxt="n", bty="n", xlab="", ylab="", pch=24, bg=pal['Superasterids'], cex=1, lwd=0.5)
points(x=subset_crown_other$age_mean[subset_crown_other$clade=="Superasterids"], y=subset_crown_other$diversity[subset_crown_other$clade=="Superasterids"], pch=21, bg=pal['Superasterids'], cex=1, lwd=0.5)

# adding axes
par(lwd=.6)
axis(side=1, at = seq(0,260,10), labels = seq(0,260,10), tick = TRUE, line = 0, lwd = .4, cex.axis=.6, tcl=NA, mgp=c(0,0.1,0),las=1, cex.lab=.4)
axis(side=2, at = 10^{0:5}, labels = 10^{0:5}, tick = TRUE, line = 0, lwd = .4, cex.axis=.6, tcl=NA, mgp=c(2.5,0.2,0),las=1)

mtext(text="Crown age (Mya)", side=1, line=1, cex=.5)
par(las=0)
mtext(text="Species richness", side=2, line=1.5, cex=.5)
par(lwd=.5)

# adding CI lines
length_x <- length(limits.df$lb[limits.df$Diversity==352000])
lines(y=limits.df$lb[limits.df$Diversity==352000], x=1:length_x)
lines(y=limits.df$ub[limits.df$Diversity==352000], x=1:length_x)
#lines(y=limits.df$lb[limits.df$Diversity==10293], x=1:length_x,lty=2)
#lines(y=limits.df$ub[limits.df$Diversity==10293], x=1:length_x,lty=2)
#lines(y=limits.df$lb[limits.df$Diversity==69335], x=1:length_x,lty=3)
#lines(y=limits.df$ub[limits.df$Diversity==69335], x=1:length_x,lty=3)
lines(y=limits.df$lb[limits.df$Diversity==87302], x=1:length_x,lty=4)
lines(y=limits.df$ub[limits.df$Diversity==87302], x=1:length_x,lty=4)
#lines(y=limits.df$lb[limits.df$Diversity==111856], x=1:length_x,lty=5)
#lines(y=limits.df$ub[limits.df$Diversity==111856], x=1:length_x,lty=5)

#legend(lty=1:5, c("All angiosperms","Magnoliids","Monocots","Superasterids","Superrosids"), x=150, y=5, bty="n", cex = .7, lwd=.8)
legend(pch=c(21,24), legend = c("other", "tendrils"), x=0, y=1000, bty="n", cex = .7, pt.cex=1.25)

dev.off()

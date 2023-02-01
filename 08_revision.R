
setwd("~/Desktop/climbers")

library(ape)
library(phangorn)
library(parallel)
library(geiger)

# Reviewer suggestion
# Try diversification analysis with different cutouts

genera75 <-readRDS("genera75.Rdata")

#####################################################################
#####################################################################
# Code for confidence interval plots (stem node)
#
#rm(list=ls())

climber_clades_stem <- readRDS("stem.bg.clades.list.Rdata")
all_climbers_stem <- climber_clades_stem[[7]]

# now let's exclude everything that is not 90% climber
genera90 <- subset(genera75, genera75$perc_NT_spp>0.9)

lsca.df.sum1 <- all_climbers_stem[1:6,]
all_climbers_stem <- subset(all_climbers_stem, all_climbers_stem$taxa %in% genera90$Genus)
all_climbers_stem <- rbind(lsca.df.sum1, all_climbers_stem)

#####################################################################
# Stem node and high extinction

# get diversification rates

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

getStemLimits <- function(time, e=epsilon, r) {
  t(sapply(time, stem.limits, epsilon=e, r=r))
}

limits <- mclapply(all_r, getStemLimits, time=2:root.age, e=epsilon, mc.cores = detectCores()-1)
limits.df <- data.frame(rbind(limits[[1]], limits[[2]],
                              limits[[3]], limits[[4]],
                              limits[[5]], limits[[6]]),
                        Diversity=factor(rep(tot.div1, each=nrow(limits[[1]]))),
                        r=rep(all_r, each=nrow(limits[[1]])))
colnames(limits.df)[1:2] <- c("lb", "ub")


################################################
subset_stem <- subset(all_climbers_stem, !is.na(all_climbers_stem$clade))
subset_stem_t <- subset(subset_stem, subset_stem$CM==2)
subset_stem_other <- subset(subset_stem, subset_stem$CM!=2)

pdf("ci_plot_stem_nodes_90cutoff.pdf",height=5,width=3)
########################################################################################
# ROSIDS
par(mar=c(3,2,1,0.5))
par(lwd=.3)
par(mfrow=c(2,1))

pal <- hcl.colors(5, palette = "Viridis", alpha = 0.7)
names(pal) <- c("Superasterids","Superrosids","Ranunculales","Monocots","Magnoliids")

plot(x=subset_stem_t$age_mean[subset_stem_t$clade=="Superrosids"], y=subset_stem_t$diversity[subset_stem_t$clade=="Superrosids"], ylim=c(1,1000), xlim=c(0,30), log="y", xaxt="n", yaxt="n", bty="n", xlab="", ylab="", pch=24, bg=pal['Superrosids'], cex=1, lwd=0.5)
points(x=subset_stem_other$age_mean[subset_stem_other$clade=="Superrosids"], y=subset_stem_other$diversity[subset_stem_other$clade=="Superrosids"], pch=21, bg=pal['Superrosids'], cex=1, lwd=0.5)

# adding axes
par(lwd=.6)
axis(side=1, at = seq(0,260,10), labels = seq(0,260,10), tick = TRUE, line = 0, lwd = .4, cex.axis=.6, tcl=NA, mgp=c(0,0.1,0),las=1, cex.lab=.4)
axis(side=2, at = 10^{0:5}, labels = 10^{0:5}, tick = TRUE, line = 0, lwd = .4, cex.axis=.6, tcl=NA, mgp=c(2.5,0.2,0),las=1)

mtext(text="Stem age (Mya)", side=1, line=1, cex=.5)
par(las=0)
mtext(text="Species richness", side=2, line=1.5, cex=.5)
par(lwd=.5)

# adding CI lines
length_x <- length(limits.df$lb[limits.df$Diversity==352000])
lines(y=limits.df$lb[limits.df$Diversity==352000], x=1:length_x)
lines(y=limits.df$ub[limits.df$Diversity==352000], x=1:length_x)
lines(y=limits.df$lb[limits.df$Diversity==111856], x=1:length_x,lty=2)
lines(y=limits.df$ub[limits.df$Diversity==111856], x=1:length_x,lty=2)

#legend(lty=1:5, c("All angiosperms","Magnoliids","Monocots","Superasterids","Superrosids"), x=150, y=5, bty="n", cex = .7, lwd=.8)
legend(pch=c(21,24), legend = c("other", "tendrils"), x=0, y=1000, bty="n", cex = .7, pt.cex=1.25)

############################################################
# ASTERIDS 
par(mar=c(3,2,1,0.5))
par(lwd=.3)

pal <- hcl.colors(5, palette = "Viridis", alpha = 0.7)
names(pal) <- c("Superasterids","Superrosids","Ranunculales","Monocots","Magnoliids")

plot(x=subset_stem_t$age_mean[subset_stem_t$clade=="Superasterids"], y=subset_stem_t$diversity[subset_stem_t$clade=="Superasterids"], ylim=c(1,1000), xlim=c(0,30), log="y", xaxt="n", yaxt="n", bty="n", xlab="", ylab="", pch=24, bg=pal['Superasterids'], cex=1, lwd=0.5)
points(x=subset_stem_other$age_mean[subset_stem_other$clade=="Superasterids"], y=subset_stem_other$diversity[subset_stem_other$clade=="Superasterids"], pch=21, bg=pal['Superasterids'], cex=1, lwd=0.5)

# adding axes
par(lwd=.6)
axis(side=1, at = seq(0,260,10), labels = seq(0,260,10), tick = TRUE, line = 0, lwd = .4, cex.axis=.6, tcl=NA, mgp=c(0,0.1,0),las=1, cex.lab=.4)
axis(side=2, at = 10^{0:5}, labels = 10^{0:5}, tick = TRUE, line = 0, lwd = .4, cex.axis=.6, tcl=NA, mgp=c(2.5,0.2,0),las=1)

mtext(text="Stem age (Mya)", side=1, line=1, cex=.5)
par(las=0)
mtext(text="Species richness", side=2, line=1.5, cex=.5)
par(lwd=.5)

# adding CI lines
length_x <- length(limits.df$lb[limits.df$Diversity==352000])
lines(y=limits.df$lb[limits.df$Diversity==352000], x=1:length_x)
lines(y=limits.df$ub[limits.df$Diversity==352000], x=1:length_x)
lines(y=limits.df$lb[limits.df$Diversity==87302], x=1:length_x,lty=2)
lines(y=limits.df$ub[limits.df$Diversity==87302], x=1:length_x,lty=2)

#legend(lty=1:5, c("All angiosperms","Magnoliids","Monocots","Superasterids","Superrosids"), x=150, y=5, bty="n", cex = .7, lwd=.8)
legend(pch=c(21,24), legend = c("other", "tendrils"), x=0, y=1000, bty="n", cex = .7, pt.cex=1.25)

dev.off()



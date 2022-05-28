# Diversification analyses
setwd("~/Desktop/WCVP_special_issue/Patricia_Climbers/climbers")
#rm(list=ls())

climber_clades <- readRDS("bg.clades.list.Rdata")

for (i in 1:length(climber_clades)) {
test <- climber_clades[[i]]

CMs <- as.numeric(test$CM[-1])
CMs[which(CMs==3)] <- 2
CMs[which(CMs==4)] <- 1
CMs[which(CMs==5)] <- 1
CMs[which(CMs==6)] <- 2
CMs[which(CMs==7)] <- 1
CMs[which(CMs==8)] <- 2

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

pdf(file=paste0(test$clade[2],".pdf"))
#dev.off()
plot(k.lower.yule.crown$time, log10(k.lower.yule.crown$k), type="l",
     col="gray", lwd=1.5, ylim = c(0, round(log10(k), 0)+1), xlim = c(0,time+time*0.1), xlab="", ylab="", yaxt = "n")
par(new=TRUE)

plot(k.upper.yule.crown$time, log10(k.upper.yule.crown$k), type="l",
     col="gray", lwd=1.5, ylim = c(0, round(log10(k), 0)+1), xlim = c(0,time+time*0.1), xlab="", ylab="", yaxt = "n")
par(new=TRUE)
plot(k.lower.high.ext.crown$time, log10(k.lower.high.ext.crown$k), type="l",lty=3,
     col="gray",lwd=1.5, ylim = c(0, round(log10(k), 0)+1), xlim = c(0,time+time*0.1), xlab="", ylab="",  yaxt = "n")
par(new=TRUE)
plot(k.upper.high.ext.crown$time, log10(k.upper.high.ext.crown$k), type="l",lty=3,
     col="gray",lwd=1.5, ylim = c(0, round(log10(k), 0)+1), xlim = c(0,time+time*0.1), xlab='Age of clade (million years)', ylab='K (expected number of species)',yaxt = "n")
par(new=TRUE)

# defining color vector
mode <- test$CM[-1]
colors_states <- c("midnightblue", "goldenrod")
mode_cols <- mode
mode_cols[mode_cols==1] <- "goldenrod"
mode_cols[mode_cols==2] <- "midnightblue"

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

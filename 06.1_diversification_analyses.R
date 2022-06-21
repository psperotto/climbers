# Diversification analyses
setwd("~/Desktop/WCVP_special_issue/Patricia_Climbers/climbers")
#rm(list=ls())
library(geiger)
library(phytools)
source("00_utility_functions.R")

#############################################################
#############################################################
climber_clades <- readRDS("bg.clades.list.Rdata")
all_climbers_stem <- climber_clades[[1]]
all_climbers_stem <- all_climbers_stem[-1,]

all_climbers_stem$div_rate_eps0 <- NA
all_climbers_stem$div_rate_eps0.5 <- NA
all_climbers_stem$div_rate_eps0.9 <- NA

for(u in 1:nrow(all_climbers_stem)) {
  all_climbers_stem$div_rate_eps0[u] <- round(bd.ms(phy=NULL, all_climbers_stem$age_mean[u], all_climbers_stem$diversity[u], crown=F, epsilon=0),3)
  all_climbers_stem$div_rate_eps0.5[u] <- round(bd.ms(phy=NULL, all_climbers_stem$age_mean[u], all_climbers_stem$diversity[u], crown=F, epsilon=0.5),3)
  all_climbers_stem$div_rate_eps0.9[u] <- round(bd.ms(phy=NULL, all_climbers_stem$age_mean[u], all_climbers_stem$diversity[u], crown=F, epsilon=0.9),3)
}

all_genera_tree <- read.tree("climber_genera_tree.tre")

#########################################################
# Comparison including both tendrils (trait 2) and adhesive roots (trait 4) as specialized mechanism
all_climbers_comparison1 = all_climbers_stem
all_climbers_comparison1$CM[which(all_climbers_comparison1$CM==1)] <- "Other"
all_climbers_comparison1$CM[which(all_climbers_comparison1$CM==2)] <- "Specialist"
all_climbers_comparison1$CM[which(all_climbers_comparison1$CM==3)] <- "Other"
all_climbers_comparison1$CM[which(all_climbers_comparison1$CM==4)] <- "Specialist" 
all_climbers_comparison1$CM[which(all_climbers_comparison1$CM==5)] <- "Other"
all_climbers_comparison1$CM[which(all_climbers_comparison1$CM==6)] <- "Other"
all_climbers_comparison1$CM[which(all_climbers_comparison1$CM==7)] <- "Other"
all_climbers_comparison1$CM[which(all_climbers_comparison1$CM==8)] <- "Other"

div_rates1 <- all_climbers_comparison1$div_rate_eps0.9
names(div_rates1) <- all_climbers_comparison1$taxa
mech <- all_climbers_comparison1$CM
names(mech) <- all_climbers_comparison1$taxa

climb_comparison1 <- phylANOVA(all_genera_tree, mech, div_rates1,nsim=10000, p.adj="bonferroni")

#########################################################
# Comparison including only tendrils (trait 2) as specialized mechanism 
all_climbers_comparison2 = all_climbers_stem
all_climbers_comparison2$CM[which(all_climbers_comparison2$CM==1)] <- "Other"
all_climbers_comparison2$CM[which(all_climbers_comparison2$CM==2)] <- "Specialist"
all_climbers_comparison2$CM[which(all_climbers_comparison2$CM==3)] <- "Other"
all_climbers_comparison2$CM[which(all_climbers_comparison2$CM==4)] <- "Other" 
all_climbers_comparison2$CM[which(all_climbers_comparison2$CM==5)] <- "Other"
all_climbers_comparison2$CM[which(all_climbers_comparison2$CM==6)] <- "Other"
all_climbers_comparison2$CM[which(all_climbers_comparison2$CM==7)] <- "Other"
all_climbers_comparison2$CM[which(all_climbers_comparison2$CM==8)] <- "Other"

div_rates1 <- all_climbers_comparison2$div_rate_eps0.9
names(div_rates1) <- all_climbers_comparison2$taxa
mech <- all_climbers_comparison2$CM
names(mech) <- all_climbers_comparison2$taxa

climb_comparison2 <- phylANOVA(all_genera_tree, mech, div_rates1,nsim=10000, p.adj="bonferroni")


# Plots
frst_comp_boxplot <- boxplot(all_climbers_comparison1$div_rate_eps0.9~all_climbers_comparison1$CM, ylim=c(0,1))
scnd_comp_boxplot <- boxplot(all_climbers_comparison2$div_rate_eps0.9~all_climbers_comparison2$CM, ylim=c(0,1))


library(ggplot2)
library(ggridges)


pdf("div_rate_comparisons.pdf",height=3,width=5.5)
plot_comparisons <- ggplot(all_climbers_comparison1, aes(x=CM, y=div_rate_eps0.9, fill=CM)) + 
  geom_boxplot(lwd=0.3, outlier.size=0.25, alpha=0.8) + 
  geom_jitter(colour = 2, alpha=0.3, size=0.9) +
  #coord_flip()  + 
  scale_fill_brewer(palette="Pastel1")+
  theme_bw(base_size = 8) + 
  theme(legend.position = "none") + 
  annotate(geom="text", x=.9, y=0.9, size=3, hjust=0, label="") + 
  xlab("") +
  ylab("Mean aridity index") +
  theme(axis.text.y = element_text(colour = 'black', size = 10),
        axis.text.x = element_text(colour = 'black', size = 10),
        axis.title.x = element_text(colour = 'black', size = 10)) 
plot_comparisons
dev.off()



#write.csv(myrtales_final, file="myrtales_div_table_full.csv", row.names = F)

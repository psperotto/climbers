# Diversification analyses
setwd("~/Desktop/WCVP_special_issue/Patricia_Climbers/climbers")
#rm(list=ls())
library(geiger)
library(phytools)
library(ggplot2)
library(ggridges)
library(gridExtra)

source("00_utility_functions.R")

#pal <- hcl.colors(5, palette = "Set3", alpha = 0.7)
#names(pal) <- c("Superasterids","Superrosids","Ranunculales","Monocots","Magnoliids")

#############################################################
#############################################################
# STEM NODE 
climber_clades_stem <- readRDS("stem.bg.clades.list.Rdata")
all_climbers_stem <- climber_clades_stem[[1]]
all_climbers_stem <- all_climbers_stem[-1,]

all_climbers_stem$div_rate_eps0 <- NA
all_climbers_stem$div_rate_eps0.5 <- NA
all_climbers_stem$div_rate_eps0.9 <- NA

for(u in 1:nrow(all_climbers_stem)) {
  all_climbers_stem$div_rate_eps0[u] <- round(bd.ms(phy=NULL, all_climbers_stem$age_mean[u], all_climbers_stem$diversity[u], crown=F, epsilon=0),3)
  all_climbers_stem$div_rate_eps0.5[u] <- round(bd.ms(phy=NULL, all_climbers_stem$age_mean[u], all_climbers_stem$diversity[u], crown=F, epsilon=0.5),3)
  all_climbers_stem$div_rate_eps0.9[u] <- round(bd.ms(phy=NULL, all_climbers_stem$age_mean[u], all_climbers_stem$diversity[u], crown=F, epsilon=0.9),3)
}

#View(all_climbers_stem)
all_genera_tree <- read.tree("climber_genera_tree.tre")

#########################################################
# Comparison including both tendrils (trait 2) and adhesive roots (trait 4) as specialized mechanism
all_climbers_comparison1 = all_climbers_stem
all_climbers_comparison1$CM[which(all_climbers_comparison1$CM==1)] <- "Other"
all_climbers_comparison1$CM[which(all_climbers_comparison1$CM==2)] <- "Specialist"
all_climbers_comparison1$CM[which(all_climbers_comparison1$CM==3)] <- "Other"
all_climbers_comparison1$CM[which(all_climbers_comparison1$CM==4)] <- "Other" 
all_climbers_comparison1$CM[which(all_climbers_comparison1$CM==5)] <- "Other"
all_climbers_comparison1$CM[which(all_climbers_comparison1$CM==6)] <- "Other"
all_climbers_comparison1$CM[which(all_climbers_comparison1$CM==7)] <- "Other"
all_climbers_comparison1$CM[which(all_climbers_comparison1$CM==8)] <- "Other"

#########################################################
# Superrosids
all_climbers_comparison1_sros <- subset(all_climbers_comparison1, all_climbers_comparison1$clade=="Superrosids")

div_rates_sros <- all_climbers_comparison1_sros$div_rate_eps0.9
names(div_rates_sros) <- all_climbers_comparison1_sros$taxa
mech_sros <- all_climbers_comparison1_sros$CM
names(mech_sros) <- all_climbers_comparison1_sros$taxa

all_genera_tree_ros <- keep.tip(all_genera_tree, all_climbers_comparison1_sros$taxa)
climb_comparison1_ros <- phylANOVA(all_genera_tree_ros, mech_sros, div_rates_sros ,nsim=1000, p.adj="bonferroni")

#########################################################
# Superasterids
all_climbers_comparison1_sast <- subset(all_climbers_comparison1, all_climbers_comparison1$clade=="Superasterids")

div_rates_sast <- all_climbers_comparison1_sast$div_rate_eps0.9
names(div_rates_sast) <- all_climbers_comparison1_sast$taxa
mech_sast <- all_climbers_comparison1_sast$CM
names(mech_sast) <- all_climbers_comparison1_sast$taxa

all_genera_tree_ast <- keep.tip(all_genera_tree, all_climbers_comparison1_sast$taxa)
climb_comparison1_ast <- phylANOVA(all_genera_tree_ast, mech_sast, div_rates_sast, nsim=1000, p.adj="bonferroni")

plot_comparisons_ros1_stem <- ggplot(all_climbers_comparison1_sros, aes(x=CM, y=div_rate_eps0.9, fill=CM)) + 
  geom_boxplot(lwd=0.3, outlier.size=0.25, alpha=0.8) + 
  #geom_jitter(colour = 2, alpha=0.3, size=0.9) +
  coord_cartesian(ylim = c(0, 1)) +
  #coord_flip()  + 
  scale_fill_manual(values = c("#00588BB3","#00588BB3")) +
  theme_bw(base_size = 8) + 
  theme(legend.position = "none") + 
  annotate(geom="text", x=.9, y=0.9, size=3, hjust=0, label="") + 
  xlab("") +
  ylab("net-diversification rates") +
  theme(axis.text.y = element_text(colour = 'black', size = 5),
        axis.text.x = element_text(colour = 'black', size = 5),
        axis.title.x = element_text(colour = 'black', size = 10)) 

plot_comparisons_ast1_stem <- ggplot(all_climbers_comparison1_sast, aes(x=CM, y=div_rate_eps0.9, fill=CM)) + 
  geom_boxplot(lwd=0.3, outlier.size=0.25, alpha=0.8) + 
  #stat_compare_means(comparisons = 3) + # Add pairwise comparisons p-value
  #geom_jitter(colour = 2, alpha=0.3, size=0.9) +
  coord_cartesian(ylim = c(0, 1)) +
  #coord_flip()  + 
  scale_fill_manual( values = c("#4B0055B3", "#4B0055B3")) +
  theme_bw(base_size = 8) + 
  theme(legend.position = "none") + 
  annotate(geom="text", x=.9, y=0.9, size=3, hjust=0, label="") + 
  xlab("") +
  ylab("net-diversification rates") +
  theme(axis.text.y = element_text(colour = 'black', size = 5),
        axis.text.x = element_text(colour = 'black', size = 5),
        axis.title.x = element_text(colour = 'black', size = 10)) 

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

#########################################################
# Superrosids
all_climbers_comparison2_sros <- subset(all_climbers_comparison2, all_climbers_comparison2$clade=="Superrosids")

div_rates_sros <- all_climbers_comparison2_sros$div_rate_eps0.9
names(div_rates_sros) <- all_climbers_comparison2_sros$taxa
mech_sros <- all_climbers_comparison2_sros$CM
names(mech_sros) <- all_climbers_comparison2_sros$taxa

all_genera_tree_ros <- keep.tip(all_genera_tree, all_climbers_comparison2_sros$taxa)
climb_comparison2 <- phylANOVA(all_genera_tree_ros, mech_sros, div_rates_sros,nsim=1000, p.adj="bonferroni")

#########################################################
# Superasterids
all_climbers_comparison2_sast <- subset(all_climbers_comparison2, all_climbers_comparison2$clade=="Superasterids")

div_rates_sast <- all_climbers_comparison2_sast$div_rate_eps0.9
names(div_rates_sast) <- all_climbers_comparison2_sast$taxa
mech_sast <- all_climbers_comparison2_sast$CM
names(mech_sast) <- all_climbers_comparison2_sast$taxa

all_genera_tree_ast <- keep.tip(all_genera_tree, all_climbers_comparison2_sast$taxa)
climb_comparison2 <- phylANOVA(all_genera_tree_ast, mech_sast, div_rates_sast,nsim=1000, p.adj="bonferroni")

# PLOTS
pal <- hcl.colors(5, palette = "Viridis", alpha = 0.7)
names(pal) <- c("Superasterids","Superrosids","Ranunculales","Monocots","Magnoliids")

plot_comparisons_ros2_stem <- ggplot(all_climbers_comparison2_sros, aes(x=CM, y=div_rate_eps0.9, fill=CM)) + 
  geom_boxplot(lwd=0.3, outlier.size=0.25, alpha=0.8) + 
  #geom_jitter(colour = 2, alpha=0.3, size=0.9) +
  coord_cartesian(ylim = c(0, 1)) +
  #coord_flip()  + 
  scale_fill_manual(values = c("#00588BB3","#00588BB3")) +
  theme_bw(base_size = 8) + 
  theme(legend.position = "none") + 
  annotate(geom="text", x=.9, y=0.9, size=3, hjust=0, label="") + 
  xlab("") +
  ylab("net-diversification rates") +
  theme(axis.text.y = element_text(colour = 'black', size = 5),
        axis.text.x = element_text(colour = 'black', size = 5),
        axis.title.x = element_text(colour = 'black', size = 10)) 

plot_comparisons_ast2_stem <- ggplot(all_climbers_comparison2_sast, aes(x=CM, y=div_rate_eps0.9, fill=CM)) + 
  geom_boxplot(lwd=0.3, outlier.size=0.25, alpha=0.8) + 
  #stat_compare_means(comparisons = 3) + # Add pairwise comparisons p-value
  #geom_jitter(colour = 2, alpha=0.3, size=0.9) +
  coord_cartesian(ylim = c(0, 1)) +
  #coord_flip()  + 
  scale_fill_manual( values = c("#4B0055B3", "#4B0055B3")) +
  theme_bw(base_size = 8) + 
  theme(legend.position = "none") + 
  annotate(geom="text", x=.9, y=0.9, size=3, hjust=0, label="") + 
  xlab("") +
  ylab("net-diversification rates") +
  theme(axis.text.y = element_text(colour = 'black', size = 5),
        axis.text.x = element_text(colour = 'black', size = 5),
        axis.title.x = element_text(colour = 'black', size = 10)) 


#############################################################
#############################################################
# CROWN NODE
climber_clades_crown <- readRDS("crown.bg.clades.list.Rdata")
all_climbers_crown <- climber_clades_crown[[1]]
all_climbers_crown <- all_climbers_crown[-1,]
all_climbers_crown <- subset(all_climbers_crown, all_climbers_crown$diversity > 1) #need to remove those with only one species

all_climbers_crown$div_rate_eps0 <- NA
all_climbers_crown$div_rate_eps0.5 <- NA
all_climbers_crown$div_rate_eps0.9 <- NA

for(u in 1:nrow(all_climbers_crown)) {
  all_climbers_crown$div_rate_eps0[u] <- round(bd.ms(phy=NULL, all_climbers_crown$age_mean[u], all_climbers_crown$diversity[u], crown=T, epsilon=0),3)
  all_climbers_crown$div_rate_eps0.5[u] <- round(bd.ms(phy=NULL, all_climbers_crown$age_mean[u], all_climbers_crown$diversity[u], crown=T, epsilon=0.5),3)
  all_climbers_crown$div_rate_eps0.9[u] <- round(bd.ms(phy=NULL, all_climbers_crown$age_mean[u], all_climbers_crown$diversity[u], crown=T, epsilon=0.9),3)
}

#View(all_climbers_crown)
all_genera_tree <- read.tree("climber_genera_tree.tre")

#########################################################
#########################################################
# Comparison including both tendrils (trait 2) and adhesive roots (trait 4) as specialized mechanism
all_climbers_comparison1 = all_climbers_crown
all_climbers_comparison1$CM[which(all_climbers_comparison1$CM==1)] <- "Other"
all_climbers_comparison1$CM[which(all_climbers_comparison1$CM==2)] <- "Specialist"
all_climbers_comparison1$CM[which(all_climbers_comparison1$CM==3)] <- "Other"
all_climbers_comparison1$CM[which(all_climbers_comparison1$CM==4)] <- "Specialist" 
all_climbers_comparison1$CM[which(all_climbers_comparison1$CM==5)] <- "Other"
all_climbers_comparison1$CM[which(all_climbers_comparison1$CM==6)] <- "Other"
all_climbers_comparison1$CM[which(all_climbers_comparison1$CM==7)] <- "Other"
all_climbers_comparison1$CM[which(all_climbers_comparison1$CM==8)] <- "Other"

#########################################################
# Superrosids
all_climbers_comparison1_sros <- subset(all_climbers_comparison1, all_climbers_comparison1$clade=="Superrosids")

div_rates_sros <- all_climbers_comparison1_sros$div_rate_eps0.9
names(div_rates_sros) <- all_climbers_comparison1_sros$taxa
mech_sros <- all_climbers_comparison1_sros$CM
names(mech_sros) <- all_climbers_comparison1_sros$taxa

all_genera_tree_ros <- keep.tip(all_genera_tree, all_climbers_comparison1_sros$taxa)
climb_comparison1_ros <- phylANOVA(all_genera_tree_ros, mech_sros, div_rates_sros ,nsim=1000, p.adj="bonferroni")

#########################################################
# Superasterids
all_climbers_comparison1_sast <- subset(all_climbers_comparison1, all_climbers_comparison1$clade=="Superasterids")

div_rates_sast <- all_climbers_comparison1_sast$div_rate_eps0.9
names(div_rates_sast) <- all_climbers_comparison1_sast$taxa
mech_sast <- all_climbers_comparison1_sast$CM
names(mech_sast) <- all_climbers_comparison1_sast$taxa

all_genera_tree_ast <- keep.tip(all_genera_tree, all_climbers_comparison1_sast$taxa)
climb_comparison1_ast <- phylANOVA(all_genera_tree_ast, mech_sast, div_rates_sast, nsim=1000, p.adj="bonferroni")


plot_comparisons_ros1_crown <- ggplot(all_climbers_comparison1_sros, aes(x=CM, y=div_rate_eps0.9, fill=CM)) + 
  geom_boxplot(lwd=0.3, outlier.size=0.25, alpha=0.8) + 
  #geom_jitter(colour = 2, alpha=0.3, size=0.9) +
  coord_cartesian(ylim = c(0, 1)) +
  #coord_flip()  + 
  scale_fill_manual(values = c("#00588BB3","#00588BB3")) +
  theme_bw(base_size = 8) + 
  theme(legend.position = "none") + 
  annotate(geom="text", x=.9, y=0.9, size=3, hjust=0, label="") + 
  xlab("") +
  ylab("net-diversification rates") +
  theme(axis.text.y = element_text(colour = 'black', size = 5),
        axis.text.x = element_text(colour = 'black', size = 5),
        axis.title.x = element_text(colour = 'black', size = 10)) 

plot_comparisons_ast1_crown <- ggplot(all_climbers_comparison1_sast, aes(x=CM, y=div_rate_eps0.9, fill=CM)) + 
  geom_boxplot(lwd=0.3, outlier.size=0.25, alpha=0.8) + 
  #stat_compare_means(comparisons = 3) + # Add pairwise comparisons p-value
  #geom_jitter(colour = 2, alpha=0.3, size=0.9) +
  coord_cartesian(ylim = c(0, 1)) +
  #coord_flip()  + 
  scale_fill_manual( values = c("#4B0055B3", "#4B0055B3")) +
  theme_bw(base_size = 8) + 
  theme(legend.position = "none") + 
  annotate(geom="text", x=.9, y=0.9, size=3, hjust=0, label="") + 
  xlab("") +
  ylab("net-diversification rates") +
  theme(axis.text.y = element_text(colour = 'black', size = 5),
        axis.text.x = element_text(colour = 'black', size = 5),
        axis.title.x = element_text(colour = 'black', size = 10)) 

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

#########################################################
# Superrosids
all_climbers_comparison2_sros <- subset(all_climbers_comparison2, all_climbers_comparison2$clade=="Superrosids")

div_rates_sros <- all_climbers_comparison2_sros$div_rate_eps0.9
names(div_rates_sros) <- all_climbers_comparison2_sros$taxa
mech_sros <- all_climbers_comparison2_sros$CM
names(mech_sros) <- all_climbers_comparison2_sros$taxa

all_genera_tree_ros <- keep.tip(all_genera_tree, all_climbers_comparison2_sros$taxa)
climb_comparison2 <- phylANOVA(all_genera_tree_ros, mech_sros, div_rates_sros,nsim=1000, p.adj="bonferroni")

#########################################################
# Superasterids
all_climbers_comparison2_sast <- subset(all_climbers_comparison2, all_climbers_comparison2$clade=="Superasterids")

div_rates_sast <- all_climbers_comparison2_sast$div_rate_eps0.9
names(div_rates_sast) <- all_climbers_comparison2_sast$taxa
mech_sast <- all_climbers_comparison2_sast$CM
names(mech_sast) <- all_climbers_comparison2_sast$taxa

all_genera_tree_ast <- keep.tip(all_genera_tree, all_climbers_comparison2_sast$taxa)
climb_comparison2 <- phylANOVA(all_genera_tree_ast, mech_sast, div_rates_sast,nsim=1000, p.adj="bonferroni")


plot_comparisons_ros2_crown <- ggplot(all_climbers_comparison2_sros, aes(x=CM, y=div_rate_eps0.9, fill=CM)) + 
  geom_boxplot(lwd=0.3, outlier.size=0.25, alpha=0.8) + 
  #geom_jitter(colour = 2, alpha=0.3, size=0.9) +
  coord_cartesian(ylim = c(0, 1)) +
  #coord_flip()  + 
  scale_fill_manual(values = c("#00588BB3","#00588BB3")) +
  theme_bw(base_size = 8) + 
  theme(legend.position = "none") + 
  annotate(geom="text", x=.9, y=0.9, size=3, hjust=0, label="") + 
  xlab("") +
  ylab("net-diversification") +
  theme(axis.text.y = element_text(colour = 'black', size = 5),
        axis.text.x = element_text(colour = 'black', size = 5),
        axis.title.x = element_text(colour = 'black', size = 10)) 

plot_comparisons_ast2_crown <- ggplot(all_climbers_comparison2_sast, aes(x=CM, y=div_rate_eps0.9, fill=CM)) + 
  geom_boxplot(lwd=0.3, outlier.size=0.25, alpha=0.8) + 
  #stat_compare_means(comparisons = 3) + # Add pairwise comparisons p-value
  #geom_jitter(colour = 2, alpha=0.3, size=0.9) +
  coord_cartesian(ylim = c(0, 1)) +
  #coord_flip()  + 
  scale_fill_manual( values = c("#4B0055B3", "#4B0055B3")) +
  theme_bw(base_size = 8) + 
  theme(legend.position = "none") + 
  annotate(geom="text", x=.9, y=0.9, size=3, hjust=0, label="") + 
  xlab("") +
  ylab("net-diversification rates") +
  theme(axis.text.y = element_text(colour = 'black', size = 5),
        axis.text.x = element_text(colour = 'black', size = 5),
        axis.title.x = element_text(colour = 'black', size = 10)) 

#########################################################
#########################################################
# Plots
#frst_comp_boxplot <- boxplot(all_climbers_comparison1$div_rate_eps0.9~all_climbers_comparison1$CM, ylim=c(0,1))
#scnd_comp_boxplot <- boxplot(all_climbers_comparison2$div_rate_eps0.9~all_climbers_comparison2$CM, ylim=c(0,1))

pdf("div_rate_comparisons_stem.pdf",height=3,width=5.5)
grid.arrange(plot_comparisons_ast1_stem, plot_comparisons_ast2_stem, plot_comparisons_ros1_stem, plot_comparisons_ros2_stem, ncol=4, nrow = 1)
dev.off()

pdf("div_rate_comparisons_crown.pdf",height=3,width=5.5)
grid.arrange(plot_comparisons_ast1_crown, plot_comparisons_ast2_crown, plot_comparisons_ros1_crown, plot_comparisons_ros2_crown, ncol=4, nrow = 1)
dev.off()



# Diversification analyses
setwd("~/Desktop/WCVP_special_issue/Patricia_Climbers/climbers")
#rm(list=ls())
library(geiger)
library(phytools)
source("00_utility_functions.R")

#--------------------------------------
climber_clades_stem <- readRDS("stem.bg.clades.list.Rdata")
climber_clades_crown <- readRDS("crown.bg.clades.list.Rdata")


for(i in 1:length(climber_clades_stem)) {
  try(results <- get.tail.probs(table=climber_clades_stem[[i]]))
  if(exists("results")) {
    saveRDS(results, file=paste0(plot_ci_div_rates,"/stem_plot_ci_",i,".Rdata"))
    rm("results")
  }
}

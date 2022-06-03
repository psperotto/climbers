
setwd("~/Desktop/WCVP_special_issue/Patricia_Climbers/climbers")

rm(list=ls())

wcvp_climbers <- readRDS("habit_climbers_wcvp.Rdata")
wcvp_climbers <- readRDS("bg.clades.list.Rdata")
wcvp_climbers <- readRDS("wcvp_nt_climbers_final.Rdata")

View(wcvp_climbers)

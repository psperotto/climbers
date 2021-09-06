setwd("G:/Meu Drive/Papers/Diversificação/climbers")

# rm(list=ls())
## Script and functions to plot Magallon and Sanderson comparison
get.tail.probs <- function (table=table) {
  if(!is.data.frame(table)){
    stop("Input table is not a data.frame.") }
  
  #if(!c("taxa", "diversity", "node", "age_mean") %in% colnames(table)){
  #  stop("Bad formatted input file. Please make sure table has the same column names as the template table. See help page for more info.") }
  # organize table
  # table <- organize.table(table, clade, diversity, node, age_mean, age_upper, age_lower) # find a better way to pass arguments to nested functions
  background = which(table$taxa =="bg_clade")
  
  if (length(background)==0) {
    stop("No clade was defined as background. Make sure one of the clades is named 'bg_clade' in the 'clade' column of the input table.")
  }
  if (length(background)>1) {
    stop("More than one clade was defined as background. Make sure just one of the clades is named 'bg_clade' in the 'clade' column of the input table.")
  }
  if(which.max(table$diversity)!=background & which.max(table$age_mean)!=background) {
    stop("Background clade is not the oldest and/or the most diverse. Make sure to ajust that in the input table.")
  }
  
  background.clade <- table[background,]
  
  results <- list()
  # define vars
  eps.yule = 0
  eps.high.ext = 0.9
  k = background.clade$diversity
  time = background.clade$age_mean
  
  stem.results <- get.tail.prob.stem(eps.yule, eps.high.ext, k, time)
  crown.results <- get.tail.prob.crown(eps.yule, eps.high.ext, k, time)
  
  results[[1]] <- stem.results
  results[[2]] <- crown.results
  names(results) <- c("stem.results", "crown.results")
  
  return(results)
}
##
plottailprobs <- function (table = example, tail.probs = results) {

  # organize table
  #table <- organize.table(table, clade, diversity, node, age_mean, age_upper, age_lower) # find a better way to pass arguments to nested functions
  table=table
  tail.probs=tail.probs

  background = which(table$taxa =="bg_clade")

  if (length(background)==0) {
    stop("No clade was defined as background. Make sure one of the clades is named 'bg_clade' in the 'clade' column of master table.")
  }
  if (length(background)>1) {
    stop("More than one clade was defined as background. Make sure just one of the clades is named 'bg_clade' in the 'clade' column of master table.")
  }
  if(which.max(table$diversity)!=background & which.max(table$age_mean)!=background) {
    stop("Background clade is not the oldest and/or the most diverse. Make sure to ajust that in the master table.")
  }

  background.clade <- table[background,]

  results <- list()
  # define vars
  k = background.clade$diversity
  time = background.clade$age_mean

  par(mfrow=c(1,2))
  ###### stem plot
  stem.results <- tail.probs[[1]]
  k.lower.yule.stem <- stem.results[[1]]
  k.upper.yule.stem <- stem.results[[2]]
  k.lower.high.ext.stem <- stem.results[[3]]
  k.upper.high.ext.stem <- stem.results[[4]]

  plot(k.lower.yule.stem$time, log10(k.lower.yule.stem$k), type="l",
       col="red", lwd=1.5, ylim = c(0, round(log10(k), 0)+1), xlim = c(0,time+time*0.1), xlab="", ylab="", yaxt = "n")
  par(new=TRUE)
  plot(k.upper.yule.stem$time, log10(k.upper.yule.stem$k), type="l",
       col="red", lwd=1.5, ylim = c(0, round(log10(k), 0)+1), xlim = c(0,time+time*0.1), xlab="", ylab="", yaxt = "n")
  par(new=TRUE)
  plot(k.lower.high.ext.stem$time, log10(k.lower.high.ext.stem$k), type="l",lty=2,
       col="red",lwd=1.5, ylim = c(0, round(log10(k), 0)+1), xlim = c(0,time+time*0.1), xlab="", ylab="",  yaxt = "n")
  par(new=TRUE)
  plot(k.upper.high.ext.stem$time, log10(k.upper.high.ext.stem$k), type="l",lty=2,
       col="red",lwd=1.5, ylim = c(0, round(log10(k), 0)+1), xlim = c(0,time+time*0.1), xlab='Age of clade (million years)', ylab='K (expected number of species)',yaxt = "n")

  points(x=table$age_mean[table$taxa=="bg_clade"], y=log10(table$diversity[table$taxa=="bg_clade"]), pch=1, col="orange", type="p")
  points(x=table$age_mean[table$taxa!="bg_clade" & table$node=="SG"], y=log10(table$diversity[table$taxa!="bg_clade" & table$node=="SG"]), pch=16, type="p")

  axis.par <- c(0:round(log10(k), 0)+1)
  axis(2, at=axis.par, labels=axis.par, las=1)
  legend(1, tail(axis.par, 1), legend=c("eps = 0", "eps = 0.9"), col=c("red", "red"), lty=1:2, cex=1)

  title(main="stem node rates")


  ###### crown plot
  crown.results <- tail.probs[[2]]
  k.lower.yule.crown <- crown.results[[1]]
  k.upper.yule.crown <- crown.results[[2]]
  k.lower.high.ext.crown <- crown.results[[3]]
  k.upper.high.ext.crown <- crown.results[[4]]

  plot(k.lower.yule.crown$time, log10(k.lower.yule.crown$k), type="l",
       col="blue", lwd=1.5, ylim = c(0, round(log10(k), 0)+1), xlim = c(0,time+time*0.1), xlab="", ylab="", yaxt = "n")
  par(new=TRUE)
  plot(k.upper.yule.crown$time, log10(k.upper.yule.crown$k), type="l",
       col="blue", lwd=1.5, ylim = c(0, round(log10(k), 0)+1), xlim = c(0,time+time*0.1), xlab="", ylab="", yaxt = "n")
  par(new=TRUE)
  plot(k.lower.high.ext.crown$time, log10(k.lower.high.ext.crown$k), type="l",lty=3,
       col="blue",lwd=1.5, ylim = c(0, round(log10(k), 0)+1), xlim = c(0,time+time*0.1), xlab="", ylab="",  yaxt = "n")
  par(new=TRUE)
  plot(k.upper.high.ext.crown$time, log10(k.upper.high.ext.crown$k), type="l",lty=3,
       col="blue",lwd=1.5, ylim = c(0, round(log10(k), 0)+1), xlim = c(0,time+time*0.1), xlab='Age of clade (million years)', ylab='K (expected number of species)',yaxt = "n")
  par(new=TRUE)

  #add points for clades
  points(x=table$age_mean[table$taxa=="bg_clade"], y=log10(table$diversity[table$taxa=="bg_clade"]), pch=1, col="orange", type="p")
  points(x=table$age_mean[table$taxa!="bg_clade" & table$node=="CG"], y=log10(table$diversity[table$taxa!="bg_clade" & table$node=="CG"]), pch=16, type="p")

  axis.par <- c(0:round(log10(k), 0)+1)
  axis(2, at=axis.par, labels=axis.par, las=1)
  legend(1, tail(axis.par, 1), legend=c("eps = 0", "eps = 0.9"), col=c("blue", "blue"), lty=c(1,3), cex=1)
  title(main="crown node rates")

}
##
get.age.nodes <- function (phy) {
  phy <- ape::ladderize(phy)
  root.node <- length(phy$tip.label)+1
  seq.nodes <- phy$edge
  dists <- phy$edge.length
  res <- numeric(max(phy$edge))
  for (i in seq_len(nrow(seq.nodes))) {
    res[seq.nodes[i, 2]] <- res[seq.nodes[i,1]] + dists[i]
  }
  return(res)
}
##
get.bg.rate.tree <- function (phy, sf=NA, eps=NA) {
  if (!inherits(phy, "phylo")) {
    stop("phy should be an object of class \"phylo\".")
  }
  if(is.na(eps)){
    eps=0
    cat("No eps specified, assuming no extinction...\n\n")
  }
  if(is.na(sf)){
    sf=1
    cat("No sf specified, assuming complete sampling...\n\n")
  }

  phy=phy
  sf=sf
  eps=eps

  diversity <- length(phy$tip.label)/sf
  time <- get.age.nodes(phy)[1]
  bg.rate <- eq6_and_7(eps, diversity, time)

  return(bg.rate)
}

####
get.template <- function(output.dir = getwd(), save.file=TRUE) {

  template <- data.frame(taxa=c("bg_clade","Calycanthales","Asterales","Myrtales","Arecaceae_c","Arecaceae_s"),
                         diversity=c(262196, 10, 25996, 10782, 2780, 2780),
                         node=c("bg","SG","CG","CG","CG","SG"),
                         age_mean=c(132, 108.8, 28.75, 88.2, 77.4, 84),
                         stringsAsFactors = FALSE)
  template$age_upper <- template$age_mean + 5
  template$age_lower <- template$age_mean - 5

  # if(save.file) {
  #    write.csv(template, file=paste0(output.dir,"/template.csv"))
  #  }

  return(template)
}

# Internal function
get.tail.prob.crown <- function (eps.yule, eps.high.ext, k, time) {
  eps.yule=eps.yule
  eps.high.ext=eps.high.ext
  k=k
  time=time

  results <- list()

  bg.rate.crown.yule <- eq6_and_7(eps.yule, k, time, crown=TRUE) # eq 7
  bg.rate.crown.high.ext <- eq6_and_7(eps.high.ext, k, time, crown=TRUE) #

  # define base vars

  birth.rate.crown.yule = bg.rate.crown.yule / (1 - eps.yule)
  birth.rate.crown.high.ext = bg.rate.crown.high.ext / (1 - eps.high.ext)

  death.rate.crown.yule = eps.yule * birth.rate.crown.yule # always 0
  death.rate.crown.high.ext = eps.high.ext * birth.rate.crown.high.ext

  # get CI for expected diversity based on background rates from the stem node (eq 10)
  prob <- 0
  time.bins <- seq(1, time, by=1)
  k.bins <- seq(1, k, by=1)

  # yule
  k.upper.yule <- matrix (nrow = 0, ncol = 2)
  colnames(k.upper.yule) <- c("k", "time")
  k.lower.yule <- matrix (nrow = 0, ncol = 2)
  colnames(k.lower.yule) <- c("k", "time")

  cat("Estimating confidence intervals of expected species diversity according to age of crown group (yule): time bins 1 to", time,  "\n")
  for (time_index in time.bins) {
    cat("\r", time_index)
    for (k_index in k.bins-1) {
      prob <- CrownNtKmore(k_index, birth.rate.crown.yule, death.rate.crown.yule, time_index) # eq 10
      # prob <- CrownNtKmore(k_index, bg.rate.crown.yule, eps.yule, time_index) # eq 10
      if (0.025 == round(prob,3)) {
        k.lower.yule <- rbind(k.lower.yule, c(k_index, time_index)) }
      else if (0.975 == round(prob,3)) {
        k.upper.yule <- rbind(k.upper.yule, c(k_index, time_index)) }
    }
  }
  cat("","\n")

  # high extinction
  k.upper.high.ext <- matrix (nrow = 0, ncol = 2)
  colnames(k.upper.high.ext) <- c("k", "time")
  k.lower.high.ext <- matrix (nrow = 0, ncol = 2)
  colnames(k.lower.high.ext) <- c("k", "time")

  cat("Estimating confidence intervals of expected species diversity according to age of crown group (high ext): time bins 1 to", time, "\n")
  for (time_index in time.bins) {
    cat("\r", time_index)
    for (k_index in k.bins-1) {
      prob <- CrownNtKmore(k_index, birth.rate.crown.high.ext, death.rate.crown.high.ext, time_index) # eq 10
      # prob <- CrownNtKmore(k_index, bg.rate.crown.high.ext, eps.high.ext, time_index) # eq 10
      if (0.025 == round(prob,3)) {
        k.lower.high.ext <- rbind(k.lower.high.ext, c(k_index, time_index)) }
      else if (0.975 == round(prob,3)) {
        k.upper.high.ext <- rbind(k.upper.high.ext, c(k_index, time_index)) }
    }
  }
  cat("","\n")

  results[[1]] <- as.data.frame(k.lower.yule)
  results[[2]] <- as.data.frame(k.upper.yule)
  results[[3]] <- as.data.frame(k.lower.high.ext)
  results[[4]] <- as.data.frame(k.upper.high.ext)
  names(results) <- c("k.lower.yule", "k.upper.yule", "k.lower.high.ext", "k.upper.high.ext")

  return(results)
}

# Internal function
get.tail.prob.stem <- function (eps.yule, eps.high.ext, k, time) {
  eps.yule=eps.yule
  eps.high.ext=eps.high.ext
  k=k
  time=time
  results <- list()

  bg.rate.stem.yule <- eq6_and_7(eps.yule, k, time, crown=FALSE) # eq 6
  bg.rate.stem.high.ext <- eq6_and_7(eps.high.ext, k, time, crown=FALSE) # eq 6

  # define base vars
  death.rate.stem.yule = eps.yule * bg.rate.stem.yule/ (1 - eps.yule) # always 0
  death.rate.stem.high.ext = eps.high.ext * bg.rate.stem.high.ext / (1 - eps.high.ext)

  birth.rate.stem.yule = bg.rate.stem.yule / (1 - eps.yule)
  birth.rate.stem.high.ext  = bg.rate.stem.high.ext / (1 - eps.high.ext)

  # get CI for expected diversity based on background rates from the stem node (eq 10)
  prob <- 0
  time.bins <- seq(1, time, by=1)
  k.bins <- seq(1, k, by=1)

  # yule
  k.upper.yule <- matrix (nrow = 0, ncol = 2)
  colnames(k.upper.yule) <- c("k", "time")
  k.lower.yule <- matrix (nrow = 0, ncol = 2)
  colnames(k.lower.yule) <- c("k", "time")

  cat("Estimating confidence intervals of expected species diversity according to the age of stem group (yule): time bins 1 to", time, "\n")
  for (time_index in time.bins) {
    cat("\r", time_index)
    for (k_index in k.bins-1) {
      prob <- StemNtKmore(k_index, birth.rate.stem.yule, death.rate.stem.yule, time_index) # eq 10
      # prob <- StemNtKmore(k_index, bg.rate.stem.yule, eps.yule, time_index) # eq 10
      if (0.025 == round(prob,3)) {
        k.lower.yule <- rbind(k.lower.yule, c(k_index, time_index)) }
      else if (0.975 == round(prob,3)) {
        k.upper.yule <- rbind(k.upper.yule, c(k_index, time_index)) }
    }
  }
  cat("","\n")

  # high extinction
  k.upper.high.ext <- matrix (nrow = 0, ncol = 2)
  colnames(k.upper.high.ext) <- c("k", "time")
  k.lower.high.ext <- matrix (nrow = 0, ncol = 2)
  colnames(k.lower.high.ext) <- c("k", "time")

  cat("Estimating confidence intervals of expected species diversity according to age of the stem group (high ext): time bins 1 to", time, "\n" )
  for (time_index in time.bins) {
    cat("\r", time_index)
    for (k_index in k.bins-1) {
      prob <- StemNtKmore(k_index, birth.rate.stem.high.ext, death.rate.stem.high.ext, time_index) # eq 10
      # prob <- StemNtKmore(k_index, bg.rate.stem.high.ext, eps.high.ext, time_index) # eq 10
      if (0.025 == round(prob,3)) {
        k.lower.high.ext <- rbind(k.lower.high.ext, c(k_index, time_index)) }
      else if (0.975 == round(prob,3)) {
        k.upper.high.ext <- rbind(k.upper.high.ext, c(k_index, time_index)) }
    }
  }
  cat("","\n")

  results[[1]] <- as.data.frame(k.lower.yule)
  results[[2]] <- as.data.frame(k.upper.yule)
  results[[3]] <- as.data.frame(k.lower.high.ext)
  results[[4]] <- as.data.frame(k.upper.high.ext)
  names(results) <- c("k.lower.yule", "k.upper.yule", "k.lower.high.ext", "k.upper.high.ext")

  return(results)
}

#-------------------------
# Equations from paper:
# equations 6 and 7
eq6_and_7 <- function(eps, diversity, time, crown=TRUE) {
  if(crown == TRUE) {
    if(eps == 0) {
      net.div.mle <- (log(diversity) - log(2))/time # eq 4
    } else {
      net.div.mle <- 1/time*(log(diversity/2*(1-eps^2)+ 2*eps+1/2*(1-eps)*sqrt(diversity*(diversity*eps^2-8*eps+2*diversity*eps+diversity)))-log(2)) # eq 7
    }
  }else{
    if(eps == 0) {
      net.div.mle <- log(diversity)/time # eq 3
    } else {
      net.div.mle  <- 1/time*log(diversity*(1-eps)+eps) # eq 6
    }
  }
  return(net.div.mle)
}

#Eq. 1a.
ProbNt0 <- function(birth.rate, death.rate, time, crown=TRUE){
  if(crown == TRUE){
    a <- 2
  }else{
    a <- 1
  }
  p0 <- AlphaT(birth.rate=birth.rate, death.rate=death.rate, time=time) ^ a
  return(p0)
}

#Eq. 2a.
AlphaT <- function(birth.rate, death.rate, time){
  beta.t <- BetaT(birth.rate=birth.rate, death.rate=death.rate, time=time)
  alpha.t <- (death.rate/birth.rate) * beta.t
  return(alpha.t)
}

#Eq. 2b.
BetaT <- function(birth.rate, death.rate, time){
  net.diver.rate <- birth.rate - death.rate
  exprt <- exp(net.diver.rate * time)
  beta.t <- (exprt - 1) / (exprt - (death.rate/birth.rate))
  return(beta.t)
}

#Eq 10a
StemNtKmore <- function(k, birth.rate, death.rate, time){
  beta.t <- BetaT(birth.rate=birth.rate, death.rate=death.rate, time=time)
  probNgeNtax <- beta.t ^ (k-1)
  return(probNgeNtax)
}

#Eq 11a
CrownNtKmore <- function(k, birth.rate, death.rate, time) {
  alpha.t <- AlphaT(birth.rate=birth.rate, death.rate=death.rate, time=time)
  beta.t <- BetaT(birth.rate=birth.rate, death.rate=death.rate, time=time)
  #probNgeNtax <- ((beta.t^(k-2))/ (1+alpha.t))* ((k * (1 - alpha.t - beta.t + (alpha.t*beta.t)) + alpha.t + (2*beta.t) - 1))
  probNgeNtax <- (beta.t^(k-2))*(k*(1 - alpha.t - beta.t + alpha.t*beta.t) + alpha.t + 2*beta.t -1)/(1+alpha.t)
  #print(probNgeNtax)
  #probNgeNtax <- (beta.t^(k-2))*(k*(1 - alpha.t - beta.t + alpha.t*beta.t) + alpha.t + 2*beta.t -1)/(1 - alpha.t + 2*alpha.t*beta.t)
  #print(probNgeNtax)
  return(probNgeNtax)
}


#####################################
#####################################
# example code

example <- get.template(save.file=F)

######## PLOT FIG 3 & 4 ######

results <- get.tail.probs(table=table)
# save(results, file="test.Rsave")

pdf(file="plot.tail.probs.pdf", height=6, width=8)
plottailprobs(example, results)
dev.off()


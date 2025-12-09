#Aug 2025, process sim results for 15 fold scen
#For 40 traps and a large regular area
#version 2 of results is after overriding the default that uses true values as start values

library(secr)
library(secrdesign)

################
#Grid designs
################

######################
#1st data, then ests
######################
load("15FoldScen/Cluster/Sims/GridResults240.RData")

#summarise a CH
#note the "group" covariate is the relevant one, not "sex"
summary(Grid.800.results[[1]]$output[[1]][[1]])

#function to extract data summaries
CH.data.summs <- function(res.obj, nreps){
  grid.data.summaries <- matrix(numeric(0), nrow = nreps, ncol = 8, dimnames = list(c(1:100),
            c("#inds 1", "#inds 2","spr 1", "spr 2", "sing 1", "sing 2", "#dets 1", "#dets 2")))
  
  for (i in 1:nreps){
    ch <- res.obj$output[[i]][[1]]
    S1 <- subset(ch, attr(ch, "covariates")$group == "S1")
    S2 <- subset(ch, attr(ch, "covariates")$group == "S2")
    n.S1 <- dim(S1)[1]
    n.S2 <- dim(S2)[1]
    dets.S1 <- sum(S1)
    dets.S2 <- sum(S2)
    ind.tots.S1 <- apply(S1, 1, sum, na.rm = TRUE)
    ind.tots.S2 <- apply(S2, 1, sum, na.rm = TRUE)
    
    recaps.S1 <- sum(apply(S1, 1, function(x){sum(x > 0)})>1)
    recaps.S2 <- sum(apply(S2, 1, function(x){sum(x > 0)})>1)
    
    singletons.S1 <- n.S1 - recaps.S1
    singletons.S2 <- n.S2 - recaps.S2

    grid.data.summaries[i,] <- c(n.S1, n.S2, recaps.S1, recaps.S2, singletons.S1, singletons.S2, dets.S1, dets.S2)
  }
  return(grid.data.summaries)
} 

#function to extract and summarise estimates for two strata
#can choose to produce boxplots, ylim.ceiling needs two values (1 for each strata)
#median also shown in blue
summ.results <- function(results.obj, par = "D", true.values, numests, ylim.ceiling, mean = TRUE, plot = TRUE, par.est = TRUE){
  S1.ests <- matrix(,nrow = numests, ncol = 2)
  S2.ests <- matrix(,nrow = numests, ncol = 2)
  colnames(S1.ests) <- c("Estimates", "Std errors")
  colnames(S2.ests) <- c("Estimates", "Std errors")
  
  for (i in 1:numests){
    if (par == "D"){
      S1.ests[i,1] <- results.obj$output[[i]][[1]][[1]][1,2]
      S2.ests[i,1] <- results.obj$output[[i]][[1]][[2]][1,2]
      S1.ests[i,2] <- results.obj$output[[i]][[1]][[1]][1,3]
      S2.ests[i,2] <- results.obj$output[[i]][[1]][[2]][1,3]
      plot.label1 <- "Strata 1 Density"
      plot.label2 <- "Strata 2 Density"
    } else {
      if (par=="Sigma"){
        S1.ests[i,1] <- results.obj$output[[i]][[1]][[1]][3,2]
        S2.ests[i,1] <- results.obj$output[[i]][[1]][[2]][3,2]
        S1.ests[i,2] <- results.obj$output[[i]][[1]][[1]][3,3]
        S2.ests[i,2] <- results.obj$output[[i]][[1]][[2]][3,3]
        plot.label1 <- "Strata 1 Sigma"
        plot.label2 <- "Strata 2 Sigma"
      } else {
        S1.ests[i,1] <- results.obj$output[[i]][[1]][[1]][2,2]
        S2.ests[i,1] <- results.obj$output[[i]][[1]][[2]][2,2]
        S1.ests[i,2] <- results.obj$output[[i]][[1]][[1]][2,3]
        S2.ests[i,2] <- results.obj$output[[i]][[1]][[2]][2,3]
        plot.label1 <- "Strata 1 Lambda0"
        plot.label2 <- "Strata 2 Lambda0"
      }
    }
  }
  
  #calculate RB and RSE
  if (mean == TRUE){
    S1.RB <- (mean(S1.ests, na.rm = T) - true.values[1]) / true.values[1] 
    S2.RB <- (mean(S2.ests, na.rm = T) - true.values[2]) / true.values[2]
  } else {
    S1.RB <- (median(S1.ests, na.rm = T) - true.values[1]) / true.values[1] 
    S2.RB <- (median(S2.ests, na.rm = T) - true.values[2]) / true.values[2]
  }
  
  S1.RSE <- sd(S1.ests, na.rm = T) / mean(S1.ests, na.rm = T) 
  S2.RSE <- sd(S2.ests, na.rm = T) / mean(S2.ests, na.rm = T)  
  
  #produce side by side boxplots
  par(mfrow=c(1,2))
  if (plot==TRUE){
    if (par.est == TRUE){
      boxplot(S1.ests[,1], main = paste(plot.label1, "estimates", sep =" "), ylim = c(0, ylim.ceiling[1]))
      abline(h = true.values[1], col = 'red')
      abline(h = median(S1.ests[,1]), col = 'blue')
      
      boxplot(S2.ests[,1], main = paste(plot.label2, "estimates", sep =" "), ylim = c(0, ylim.ceiling[2]))
      abline(h = true.values[2], col = 'red')
      abline(h = median(S2.ests[,1]), col = 'blue')
    } else {
      boxplot(S1.ests[,2], main = paste(plot.label1, "std errors", sep =" "), ylim = c(0, ylim.ceiling[1]))
      boxplot(S2.ests[,2], main = paste(plot.label2, "std errors", sep =" "), ylim = c(0, ylim.ceiling[2]))
    }
  }
  
  res <- list("Rel bias" = c(S1.RB, S2.RB), "Rel SE" = c(S1.RSE, S2.RSE), "Strata 1 ests" = S1.ests, "Strata 2 ests" = S2.ests)
  return(res)
}

#another version that extracts all estimates together for each strata
#will plot and summ in another function
extract.results <- function(results.obj, numests){
  S1.ests <- matrix(,nrow = numests, ncol = 6)
  S2.ests <- matrix(,nrow = numests, ncol = 6)
  colnames(S1.ests) <- c("D1", "D1.se", "L01", "L01.se", "Sigma1", "S1.se")
  colnames(S2.ests) <- c("D2", "D2.se", "L02", "L02.se", "Sigma2", "S2.se")
  
  for (i in 1:numests){
      S1.ests[i,] <- as.vector(t(as.matrix(results.obj$output[[i]][[1]][[1]][,2:3])))
      S2.ests[i,] <- as.vector(t(as.matrix(results.obj$output[[i]][[1]][[2]][,2:3])))
  }
  return(list("S1 ests" = S1.ests, "S2 ests" = S2.ests))
}

#new fn to plot and summarise
#obj to pass is a set of ests and ses for one parameter
summ.plot <- function(results.df, plot.label = "Test", true.value, numests, ylim.ceiling, median = FALSE, plot = TRUE, par.est = TRUE){
  RB <- NULL ; RSE <- NULL ; mean <- NULL   
  if (par.est == TRUE){
    mean <- mean(results.df[,1], na.rm = T)
    RSE <- sd(results.df[,1], na.rm = T) / mean(results.df[,1], na.rm = T)
    if (median == FALSE){
      RB <- (mean(results.df[,1], na.rm = T) - true.value) / true.value
    } else {
      RB <- (median(results.df[,1], na.rm = T) - true.value) / true.value
    }
  } else {
    mean <- mean(results.df[,2], na.rm = T)
    RSE <- sd(results.df[,2], na.rm = T) / mean(results.df[,2], na.rm = T)
  }

  #produce side by side boxplots
  if (plot==TRUE){
    if (par.est == TRUE){
      boxplot(results.df[,1], main = paste(plot.label, "estimates", sep =" "), ylim = c(0, ylim.ceiling))
      abline(h = true.value, col = 'red')
      abline(h = median(results.df[,1]), col = 'blue')
    } else {
      boxplot(results.df[,2], main = paste(plot.label, "std errors", sep =" "), ylim = c(0, ylim.ceiling))
    }
  }
  
  res <- list("Mean" = mean,"Rel bias" = RB, "Rel SE" = RSE)
  return(res)
}

#function to find rogue values for a strata
#true values provided as D, L0, sigma
#detects NA and infinite values, 
#and also checks for any estimates greater or less than mag fold 
find.rogue <- function(df, mag = 10, true){
  bad.rows <- NULL
  df <- as.data.frame(df)
  colnames(df) <- c("D", "D.se", "L0", "L0.se", "Sigma", "Sigma.se")
  
  ##identify rows with NA or infinite
  rogue.rows <- which(is.na(rowSums(df))|is.infinite(rowSums(df)))
  if (length(rogue.rows)>0) df <- df[-rogue.rows,]
  
  #identify rows with wild values for any estimates (mag fold)
  D.rogue <- which(df$D > true[1]*mag | df$D < true[1]/mag)
  L0.rogue <- which(df$L0 > true[2]*mag | df$L0 < true[2]/mag)
  Sigma.rogue <- which(df$Sigma > true[3]*mag | df$Sigma < true[3]/mag)
  
  if (length(D.rogue)>0) bad.rows <- c(bad.rows, D.rogue)
  if (length(L0.rogue)>0) bad.rows <- c(bad.rows, L0.rogue)
  if (length(Sigma.rogue)>0) bad.rows <- c(bad.rows, Sigma.rogue)
  if (!is.null(bad.rows)) df <- df[-bad.rows,]
  
  return("Clean data" = df)
}

Grid.800.data <- CH.data.summs(Grid.800.results[[1]], nreps = 100)
Grid.1600.data <- CH.data.summs(Grid.1600.results[[1]], nreps = 100)

pdf("15FoldScen/GridDataSumms.pdf", height = 8, width = 10, pointsize = 11)
par(mfrow = c(2,1))
boxplot(Grid.800.data)
boxplot(Grid.1600.data)
dev.off()

########################################
#Grid estimates
#the results obj has both 800m and 1600m spacing
########################################

Grid.800.Ests <- extract.results(Grid.800.results[[2]], numests = 100)

#filter data, NAs and infinite, and wildly out values (20 fold)
Grid.800.S1.red <- find.rogue(Grid.800.Ests$`S1 ests`, mag = 10, true = c(0.05, 2, 200))
Grid.800.S2.red <- find.rogue(Grid.800.Ests$`S2 ests`, mag = 10, true = c(0.05/15, 2/15, 3000))

summary(Grid.800.S1.red)
summary(Grid.800.S2.red)

#there is one row (50) with a very large se, the D estimate is just within the 10 fold diff
Grid.800.S2.red <- Grid.800.S2.red[-50,]

summary(Grid.800.S2.red)

#plot and summarise
#estimates first
summ.plot(Grid.800.S1.red[,1:2], plot.label = "Strata 1 Density", true.value = 0.05, numests = 98, 
          ylim.ceiling = 0.15, median = FALSE, plot = TRUE, par.est = TRUE)
summ.plot(Grid.800.S2.red[,1:2], plot.label = "Strata 2 Density", true.value = 0.05/15, numests = 93, 
          ylim.ceiling = 0.01, median = FALSE, plot = TRUE, par.est = TRUE)

summ.plot(Grid.800.S1.red[,3:4], plot.label = "Strata 1 L0", true.value = 2, numests = 98, 
          ylim.ceiling = 3.2, median = FALSE, plot = TRUE, par.est = TRUE)
summ.plot(Grid.800.S2.red[,3:4], plot.label = "Strata 2 L0", true.value = 2/15, numests = 93, 
          ylim.ceiling = 0.21, median = FALSE, plot = TRUE, par.est = TRUE)

summ.plot(Grid.800.S1.red[,5:6], plot.label = "Strata 1 Sigma", true.value = 200, numests = 98, 
          ylim.ceiling = 250, median = FALSE, plot = TRUE, par.est = TRUE)
summ.plot(Grid.800.S2.red[,5:6], plot.label = "Strata 2 Sigma", true.value = 3000, numests = 93, 
          ylim.ceiling = 30000, median = FALSE, plot = TRUE, par.est = TRUE)

#precision
summ.plot(Grid.800.S1.red[,1:2], plot.label = "Strata 1 D", numests = 98, 
          ylim.ceiling = 0.1, median = FALSE, plot = TRUE, par.est = FALSE)
summ.plot(Grid.800.S2.red[,1:2], plot.label = "Strata 2 Density", numests = 93, 
          ylim.ceiling = 0.01, median = FALSE, plot = TRUE, par.est = FALSE)

summ.plot(Grid.800.S1.red[,3:4], plot.label = "Strata 1 L0", numests = 98, 
          ylim.ceiling = 0.65, median = FALSE, plot = TRUE, par.est = FALSE)
summ.plot(Grid.800.S2.red[,3:4], plot.label = "Strata 2 L0", numests = 93, 
          ylim.ceiling = 0.05, median = FALSE, plot = TRUE, par.est = FALSE)

summ.plot(Grid.800.S1.red[,5:6], plot.label = "Strata 1 Sigma", numests = 98, 
          ylim.ceiling = 45, median = FALSE, plot = TRUE, par.est = FALSE)
summ.plot(Grid.800.S2.red[,5:6], plot.label = "Strata 2 Sigma", numests = 93, 
          ylim.ceiling = 15000, median = FALSE, plot = TRUE, par.est = FALSE)


################
#1600 m spacing
################

Grid.1600.Ests <- extract.results(Grid.1600.results[[2]], numests = 100)

#filter data, NAs and infinite, and wildly out values (20 fold)
Grid.1600.S1.red <- find.rogue(Grid.1600.Ests$`S1 ests`, mag = 10, true = c(0.05, 2, 200))
Grid.1600.S2.red <- find.rogue(Grid.1600.Ests$`S2 ests`, mag = 10, true = c(0.05/15, 2/15, 3000))

summary(Grid.1600.S1.red)
summary(Grid.1600.S2.red)

#there is one row (1) with a large se, the D estimate is not that far off
Grid.1600.S1.red <- Grid.1600.S1.red[-1,]

summary(Grid.1600.S1.red)

#plot and summarise
#estimates first
summ.plot(Grid.1600.S1.red[,1:2], plot.label = "Strata 1 Density", true.value = 0.05, numests = 96, 
          ylim.ceiling = 0.20, median = FALSE, plot = TRUE, par.est = TRUE)
summ.plot(Grid.1600.S2.red[,1:2], plot.label = "Strata 2 Density", true.value = 0.05/15, numests = 100, 
          ylim.ceiling = 0.01, median = FALSE, plot = TRUE, par.est = TRUE)

summ.plot(Grid.1600.S1.red[,3:4], plot.label = "Strata 1 L0", true.value = 2, numests = 96, 
          ylim.ceiling = 3.3, median = FALSE, plot = TRUE, par.est = TRUE)
summ.plot(Grid.1600.S2.red[,3:4], plot.label = "Strata 2 L0", true.value = 2/15, numests = 100, 
          ylim.ceiling = 0.25, median = FALSE, plot = TRUE, par.est = TRUE)

summ.plot(Grid.1600.S1.red[,5:6], plot.label = "Strata 1 Sigma", true.value = 200, numests = 96, 
          ylim.ceiling = 260, median = FALSE, plot = TRUE, par.est = TRUE)
summ.plot(Grid.1600.S2.red[,5:6], plot.label = "Strata 2 Sigma", true.value = 3000, numests = 100, 
          ylim.ceiling = 4500, median = FALSE, plot = TRUE, par.est = TRUE)

#precision
summ.plot(Grid.1600.S1.red[,1:2], plot.label = "Strata 1 D", numests = 96, 
          ylim.ceiling = 15, median = FALSE, plot = TRUE, par.est = FALSE)
summ.plot(Grid.1600.S2.red[,1:2], plot.label = "Strata 2 Density", numests = 100, 
          ylim.ceiling = 0.0015, median = FALSE, plot = TRUE, par.est = FALSE)

summ.plot(Grid.1600.S1.red[,3:4], plot.label = "Strata 1 L0", numests = 96, 
          ylim.ceiling = 2.6, median = FALSE, plot = TRUE, par.est = FALSE)
summ.plot(Grid.1600.S2.red[,3:4], plot.label = "Strata 2 L0", numests = 100, 
          ylim.ceiling = 0.05, median = FALSE, plot = TRUE, par.est = FALSE)

summ.plot(Grid.1600.S1.red[,5:6], plot.label = "Strata 1 Sigma", numests = 96, 
          ylim.ceiling = 1000, median = FALSE, plot = TRUE, par.est = FALSE)
summ.plot(Grid.1600.S2.red[,5:6], plot.label = "Strata 2 Sigma", numests = 100, 
          ylim.ceiling = 1200, median = FALSE, plot = TRUE, par.est = FALSE)

#######################################################################



################
#GA4 designs
################

######################
#GA4 with F sigma
######################

######################
#GA4 with M sigma
######################

######################
#GA4 with Avg sigma
######################

######################
#GA4 with both
######################

######################
#GA4 with both (max = T)
######################

###############################################################
################
#GA5 designs
################

######################
#GA5 with F sigma
######################

######################
#GA5 with M sigma
######################

######################
#GA5 with Avg sigma
######################

######################
#GA5 with both
######################

#use function to extract, summ and plot

######################
#GA5 with both (max = T)
######################

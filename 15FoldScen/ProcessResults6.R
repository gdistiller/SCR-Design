#June 2026, process sim results for 15 fold scen
#modified June to save actual results not just summaries
#For 40 traps and 120 traps and a large regular area
#this version to process results from fitting models outside secrdesign 
#uses 500 reps
#also filters out reps where se/est > 100

library(secr)
library(secrdesign)

################
#Functions
################

#function to extract data summaries
#input is a list of CHs
CH.data.summs <- function(res.list){
  nreps <- length(res.list)
  grid.data.summaries <- matrix(numeric(0), nrow = nreps, ncol = 8, dimnames = list(c(1:nreps),
            c("#inds 1", "#inds 2","spr 1", "spr 2", "sing 1", "sing 2", "#dets 1", "#dets 2")))
  
  for (i in 1:nreps){
    ch <- res.list[[i]][[1]]
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

#function to summarise estimates for two strata
#can choose to produce boxplots, 
#ylim.ceiling needs two values (1 for each strata)
#true values also needs the true values for the chosen parameter for each strata
#now stores all results
summ.results <- function(results.df, par = "D", true.values, ylim.ceiling, plot = TRUE, se = FALSE, rel.se = TRUE){
    if (par == "D"){
      index <- 1
      plot.label1 <- "Strata 1 Density"
      plot.label2 <- "Strata 2 Density"
    } else {
      if (par=="Sigma"){
        index <- 5
        plot.label1 <- "Strata 1 Sigma"
        plot.label2 <- "Strata 2 Sigma"
      } else {
        index <- 3
        plot.label1 <- "Strata 1 Lambda0"
        plot.label2 <- "Strata 2 Lambda0"
      }
    }
  index <- index + 1 #this is needed because meging the dfs leads to the 1st column being the row number / index
  #calculate RB for est before checking se argument
  S1.RB.mean <- (mean(results.df[,index], na.rm = T) - true.values[1]) / true.values[1] 
  S2.RB.mean <- (mean(results.df[,index+6], na.rm = T) - true.values[2]) / true.values[2]
  S1.RB.median <- (median(results.df[,index], na.rm = T) - true.values[1]) / true.values[1] 
  S2.RB.median <- (median(results.df[,index+6], na.rm = T) - true.values[2]) / true.values[2]

  #produce side by side boxplots
  par(mfrow=c(1,2))
  if (plot==TRUE){
    if (se == FALSE){
      boxplot(results.df[,index], main = paste(plot.label1, "estimates", sep =" "), ylim = c(0, ylim.ceiling[1]))
      abline(h = true.values[1], col = 'red')
      boxplot(results.df[,index+6], main = paste(plot.label2, "estimates", sep =" "), ylim = c(0, ylim.ceiling[2]))
      abline(h = true.values[2], col = 'red')
    } else {
      index <- index + 1
      if (rel.se == TRUE){
        boxplot(results.df[,index] / mean(results.df[,index], na.rm = T), main = paste(plot.label1, "Rel SE", sep = " "), ylim = c(0, ylim.ceiling[1]))
        boxplot(results.df[,index+6] / mean(results.df[,index+6], na.rm = T), main = paste(plot.label2, "Rel SE", sep = " "), ylim = c(0, ylim.ceiling[2]))
      } else {
        boxplot(results.df[,index], main = paste(plot.label1, "SE", sep = " "), ylim = c(0, ylim.ceiling[1]))
        boxplot(results.df[,index+6], main = paste(plot.label2, "SE", sep = " "), ylim = c(0, ylim.ceiling[2]))
      }
    }
  }
  
  #Calc RSE
  S1.RSE <- sd(results.df[,index], na.rm = T) / mean(results.df[,index], na.rm = T) 
  S2.RSE <- sd(results.df[,index+6], na.rm = T) / mean(results.df[,index+6], na.rm = T) 
  
  res <- list("Results" = results.df,"Rel bias (mean)" = c(S1.RB.mean, S2.RB.mean), "Rel bias (median)" = c(S1.RB.median, S2.RB.median), "Rel SE" = c(S1.RSE, S2.RSE))
  return(res)
}

#reduced function to just plot from the df of results
#merging the dfs leads to the 1st column being the row number so need + 1 to index
#true values needed for plots of estimates
plot.results <- function(results.df, par, se = T, rel.se = T, ylim.ceiling, true.values = NULL){
  if (par == "D"){
    index <- 2
    plot.label1 <- "Strata 1 Density"
    plot.label2 <- "Strata 2 Density"
  } else {
    if (par=="Sigma"){
      index <- 6
      plot.label1 <- "Strata 1 Sigma"
      plot.label2 <- "Strata 2 Sigma"
    } else {
      index <- 4
      plot.label1 <- "Strata 1 Lambda0"
      plot.label2 <- "Strata 2 Lambda0"
    }
  }

  #produce side by side boxplots
  #if rel.se = T it divides the estimates by the mean
  par(mfrow=c(1,2))
  if (se == FALSE){
    boxplot(results.df[,index], main = paste(plot.label1, "estimates", sep =" "), ylim = c(0, ylim.ceiling[1]))
    abline(h = true.values[1], col = 'red')
    boxplot(results.df[,index+6], main = paste(plot.label2, "estimates", sep =" "), ylim = c(0, ylim.ceiling[2]))
    abline(h = true.values[2], col = 'red')
  } else {
    index <- index + 1
    if (rel.se == TRUE){
      boxplot(results.df[,index] / mean(results.df[,index], na.rm = T), main = paste(plot.label1, "Rel SE", sep = " "), ylim = c(0, ylim.ceiling[1]))
      boxplot(results.df[,index+6] / mean(results.df[,index+6], na.rm = T), main = paste(plot.label2, "Rel SE", sep = " "), ylim = c(0, ylim.ceiling[2]))
    } else {
      boxplot(results.df[,index], main = paste(plot.label1, "SE", sep = " "), ylim = c(0, ylim.ceiling[1]))
      boxplot(results.df[,index+6], main = paste(plot.label2, "SE", sep = " "), ylim = c(0, ylim.ceiling[2]))
    }
  }
}


#function to find rogue values for a single strata
#cannot do both together as diff strata cant be in the same row
#true values provided as D, L0, sigma
#detects NA and infinite values, 
#and also checks for any estimates greater or less than mag fold 
#diff version of the fn, more efficient
find.rogue <- function(df, mag = 10, true){
  bad.rows <- NULL
  df <- as.data.frame(df)
  colnames(df) <- c("D", "D.se", "L0", "L0.se", "Sig", "Sig.se")
  
  ##identify rows with NA or infinite
  inds1 <- which(is.na(rowSums(df))|is.infinite(rowSums(df)))
  
  ##identify rows where the relative se is v large relative to the est
  RSE.D <- df[,2] / df[,1] 
  RSE.L0 <- df[,4] / df[,3] 
  RSE.Sig <- df[,6] / df[,5]
  
  inds2 <- which(RSE.D > 100 | RSE.L0 > 100 | RSE.Sig > 100)
  
  #identify rows with wild values for any estimates (mag fold)
  inds3 <- which(df$D > true[1]*mag | df$D < true[1]/mag | df$L0 > true[2]*mag | df$L0 < true[2]/mag | df$Sig > true[3]*mag | df$Sig < true[3]/mag)  
  
  bad.rows <- c(inds1, inds2, inds3)
  
  if (length(bad.rows)==0) bad.rows <- NULL
  
  if (!is.null(bad.rows)) {
    bad.rows <- sort(unique(bad.rows)) 
    df <- df[-bad.rows,]
  }
  
  df$row <- row.names(df)
  return("Clean data" = df)
}

####################################################

mag.factor <- 5 
D1.ylim <- 0.25 ; D2.ylim <- 0.015
L1.ylim <- 3.5 ; L2.ylim <- 0.25
Sig1.ylim <- 350 ; Sig2.ylim <- 9000

################
#Grid designs
#now includes grids with optimal spacing 2G, and clustered grids
#first for 40 traps
################

load("15FoldScen/Cluster/Sims/40Traps/GridsdResults40.RData")

###############
#800 m spacing
###############

##first data summary, then estimates
Grid40.800.data.summ <- CH.data.summs(Grid.800.Data$output)

pdf("15FoldScen/Grid800DataSumms.pdf", height = 8, width = 10, pointsize = 11)
boxplot(Grid40.800.data.summ)

#filter data, NAs and infinite, and values out by at least 5 fold
#one strata at a time
Grid40.800.red1 <- find.rogue(Grid.800.results[,c(1:6)], mag = mag.factor, true = c(0.05, 2, 200))
Grid40.800.red2 <- find.rogue(Grid.800.results[,c(7:12)], mag = mag.factor, true = c(0.05/15, 2/15, 3000))

Grid40.800.clean <- merge(Grid40.800.red1, Grid40.800.red2, by = "row", all.x = TRUE, all.y = TRUE)
summary(Grid.800.clean)

#plot and summarise
#estimates first
Grid40.800.D.est <- summ.results(Grid40.800.clean, par = "D", true.value = c(0.05,0.05/15) , ylim.ceiling = c(D1.ylim, D2.ylim), plot = TRUE, se = FALSE)
Grid40.800.L0.est <- summ.results(Grid40.800.clean, par = "L0", true.value = c(2, 2/15) , ylim.ceiling = c(L1.ylim, L2.ylim), plot = TRUE, se = FALSE)
Grid40.800.Sig.est <- summ.results(Grid40.800.clean, par = "Sigma", true.value = c(200, 3000) , ylim.ceiling = c(Sig1.ylim, Sig2.ylim), plot = TRUE, se = FALSE)

#relative standard errors
Grid40.800.D.se <- summ.results(Grid40.800.clean, par = "D", true.value = c(0.05,0.05/15) , ylim.ceiling = c(10,10), plot = TRUE, se = TRUE, rel.se = TRUE)
Grid40.800.L0.se <- summ.results(Grid40.800.clean, par = "L0", true.value = c(2, 2/15) , ylim.ceiling = c(10,10), plot = TRUE, se = TRUE, rel.se = TRUE)
Grid40.800.Sig.se <- summ.results(Grid40.800.clean, par = "Sigma", true.value = c(200, 3000) , ylim.ceiling = c(10,10), plot = TRUE, se = TRUE, rel.se = TRUE)

################
#Optimal Spacing
#using OS of 700m
#was initially run with 600 so redone -> OSGridsResults40b
################

load("15FoldScen/Cluster/Sims/40Traps/OSGridsResults40b.RData")

##first data summary, then estimates
Grid40.700.data.summ <- CH.data.summs(Grid.700.Data$output)

pdf("15FoldScen/Grid700DataSumms.pdf", height = 8, width = 10, pointsize = 11)
boxplot(Grid40.700.data.summ)
dev.off()

#filter data, NAs and infinite, and wildly out values (20 fold)
#one strata at a time
Grid40.700.red1 <- find.rogue(Grid.700.results[,c(1:6)], mag = mag.factor, true = c(0.05, 2, 200))
Grid40.700.red2 <- find.rogue(Grid.700.results[,c(7:12)], mag = mag.factor, true = c(0.05/15, 2/15, 3000))

Grid40.700.clean <- merge(Grid40.700.red1, Grid40.700.red2, by = "row", all.x = TRUE, all.y = TRUE)
summary(Grid40.700.clean)

#plot and summarise
#estimates first
Grid40.700.D.est <- summ.results(Grid40.700.clean, par = "D", true.values = c(0.05,0.05/15) , ylim.ceiling = c(D1.ylim, D2.ylim), plot = TRUE, se = FALSE)
Grid40.700.L0.est <- summ.results(Grid40.700.clean, par = "L0", true.value = c(2, 2/15) , ylim.ceiling = c(L1.ylim, L2.ylim), plot = TRUE, se = FALSE)
Grid40.700.Sig.est <- summ.results(Grid40.700.clean, par = "Sigma", true.value = c(200, 3000) , ylim.ceiling = c(Sig1.ylim, Sig2.ylim), plot = TRUE, se = FALSE)

#rel standard errors
Grid40.700.D.se <- summ.results(Grid40.700.clean, par = "D", true.value = c(0.05,0.05/15) , ylim.ceiling = c(10, 10), plot = TRUE, se = TRUE, rel.se = TRUE)
Grid40.700.L0.se <- summ.results(Grid40.700.clean, par = "L0", true.value = c(2, 2/15) , ylim.ceiling = c(10, 10), plot = TRUE, se = TRUE, rel.se = TRUE)
Grid40.700.Sig.se <- summ.results(Grid40.700.clean, par = "Sigma", true.value = c(200, 3000) , ylim.ceiling = c(10, 10), plot = TRUE, se = TRUE, rel.se = TRUE)

################
#Clustered designs
#three configs tried
#using 2 configs: OS spacings and 2 * sigma
################

#note that used to load Clusters1.RData but reran the OS spacings in May (to ensure spacing of 700) -> Clusters1b.RData
load("15FoldScen/Cluster/Sims/40Traps/Clusters1b.RData")

#first trying OS that gave 600 / 1200
##first data summary, then estimates
Cluster40.os.data.summ <- CH.data.summs(Grid.os.Data$output)

boxplot(Cluster40.os.data.summ)

#filter data, NAs and infinite, and wildly out values (20 fold)
#one strata at a time
Cluster40.os.red1 <- find.rogue(Grid.os.results[,c(1:6)], mag = mag.factor, true = c(0.05, 2, 200))
Cluster40.os.red2 <- find.rogue(Grid.os.results[,c(7:12)], mag = mag.factor, true = c(0.05/15, 2/15, 3000))

Cluster40.os.clean <- merge(Cluster40.os.red1, Cluster40.os.red2, by = "row", all.x = TRUE, all.y = TRUE)
summary(Cluster40.os.clean)

#plot and summarise
#estimates first
Cluster40.os.D.est <- summ.results(Cluster40.os.clean, par = "D", true.value = c(0.05,0.05/15) , ylim.ceiling = c(D1.ylim, D2.ylim), plot = TRUE, se = FALSE)
Cluster40.os.L0.est <- summ.results(Cluster40.os.clean, par = "L0", true.value = c(2, 2/15) , ylim.ceiling = c(L1.ylim, L2.ylim), plot = TRUE, se = FALSE)
Cluster40.os.Sig.est <- summ.results(Cluster40.os.clean, par = "Sigma", true.value = c(200, 3000) , ylim.ceiling = c(Sig1.ylim, Sig2.ylim), plot = TRUE, se = FALSE)

#rel standard errors
Cluster40.os.D.se <- summ.results(Cluster40.os.clean, par = "D", true.value = c(0.05,0.05/15) , ylim.ceiling = c(10, 10), plot = TRUE, se = TRUE, rel.se = TRUE)
Cluster40.os.L0.se <- summ.results(Cluster40.os.clean, par = "L0", true.value = c(2, 2/15) , ylim.ceiling = c(10, 10), plot = TRUE, se = TRUE, rel.se = TRUE)
Cluster40.os.Sig.se <- summ.results(Cluster40.os.clean, par = "Sigma", true.value = c(200, 3000) , ylim.ceiling = c(10, 10), plot = TRUE, se = TRUE, rel.se = TRUE)

#now trying 2 sigma spacing
load("15FoldScen/Cluster/Sims/40Traps/Clusters2.RData")

##first data summary, then estimates
Cluster40.2sig.data.summ <- CH.data.summs(Grid.2sig.Data$output)

boxplot(Cluster40.2sig.data.summ)

#filter data, NAs and infinite, and wildly out values (20 fold)
#one strata at a time
Cluster40.2sig.red1 <- find.rogue(Grid.2sig.results[,c(1:6)], mag = mag.factor, true = c(0.05, 2, 200))
Cluster40.2sig.red2 <- find.rogue(Grid.2sig.results[,c(7:12)], mag = mag.factor, true = c(0.05/15, 2/15, 3000))

Cluster40.2sig.clean <- merge(Cluster40.2sig.red1, Cluster40.2sig.red2, by = "row", all.x = TRUE, all.y = TRUE)
summary(Cluster40.2sig.clean)

#plot and summarise
#estimates first
Cluster40.2sig.D.est <- summ.results(Cluster40.2sig.clean, par = "D", true.value = c(0.05,0.05/15) , ylim.ceiling = c(D1.ylim, D2.ylim), plot = TRUE, se = FALSE)
Cluster40.2sig.L0.est <- summ.results(Cluster40.2sig.clean, par = "L0", true.value = c(2, 2/15) , ylim.ceiling = c(L1.ylim, L2.ylim), plot = TRUE, se = FALSE)
Cluster40.2sig.Sig.est <- summ.results(Cluster40.2sig.clean, par = "Sigma", true.value = c(200, 3000) , ylim.ceiling = c(Sig1.ylim, Sig2.ylim), plot = TRUE, se = FALSE)

#standard errors
Cluster40.2sig.D.se <- summ.results(Cluster40.2sig.clean, par = "D", true.value = c(0.05,0.05/15) , ylim.ceiling = c(10, 10), plot = TRUE, se = TRUE, rel.se = TRUE)
Cluster40.2sig.L0.se <- summ.results(Cluster40.2sig.clean, par = "L0", true.value = c(2, 2/15) , ylim.ceiling = c(10, 10), plot = TRUE, se = TRUE, rel.se = TRUE)
Cluster40.2sig.Sig.se <- summ.results(Cluster40.2sig.clean, par = "Sigma", true.value = c(200, 3000) , ylim.ceiling = c(10, 10), plot = TRUE, se = TRUE, rel.se = TRUE)

#######################################################################

################
#GA4 designs
################

#extract and collate results

load("15FoldScen/Cluster/Sims/40Traps/GA4dResults40.RData")

######################
#GA4 with S1 pars
######################

GA4.40.S1.data.summ <- CH.data.summs(G4.S1.Data$output)
boxplot(GA4.40.S1.data.summ)

#filter estimates, NAs and infinite, and wildly out values (10 fold)
#one strata at a time
GA4.40.S1.red1 <- find.rogue(G4.S1.results[,c(1:6)], mag = mag.factor, true = c(0.05, 2, 200))
GA4.40.S1.red2 <- find.rogue(G4.S1.results[,c(7:12)], mag = mag.factor, true = c(0.05/15, 2/15, 3000))

GA4.40.S1.clean <- merge(GA4.40.S1.red1, GA4.40.S1.red2, by = "row", all.x = TRUE, all.y = TRUE)
summary(GA4.S1.clean)

#plot and summarise
#estimates first
GA4.40.S1.D.est <- summ.results(GA4.40.S1.clean, par = "D", true.value = c(0.05,0.05/15) , ylim.ceiling = c(D1.ylim, D2.ylim), plot = TRUE, se = FALSE)
GA4.40.S1.L0.est <- summ.results(GA4.40.S1.clean, par = "L0", true.value = c(2, 2/15) , ylim.ceiling = c(L1.ylim, L2.ylim), plot = TRUE, se = FALSE)
GA4.40.S1.Sig.est <- summ.results(GA4.40.S1.clean, par = "Sigma", true.value = c(200, 3000) , ylim.ceiling = c(Sig1.ylim, Sig2.ylim), plot = TRUE, se = FALSE)

#standard errors
GA4.40.S1.D.se <- summ.results(GA4.40.S1.clean, par = "D", true.value = c(0.05,0.05/15) , ylim.ceiling = c(10, 10), plot = TRUE, se = TRUE, rel.se = TRUE)
GA4.40.S1.L0.se <- summ.results(GA4.40.S1.clean, par = "L0", true.value = c(2, 2/15) , ylim.ceiling = c(10, 10), plot = TRUE, se = TRUE, rel.se = TRUE)
GA4.40.S1.Sig.se <- summ.results(GA4.40.S1.clean, par = "Sigma", true.value = c(200, 3000) , ylim.ceiling = c(10, 10), plot = TRUE, se = TRUE, rel.se = TRUE)

######################
#GA4 with S2 pars
######################

GA4.40.S2.data.summ <- CH.data.summs(G4.S2.Data$output)
boxplot(GA4.40.S2.data.summ)

#filter estimates, NAs and infinite, and wildly out values (10 fold)
#one strata at a time
GA4.40.S2.red1 <- find.rogue(G4.S2.results[,c(1:6)], mag = mag.factor, true = c(0.05, 2, 200))
GA4.40.S2.red2 <- find.rogue(G4.S2.results[,c(7:12)], mag = mag.factor, true = c(0.05/15, 2/15, 3000))

GA4.40.S2.clean <- merge(GA4.40.S2.red1, GA4.40.S2.red2, by = "row", all.x = TRUE, all.y = TRUE)
summary(GA4.S2.clean)

#plot and summarise
#estimates first
GA4.40.S2.D.est <- summ.results(GA4.40.S2.clean, par = "D", true.value = c(0.05,0.05/15) , ylim.ceiling = c(D1.ylim, D2.ylim), plot = TRUE, se = FALSE)
GA4.40.S2.L0.est <- summ.results(GA4.40.S2.clean, par = "L0", true.value = c(2, 2/15) , ylim.ceiling = c(L1.ylim, L2.ylim), plot = TRUE, se = FALSE)
GA4.40.S2.Sig.est <- summ.results(GA4.40.S2.clean, par = "Sigma", true.value = c(200, 3000) , ylim.ceiling = c(Sig1.ylim, Sig2.ylim), plot = TRUE, se = FALSE)

#standard errors
GA4.40.S2.D.se <- summ.results(GA4.40.S2.clean, par = "D", true.value = c(0.05,0.05/15) , ylim.ceiling = c(10, 10), plot = TRUE, se = TRUE, rel.se = TRUE)
GA4.40.S2.L0.se <- summ.results(GA4.40.S2.clean, par = "L0", true.value = c(2, 2/15) , ylim.ceiling = c(10, 10), plot = TRUE, se = TRUE, rel.se = TRUE)
GA4.40.S2.Sig.se <- summ.results(GA4.40.S2.clean, par = "Sigma", true.value = c(200, 3000) , ylim.ceiling = c(10, 10), plot = TRUE, se = TRUE, rel.se = TRUE)

######################
#GA4 with avg pars
######################

GA4.40.Avg.data.summ <- CH.data.summs(G4.Avg.Data$output)
boxplot(GA4.Avg.data.summ)

#filter estimates, NAs and infinite, and wildly out values (10 fold)
#one strata at a time
GA4.40.Avg.red1 <- find.rogue(G4.Avg.results[,c(1:6)], mag = mag.factor, true = c(0.05, 2, 200))
GA4.40.Avg.red2 <- find.rogue(G4.Avg.results[,c(7:12)], mag = mag.factor, true = c(0.05/15, 2/15, 3000))

GA4.40.Avg.clean <- merge(GA4.40.Avg.red1, GA4.40.Avg.red2, by = "row", all.x = TRUE, all.y = TRUE)
summary(GA4.Avg.clean)

#plot and summarise
#estimates first
GA4.40.Avg.D.est <- summ.results(GA4.40.Avg.clean, par = "D", true.value = c(0.05,0.05/15) , ylim.ceiling = c(D1.ylim, D2.ylim), plot = TRUE, se = FALSE)
GA4.40.Avg.L0.est <- summ.results(GA4.40.Avg.clean, par = "L0", true.value = c(2, 2/15) , ylim.ceiling = c(L1.ylim, L2.ylim), plot = TRUE, se = FALSE)
GA4.40.Avg.Sig.est <- summ.results(GA4.40.Avg.clean, par = "Sigma", true.value = c(200, 3000) , ylim.ceiling = c(Sig1.ylim, Sig2.ylim), plot = TRUE, se = FALSE)

#standard errors
GA4.40.Avg.D.se <- summ.results(GA4.40.Avg.clean, par = "D", true.value = c(0.05,0.05/15) , ylim.ceiling = c(10, 10), plot = TRUE, se = TRUE, rel.se = TRUE)
GA4.40.Avg.L0.se <- summ.results(GA4.40.Avg.clean, par = "L0", true.value = c(2, 2/15) , ylim.ceiling = c(10, 10), plot = TRUE, se = TRUE, rel.se = TRUE)
GA4.40.Avg.Sig.se <- summ.results(GA4.40.Avg.clean, par = "Sigma", true.value = c(200, 3000) , ylim.ceiling = c(10, 10), plot = TRUE, se = TRUE, rel.se = TRUE)

######################
#GA4 with both
######################

GA4.40.Both.data.summ <- CH.data.summs(G4.Both.Data$output)
boxplot(GA4.Both.data.summ)

#filter estimates, NAs and infinite, and wildly out values (10 fold)
#one strata at a time
GA4.40.Both.red1 <- find.rogue(G4.Both.results[,c(1:6)], mag = mag.factor, true = c(0.05, 2, 200))
GA4.40.Both.red2 <- find.rogue(G4.Both.results[,c(7:12)], mag = mag.factor, true = c(0.05/15, 2/15, 3000))

GA4.40.Both.clean <- merge(GA4.40.Both.red1, GA4.40.Both.red2, by = "row", all.x = TRUE, all.y = TRUE)
summary(GA4.Both.clean)

#plot and summarise
#estimates first
GA4.40.Both.D.est <- summ.results(GA4.40.Both.clean, par = "D", true.value = c(0.05,0.05/15) , ylim.ceiling = c(D1.ylim, D2.ylim), plot = TRUE, se = FALSE)
GA4.40.Both.L0.est <- summ.results(GA4.40.Both.clean, par = "L0", true.value = c(2, 2/15) , ylim.ceiling = c(L1.ylim, L2.ylim), plot = TRUE, se = FALSE)
GA4.40.Both.Sig.est <- summ.results(GA4.40.Both.clean, par = "Sigma", true.value = c(200, 3000) , ylim.ceiling = c(Sig1.ylim, Sig2.ylim), plot = TRUE, se = FALSE)

#standard errors
GA4.40.Both.D.se <- summ.results(GA4.40.Both.clean, par = "D", true.value = c(0.05,0.05/15) , ylim.ceiling = c(10, 10), plot = TRUE, se = TRUE, rel.se = TRUE)
GA4.40.Both.L0.se <- summ.results(GA4.40.Both.clean, par = "L0", true.value = c(2, 2/15) , ylim.ceiling = c(10, 10), plot = TRUE, se = TRUE, rel.se = TRUE)
GA4.40.Both.Sig.se <- summ.results(GA4.40.Both.clean, par = "Sigma", true.value = c(200, 3000) , ylim.ceiling = c(10, 10), plot = TRUE, se = TRUE, rel.se = TRUE)

######################
#GA4 with both (max = T)
######################

GA4.40.BothMax.data.summ <- CH.data.summs(G4.BothMax.Data$output)
boxplot(GA4.BothMax.data.summ)

#filter estimates, NAs and infinite, and wildly out values (10 fold)
#one strata at a time
GA4.40.BothMax.red1 <- find.rogue(G4.BothMax.results[,c(1:6)], mag = mag.factor, true = c(0.05, 2, 200))
GA4.40.BothMax.red2 <- find.rogue(G4.BothMax.results[,c(7:12)], mag = mag.factor, true = c(0.05/15, 2/15, 3000))

GA4.40.BothMax.clean <- merge(GA4.40.BothMax.red1, GA4.40.BothMax.red2, by = "row", all.x = TRUE, all.y = TRUE)
summary(GA4.BothMax.clean)

#plot and summarise
#estimates first
GA4.40.BothMax.D.est <- summ.results(GA4.40.BothMax.clean, par = "D", true.value = c(0.05,0.05/15) , ylim.ceiling = c(D1.ylim, D2.ylim), plot = TRUE, se = FALSE)
GA4.40.BothMax.L0.est <- summ.results(GA4.40.BothMax.clean, par = "L0", true.value = c(2, 2/15) , ylim.ceiling = c(L1.ylim, L2.ylim), plot = TRUE, se = FALSE)
GA4.40.BothMax.Sig.est <- summ.results(GA4.40.BothMax.clean, par = "Sigma", true.value = c(200, 3000) , ylim.ceiling = c(Sig1.ylim, Sig2.ylim), plot = TRUE, se = FALSE)

#standard errors
GA4.40.BothMax.D.se <- summ.results(GA4.40.BothMax.clean, par = "D", true.value = c(0.05,0.05/15) , ylim.ceiling = c(10, 10), plot = TRUE, se = TRUE, rel.se = TRUE)
GA4.40.BothMax.L0.se <- summ.results(GA4.40.BothMax.clean, par = "L0", true.value = c(2, 2/15) , ylim.ceiling = c(10, 10), plot = TRUE, se = TRUE, rel.se = TRUE)
GA4.40.BothMax.Sig.se <- summ.results(GA4.40.BothMax.clean, par = "Sigma", true.value = c(200, 3000) , ylim.ceiling = c(10, 10), plot = TRUE, se = TRUE, rel.se = TRUE)

################################
#GA5 designs
################################

load("15FoldScen/Cluster/Sims/40Traps/GA5dResults40.RData")

######################
#GA5 with S1 pars
######################

GA5.40.S1.data.summ <- CH.data.summs(G5.S1.Data$output)
boxplot(GA5.40.S1.data.summ)

#filter estimates, NAs and infinite, and wildly out values (10 fold)
#one strata at a time
GA5.40.S1.red1 <- find.rogue(G5.S1.results[,c(1:6)], mag = mag.factor, true = c(0.05, 2, 200))
GA5.40.S1.red2 <- find.rogue(G5.S1.results[,c(7:12)], mag = mag.factor, true = c(0.05/15, 2/15, 3000))

GA5.40.S1.clean <- merge(GA5.40.S1.red1, GA5.40.S1.red2, by = "row", all.x = TRUE, all.y = TRUE)
summary(GA5.40.S1.clean)

#plot and summarise
#estimates first
GA5.40.S1.D.est <- summ.results(GA5.40.S1.clean, par = "D", true.value = c(0.05,0.05/15) , ylim.ceiling = c(D1.ylim, D2.ylim), plot = TRUE, se = FALSE)
GA5.40.S1.L0.est <- summ.results(GA5.40.S1.clean, par = "L0", true.value = c(2, 2/15) , ylim.ceiling = c(L1.ylim, L2.ylim), plot = TRUE, se = FALSE)
GA5.40.S1.Sig.est <- summ.results(GA5.40.S1.clean, par = "Sigma", true.value = c(200, 3000) , ylim.ceiling = c(Sig1.ylim, Sig2.ylim), plot = TRUE, se = FALSE)

#rel standard errors
GA5.40.S1.D.se <- summ.results(GA5.40.S1.clean, par = "D", true.value = c(0.05,0.05/15) , ylim.ceiling = c(10, 10), plot = TRUE, se = TRUE, rel.se = TRUE)
GA5.40.S1.L0.se <- summ.results(GA5.40.S1.clean, par = "L0", true.value = c(2, 2/15) , ylim.ceiling = c(10, 10), plot = TRUE, se = TRUE, rel.se = TRUE)
GA5.40.S1.Sig.se <- summ.results(GA5.40.S1.clean, par = "Sigma", true.value = c(200, 3000) , ylim.ceiling = c(10, 10), plot = TRUE, se = TRUE, rel.se = TRUE)

######################
#GA5 with S2 pars
######################

GA5.40.S2.data.summ <- CH.data.summs(G5.S2.Data$output)
boxplot(GA5.40.S2.data.summ)

#filter estimates, NAs and infinite, and wildly out values (10 fold)
#one strata at a time
GA5.40.S2.red1 <- find.rogue(G5.S2.results[,c(1:6)], mag = mag.factor, true = c(0.05, 2, 200))
GA5.40.S2.red2 <- find.rogue(G5.S2.results[,c(7:12)], mag = mag.factor, true = c(0.05/15, 2/15, 3000))

GA5.40.S2.clean <- merge(GA5.40.S2.red1, GA5.40.S2.red2, by = "row", all.x = TRUE, all.y = TRUE)
summary(GA5.40.S2.clean)

#plot and summarise
#estimates first
GA5.40.S2.D.est  <- summ.results(GA5.40.S2.clean, par = "D", true.value = c(0.05,0.05/15) , ylim.ceiling = c(D1.ylim, D2.ylim), plot = TRUE, se = FALSE)
GA5.40.S2.L0.est  <- summ.results(GA5.40.S2.clean, par = "L0", true.value = c(2, 2/15) , ylim.ceiling = c(L1.ylim, L2.ylim), plot = TRUE, se = FALSE)
GA5.40.S2.Sig.est  <- summ.results(GA5.40.S2.clean, par = "Sigma", true.value = c(200, 3000) , ylim.ceiling = c(Sig1.ylim, Sig2.ylim), plot = TRUE, se = FALSE)

#rel standard errors
GA5.40.S2.D.se  <- summ.results(GA5.40.S2.clean, par = "D", true.value = c(0.05,0.05/15) , ylim.ceiling = c(10, 10), plot = TRUE, se = TRUE, rel.se = TRUE)
GA5.40.S2.L0.se  <- summ.results(GA5.40.S2.clean, par = "L0", true.value = c(2, 2/15) , ylim.ceiling = c(10, 10), plot = TRUE, se = TRUE, rel.se = TRUE)
GA5.40.S2.Sig.se  <- summ.results(GA5.40.S2.clean, par = "Sigma", true.value = c(200, 3000) , ylim.ceiling = c(10, 10), plot = TRUE, se = TRUE, rel.se = TRUE)

######################
#GA5 with avg pars
######################

GA5.40.Avg.data.summ <- CH.data.summs(G5.Avg.Data$output)
boxplot(GA5.40.40.Avg.data.summ)

#filter estimates, NAs and infinite, and wildly out values (10 fold)
#one strata at a time
GA5.40.Avg.red1 <- find.rogue(G5.Avg.results[,c(1:6)], mag = mag.factor, true = c(0.05, 2, 200))
GA5.40.Avg.red2 <- find.rogue(G5.Avg.results[,c(7:12)], mag = mag.factor, true = c(0.05/15, 2/15, 3000))

GA5.40.Avg.clean <- merge(GA5.40.Avg.red1, GA5.40.Avg.red2, by = "row", all.x = TRUE, all.y = TRUE)
summary(GA5.40.Avg.clean)

#plot and summarise
#estimates first
GA5.40.Avg.D.est  <- summ.results(GA5.40.Avg.clean, par = "D", true.value = c(0.05,0.05/15) , ylim.ceiling = c(D1.ylim, D2.ylim), plot = TRUE, se = FALSE)
GA5.40.Avg.L0.est  <- summ.results(GA5.40.Avg.clean, par = "L0", true.value = c(2, 2/15) , ylim.ceiling = c(L1.ylim, L2.ylim), plot = TRUE, se = FALSE)
GA5.40.Avg.Sig.est  <- summ.results(GA5.40.Avg.clean, par = "Sigma", true.value = c(200, 3000) , ylim.ceiling = c(Sig1.ylim, Sig2.ylim), plot = TRUE, se = FALSE)

#rel standard errors
GA5.40.Avg.D.se <- summ.results(GA5.40.Avg.clean, par = "D", true.value = c(0.05,0.05/15) , ylim.ceiling = c(10, 10), plot = TRUE, se = TRUE, rel.se = TRUE)
GA5.40.Avg.L0.se <- summ.results(GA5.40.Avg.clean, par = "L0", true.value = c(2, 2/15) , ylim.ceiling = c(10, 10), plot = TRUE, se = TRUE, rel.se = TRUE)
GA5.40.Avg.Sig.se <- summ.results(GA5.40.Avg.clean, par = "Sigma", true.value = c(200, 3000) , ylim.ceiling = c(10, 10), plot = TRUE, se = TRUE, rel.se = TRUE)

######################
#GA5 with both
######################

GA5.40.Both.data.summ <- CH.data.summs(G5.Both.Data$output)
boxplot(GA5.40.Both.data.summ)

#filter estimates, NAs and infinite, and wildly out values (10 fold)
#one strata at a time
GA5.40.Both.red1 <- find.rogue(G5.Both.results[,c(1:6)], mag = mag.factor, true = c(0.05, 2, 200))
GA5.40.Both.red2 <- find.rogue(G5.Both.results[,c(7:12)], mag = mag.factor, true = c(0.05/15, 2/15, 3000))

GA5.40.Both.clean <- merge(GA5.40.Both.red1, GA5.40.Both.red2, by = "row", all.x = TRUE, all.y = TRUE)
summary(GA5.40.Both.clean)

#plot and summarise
#estimates first
GA5.40.Both.D.est <- summ.results(GA5.40.Both.clean, par = "D", true.value = c(0.05,0.05/15) , ylim.ceiling = c(D1.ylim, D2.ylim), plot = TRUE, se = FALSE)
GA5.40.Both.L0.est <- summ.results(GA5.40.Both.clean, par = "L0", true.value = c(2, 2/15) , ylim.ceiling = c(L1.ylim, L2.ylim), plot = TRUE, se = FALSE)
GA5.40.Both.Sig.est <- summ.results(GA5.40.Both.clean, par = "Sigma", true.value = c(200, 3000) , ylim.ceiling = c(Sig1.ylim, Sig2.ylim), plot = TRUE, se = FALSE)

#rel standard errors
GA5.40.Both.D.se <- summ.results(GA5.40.Both.clean, par = "D", true.value = c(0.05,0.05/15) , ylim.ceiling = c(10, 10), plot = TRUE, se = TRUE, rel.se = TRUE)
GA5.40.Both.L0.se <- summ.results(GA5.40.Both.clean, par = "L0", true.value = c(2, 2/15) , ylim.ceiling = c(10, 10), plot = TRUE, se = TRUE, rel.se = TRUE)
GA5.40.Both.Sig.se <- summ.results(GA5.40.Both.clean, par = "Sigma", true.value = c(200, 3000) , ylim.ceiling = c(10, 10), plot = TRUE, se = TRUE, rel.se = TRUE)

######################
#GA5 with both (max = T)
######################

GA5.40.BothMax.data.summ <- CH.data.summs(G5.BothMax.Data$output)
boxplot(GA5.40.BothMax.data.summ)

#filter estimates, NAs and infinite, and wildly out values (10 fold)
#one strata at a time
GA5.40.BothMax.red1 <- find.rogue(G5.BothMax.results[,c(1:6)], mag = mag.factor, true = c(0.05, 2, 200))
GA5.40.BothMax.red2 <- find.rogue(G5.BothMax.results[,c(7:12)], mag = mag.factor, true = c(0.05/15, 2/15, 3000))

GA5.40.BothMax.clean <- merge(GA5.40.BothMax.red1, GA5.40.BothMax.red2, by = "row", all.x = TRUE, all.y = TRUE)
summary(GA5.40.BothMax.clean)

#plot and summarise
#estimates first
GA5.40.BothMax.D.est <- summ.results(GA5.40.BothMax.clean, par = "D", true.value = c(0.05,0.05/15) , ylim.ceiling = c(D1.ylim, D2.ylim), plot = TRUE, se = FALSE)
GA5.40.BothMax.L0.est <- summ.results(GA5.40.BothMax.clean, par = "L0", true.value = c(2, 2/15) , ylim.ceiling = c(L1.ylim, L2.ylim), plot = TRUE, se = FALSE)
GA5.40.BothMax.Sig.est <- summ.results(GA5.40.BothMax.clean, par = "Sigma", true.value = c(200, 3000) , ylim.ceiling = c(Sig1.ylim, Sig2.ylim), plot = TRUE, se = FALSE)

#rel standard errors
GA5.40.BothMax.D.se <- summ.results(GA5.40.BothMax.clean, par = "D", true.value = c(0.05,0.05/15) , ylim.ceiling = c(10, 10), plot = TRUE, se = TRUE, rel.se = TRUE)
GA5.40.BothMax.L0.se <- summ.results(GA5.40.BothMax.clean, par = "L0", true.value = c(2, 2/15) , ylim.ceiling = c(10, 10), plot = TRUE, se = TRUE, rel.se = TRUE)
GA5.40.BothMax.Sig.se <- summ.results(GA5.40.BothMax.clean, par = "Sigma", true.value = c(200, 3000) , ylim.ceiling = c(10, 10), plot = TRUE, se = TRUE, rel.se = TRUE)

############################################################################################
#Lacework

load("15FoldScen/Cluster/Sims/40Traps/LWResults40.RData")

LW.40.Data.reformatted <- lapply(LW.40.Data$output$`1`, function(x) list(x))

LW40.data.summ <- CH.data.summs(LW.40.Data.reformatted)
boxplot(LW40.data.summ)

#filter estimates, NAs and infinite, and wildly out values (10 fold)
#one strata at a time
LW40.red1 <- find.rogue(LW.40.results[,c(1:6)], mag = mag.factor, true = c(0.05, 2, 200))
LW40.red2 <- find.rogue(LW.40.results[,c(7:12)], mag = mag.factor, true = c(0.05/15, 2/15, 3000))

LW40.clean <- merge(LW40.red1, LW40.red2, by = "row", all.x = TRUE, all.y = TRUE)
summary(LW40.clean)

#plot and summarise
#estimates first
LW40.D.est <- summ.results(LW40.clean, par = "D", true.value = c(0.05,0.05/15) , ylim.ceiling = c(D1.ylim, D2.ylim), plot = TRUE, se = FALSE)
LW40.L0.est <- summ.results(LW40.clean, par = "L0", true.value = c(2, 2/15) , ylim.ceiling = c(L1.ylim, L2.ylim), plot = TRUE, se = FALSE)
LW40.Sig.est <- summ.results(LW40.clean, par = "Sigma", true.value = c(200, 3000) , ylim.ceiling = c(Sig1.ylim, Sig2.ylim), plot = TRUE, se = FALSE)

#rel standard errors
LW40.D.se <- summ.results(LW40.clean, par = "D", true.value = c(0.05,0.05/15) , ylim.ceiling = c(10, 10), plot = TRUE, se = TRUE, rel.se = TRUE)
LW40.L0.se <- summ.results(LW40.clean, par = "L0", true.value = c(2, 2/15) , ylim.ceiling = c(10, 10), plot = TRUE, se = TRUE, rel.se = TRUE)
LW40.Sig.se <- summ.results(LW40.clean, par = "Sigma", true.value = c(200, 3000) , ylim.ceiling = c(10, 10), plot = TRUE, se = TRUE, rel.se = TRUE)

##########################################################
#save all results objects into a list for easy loading
##########################################################

datasumms.40 <- list("Grid40.700" = Grid40.700.data.summ, "Grid40.800" = Grid40.800.data.summ,
                  "Clust40.OS" = Cluster40.os.data.summ, "Clust40.2sig" = Cluster40.2sig.data.summ, "G4.40.S1" = GA4.40.S1.data.summ, 
                  "G4.40.S2" = GA4.40.S2.data.summ, "G4.40.Avg" = GA4.40.Avg.data.summ, "G4.40.Both" = GA4.40.Both.data.summ, "G4.40.BothMax" = GA4.40.BothMax.data.summ,
                  "G5.40.S1" = GA5.40.S1.data.summ, "G5.40.S2" = GA5.40.S2.data.summ, "G5.40.Avg" = GA5.40.Avg.data.summ, "G5.40.Both" = GA5.40.Both.data.summ, 
                  "G5.40.BothMax" = GA5.40.BothMax.data.summ, "LW40" = LW40.data.summ)
save(datasumms.40, file = "15FoldScen/DataSumms40Traps500.RData")

simresults.40 <- list("Grid40.700.D" = Grid40.700.D.est, "Grid40.700.L0" = Grid40.700.L0.est, "Grid40.700.Sig" = Grid40.700.Sig.est,
                      "Grid40.800.D" = Grid40.800.D.est, "Grid40.800.L0" = Grid40.800.L0.est, "Grid40.800.Sig" = Grid40.800.Sig.est,
                      "Clust40.OS.D" = Cluster40.os.D.est, "Clust40.OS.L0" = Cluster40.os.L0.est, "Clust40.OS.Sig" = Cluster40.os.Sig.est,
                      "Clust40.2sig.D" = Cluster40.2sig.D.est, "Clust40.2sig.L0" = Cluster40.2sig.L0.est, "Clust40.2sig.Sig" = Cluster40.2sig.Sig.est,
                      "GA4.40.S1.D" = GA4.40.S1.D.est, "GA4.40.S1.L0" = GA4.40.S1.L0.est, "GA4.40.S1.Sig" = GA4.40.S1.Sig.est,
                      "GA4.40.S2.D" = GA4.40.S2.D.est, "GA4.40.S2.L0" = GA4.40.S2.L0.est, "GA4.40.S2.Sig" = GA4.40.S2.Sig.est,
                      "GA4.40.Avg.D" = GA4.40.Avg.D.est, "GA4.40.Avg.L0" = GA4.40.Avg.L0.est, "GA4.40.Avg.Sig" = GA4.40.Avg.Sig.est,
                      "GA4.40.Both.D" = GA4.40.Both.D.est, "GA4.40.Both.L0" = GA4.40.Both.L0.est, "GA4.40.Both.Sig" = GA4.40.Both.Sig.est,
                      "GA4.40.BothMax.D" = GA4.40.BothMax.D.est, "GA4.40.BothMax.L0" = GA4.40.BothMax.L0.est, "GA4.40.BothMax.Sig" = GA4.40.BothMax.Sig.est,
                      "GA5.40.S1.D" = GA5.40.S1.D.est, "GA5.40.S1.L0" = GA5.40.S1.L0.est, "GA5.40.S1.Sig" = GA5.40.S1.Sig.est,
                      "GA5.40.S2.D" = GA5.40.S2.D.est, "GA5.40.S2.L0" = GA5.40.S2.L0.est, "GA5.40.S2.Sig" = GA5.40.S2.Sig.est,
                      "GA5.40.Avg.D" = GA5.40.Avg.D.est, "GA5.40.Avg.L0" = GA5.40.Avg.L0.est, "GA5.40.Avg.Sig" = GA5.40.Avg.Sig.est,
                      "GA5.40.Both.D" = GA5.40.Both.D.est, "GA5.40.Both.L0" = GA5.40.Both.L0.est, "GA5.40.Both.Sig" = GA5.40.Both.Sig.est,
                      "GA5.40.BothMax.D" = GA5.40.BothMax.D.est, "GA5.40.BothMax.L0" = GA5.40.BothMax.L0.est, "GA5.40.BothMax.Sig" = GA5.40.BothMax.Sig.est,
                      "LW40.D" = LW40.D.est, "LW40.L0" = LW40.L0.est, "LW40.Sig" = LW40.Sig.est)
save(simresults.40, file = "15FoldScen/SimResults40Traps500.RData")

seresults.40 <- list("Grid40.700.D" = Grid40.700.D.se, "Grid40.700.L0" = Grid40.700.L0.se, "Grid40.700.Sig" = Grid40.700.Sig.se,
                     "Grid40.800.D" = Grid40.800.D.se, "Grid40.800.L0" = Grid40.800.L0.se, "Grid40.800.Sig" = Grid40.800.Sig.se,
                     "Clust40.OS.D" = Cluster40.os.D.se, "Clust40.OS.L0" = Cluster40.os.L0.se, "Clust40.OS.Sig" = Cluster40.os.Sig.se,
                     "Clust40.2sig.D" = Cluster40.2sig.D.se, "Clust40.2sig.L0" = Cluster40.2sig.L0.se, "Clust40.2sig.Sig" = Cluster40.2sig.Sig.se,
                     "GA4.40.S1.D" = GA4.40.S1.D.se, "GA4.40.S1.L0" = GA4.40.S1.L0.se, "GA4.40.S1.Sig" = GA4.40.S1.Sig.se,
                     "GA4.40.S2.D" = GA4.40.S2.D.se, "GA4.40.S2.L0" = GA4.40.S2.L0.se, "GA4.40.S2.Sig" = GA4.40.S2.Sig.se,
                     "GA4.40.Avg.D" = GA4.40.Avg.D.se, "GA4.40.Avg.L0" = GA4.40.Avg.L0.se, "GA4.40.Avg.Sig" = GA4.40.Avg.Sig.se,
                     "GA4.40.Both.D" = GA4.40.Both.D.se, "GA4.40.Both.L0" = GA4.40.Both.L0.se, "GA4.40.Both.Sig" = GA4.40.Both.Sig.se,
                     "GA4.40.BothMax.D" = GA4.40.BothMax.D.se, "GA4.40.BothMax.L0" = GA4.40.BothMax.L0.se, "GA4.40.BothMax.Sig" = GA4.40.BothMax.Sig.se,
                     "GA5.40.S1.D" = GA5.40.S1.D.se, "GA5.40.S1.L0" = GA5.40.S1.L0.se, "GA5.40.S1.Sig" = GA5.40.S1.Sig.se,
                     "GA5.40.S2.D" = GA5.40.S2.D.se, "GA5.40.S2.L0" = GA5.40.S2.L0.se, "GA5.40.S2.Sig" = GA5.40.S2.Sig.se,
                     "GA5.40.Avg.D" = GA5.40.Avg.D.se, "GA5.40.Avg.L0" = GA5.40.Avg.L0.se, "GA5.40.Avg.Sig" = GA5.40.Avg.Sig.se,
                     "GA5.40.Both.D" = GA5.40.Both.D.se, "GA5.40.Both.L0" = GA5.40.Both.L0.se, "GA5.40.Both.Sig" = GA5.40.Both.Sig.se,
                     "GA5.40.BothMax.D" = GA5.40.BothMax.D.se, "GA5.40.BothMax.L0" = GA5.40.BothMax.L0.se, "GA5.40.BothMax.Sig" = GA5.40.BothMax.Sig.se,
                     "LW40.D" = LW40.D.se, "LW40.L0" = LW40.L0.se, "LW40.Sig" = LW40.Sig.se)

save(seresults.40, file = "15FoldScen/seResults40Traps500.RData")

excluded.40 <- list("Grid40.700.S1" = Grid40.700.red1, "Grid40.700.S2" = Grid40.700.red2,
                 "Grid40.800.S1" = Grid40.800.red1, "Grid40.800.S2" = Grid40.800.red2,
                 "Clust40.os.S1" = Cluster40.os.red1, "Clust40.os.S2" = Cluster40.os.red2,
                 "Clust40.2sig.S1" = Cluster40.2sig.red1, "Clust40.2sig.S2" = Cluster40.2sig.red2,
                 "GA4.40.S1.S1" = GA4.40.S1.red1, "GA4.40.S1.S2" = GA4.40.S1.red2,
                 "GA4.40.S2.S1" = GA4.40.S2.red1, "GA4.40.S2.S2" = GA4.40.S2.red2,
                 "GA4.40.Avg.S1" = GA4.40.Avg.red1, "GA4.40.Avg.S2" = GA4.40.Avg.red2,
                 "GA4.40.Both.S1" = GA4.40.Both.red1, "GA4.40.Both.S2" = GA4.40.Both.red2,
                 "GA4.40.BothMax.S1" = GA4.40.BothMax.red1, "GA4.40.BothMax.S2" = GA4.40.BothMax.red2,
                 "GA5.40.S1.S1" = GA5.40.S1.red1, "GA5.40.S1.S2" = GA5.40.S1.red2,
                 "GA5.40.S2.S1" = GA5.40.S2.red1, "GA5.40.S2.S2" = GA5.40.S2.red2,
                 "GA5.40.Avg.S1" = GA5.40.Avg.red1, "GA5.40.Avg.S2" = GA5.40.Avg.red2,
                 "GA5.40.Both.S1" = GA5.40.Both.red1, "GA5.40.Both.S2" = GA5.40.Both.red2,
                 "GA5.40.BothMax.S1" = GA5.40.BothMax.red1, "GA5.40.BothMax.S2" = GA5.40.BothMax.red2,
                 "LW40.S1" = LW40.red1, "LW40.S2" = LW40.red2)
save(excluded.40, file = "15FoldScen/Excluded40Traps500.RData")

##############################################################################################################
################
#120 traps
################

################
#Grid designs
################

load("15FoldScen/Cluster/Sims/120Traps/Grids800Results120.RData")

###############
#800 m spacing
###############

##first data summary, then estimates
Grid120.800.data.summ <- CH.data.summs(Grid.800.Data$output)

pdf("15FoldScen/Grid800.120.DataSumms.pdf", height = 8, width = 10, pointsize = 11)
boxplot(Grid120.800.data.summ)

#filter data, NAs and infinite, and wildly out values (20 fold)
#one strata at a time
Grid120.800.red1 <- find.rogue(Grid.800.results[,c(1:6)], mag = mag.factor, true = c(0.05, 2, 200))
Grid120.800.red2 <- find.rogue(Grid.800.results[,c(7:12)], mag = mag.factor, true = c(0.05/15, 2/15, 3000))

Grid120.800.clean <- merge(Grid120.800.red1, Grid120.800.red2, by = "row", all.x = TRUE, all.y = TRUE)
summary(Grid120.800.clean)

#plot and summarise
#estimates first
Grid120.800.D.est <- summ.results(Grid120.800.clean, par = "D", true.value = c(0.05,0.05/15) , ylim.ceiling = c(D1.ylim, D2.ylim), plot = TRUE, se = FALSE)
Grid120.800.L0.est <- summ.results(Grid120.800.clean, par = "L0", true.value = c(2, 2/15) , ylim.ceiling = c(L1.ylim, L2.ylim), plot = TRUE, se = FALSE)
Grid120.800.Sig.est <- summ.results(Grid120.800.clean, par = "Sigma", true.value = c(200, 3000) , ylim.ceiling = c(Sig1.ylim, Sig2.ylim), plot = TRUE, se = FALSE)

#relative standard errors
Grid120.800.D.se <- summ.results(Grid120.800.clean, par = "D", true.value = c(0.05,0.05/15) , ylim.ceiling = c(10,10), plot = TRUE, se = TRUE, rel.se = TRUE)
Grid120.800.L0.se <- summ.results(Grid120.800.clean, par = "L0", true.value = c(2, 2/15) , ylim.ceiling = c(10,10), plot = TRUE, se = TRUE, rel.se = TRUE)
Grid120.800.Sig.se <- summ.results(Grid120.800.clean, par = "Sigma", true.value = c(200, 3000) , ylim.ceiling = c(10,10), plot = TRUE, se = TRUE, rel.se = TRUE)

################
#Optimal Spacing
#using OS of 700m
################

load("15FoldScen/Cluster/Sims/120Traps/Grids700Results120.RData")

##first data summary, then estimates
Grid120.700.data.summ <- CH.data.summs(Grid.700.Data$output)

pdf("15FoldScen/Grid700.120.DataSumms.pdf", height = 8, width = 10, pointsize = 11)
boxplot(Grid120.700.data.summ)
dev.off()

#filter data, NAs and infinite, and wildly out values (20 fold)
#one strata at a time
Grid120.700.red1 <- find.rogue(Grid.700.results[,c(1:6)], mag = mag.factor, true = c(0.05, 2, 200))
Grid120.700.red2 <- find.rogue(Grid.700.results[,c(7:12)], mag = mag.factor, true = c(0.05/15, 2/15, 3000))

Grid120.700.clean <- merge(Grid120.700.red1, Grid120.700.red2, by = "row", all.x = TRUE, all.y = TRUE)
summary(Grid120.700.clean)

#plot and summarise
#estimates first
Grid120.700.D.est <- summ.results(Grid120.700.clean, par = "D", true.values = c(0.05,0.05/15) , ylim.ceiling = c(D1.ylim, D2.ylim), plot = TRUE, se = FALSE)
Grid120.700.L0.est <- summ.results(Grid120.700.clean, par = "L0", true.value = c(2, 2/15) , ylim.ceiling = c(L1.ylim, L2.ylim), plot = TRUE, se = FALSE)
Grid120.700.Sig.est <- summ.results(Grid120.700.clean, par = "Sigma", true.value = c(200, 3000) , ylim.ceiling = c(Sig1.ylim, Sig2.ylim), plot = TRUE, se = FALSE)

#rel standard errors
Grid120.700.D.se <- summ.results(Grid120.700.clean, par = "D", true.value = c(0.05,0.05/15) , ylim.ceiling = c(10, 10), plot = TRUE, se = TRUE, rel.se = TRUE)
Grid120.700.L0.se <- summ.results(Grid120.700.clean, par = "L0", true.value = c(2, 2/15) , ylim.ceiling = c(10, 10), plot = TRUE, se = TRUE, rel.se = TRUE)
Grid120.700.Sig.se <- summ.results(Grid120.700.clean, par = "Sigma", true.value = c(200, 3000) , ylim.ceiling = c(10, 10), plot = TRUE, se = TRUE, rel.se = TRUE)

################
#Clustered designs
#three configs tried
#using 2 configs: OS spacings and 2 * sigma
################

#Clustersa uses OS and Clustersb uses 2 sigma spacing
#500 reps 
load("15FoldScen/Cluster/Sims/120Traps/Clusters120a.RData")

#OS of 500 / 1200, objs named: Grid.os.xx
##first data summary, then estimates
Cluster120.os.data.summ <- CH.data.summs(Grid.os.Data$output)

boxplot(Cluster120.os.data.summ)

#filter data, NAs and infinite, and wildly out values (20 fold)
#one strata at a time
Cluster120.os.red1 <- find.rogue(Grid.os.results[,c(1:6)], mag = mag.factor, true = c(0.05, 2, 200))
Cluster120.os.red2 <- find.rogue(Grid.os.results[,c(7:12)], mag = mag.factor, true = c(0.05/15, 2/15, 3000))

Cluster120.os.clean <- merge(Cluster120.os.red1, Cluster120.os.red2, by = "row", all.x = TRUE, all.y = TRUE)
summary(Cluster120.os.clean)

#plot and summarise
#estimates first
Cluster120.os.D.est <- summ.results(Cluster120.os.clean, par = "D", true.value = c(0.05,0.05/15) , ylim.ceiling = c(D1.ylim, D2.ylim), plot = TRUE, se = FALSE)
Cluster120.os.L0.est <- summ.results(Cluster120.os.clean, par = "L0", true.value = c(2, 2/15) , ylim.ceiling = c(L1.ylim, L2.ylim), plot = TRUE, se = FALSE)
Cluster120.os.Sig.est <- summ.results(Cluster120.os.clean, par = "Sigma", true.value = c(200, 3000) , ylim.ceiling = c(Sig1.ylim, Sig2.ylim), plot = TRUE, se = FALSE)

#rel standard errors
Cluster120.os.D.se <- summ.results(Cluster120.os.clean, par = "D", true.value = c(0.05,0.05/15) , ylim.ceiling = c(10, 10), plot = TRUE, se = TRUE, rel.se = TRUE)
Cluster120.os.L0.se <- summ.results(Cluster120.os.clean, par = "L0", true.value = c(2, 2/15) , ylim.ceiling = c(10, 10), plot = TRUE, se = TRUE, rel.se = TRUE)
Cluster120.os.Sig.se <- summ.results(Cluster120.os.clean, par = "Sigma", true.value = c(200, 3000) , ylim.ceiling = c(10, 10), plot = TRUE, se = TRUE, rel.se = TRUE)

#now trying 2 sigma spacing
load("15FoldScen/Cluster/Sims/120Traps/Clusters120b.RData")

##first data summary, then estimates
Cluster120.2sig.data.summ <- CH.data.summs(Grid.2sig.Data$output)

boxplot(Cluster120.2sig.data.summ)

#filter data, NAs and infinite, and wildly out values
#one strata at a time
Cluster120.2sig.red1 <- find.rogue(Grid.2sig.results[,c(1:6)], mag = mag.factor, true = c(0.05, 2, 200))
Cluster120.2sig.red2 <- find.rogue(Grid.2sig.results[,c(7:12)], mag = mag.factor, true = c(0.05/15, 2/15, 3000))

Cluster120.2sig.clean <- merge(Cluster120.2sig.red1, Cluster120.2sig.red2, by = "row", all.x = TRUE, all.y = TRUE)
summary(Cluster120.2sig.clean)

#plot and summarise
#estimates first
Cluster120.2sig.D.est <- summ.results(Cluster120.2sig.clean, par = "D", true.value = c(0.05,0.05/15) , ylim.ceiling = c(D1.ylim, D2.ylim), plot = TRUE, se = FALSE)
Cluster120.2sig.L0.est <- summ.results(Cluster120.2sig.clean, par = "L0", true.value = c(2, 2/15) , ylim.ceiling = c(L1.ylim, L2.ylim), plot = TRUE, se = FALSE)
Cluster120.2sig.Sig.est <- summ.results(Cluster120.2sig.clean, par = "Sigma", true.value = c(200, 3000) , ylim.ceiling = c(Sig1.ylim, Sig2.ylim), plot = TRUE, se = FALSE)

#standard errors
Cluster120.2sig.D.se <- summ.results(Cluster120.2sig.clean, par = "D", true.value = c(0.05,0.05/15) , ylim.ceiling = c(10, 10), plot = TRUE, se = TRUE, rel.se = TRUE)
Cluster120.2sig.L0.se <- summ.results(Cluster120.2sig.clean, par = "L0", true.value = c(2, 2/15) , ylim.ceiling = c(10, 10), plot = TRUE, se = TRUE, rel.se = TRUE)
Cluster120.2sig.Sig.se <- summ.results(Cluster120.2sig.clean, par = "Sigma", true.value = c(200, 3000) , ylim.ceiling = c(10, 10), plot = TRUE, se = TRUE, rel.se = TRUE)

################
#GA4 designs
################

#extract and collate results
#GA 500 reps in 2 parts eg: a (100) and a2 (400) ; S1 and S2 exceptions
#letters refer to diff GA types

######################
#GA4 with S1 pars (a)
######################

load("15FoldScen/Cluster/Sims/120Traps/GA4Results120a.RData")
load("15FoldScen/Cluster/Sims/120Traps/GA4Results120a2.RData")

G4.S1.Data <- list(
  output = c(G4.S1.Dataa$output, G4.S1.Dataa2$output)
)

G4.S1.results <- rbind(G4.S1.resultsa, G4.S1.resultsa2)

GA4.120.S1.data.summ <- CH.data.summs(G4.S1.Data$output)
boxplot(GA4.120.S1.data.summ)

#filter estimates, NAs and infinite, and wildly out values (10 fold)
#one strata at a time
GA4.120.S1.red1 <- find.rogue(G4.S1.results[,c(1:6)], mag = mag.factor, true = c(0.05, 2, 200))
GA4.120.S1.red2 <- find.rogue(G4.S1.results[,c(7:12)], mag = mag.factor, true = c(0.05/15, 2/15, 3000))

GA4.120.S1.clean <- merge(GA4.120.S1.red1, GA4.120.S1.red2, by = "row", all.x = TRUE, all.y = TRUE)
summary(GA4.120.S1.clean)

#plot and summarise
#estimates first
GA4.120.S1.D.est <- summ.results(GA4.120.S1.clean, par = "D", true.value = c(0.05,0.05/15) , ylim.ceiling = c(D1.ylim, D2.ylim), plot = TRUE, se = FALSE)
GA4.120.S1.L0.est <- summ.results(GA4.120.S1.clean, par = "L0", true.value = c(2, 2/15) , ylim.ceiling = c(L1.ylim, L2.ylim), plot = TRUE, se = FALSE)
GA4.120.S1.Sig.est <- summ.results(GA4.120.S1.clean, par = "Sigma", true.value = c(200, 3000) , ylim.ceiling = c(Sig1.ylim, Sig2.ylim), plot = TRUE, se = FALSE)

#standard errors
GA4.120.S1.D.se <- summ.results(GA4.120.S1.clean, par = "D", true.value = c(0.05,0.05/15) , ylim.ceiling = c(10, 10), plot = TRUE, se = TRUE, rel.se = TRUE)
GA4.120.S1.L0.se <- summ.results(GA4.120.S1.clean, par = "L0", true.value = c(2, 2/15) , ylim.ceiling = c(10, 10), plot = TRUE, se = TRUE, rel.se = TRUE)
GA4.120.S1.Sig.se <- summ.results(GA4.120.S1.clean, par = "Sigma", true.value = c(200, 3000) , ylim.ceiling = c(10, 10), plot = TRUE, se = TRUE, rel.se = TRUE)

######################
#GA4 with S2 pars (b)
######################

load("15FoldScen/Cluster/Sims/120Traps/GA4Results120b.RData")
load("15FoldScen/Cluster/Sims/120Traps/GA4Results120b2.RData")

G4.S2.Data <- list(
  output = c(G4.S2.Datab$output, G4.S2.Datab2$output)
)

G4.S2.results <- rbind(G4.S2.resultsb, G4.S2.resultsb2)

GA4.120.S2.data.summ <- CH.data.summs(G4.S2.Data$output)
boxplot(GA4.120.S2.data.summ)

#filter estimates, NAs and infinite, and wildly out values (10 fold)
#one strata at a time
GA4.120.S2.red1 <- find.rogue(G4.S2.results[,c(1:6)], mag = mag.factor, true = c(0.05, 2, 200))
GA4.120.S2.red2 <- find.rogue(G4.S2.results[,c(7:12)], mag = mag.factor, true = c(0.05/15, 2/15, 3000))

GA4.120.S2.clean <- merge(GA4.120.S2.red1, GA4.120.S2.red2, by = "row", all.x = TRUE, all.y = TRUE)
summary(GA4.120.S2.clean)

#plot and summarise
#estimates first
GA4.120.S2.D.est <- summ.results(GA4.120.S2.clean, par = "D", true.value = c(0.05,0.05/15) , ylim.ceiling = c(D1.ylim, D2.ylim), plot = TRUE, se = FALSE)
GA4.120.S2.L0.est <- summ.results(GA4.120.S2.clean, par = "L0", true.value = c(2, 2/15) , ylim.ceiling = c(L1.ylim, L2.ylim), plot = TRUE, se = FALSE)
GA4.120.S2.Sig.est <- summ.results(GA4.120.S2.clean, par = "Sigma", true.value = c(200, 3000) , ylim.ceiling = c(Sig1.ylim, Sig2.ylim), plot = TRUE, se = FALSE)

#standard errors
GA4.120.S2.D.se <- summ.results(GA4.120.S2.clean, par = "D", true.value = c(0.05,0.05/15) , ylim.ceiling = c(10, 10), plot = TRUE, se = TRUE, rel.se = TRUE)
GA4.120.S2.L0.se <- summ.results(GA4.120.S2.clean, par = "L0", true.value = c(2, 2/15) , ylim.ceiling = c(10, 10), plot = TRUE, se = TRUE, rel.se = TRUE)
GA4.120.S2.Sig.se <- summ.results(GA4.120.S2.clean, par = "Sigma", true.value = c(200, 3000) , ylim.ceiling = c(10, 10), plot = TRUE, se = TRUE, rel.se = TRUE)

######################
#GA4 with avg pars (c)
######################

load("15FoldScen/Cluster/Sims/120Traps/GA4Results120C.RData")

GA4.120.Avg.data.summ <- CH.data.summs(G4.Avg.Datac$output)
boxplot(GA4.120.Avg.data.summ)

#filter estimates, NAs and infinite, and wildly out values (10 fold)
#one strata at a time
GA4.120.Avg.red1 <- find.rogue(G4.Avg.resultsc[,c(1:6)], mag = mag.factor, true = c(0.05, 2, 200))
GA4.120.Avg.red2 <- find.rogue(G4.Avg.resultsc[,c(7:12)], mag = mag.factor, true = c(0.05/15, 2/15, 3000))

GA4.120.Avg.clean <- merge(GA4.120.Avg.red1, GA4.120.Avg.red2, by = "row", all.x = TRUE, all.y = TRUE)
summary(GA4.120.Avg.clean)

#plot and summarise
#estimates first
GA4.120.Avg.D.est <- summ.results(GA4.120.Avg.clean, par = "D", true.value = c(0.05,0.05/15) , ylim.ceiling = c(D1.ylim, D2.ylim), plot = TRUE, se = FALSE)
GA4.120.Avg.L0.est <- summ.results(GA4.120.Avg.clean, par = "L0", true.value = c(2, 2/15) , ylim.ceiling = c(L1.ylim, L2.ylim), plot = TRUE, se = FALSE)
GA4.120.Avg.Sig.est <- summ.results(GA4.120.Avg.clean, par = "Sigma", true.value = c(200, 3000) , ylim.ceiling = c(Sig1.ylim, Sig2.ylim), plot = TRUE, se = FALSE)

#standard errors
GA4.120.Avg.D.se <- summ.results(GA4.120.Avg.clean, par = "D", true.value = c(0.05,0.05/15) , ylim.ceiling = c(10, 10), plot = TRUE, se = TRUE, rel.se = TRUE)
GA4.120.Avg.L0.se <- summ.results(GA4.120.Avg.clean, par = "L0", true.value = c(2, 2/15) , ylim.ceiling = c(10, 10), plot = TRUE, se = TRUE, rel.se = TRUE)
GA4.120.Avg.Sig.se <- summ.results(GA4.120.Avg.clean, par = "Sigma", true.value = c(200, 3000) , ylim.ceiling = c(10, 10), plot = TRUE, se = TRUE, rel.se = TRUE)

######################
#GA4 with both (d)
######################

load("15FoldScen/Cluster/Sims/120Traps/GA4Results120d.RData")

GA4.120.Both.data.summ <- CH.data.summs(G4.Both.Datad$output)
boxplot(GA4.120.Both.data.summ)

#filter estimates, NAs and infinite, and wildly out values (10 fold)
#one strata at a time
GA4.120.Both.red1 <- find.rogue(G4.Both.resultsd[,c(1:6)], mag = mag.factor, true = c(0.05, 2, 200))
GA4.120.Both.red2 <- find.rogue(G4.Both.resultsd[,c(7:12)], mag = mag.factor, true = c(0.05/15, 2/15, 3000))

GA4.120.Both.clean <- merge(GA4.120.Both.red1, GA4.120.Both.red2, by = "row", all.x = TRUE, all.y = TRUE)
summary(GA4.120.Both.clean)

#plot and summarise
#estimates first
GA4.120.Both.D.est <- summ.results(GA4.120.Both.clean, par = "D", true.value = c(0.05,0.05/15) , ylim.ceiling = c(D1.ylim, D2.ylim), plot = TRUE, se = FALSE)
GA4.120.Both.L0.est <- summ.results(GA4.120.Both.clean, par = "L0", true.value = c(2, 2/15) , ylim.ceiling = c(L1.ylim, L2.ylim), plot = TRUE, se = FALSE)
GA4.120.Both.Sig.est <- summ.results(GA4.120.Both.clean, par = "Sigma", true.value = c(200, 3000) , ylim.ceiling = c(Sig1.ylim, Sig2.ylim), plot = TRUE, se = FALSE)

#standard errors
GA4.120.Both.D.se <- summ.results(GA4.120.Both.clean, par = "D", true.value = c(0.05,0.05/15) , ylim.ceiling = c(10, 10), plot = TRUE, se = TRUE, rel.se = TRUE)
GA4.120.Both.L0.se <- summ.results(GA4.120.Both.clean, par = "L0", true.value = c(2, 2/15) , ylim.ceiling = c(10, 10), plot = TRUE, se = TRUE, rel.se = TRUE)
GA4.120.Both.Sig.se <- summ.results(GA4.120.Both.clean, par = "Sigma", true.value = c(200, 3000) , ylim.ceiling = c(10, 10), plot = TRUE, se = TRUE, rel.se = TRUE)

############################
#GA4 with both (max = T) (e)
############################

load("15FoldScen/Cluster/Sims/120Traps/GA4Results120e.RData")

GA4.120.BothMax.data.summ <- CH.data.summs(G4.BothMax.Datae$output)
boxplot(GA4.BothMax.data.summ)

#filter estimates, NAs and infinite, and wildly out values (10 fold)
#one strata at a time
GA4.120.BothMax.red1 <- find.rogue(G4.BothMax.resultse[,c(1:6)], mag = mag.factor, true = c(0.05, 2, 200))
GA4.120.BothMax.red2 <- find.rogue(G4.BothMax.resultse[,c(7:12)], mag = mag.factor, true = c(0.05/15, 2/15, 3000))

GA4.120.BothMax.clean <- merge(GA4.120.BothMax.red1, GA4.120.BothMax.red2, by = "row", all.x = TRUE, all.y = TRUE)
summary(GA4.120.BothMax.clean)

#plot and summarise
#estimates first
GA4.120.BothMax.D.est <- summ.results(GA4.120.BothMax.clean, par = "D", true.value = c(0.05,0.05/15) , ylim.ceiling = c(D1.ylim, D2.ylim), plot = TRUE, se = FALSE)
GA4.120.BothMax.L0.est <- summ.results(GA4.120.BothMax.clean, par = "L0", true.value = c(2, 2/15) , ylim.ceiling = c(L1.ylim, L2.ylim), plot = TRUE, se = FALSE)
GA4.120.BothMax.Sig.est <- summ.results(GA4.120.BothMax.clean, par = "Sigma", true.value = c(200, 3000) , ylim.ceiling = c(Sig1.ylim, Sig2.ylim), plot = TRUE, se = FALSE)

#standard errors
GA4.120.BothMax.D.se <- summ.results(GA4.120.BothMax.clean, par = "D", true.value = c(0.05,0.05/15) , ylim.ceiling = c(10, 10), plot = TRUE, se = TRUE, rel.se = TRUE)
GA4.120.BothMax.L0.se <- summ.results(GA4.120.BothMax.clean, par = "L0", true.value = c(2, 2/15) , ylim.ceiling = c(10, 10), plot = TRUE, se = TRUE, rel.se = TRUE)
GA4.120.BothMax.Sig.se <- summ.results(GA4.120.BothMax.clean, par = "Sigma", true.value = c(200, 3000) , ylim.ceiling = c(10, 10), plot = TRUE, se = TRUE, rel.se = TRUE)

################################
#GA5 designs
################################

######################
#GA5 with S1 pars (a)
######################

load("15FoldScen/Cluster/Sims/120Traps/GA5Results120a.RData")
load("15FoldScen/Cluster/Sims/120Traps/GA5Results120a2.RData")

G5.S1.Data <- list(
  output = c(G5.S1.Dataa$output, G5.S1.Dataa2$output)
)

G5.S1.results <- rbind(G5.S1.resultsa, G5.S1.resultsa2)

GA5.120.S1.data.summ <- CH.data.summs(G5.S1.Data$output)
boxplot(GA5.120.S1.data.summ)

#filter estimates, NAs and infinite, and wildly out values (10 fold)
#one strata at a time
GA5.120.S1.red1 <- find.rogue(G5.S1.results[,c(1:6)], mag = mag.factor, true = c(0.05, 2, 200))
GA5.120.S1.red2 <- find.rogue(G5.S1.results[,c(7:12)], mag = mag.factor, true = c(0.05/15, 2/15, 3000))

GA5.120.S1.clean <- merge(GA5.120.S1.red1, GA5.120.S1.red2, by = "row", all.x = TRUE, all.y = TRUE)
summary(GA5.120.S1.clean)

#plot and summarise
#estimates first
GA5.120.S1.D.est <- summ.results(GA5.120.S1.clean, par = "D", true.value = c(0.05,0.05/15) , ylim.ceiling = c(D1.ylim, D2.ylim), plot = TRUE, se = FALSE)
GA5.120.S1.L0.est <- summ.results(GA5.120.S1.clean, par = "L0", true.value = c(2, 2/15) , ylim.ceiling = c(L1.ylim, L2.ylim), plot = TRUE, se = FALSE)
GA5.120.S1.Sig.est <- summ.results(GA5.120.S1.clean, par = "Sigma", true.value = c(200, 3000) , ylim.ceiling = c(Sig1.ylim, Sig2.ylim), plot = TRUE, se = FALSE)

#rel standard errors
GA5.120.S1.D.se <- summ.results(GA5.120.S1.clean, par = "D", true.value = c(0.05,0.05/15) , ylim.ceiling = c(10, 10), plot = TRUE, se = TRUE, rel.se = TRUE)
GA5.120.S1.L0.se <- summ.results(GA5.120.S1.clean, par = "L0", true.value = c(2, 2/15) , ylim.ceiling = c(10, 10), plot = TRUE, se = TRUE, rel.se = TRUE)
GA5.120.S1.Sig.se <- summ.results(GA5.120.S1.clean, par = "Sigma", true.value = c(200, 3000) , ylim.ceiling = c(10, 10), plot = TRUE, se = TRUE, rel.se = TRUE)

######################
#GA5 with S2 pars (b)
######################

load("15FoldScen/Cluster/Sims/120Traps/GA5Results120b.RData")
load("15FoldScen/Cluster/Sims/120Traps/GA5Results120b2.RData")

G5.S2.Data <- list(
  output = c(G5.S2.Datab$output, G5.S2.Datab2$output)
)

G5.S2.results <- rbind(G5.S2.resultsb, G5.S2.resultsb2)

GA5.120.S2.data.summ <- CH.data.summs(G5.S2.Data$output)
boxplot(GA5.120.S2.data.summ)

#filter estimates, NAs and infinite, and wildly out values (10 fold)
#one strata at a time
GA5.120.S2.red1 <- find.rogue(G5.S2.results[,c(1:6)], mag = mag.factor, true = c(0.05, 2, 200))
GA5.120.S2.red2 <- find.rogue(G5.S2.results[,c(7:12)], mag = mag.factor, true = c(0.05/15, 2/15, 3000))

GA5.120.S2.clean <- merge(GA5.120.S2.red1, GA5.120.S2.red2, by = "row", all.x = TRUE, all.y = TRUE)
summary(GA5.120.S2.clean)

#plot and summarise
#estimates first
GA5.120.S2.D.est  <- summ.results(GA5.120.S2.clean, par = "D", true.value = c(0.05,0.05/15) , ylim.ceiling = c(D1.ylim, D2.ylim), plot = TRUE, se = FALSE)
GA5.120.S2.L0.est  <- summ.results(GA5.120.S2.clean, par = "L0", true.value = c(2, 2/15) , ylim.ceiling = c(L1.ylim, L2.ylim), plot = TRUE, se = FALSE)
GA5.120.S2.Sig.est  <- summ.results(GA5.120.S2.clean, par = "Sigma", true.value = c(200, 3000) , ylim.ceiling = c(Sig1.ylim, Sig2.ylim), plot = TRUE, se = FALSE)

#rel standard errors
GA5.120.S2.D.se  <- summ.results(GA5.120.S2.clean, par = "D", true.value = c(0.05,0.05/15) , ylim.ceiling = c(10, 10), plot = TRUE, se = TRUE, rel.se = TRUE)
GA5.120.S2.L0.se  <- summ.results(GA5.120.S2.clean, par = "L0", true.value = c(2, 2/15) , ylim.ceiling = c(10, 10), plot = TRUE, se = TRUE, rel.se = TRUE)
GA5.120.S2.Sig.se  <- summ.results(GA5.120.S2.clean, par = "Sigma", true.value = c(200, 3000) , ylim.ceiling = c(10, 10), plot = TRUE, se = TRUE, rel.se = TRUE)

######################
#GA5 with avg pars (c)
######################

load("15FoldScen/Cluster/Sims/120Traps/GA5Results120c.RData")

GA5.120.Avg.data.summ <- CH.data.summs(G5.Avg.Datac$output)
boxplot(GA5.120.Avg.data.summ)

#filter estimates, NAs and infinite, and wildly out values (10 fold)
#one strata at a time
GA5.120.Avg.red1 <- find.rogue(G5.Avg.resultsc[,c(1:6)], mag = mag.factor, true = c(0.05, 2, 200))
GA5.120.Avg.red2 <- find.rogue(G5.Avg.resultsc[,c(7:12)], mag = mag.factor, true = c(0.05/15, 2/15, 3000))

GA5.120.Avg.clean <- merge(GA5.120.Avg.red1, GA5.120.Avg.red2, by = "row", all.x = TRUE, all.y = TRUE)
summary(GA5.120.Avg.clean)

#plot and summarise
#estimates first
GA5.120.Avg.D.est  <- summ.results(GA5.120.Avg.clean, par = "D", true.value = c(0.05,0.05/15) , ylim.ceiling = c(D1.ylim, D2.ylim), plot = TRUE, se = FALSE)
GA5.120.Avg.L0.est  <- summ.results(GA5.120.Avg.clean, par = "L0", true.value = c(2, 2/15) , ylim.ceiling = c(L1.ylim, L2.ylim), plot = TRUE, se = FALSE)
GA5.120.Avg.Sig.est  <- summ.results(GA5.120.Avg.clean, par = "Sigma", true.value = c(200, 3000) , ylim.ceiling = c(Sig1.ylim, Sig2.ylim), plot = TRUE, se = FALSE)

#rel standard errors
GA5.120.Avg.D.se <- summ.results(GA5.120.Avg.clean, par = "D", true.value = c(0.05,0.05/15) , ylim.ceiling = c(10, 10), plot = TRUE, se = TRUE, rel.se = TRUE)
GA5.120.Avg.L0.se <- summ.results(GA5.120.Avg.clean, par = "L0", true.value = c(2, 2/15) , ylim.ceiling = c(10, 10), plot = TRUE, se = TRUE, rel.se = TRUE)
GA5.120.Avg.Sig.se <- summ.results(GA5.120.Avg.clean, par = "Sigma", true.value = c(200, 3000) , ylim.ceiling = c(10, 10), plot = TRUE, se = TRUE, rel.se = TRUE)

######################
#GA5 with both (d)
######################

load("15FoldScen/Cluster/Sims/120Traps/GA5Results120d.RData")

GA5.120.Both.data.summ <- CH.data.summs(G5.Both.Datad$output)
boxplot(GA5.120.Both.data.summ)

#filter estimates, NAs and infinite, and wildly out values (10 fold)
#one strata at a time
GA5.120.Both.red1 <- find.rogue(G5.Both.resultsd[,c(1:6)], mag = mag.factor, true = c(0.05, 2, 200))
GA5.120.Both.red2 <- find.rogue(G5.Both.resultsd[,c(7:12)], mag = mag.factor, true = c(0.05/15, 2/15, 3000))

GA5.120.Both.clean <- merge(GA5.120.Both.red1, GA5.120.Both.red2, by = "row", all.x = TRUE, all.y = TRUE)
summary(GA5.120.Both.clean)

#plot and summarise
#estimates first
GA5.120.Both.D.est <- summ.results(GA5.120.Both.clean, par = "D", true.value = c(0.05,0.05/15) , ylim.ceiling = c(D1.ylim, D2.ylim), plot = TRUE, se = FALSE)
GA5.120.Both.L0.est <- summ.results(GA5.120.Both.clean, par = "L0", true.value = c(2, 2/15) , ylim.ceiling = c(L1.ylim, L2.ylim), plot = TRUE, se = FALSE)
GA5.120.Both.Sig.est <- summ.results(GA5.120.Both.clean, par = "Sigma", true.value = c(200, 3000) , ylim.ceiling = c(Sig1.ylim, Sig2.ylim), plot = TRUE, se = FALSE)

#rel standard errors
GA5.120.Both.D.se <- summ.results(GA5.120.Both.clean, par = "D", true.value = c(0.05,0.05/15) , ylim.ceiling = c(10, 10), plot = TRUE, se = TRUE, rel.se = TRUE)
GA5.120.Both.L0.se <- summ.results(GA5.120.Both.clean, par = "L0", true.value = c(2, 2/15) , ylim.ceiling = c(10, 10), plot = TRUE, se = TRUE, rel.se = TRUE)
GA5.120.Both.Sig.se <- summ.results(GA5.120.Both.clean, par = "Sigma", true.value = c(200, 3000) , ylim.ceiling = c(10, 10), plot = TRUE, se = TRUE, rel.se = TRUE)

######################
#GA5 with both (max = T) (e)
######################

load("15FoldScen/Cluster/Sims/120Traps/GA5Results120e.RData")

GA5.120.BothMax.data.summ <- CH.data.summs(G5.BothMax.Datae$output)
boxplot(GA5.120.BothMax.data.summ)

#filter estimates, NAs and infinite, and wildly out values (10 fold)
#one strata at a time
GA5.120.BothMax.red1 <- find.rogue(G5.BothMax.resultse[,c(1:6)], mag = mag.factor, true = c(0.05, 2, 200))
GA5.120.BothMax.red2 <- find.rogue(G5.BothMax.resultse[,c(7:12)], mag = mag.factor, true = c(0.05/15, 2/15, 3000))

GA5.120.BothMax.clean <- merge(GA5.120.BothMax.red1, GA5.120.BothMax.red2, by = "row", all.x = TRUE, all.y = TRUE)
summary(GA5.120.BothMax.clean)

#plot and summarise
#estimates first
GA5.120.BothMax.D.est <- summ.results(GA5.120.BothMax.clean, par = "D", true.value = c(0.05,0.05/15) , ylim.ceiling = c(D1.ylim, D2.ylim), plot = TRUE, se = FALSE)
GA5.120.BothMax.L0.est <- summ.results(GA5.120.BothMax.clean, par = "L0", true.value = c(2, 2/15) , ylim.ceiling = c(L1.ylim, L2.ylim), plot = TRUE, se = FALSE)
GA5.120.BothMax.Sig.est <- summ.results(GA5.120.BothMax.clean, par = "Sigma", true.value = c(200, 3000) , ylim.ceiling = c(Sig1.ylim, Sig2.ylim), plot = TRUE, se = FALSE)

#rel standard errors
GA5.120.BothMax.D.se <- summ.results(GA5.120.BothMax.clean, par = "D", true.value = c(0.05,0.05/15) , ylim.ceiling = c(10, 10), plot = TRUE, se = TRUE, rel.se = TRUE)
GA5.120.BothMax.L0.se <- summ.results(GA5.120.BothMax.clean, par = "L0", true.value = c(2, 2/15) , ylim.ceiling = c(10, 10), plot = TRUE, se = TRUE, rel.se = TRUE)
GA5.120.BothMax.Sig.se <- summ.results(GA5.120.BothMax.clean, par = "Sigma", true.value = c(200, 3000) , ylim.ceiling = c(10, 10), plot = TRUE, se = TRUE, rel.se = TRUE)

############################################################################################
#Lacework

load("15FoldScen/Cluster/Sims/120Traps/LWResults120.RData")

LW.120.Data.reformatted <- lapply(LW.120.Data$output$`1`, function(x) list(x))

LW120.data.summ <- CH.data.summs(LW.120.Data.reformatted)
boxplot(LW120.data.summ)

#filter estimates, NAs and infinite, and wildly out values (10 fold)
#one strata at a time
LW120.red1 <- find.rogue(LW.120.results[,c(1:6)], mag = mag.factor, true = c(0.05, 2, 200))
LW120.red2 <- find.rogue(LW.120.results[,c(7:12)], mag = mag.factor, true = c(0.05/15, 2/15, 3000))

LW120.clean <- merge(LW120.red1, LW120.red2, by = "row", all.x = TRUE, all.y = TRUE)
summary(LW120.clean)

#plot and summarise
#estimates first
LW120.D.est <- summ.results(LW120.clean, par = "D", true.value = c(0.05,0.05/15) , ylim.ceiling = c(D1.ylim, D2.ylim), plot = TRUE, se = FALSE)
LW120.L0.est <- summ.results(LW120.clean, par = "L0", true.value = c(2, 2/15) , ylim.ceiling = c(L1.ylim, L2.ylim), plot = TRUE, se = FALSE)
LW120.Sig.est <- summ.results(LW120.clean, par = "Sigma", true.value = c(200, 3000) , ylim.ceiling = c(Sig1.ylim, Sig2.ylim), plot = TRUE, se = FALSE)

#rel standard errors
LW120.D.se <- summ.results(LW120.clean, par = "D", true.value = c(0.05,0.05/15) , ylim.ceiling = c(10, 10), plot = TRUE, se = TRUE, rel.se = TRUE)
LW120.L0.se <- summ.results(LW120.clean, par = "L0", true.value = c(2, 2/15) , ylim.ceiling = c(10, 10), plot = TRUE, se = TRUE, rel.se = TRUE)
LW120.Sig.se <- summ.results(LW120.clean, par = "Sigma", true.value = c(200, 3000) , ylim.ceiling = c(10, 10), plot = TRUE, se = TRUE, rel.se = TRUE)

############################################################################################
#Two stage - busy running

##########################################################
#save all results objects into a list for easy loading
datasumms.120 <- list("Grid120.700" = Grid120.700.data.summ, "Grid120.800" = Grid120.800.data.summ,
                      "Clust120.OS" = Cluster120.os.data.summ, "Clust120.2sig" = Cluster120.2sig.data.summ, "G4.120.S1" = GA4.120.S1.data.summ, 
                      "G4.120.S2" = GA4.120.S2.data.summ, "G4.120.Avg" = GA4.120.Avg.data.summ, "G4.120.Both" = GA4.120.Both.data.summ, "G4.120.BothMax" = GA4.120.BothMax.data.summ,
                      "G5.120.S1" = GA5.120.S1.data.summ, "G5.120.S2" = GA5.120.S2.data.summ, "G5.120.Avg" = GA5.120.Avg.data.summ, "G5.120.Both" = GA5.120.Both.data.summ, 
                      "G5.120.BothMax" = GA5.120.BothMax.data.summ, "LW40" = LW120.data.summ)
save(datasumms.120, file = "15FoldScen/DataSumms120Traps500.RData")

simresults.120 <- list("Grid120.700.D" = Grid120.700.D.est, "Grid120.700.L0" = Grid120.700.L0.est, "Grid120.700.Sig" = Grid120.700.Sig.est,
                       "Grid120.800.D" = Grid120.800.D.est, "Grid120.800.L0" = Grid120.800.L0.est, "Grid120.800.Sig" = Grid120.800.Sig.est,
                       "Clust120.OS.D" = Cluster120.os.D.est, "Clust120.OS.L0" = Cluster120.os.L0.est, "Clust120.OS.Sig" = Cluster120.os.Sig.est,
                       "Clust120.2sig.D" = Cluster120.2sig.D.est, "Clust120.2sig.L0" = Cluster120.2sig.L0.est, "Clust120.2sig.Sig" = Cluster120.2sig.Sig.est,
                       "GA4.120.S1.D" = GA4.120.S1.D.est, "GA4.120.S1.L0" = GA4.120.S1.L0.est, "GA4.120.S1.Sig" = GA4.120.S1.Sig.est,
                       "GA4.120.S2.D" = GA4.120.S2.D.est, "GA4.120.S2.L0" = GA4.120.S2.L0.est, "GA4.120.S2.Sig" = GA4.120.S2.Sig.est,
                       "GA4.120.Avg.D" = GA4.120.Avg.D.est, "GA4.120.Avg.L0" = GA4.120.Avg.L0.est, "GA4.120.Avg.Sig" = GA4.120.Avg.Sig.est,
                       "GA4.120.Both.D" = GA4.120.Both.D.est, "GA4.120.Both.L0" = GA4.120.Both.L0.est, "GA4.120.Both.Sig" = GA4.120.Both.Sig.est,
                       "GA4.120.BothMax.D" = GA4.120.BothMax.D.est, "GA4.120.BothMax.L0" = GA4.120.BothMax.L0.est, "GA4.120.BothMax.Sig" = GA4.120.BothMax.Sig.est,
                       "GA5.120.S1.D" = GA5.120.S1.D.est, "GA5.120.S1.L0" = GA5.120.S1.L0.est, "GA5.120.S1.Sig" = GA5.120.S1.Sig.est,
                       "GA5.120.S2.D" = GA5.120.S2.D.est, "GA5.120.S2.L0" = GA5.120.S2.L0.est, "GA5.120.S2.Sig" = GA5.120.S2.Sig.est,
                       "GA5.120.Avg.D" = GA5.120.Avg.D.est, "GA5.120.Avg.L0" = GA5.120.Avg.L0.est, "GA5.120.Avg.Sig" = GA5.120.Avg.Sig.est,
                       "GA5.120.Both.D" = GA5.120.Both.D.est, "GA5.120.Both.L0" = GA5.120.Both.L0.est, "GA5.120.Both.Sig" = GA5.120.Both.Sig.est,
                       "GA5.120.BothMax.D" = GA5.120.BothMax.D.est, "GA5.120.BothMax.L0" = GA5.120.BothMax.L0.est, "GA5.120.BothMax.Sig" = GA5.120.BothMax.Sig.est,
                       "LW120.D" = LW120.D.est, "LW120.L0" = LW120.L0.est, "LW120.Sig" = LW120.Sig.est)
save(simresults.120, file = "15FoldScen/SimResults120Traps500.RData")

seresults.120 <- list("Grid120.700.D" = Grid120.700.D.se, "Grid120.700.L0" = Grid120.700.L0.se, "Grid120.700.Sig" = Grid120.700.Sig.se,
                      "Grid120.800.D" = Grid120.800.D.se, "Grid120.800.L0" = Grid120.800.L0.se, "Grid120.800.Sig" = Grid120.800.Sig.se,
                      "Clust120.OS.D" = Cluster120.os.D.se, "Clust120.OS.L0" = Cluster120.os.L0.se, "Clust120.OS.Sig" = Cluster120.os.Sig.se,
                      "Clust120.2sig.D" = Cluster120.2sig.D.se, "Clust120.2sig.L0" = Cluster120.2sig.L0.se, "Clust120.2sig.Sig" = Cluster120.2sig.Sig.se,
                      "GA4.120.S1.D" = GA4.120.S1.D.se, "GA4.120.S1.L0" = GA4.120.S1.L0.se, "GA4.120.S1.Sig" = GA4.120.S1.Sig.se,
                      "GA4.120.S2.D" = GA4.120.S2.D.se, "GA4.120.S2.L0" = GA4.120.S2.L0.se, "GA4.120.S2.Sig" = GA4.120.S2.Sig.se,
                      "GA4.120.Avg.D" = GA4.120.Avg.D.se, "GA4.120.Avg.L0" = GA4.120.Avg.L0.se, "GA4.120.Avg.Sig" = GA4.120.Avg.Sig.se,
                      "GA4.120.Both.D" = GA4.120.Both.D.se, "GA4.120.Both.L0" = GA4.120.Both.L0.se, "GA4.120.Both.Sig" = GA4.120.Both.Sig.se,
                      "GA4.120.BothMax.D" = GA4.120.BothMax.D.se, "GA4.120.BothMax.L0" = GA4.120.BothMax.L0.se, "GA4.120.BothMax.Sig" = GA4.120.BothMax.Sig.se,
                      "GA5.120.S1.D" = GA5.120.S1.D.se, "GA5.120.S1.L0" = GA5.120.S1.L0.se, "GA5.120.S1.Sig" = GA5.120.S1.Sig.se,
                      "GA5.120.S2.D" = GA5.120.S2.D.se, "GA5.120.S2.L0" = GA5.120.S2.L0.se, "GA5.120.S2.Sig" = GA5.120.S2.Sig.se,
                      "GA5.120.Avg.D" = GA5.120.Avg.D.se, "GA5.120.Avg.L0" = GA5.120.Avg.L0.se, "GA5.120.Avg.Sig" = GA5.120.Avg.Sig.se,
                      "GA5.120.Both.D" = GA5.120.Both.D.se, "GA5.120.Both.L0" = GA5.120.Both.L0.se, "GA5.120.Both.Sig" = GA5.120.Both.Sig.se,
                      "GA5.120.BothMax.D" = GA5.120.BothMax.D.se, "GA5.120.BothMax.L0" = GA5.120.BothMax.L0.se, "GA5.120.BothMax.Sig" = GA5.120.BothMax.Sig.se,
                      "LW120.D" = LW120.D.se, "LW120.L0" = LW120.L0.se, "LW120.Sig" = LW120.Sig.se)

save(seresults.120, file = "15FoldScen/seResults120Traps500.RData")

excluded.120 <- list("Grid120.700.S1" = Grid120.700.red1, "Grid120.700.S2" = Grid120.700.red2,
                     "Grid120.800.S1" = Grid120.800.red1, "Grid120.800.S2" = Grid120.800.red2,
                     "Clust120.os.S1" = Cluster120.os.red1, "Clust120.os.S2" = Cluster120.os.red2,
                     "Clust120.2sig.S1" = Cluster120.2sig.red1, "Clust120.2sig.S2" = Cluster120.2sig.red2,
                     "GA4.120.S1.S1" = GA4.120.S1.red1, "GA4.120.S1.S2" = GA4.120.S1.red2,
                     "GA4.120.S2.S1" = GA4.120.S2.red1, "GA4.120.S2.S2" = GA4.120.S2.red2,
                     "GA4.120.Avg.S1" = GA4.120.Avg.red1, "GA4.120.Avg.S2" = GA4.120.Avg.red2,
                     "GA4.120.Both.S1" = GA4.120.Both.red1, "GA4.120.Both.S2" = GA4.120.Both.red2,
                     "GA4.120.BothMax.S1" = GA4.120.BothMax.red1, "GA4.120.BothMax.S2" = GA4.120.BothMax.red2,
                     "GA5.120.S1.S1" = GA5.120.S1.red1, "GA5.120.S1.S2" = GA5.120.S1.red2,
                     "GA5.120.S2.S1" = GA5.120.S2.red1, "GA5.120.S2.S2" = GA5.120.S2.red2,
                     "GA5.120.Avg.S1" = GA5.120.Avg.red1, "GA5.120.Avg.S2" = GA5.120.Avg.red2,
                     "GA5.120.Both.S1" = GA5.120.Both.red1, "GA5.120.Both.S2" = GA5.120.Both.red2,
                     "GA5.120.BothMax.S1" = GA5.120.BothMax.red1, "GA5.120.BothMax.S2" = GA5.120.BothMax.red2,
                     "LW120.S1" = LW120.red1, "LW120.S2" = LW120.red2)
save(excluded.120, file = "15FoldScen/Excluded120Traps500.RData")

########################################################################

#kable to build tables
library(knitr)
library(kableExtra)

####################
#Grid
####################
grids <- data.frame(
  C1 = c("S1", "S2"),
  C2 = c(round(Grid800.D.est$`Rel bias (mean)`*100,1)),
  C3 = c(round(Grid800.D.est$`Rel SE`,3)),
  C4 = c(round(Grid800.L0.est$`Rel bias (mean)`*100,1)),
  C5 = c(round(Grid800.L0.est$`Rel SE`,3)),
  C6 = c(round(Grid800.Sig.est$`Rel bias (mean)`*100,1)),
  C7 = c(round(Grid800.Sig.est$`Rel SE`,3)),
  C8 = c(dim(Grid.800.red1)[1], dim(Grid.800.red2)[1])
)
colnames(grids) <- c("Strata","RB", "RSE", "RB", "RSE","RB", "RSE", "Reps") 

grids %>%
  kbl(format = "html", booktabs = TRUE, caption = "Sim Results from 800 m Grids", align = c('l', rep('c',6))) %>%
  kable_styling(font_size = 24) %>%
  kable_classic_2(full_width = F) %>%
  column_spec(c(1,3,5,7), border_right = T) %>%
  row_spec(0, bold = T, extra_css = "border-bottom: 1px solid") %>%
  row_spec(2,  extra_css = "border-bottom: 1px solid") %>%
  add_header_above(header = c(" " = 1, "D" = 2, "L0" = 2, "Sigma" = 2, " " = 1), border_right = TRUE, bold = TRUE) 

grids2 <- data.frame(
  C1 = c("S1", "S2"),
  C2 = c(round(Grid1600.D.est$`Rel bias (mean)`*100,1)),
  C3 = c(round(Grid1600.D.est$`Rel SE`,3)),
  C4 = c(round(Grid1600.L0.est$`Rel bias (mean)`*100,1)),
  C5 = c(round(Grid1600.L0.est$`Rel SE`,3)),
  C6 = c(round(Grid1600.Sig.est$`Rel bias (mean)`*100,1)),
  C7 = c(round(Grid1600.Sig.est$`Rel SE`,3)),
  C8 = c(dim(Grid.1600.red1)[1], dim(Grid.1600.red2)[1])
)
colnames(grids2) <- c("Strata","RB", "RSE", "RB", "RSE","RB", "RSE", "Reps") 

grids2 %>%
  kbl(format = "html", booktabs = TRUE, caption = "Sim Results from 1600 m Grids", align = c('l', rep('c',6))) %>%
  kable_styling(font_size = 24) %>%
  kable_classic_2(full_width = F) %>%
  column_spec(c(1,3,5,7), border_right = T) %>%
  row_spec(0, bold = T, extra_css = "border-bottom: 1px solid") %>%
  row_spec(2,  extra_css = "border-bottom: 1px solid") %>%
  add_header_above(header = c(" " = 1, "D" = 2, "L0" = 2, "Sigma" = 2, " " = 1), border_right = TRUE, bold = TRUE) 


####################
#GA4
####################

GA4 <- data.frame(
  C1 = c(rep("S1",2), rep("S2",2), rep("Avg",2), rep("Both",2), rep("Both (Max)",2)),
  C2 = c(rep(c("S1","S2"),5)),
  C3 = c(round(GA4.S1.D.est$`Rel bias (mean)`*100,1), round(GA4.S2.D.est$`Rel bias (mean)`*100,1), 
         round(GA4.Avg.D.est$`Rel bias (mean)`*100,1), round(GA4.Both.D.est$`Rel bias (mean)`*100,1), 
         round(GA4.BothMax.D.est$`Rel bias (mean)`*100,1)),
  C4 = c(round(GA4.S1.D.est$`Rel SE`,3), round(GA4.S2.D.est$`Rel SE`,3), round(GA4.Avg.D.est$`Rel SE`,3),
         round(GA4.Both.D.est$`Rel SE`,3), round(GA4.BothMax.D.est$`Rel SE`,3)),
  C5 = c(round(GA4.S1.L0.est$`Rel bias (mean)`*100,1), round(GA4.S2.L0.est$`Rel bias (mean)`*100,1), 
         round(GA4.Avg.L0.est$`Rel bias (mean)`*100,1), round(GA4.Both.L0.est$`Rel bias (mean)`*100,1), 
         round(GA4.BothMax.L0.est$`Rel bias (mean)`*100,1)),
  C6 = c(round(GA4.S1.L0.est$`Rel SE`,3), round(GA4.S2.L0.est$`Rel SE`,3), round(GA4.Avg.L0.est$`Rel SE`,3),
         round(GA4.Both.L0.est$`Rel SE`,3), round(GA4.BothMax.L0.est$`Rel SE`,3)),
  C7 = c(round(GA4.S1.Sig.est$`Rel bias (mean)`*100,1), round(GA4.S2.Sig.est$`Rel bias (mean)`*100,1), 
         round(GA4.Avg.Sig.est$`Rel bias (mean)`*100,1), round(GA4.Both.Sig.est$`Rel bias (mean)`*100,1), 
         round(GA4.Both.Sig.est$`Rel bias (mean)`*100,1)),
  C8 = c(round(GA4.S1.Sig.est$`Rel SE`,3), round(GA4.S2.Sig.est$`Rel SE`,3), round(GA4.Avg.Sig.est$`Rel SE`,3),
         round(GA4.Both.Sig.est$`Rel SE`,3), round(GA4.BothMax.Sig.est$`Rel SE`,3)),
  C9 = c(dim(GA4.S1.red1)[1], dim(GA4.S1.red2)[1], dim(GA4.S2.red1)[1], dim(GA4.S2.red2)[1], dim(GA4.Avg.red1)[1],
         dim(GA4.Avg.red2)[1], dim(GA4.Both.red1)[1], dim(GA4.Both.red2)[1], dim(GA4.BothMax.red1)[1], dim(GA4.BothMax.red2)[1])
)
colnames(GA4) <- c("Values","Strata","RB", "RSE", "RB", "RSE","RB", "RSE", "Reps")

GA4 %>%
  kbl(format = "html", booktabs = TRUE, caption = "Sim Results from GA5 designs", align = c(rep('l',2), rep('c',7))) %>%
  kable_styling(font_size = 24) %>%
  kable_classic_2(full_width = F) %>%
  column_spec(c(1,2,4,6,8), border_right = T) %>%
  collapse_rows(columns = 1) %>%
  row_spec(0, bold = T, extra_css = "border-bottom: 1px solid") %>%
  row_spec(c(2,4,6,8,10),  extra_css = "border-bottom: 1px solid") %>%
  add_header_above(header = c(" " = 2, "D" = 2, "L0" = 2, "Sigma" = 2, " " = 1), border_right = TRUE, bold = TRUE) 

####################
#GA5
####################

GA5 <- data.frame(
  C1 = c(rep("S1",2), rep("S2",2), rep("Avg",2), rep("Both",2), rep("Both (Max)",2)),
  C2 = c(rep(c("S1","S2"),5)),
  C3 = c(round(GA5.S1.D.est$`Rel bias (mean)`*100,1), round(GA5.S2.D.est$`Rel bias (mean)`*100,1), 
         round(GA5.Avg.D.est$`Rel bias (mean)`*100,1), round(GA5.Both.D.est$`Rel bias (mean)`*100,1), 
         round(GA5.BothMax.D.est$`Rel bias (mean)`*100,1)),
  C4 = c(round(GA5.S1.D.est$`Rel SE`,3), round(GA5.S2.D.est$`Rel SE`,3), round(GA5.Avg.D.est$`Rel SE`,3),
         round(GA5.Both.D.est$`Rel SE`,3), round(GA5.BothMax.D.est$`Rel SE`,3)),
  C5 = c(round(GA5.S1.L0.est$`Rel bias (mean)`*100,1), round(GA5.S2.L0.est$`Rel bias (mean)`*100,1), 
         round(GA5.Avg.L0.est$`Rel bias (mean)`*100,1), round(GA5.Both.L0.est$`Rel bias (mean)`*100,1), 
         round(GA5.BothMax.L0.est$`Rel bias (mean)`*100,1)),
  C6 = c(round(GA5.S1.L0.est$`Rel SE`,3), round(GA5.S2.L0.est$`Rel SE`,3), round(GA5.Avg.L0.est$`Rel SE`,3),
         round(GA5.Both.L0.est$`Rel SE`,3), round(GA5.BothMax.L0.est$`Rel SE`,3)),
  C7 = c(round(GA5.S1.Sig.est$`Rel bias (mean)`*100,1), round(GA5.S2.Sig.est$`Rel bias (mean)`*100,1), 
         round(GA5.Avg.Sig.est$`Rel bias (mean)`*100,1), round(GA5.Both.Sig.est$`Rel bias (mean)`*100,1), 
         round(GA5.Both.Sig.est$`Rel bias (mean)`*100,1)),
  C8 = c(round(GA5.S1.Sig.est$`Rel SE`,3), round(GA5.S2.Sig.est$`Rel SE`,3), round(GA5.Avg.Sig.est$`Rel SE`,3),
         round(GA5.Both.Sig.est$`Rel SE`,3), round(GA5.BothMax.Sig.est$`Rel SE`,3)),
  C9 = c(dim(GA5.S1.red1)[1], dim(GA5.S1.red2)[1], dim(GA5.S2.red1)[1], dim(GA5.S2.red2)[1], dim(GA5.Avg.red1)[1],
         dim(GA5.Avg.red2)[1], dim(GA5.Both.red1)[1], dim(GA5.Both.red2)[1], dim(GA5.BothMax.red1)[1], dim(GA5.BothMax.red2)[1])
)
colnames(GA5) <- c("Values","Strata","RB", "RSE", "RB", "RSE","RB", "RSE", "Reps") 

GA5 %>%
  kbl(format = "html", booktabs = TRUE, caption = "Sim Results from GA5 designs", align = c(rep('l',2), rep('c',7))) %>%
  kable_styling(font_size = 24) %>%
  kable_classic_2(full_width = F) %>%
  column_spec(c(1,2,4,6,8), border_right = T) %>%
  collapse_rows(columns = 1) %>%
  row_spec(0, bold = T, extra_css = "border-bottom: 1px solid") %>%
  row_spec(c(2,4,6,8,10),  extra_css = "border-bottom: 1px solid") %>%
  add_header_above(header = c(" " = 2, "D" = 2, "L0" = 2, "Sigma" = 2, " " = 1), border_right = TRUE, bold = TRUE) 



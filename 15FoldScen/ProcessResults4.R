#Oct 2025, process sim results for 15 fold scen
#For 40 traps and a large regular area
#this version to process results from fitting models outside secrdesign 
#uses 2nd set of 100 reps
#also filters out reps where se/est > 1000

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
summ.results <- function(results.df, par = "D", true.values, ylim.ceiling, plot = TRUE, se = FALSE){
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
      boxplot(results.df[,index], main = paste(plot.label1, "SE", sep = " "), ylim = c(0, ylim.ceiling[1]))
      boxplot(results.df[,index+6], main = paste(plot.label2, "SE", sep = " "), ylim = c(0, ylim.ceiling[2]))
    }
  }
  
  #Calc RSE
  S1.RSE <- sd(results.df[,index], na.rm = T) / mean(results.df[,index], na.rm = T) 
  S2.RSE <- sd(results.df[,index+6], na.rm = T) / mean(results.df[,index+6], na.rm = T) 
  
  res <- list("Rel bias (mean)" = c(S1.RB.mean, S2.RB.mean), "Rel bias (median)" = c(S1.RB.median, S2.RB.median), "Rel SE" = c(S1.RSE, S2.RSE))
  return(res)
}

#function to find rogue values for a single strata
#cannot do both together as diff strata cant be in the same row
#true values provided as D, L0, sigma
#detects NA and infinite values, 
#and also checks for any estimates greater or less than mag fold 
find.rogue <- function(df, mag = 10, true){
  bad.rows <- NULL
  df <- as.data.frame(df)
  colnames(df) <- c("D", "D.se", "L0", "L0.se", "Sig", "Sig.se")
  
  ##identify rows with NA or infinite
  rogue.rows <- which(is.na(rowSums(df))|is.infinite(rowSums(df)))
  if (length(rogue.rows)==0) rogue.rows <- NULL
  bad.rows <- c(bad.rows,rogue.rows)
  
  ##identify rows where the estimated se is v large relative to the est
  RSE.D <- df[,2] / df[,1] 
  RSE.L0 <- df[,4] / df[,3] 
  RSE.Sig <- df[,6] / df[,5]
  rogue.seD <- which(RSE.D > 100)
  rogue.seL0 <- which(RSE.L0 > 100)
  rogue.seSig <- which(RSE.Sig > 100)
  
  if (length(rogue.seD)>0) bad.rows <- c(bad.rows, rogue.seD)
  if (length(rogue.seL0)>0) bad.rows <- c(bad.rows, rogue.seL0)
  if (length(rogue.seSig)>0) bad.rows <- c(bad.rows, rogue.seSig)
  
#  if (length(rogue.rows)>0) df <- df[-rogue.rows,]
  
  #identify rows with wild values for any estimates (mag fold)
  D.rogue <- which(df$D > true[1]*mag | df$D < true[1]/mag)
  L0.rogue <- which(df$L0 > true[2]*mag | df$L0 < true[2]/mag)
  Sig.rogue <- which(df$Sig > true[3]*mag | df$Sig < true[3]/mag)
  
  if (length(D.rogue)>0) bad.rows <- c(bad.rows, D.rogue)
  if (length(L0.rogue)>0) bad.rows <- c(bad.rows, L0.rogue)
  if (length(Sig.rogue)>0) bad.rows <- c(bad.rows, Sig.rogue)

  if (!is.null(bad.rows)) {
    bad.rows <- sort(unique(bad.rows)) 
    df <- df[-bad.rows,]
  }
  
  df$row <- row.names(df)
  return("Clean data" = df)
}

#diff version of the fn, more efficient
find.rogue <- function(df, mag = 10, true){
  bad.rows <- NULL
  df <- as.data.frame(df)
  colnames(df) <- c("D", "D.se", "L0", "L0.se", "Sig", "Sig.se")
  
  ##identify rows with NA or infinite
  inds1 <- which(is.na(rowSums(df))|is.infinite(rowSums(df)))
  
  ##identify rows where the estimated se is v large relative to the est
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
Sig1.ylim <- 350 ; Sig2.ylim <- 7500

################
#Grid designs
################

Grid.800.Data200 <- vector(mode = "list", length = 200)
Grid.1600.Data100 <- vector(mode = "list", length = 100)

load("15FoldScen/Cluster/Sims/GridsbResults40.RData")

Grid.800.Data200[1:100] <- Grid.800.Data$output
Grid.800.Results.temp <- as.data.frame(Grid.800.results)

load("15FoldScen/Cluster/Sims/GridscResults40.RData")

Grid.800.Data200[101:200] <- Grid.800.Data$output
Grid.800.Results200 <- rbind(Grid.800.Results.temp, as.data.frame(Grid.800.results))

Grid.1600.Data100 <- Grid.1600.Data$output
Grid.1600.Results100 <- as.data.frame(Grid.1600.results)

###############
#800 m spacing
###############

##first data summary, then estimates
Grid.800.data.summ <- CH.data.summs(Grid.800.Data200)

pdf("15FoldScen/Grid800DataSumms.pdf", height = 8, width = 10, pointsize = 11)
boxplot(Grid.800.data.summ)

#filter data, NAs and infinite, and wildly out values (20 fold)
#one strata at a time
Grid.800.red1 <- find.rogue(Grid.800.Results200[,c(1:6)], mag = mag.factor, true = c(0.05, 2, 200))
Grid.800.red2 <- find.rogue(Grid.800.Results200[,c(7:12)], mag = mag.factor, true = c(0.05/15, 2/15, 3000))

Grid.800.clean <- merge(Grid.800.red1, Grid.800.red2, by = "row", all.x = TRUE, all.y = TRUE)
summary(Grid.800.clean)

#plot and summarise
#estimates first
Grid800.D.est <- summ.results(Grid.800.clean, par = "D", true.value = c(0.05,0.05/15) , ylim.ceiling = c(D1.ylim, D2.ylim), plot = TRUE, se = FALSE)
Grid800.L0.est <- summ.results(Grid.800.clean, par = "L0", true.value = c(2, 2/15) , ylim.ceiling = c(L1.ylim, L2.ylim), plot = TRUE, se = FALSE)
Grid800.Sig.est <- summ.results(Grid.800.clean, par = "Sigma", true.value = c(200, 3000) , ylim.ceiling = c(Sig1.ylim, Sig2.ylim), plot = TRUE, se = FALSE)

#standard errors
Grid800.D.se <- summ.results(Grid.800.clean, par = "D", true.value = c(0.05,0.05/15) , ylim.ceiling = c(0.1,0.005), plot = TRUE, se = TRUE)
Grid800.L0.se <- summ.results(Grid.800.clean, par = "L0", true.value = c(2, 2/15) , ylim.ceiling = c(0.75,0.05), plot = TRUE, se = TRUE)
Grid800.Sig.se <- summ.results(Grid.800.clean, par = "Sigma", true.value = c(200, 3000) , ylim.ceiling = c(50,7500), plot = TRUE, se = TRUE)

#############
#1600
#############

##first data summary, then estimates
Grid.1600.data.summ <- CH.data.summs(Grid.1600.Data100)

pdf("15FoldScen/Grid1600DataSumms.pdf", height = 8, width = 10, pointsize = 11)
boxplot(Grid.1600.data.summ)

#filter data, NAs and infinite, and wildly out values (20 fold)
#one strata at a time
Grid.1600.red1 <- find.rogue(Grid.1600.Results100[,c(1:6)], mag = mag.factor, true = c(0.05, 2, 200))
Grid.1600.red2 <- find.rogue(Grid.1600.Results100[,c(7:12)], mag = mag.factor, true = c(0.05/15, 2/15, 3000))

Grid.1600.clean <- merge(Grid.1600.red1, Grid.1600.red2, by = "row", all.x = TRUE, all.y = TRUE)
summary(Grid.1600.clean)

#plot and summarise
#estimates first
Grid1600.D.est <- summ.results(Grid.1600.clean, par = "D", true.value = c(0.05,0.05/15) , ylim.ceiling = c(D1.ylim, D2.ylim), plot = TRUE, se = FALSE)
Grid1600.L0.est <- summ.results(Grid.1600.clean, par = "L0", true.value = c(2, 2/15) , ylim.ceiling = c(L1.ylim, L2.ylim), plot = TRUE, se = FALSE)
Grid1600.Sig.est <- summ.results(Grid.1600.clean, par = "Sigma", true.value = c(200, 3000) , ylim.ceiling = c(Sig1.ylim, Sig2.ylim), plot = TRUE, se = FALSE)

#standard errors
Grid1600.D.se <- summ.results(Grid.1600.clean, par = "D", true.value = c(0.05,0.05/15) , ylim.ceiling = c(210,0.005), plot = TRUE, se = TRUE)
Grid1600.L0.se <- summ.results(Grid.1600.clean, par = "L0", true.value = c(2, 2/15) , ylim.ceiling = c(0.75,0.05), plot = TRUE, se = TRUE)
Grid1600.Sig.se <- summ.results(Grid.1600.clean, par = "Sigma", true.value = c(200, 3000) , ylim.ceiling = c(1750,1300), plot = TRUE, se = TRUE)

#######################################################################

################
#GA4 designs
################

#extract and collate results from two objs with 100 each 

G4.S1.Data200 <- vector(mode = "list", length = 200)
G4.S2.Data200 <- vector(mode = "list", length = 200)
G4.Avg.Data200 <- vector(mode = "list", length = 200)
G4.Both.Data200 <- vector(mode = "list", length = 200)
G4.BothMax.Data200 <- vector(mode = "list", length = 200)

load("15FoldScen/Cluster/Sims/GA4bResults40.RData")

G4.S1.Data200[1:100] <- G4.S1.Data$output
G4.S1.Results.temp <- as.data.frame(G4.S1.results)

G4.S2.Data200[1:100] <- G4.S2.Data$output
G4.S2.Results.temp <- as.data.frame(G4.S2.results)

G4.Avg.Data200[1:100] <- G4.Avg.Data$output
G4.Avg.Results.temp <- as.data.frame(G4.Avg.results)

G4.Both.Data200[1:100] <- G4.Both.Data$output
G4.Both.Results.temp <- as.data.frame(G4.Both.results)

G4.BothMax.Data200[1:100] <- G4.BothMax.Data$output
G4.BothMax.Results.temp <- as.data.frame(G4.BothMax.results)

load("15FoldScen/Cluster/Sims/GA4cResults40.RData")

G4.S1.Data200[101:200] <- G4.S1.Data$output
G4.S1.Results200 <- rbind(G4.S1.Results.temp, as.data.frame(G4.S1.results))

G4.S2.Data200[101:200] <- G4.S2.Data$output
G4.S2.Results200 <- rbind(G4.S2.Results.temp, as.data.frame(G4.S2.results))

G4.Avg.Data200[101:200] <- G4.Avg.Data$output
G4.Avg.Results200 <- rbind(G4.Avg.Results.temp, as.data.frame(G4.Avg.results))

G4.Both.Data200[101:200] <- G4.Both.Data$output
G4.Both.Results200 <- rbind(G4.Both.Results.temp, as.data.frame(G4.Both.results))

G4.BothMax.Data200[101:200] <- G4.BothMax.Data$output
G4.BothMax.Results200 <- rbind(G4.BothMax.Results.temp, as.data.frame(G4.BothMax.results))

######################
#GA4 with S1 pars
######################

GA4.S1.data.summ <- CH.data.summs(G4.S1.Data200)
boxplot(GA4.S1.data.summ)

#filter estimates, NAs and infinite, and wildly out values (10 fold)
#one strata at a time
GA4.S1.red1 <- find.rogue(G4.S1.Results200[,c(1:6)], mag = mag.factor, true = c(0.05, 2, 200))
GA4.S1.red2 <- find.rogue(G4.S1.Results200[,c(7:12)], mag = mag.factor, true = c(0.05/15, 2/15, 3000))

GA4.S1.clean <- merge(GA4.S1.red1, GA4.S1.red2, by = "row", all.x = TRUE, all.y = TRUE)
summary(GA4.S1.clean)

#plot and summarise
#estimates first
GA4.S1.D.est <- summ.results(GA4.S1.clean, par = "D", true.value = c(0.05,0.05/15) , ylim.ceiling = c(D1.ylim, D2.ylim), plot = TRUE, se = FALSE)
GA4.S1.L0.est <- summ.results(GA4.S1.clean, par = "L0", true.value = c(2, 2/15) , ylim.ceiling = c(L1.ylim, L2.ylim), plot = TRUE, se = FALSE)
GA4.S1.Sig.est <- summ.results(GA4.S1.clean, par = "Sigma", true.value = c(200, 3000) , ylim.ceiling = c(Sig1.ylim, Sig2.ylim), plot = TRUE, se = FALSE)

#standard errors
GA4.S1.D.se <- summ.results(GA4.S1.clean, par = "D", true.value = c(0.05,0.05/15) , ylim.ceiling = c(0.05,0.005), plot = TRUE, se = TRUE)
GA4.S1.L0.se <- summ.results(GA4.S1.clean, par = "L0", true.value = c(2, 2/15) , ylim.ceiling = c(0.5,0.05), plot = TRUE, se = TRUE)
GA4.S1.Sig.se <- summ.results(GA4.S1.clean, par = "Sigma", true.value = c(200, 3000) , ylim.ceiling = c(60,600), plot = TRUE, se = TRUE)

######################
#GA4 with S2 pars
######################

GA4.S2.data.summ <- CH.data.summs(G4.S2.Data200)
boxplot(GA4.S2.data.summ)

#filter estimates, NAs and infinite, and wildly out values (10 fold)
#one strata at a time
GA4.S2.red1 <- find.rogue(G4.S2.Results200[,c(1:6)], mag = mag.factor, true = c(0.05, 2, 200))
GA4.S2.red2 <- find.rogue(G4.S2.Results200[,c(7:12)], mag = mag.factor, true = c(0.05/15, 2/15, 3000))

GA4.S2.clean <- merge(GA4.S2.red1, GA4.S2.red2, by = "row", all.x = TRUE, all.y = TRUE)
summary(GA4.S2.clean)

#plot and summarise
#estimates first
GA4.S2.D.est <- summ.results(GA4.S2.clean, par = "D", true.value = c(0.05,0.05/15) , ylim.ceiling = c(D1.ylim, D2.ylim), plot = TRUE, se = FALSE)
GA4.S2.L0.est <- summ.results(GA4.S2.clean, par = "L0", true.value = c(2, 2/15) , ylim.ceiling = c(L1.ylim, L2.ylim), plot = TRUE, se = FALSE)
GA4.S2.Sig.est <- summ.results(GA4.S2.clean, par = "Sigma", true.value = c(200, 3000) , ylim.ceiling = c(Sig1.ylim, Sig2.ylim), plot = TRUE, se = FALSE)

#standard errors
GA4.S2.D.se <- summ.results(GA4.S2.clean, par = "D", true.value = c(0.05,0.05/15) , ylim.ceiling = c(0.05,0.005), plot = TRUE, se = TRUE)
GA4.S2.L0.se <- summ.results(GA4.S2.clean, par = "L0", true.value = c(2, 2/15) , ylim.ceiling = c(0.5,0.05), plot = TRUE, se = TRUE)
GA4.S2.Sig.se <- summ.results(GA4.S2.clean, par = "Sigma", true.value = c(200, 3000) , ylim.ceiling = c(35,2600), plot = TRUE, se = TRUE)

######################
#GA4 with avg pars
######################

GA4.Avg.data.summ <- CH.data.summs(G4.Avg.Data200)
boxplot(GA4.Avg.data.summ)

#filter estimates, NAs and infinite, and wildly out values (10 fold)
#one strata at a time
GA4.Avg.red1 <- find.rogue(G4.Avg.Results200[,c(1:6)], mag = mag.factor, true = c(0.05, 2, 200))
GA4.Avg.red2 <- find.rogue(G4.Avg.Results200[,c(7:12)], mag = mag.factor, true = c(0.05/15, 2/15, 3000))

GA4.Avg.clean <- merge(GA4.Avg.red1, GA4.Avg.red2, by = "row", all.x = TRUE, all.y = TRUE)
summary(GA4.Avg.clean)

#plot and summarise
#estimates first
GA4.Avg.D.est <- summ.results(GA4.Avg.clean, par = "D", true.value = c(0.05,0.05/15) , ylim.ceiling = c(D1.ylim, D2.ylim), plot = TRUE, se = FALSE)
GA4.Avg.L0.est <- summ.results(GA4.Avg.clean, par = "L0", true.value = c(2, 2/15) , ylim.ceiling = c(L1.ylim, L2.ylim), plot = TRUE, se = FALSE)
GA4.Avg.Sig.est <- summ.results(GA4.Avg.clean, par = "Sigma", true.value = c(200, 3000) , ylim.ceiling = c(Sig1.ylim, Sig2.ylim), plot = TRUE, se = FALSE)

#standard errors
GA4.Avg.D.se <- summ.results(GA4.Avg.clean, par = "D", true.value = c(0.05,0.05/15) , ylim.ceiling = c(0.15,0.005), plot = TRUE, se = TRUE)
GA4.Avg.L0.se <- summ.results(GA4.Avg.clean, par = "L0", true.value = c(2, 2/15) , ylim.ceiling = c(0.5,0.05), plot = TRUE, se = TRUE)
GA4.Avg.Sig.se <- summ.results(GA4.Avg.clean, par = "Sigma", true.value = c(200, 3000) , ylim.ceiling = c(150,650), plot = TRUE, se = TRUE)

######################
#GA4 with both
######################

GA4.Both.data.summ <- CH.data.summs(G4.Both.Data200)
boxplot(GA4.Both.data.summ)

#filter estimates, NAs and infinite, and wildly out values (10 fold)
#one strata at a time
GA4.Both.red1 <- find.rogue(G4.Both.Results200[,c(1:6)], mag = mag.factor, true = c(0.05, 2, 200))
GA4.Both.red2 <- find.rogue(G4.Both.Results200[,c(7:12)], mag = mag.factor, true = c(0.05/15, 2/15, 3000))

GA4.Both.clean <- merge(GA4.Both.red1, GA4.Both.red2, by = "row", all.x = TRUE, all.y = TRUE)
summary(GA4.Both.clean)

#plot and summarise
#estimates first
GA4.Both.D.est <- summ.results(GA4.Both.clean, par = "D", true.value = c(0.05,0.05/15) , ylim.ceiling = c(D1.ylim, D2.ylim), plot = TRUE, se = FALSE)
GA4.Both.L0.est <- summ.results(GA4.Both.clean, par = "L0", true.value = c(2, 2/15) , ylim.ceiling = c(L1.ylim, L2.ylim), plot = TRUE, se = FALSE)
GA4.Both.Sig.est <- summ.results(GA4.Both.clean, par = "Sigma", true.value = c(200, 3000) , ylim.ceiling = c(Sig1.ylim, Sig2.ylim), plot = TRUE, se = FALSE)

#standard errors
GA4.Both.D.se <- summ.results(GA4.Both.clean, par = "D", true.value = c(0.05,0.05/15) , ylim.ceiling = c(0.05,0.005), plot = TRUE, se = TRUE)
GA4.Both.L0.se <- summ.results(GA4.Both.clean, par = "L0", true.value = c(2, 2/15) , ylim.ceiling = c(0.55,0.05), plot = TRUE, se = TRUE)
GA4.Both.Sig.se <- summ.results(GA4.Both.clean, par = "Sigma", true.value = c(200, 3000) , ylim.ceiling = c(30,750), plot = TRUE, se = TRUE)

######################
#GA4 with both (max = T)
######################

GA4.BothMax.data.summ <- CH.data.summs(G4.BothMax.Data200)
boxplot(GA4.BothMax.data.summ)

#filter estimates, NAs and infinite, and wildly out values (10 fold)
#one strata at a time
GA4.BothMax.red1 <- find.rogue(G4.BothMax.Results200[,c(1:6)], mag = mag.factor, true = c(0.05, 2, 200))
GA4.BothMax.red2 <- find.rogue(G4.BothMax.Results200[,c(7:12)], mag = mag.factor, true = c(0.05/15, 2/15, 3000))

GA4.BothMax.clean <- merge(GA4.BothMax.red1, GA4.BothMax.red2, by = "row", all.x = TRUE, all.y = TRUE)
summary(GA4.BothMax.clean)

#plot and summarise
#estimates first
GA4.BothMax.D.est <- summ.results(GA4.BothMax.clean, par = "D", true.value = c(0.05,0.05/15) , ylim.ceiling = c(D1.ylim, D2.ylim), plot = TRUE, se = FALSE)
GA4.BothMax.L0.est <- summ.results(GA4.BothMax.clean, par = "L0", true.value = c(2, 2/15) , ylim.ceiling = c(L1.ylim, L2.ylim), plot = TRUE, se = FALSE)
GA4.BothMax.Sig.est <- summ.results(GA4.BothMax.clean, par = "Sigma", true.value = c(200, 3000) , ylim.ceiling = c(Sig1.ylim, Sig2.ylim), plot = TRUE, se = FALSE)

#standard errors
GA4.BothMax.D.se <- summ.results(GA4.BothMax.clean, par = "D", true.value = c(0.05,0.05/15) , ylim.ceiling = c(0.05,0.005), plot = TRUE, se = TRUE)
GA4.BothMax.L0.se <- summ.results(GA4.BothMax.clean, par = "L0", true.value = c(2, 2/15) , ylim.ceiling = c(0.5,0.05), plot = TRUE, se = TRUE)
GA4.BothMax.Sig.se <- summ.results(GA4.BothMax.clean, par = "Sigma", true.value = c(200, 3000) , ylim.ceiling = c(30,10000), plot = TRUE, se = TRUE)

################################
#GA5 designs
################################

#extract and collate results from two objs with 100 each 

G5.S1.Data200 <- vector(mode = "list", length = 200)
G5.S2.Data200 <- vector(mode = "list", length = 200)
G5.Avg.Data200 <- vector(mode = "list", length = 200)
G5.Both.Data200 <- vector(mode = "list", length = 200)
G5.BothMax.Data200 <- vector(mode = "list", length = 200)

load("15FoldScen/Cluster/Sims/GA5bResults40.RData")

G5.S1.Data200[1:100] <- G5.S1.Data$output
G5.S1.Results.temp <- as.data.frame(G5.S1.results)

G5.S2.Data200[1:100] <- G5.S2.Data$output
G5.S2.Results.temp <- as.data.frame(G5.S2.results)

G5.Avg.Data200[1:100] <- G5.Avg.Data$output
G5.Avg.Results.temp <- as.data.frame(G5.Avg.results)

G5.Both.Data200[1:100] <- G5.Both.Data$output
G5.Both.Results.temp <- as.data.frame(G5.Both.results)

G5.BothMax.Data200[1:100] <- G5.BothMax.Data$output
G5.BothMax.Results.temp <- as.data.frame(G5.BothMax.results)

load("15FoldScen/Cluster/Sims/GA5cResults40.RData")

G5.S1.Data200[101:200] <- G5.S1.Data$output
G5.S1.Results200 <- rbind(G5.S1.Results.temp, as.data.frame(G5.S1.results))

G5.S2.Data200[101:200] <- G5.S2.Data$output
G5.S2.Results200 <- rbind(G5.S2.Results.temp, as.data.frame(G5.S2.results))

G5.Avg.Data200[101:200] <- G5.Avg.Data$output
G5.Avg.Results200 <- rbind(G5.Avg.Results.temp, as.data.frame(G5.Avg.results))

G5.Both.Data200[101:200] <- G5.Both.Data$output
G5.Both.Results200 <- rbind(G5.Both.Results.temp, as.data.frame(G5.Both.results))

G5.BothMax.Data200[101:200] <- G5.BothMax.Data$output
G5.BothMax.Results200 <- rbind(G5.BothMax.Results.temp, as.data.frame(G5.BothMax.results))

######################
#GA5 with S1 pars
######################

GA5.S1.data.summ <- CH.data.summs(G5.S1.Data200)
boxplot(GA5.S1.data.summ)

#filter estimates, NAs and infinite, and wildly out values (10 fold)
#one strata at a time
GA5.S1.red1 <- find.rogue(G5.S1.Results200[,c(1:6)], mag = mag.factor, true = c(0.05, 2, 200))
GA5.S1.red2 <- find.rogue(G5.S1.Results200[,c(7:12)], mag = mag.factor, true = c(0.05/15, 2/15, 3000))

GA5.S1.clean <- merge(GA5.S1.red1, GA5.S1.red2, by = "row", all.x = TRUE, all.y = TRUE)
summary(GA5.S1.clean)

#plot and summarise
#estimates first
GA5.S1.D.est <- summ.results(GA5.S1.clean, par = "D", true.value = c(0.05,0.05/15) , ylim.ceiling = c(D1.ylim, D2.ylim), plot = TRUE, se = FALSE)
GA5.S1.L0.est <- summ.results(GA5.S1.clean, par = "L0", true.value = c(2, 2/15) , ylim.ceiling = c(L1.ylim, L2.ylim), plot = TRUE, se = FALSE)
GA5.S1.Sig.est <- summ.results(GA5.S1.clean, par = "Sigma", true.value = c(200, 3000) , ylim.ceiling = c(Sig1.ylim, Sig2.ylim), plot = TRUE, se = FALSE)

#standard errors
GA5.S1.D.se <- summ.results(GA5.S1.clean, par = "D", true.value = c(0.05,0.05/15) , ylim.ceiling = c(0.05,0.005), plot = TRUE, se = TRUE)
GA5.S1.L0.se <- summ.results(GA5.S1.clean, par = "L0", true.value = c(2, 2/15) , ylim.ceiling = c(0.6,0.055), plot = TRUE, se = TRUE)
GA5.S1.Sig.se <- summ.results(GA5.S1.clean, par = "Sigma", true.value = c(200, 3000) , ylim.ceiling = c(100,900), plot = TRUE, se = TRUE)

######################
#GA5 with S2 pars
######################

GA5.S2.data.summ <- CH.data.summs(G5.S2.Data200)
boxplot(GA5.S2.data.summ)

#filter estimates, NAs and infinite, and wildly out values (10 fold)
#one strata at a time
GA5.S2.red1 <- find.rogue(G5.S2.Results200[,c(1:6)], mag = mag.factor, true = c(0.05, 2, 200))
GA5.S2.red2 <- find.rogue(G5.S2.Results200[,c(7:12)], mag = mag.factor, true = c(0.05/15, 2/15, 3000))

GA5.S2.clean <- merge(GA5.S2.red1, GA5.S2.red2, by = "row", all.x = TRUE, all.y = TRUE)
summary(GA5.S2.clean)

#plot and summarise
#estimates first
GA5.S2.D.est  <- summ.results(GA5.S2.clean, par = "D", true.value = c(0.05,0.05/15) , ylim.ceiling = c(D1.ylim, D2.ylim), plot = TRUE, se = FALSE)
GA5.S2.L0.est  <- summ.results(GA5.S2.clean, par = "L0", true.value = c(2, 2/15) , ylim.ceiling = c(L1.ylim, L2.ylim), plot = TRUE, se = FALSE)
GA5.S2.Sig.est  <- summ.results(GA5.S2.clean, par = "Sigma", true.value = c(200, 3000) , ylim.ceiling = c(Sig1.ylim, Sig2.ylim), plot = TRUE, se = FALSE)

#standard errors
GA5.S2.D.se  <- summ.results(GA5.S2.clean, par = "D", true.value = c(0.05,0.05/15) , ylim.ceiling = c(0.05,0.005), plot = TRUE, se = TRUE)
GA5.S2.L0.se  <- summ.results(GA5.S2.clean, par = "L0", true.value = c(2, 2/15) , ylim.ceiling = c(0.8,0.05), plot = TRUE, se = TRUE)
GA5.S2.Sig.se  <- summ.results(GA5.S2.clean, par = "Sigma", true.value = c(200, 3000) , ylim.ceiling = c(25,1400), plot = TRUE, se = TRUE)

######################
#GA5 with avg pars
######################

GA5.Avg.data.summ <- CH.data.summs(G5.Avg.Data200)
boxplot(GA5.Avg.data.summ)

#filter estimates, NAs and infinite, and wildly out values (10 fold)
#one strata at a time
GA5.Avg.red1 <- find.rogue(G5.Avg.Results200[,c(1:6)], mag = mag.factor, true = c(0.05, 2, 200))
GA5.Avg.red2 <- find.rogue(G5.Avg.Results200[,c(7:12)], mag = mag.factor, true = c(0.05/15, 2/15, 3000))

GA5.Avg.clean <- merge(GA5.Avg.red1, GA5.Avg.red2, by = "row", all.x = TRUE, all.y = TRUE)
summary(GA5.Avg.clean)

#plot and summarise
#estimates first
GA5.Avg.D.est  <- summ.results(GA5.Avg.clean, par = "D", true.value = c(0.05,0.05/15) , ylim.ceiling = c(D1.ylim, D2.ylim), plot = TRUE, se = FALSE)
GA5.Avg.L0.est  <- summ.results(GA5.Avg.clean, par = "L0", true.value = c(2, 2/15) , ylim.ceiling = c(L1.ylim, L2.ylim), plot = TRUE, se = FALSE)
GA5.Avg.Sig.est  <- summ.results(GA5.Avg.clean, par = "Sigma", true.value = c(200, 3000) , ylim.ceiling = c(Sig1.ylim, Sig2.ylim), plot = TRUE, se = FALSE)

#standard errors
GA5.Avg.D.se <- summ.results(GA5.Avg.clean, par = "D", true.value = c(0.05,0.05/15) , ylim.ceiling = c(0.15,0.005), plot = TRUE, se = TRUE)
GA5.Avg.L0.se <- summ.results(GA5.Avg.clean, par = "L0", true.value = c(2, 2/15) , ylim.ceiling = c(0.5,0.05), plot = TRUE, se = TRUE)
GA5.Avg.Sig.se <- summ.results(GA5.Avg.clean, par = "Sigma", true.value = c(200, 3000) , ylim.ceiling = c(150,650), plot = TRUE, se = TRUE)

######################
#GA5 with both
######################

GA5.Both.data.summ <- CH.data.summs(G5.Both.Data200)
boxplot(GA5.Both.data.summ)

#filter estimates, NAs and infinite, and wildly out values (10 fold)
#one strata at a time
GA5.Both.red1 <- find.rogue(G5.Both.Results200[,c(1:6)], mag = mag.factor, true = c(0.05, 2, 200))
GA5.Both.red2 <- find.rogue(G5.Both.Results200[,c(7:12)], mag = mag.factor, true = c(0.05/15, 2/15, 3000))

GA5.Both.clean <- merge(GA5.Both.red1, GA5.Both.red2, by = "row", all.x = TRUE, all.y = TRUE)
summary(GA5.Both.clean)

#plot and summarise
#estimates first
GA5.Both.D.est <- summ.results(GA5.Both.clean, par = "D", true.value = c(0.05,0.05/15) , ylim.ceiling = c(D1.ylim, D2.ylim), plot = TRUE, se = FALSE)
GA5.Both.L0.est <- summ.results(GA5.Both.clean, par = "L0", true.value = c(2, 2/15) , ylim.ceiling = c(L1.ylim, L2.ylim), plot = TRUE, se = FALSE)
GA5.Both.Sig.est <- summ.results(GA5.Both.clean, par = "Sigma", true.value = c(200, 3000) , ylim.ceiling = c(Sig1.ylim, Sig2.ylim), plot = TRUE, se = FALSE)

#standard errors
GA5.Both.D.se <- summ.results(GA5.Both.clean, par = "D", true.value = c(0.05,0.05/15) , ylim.ceiling = c(0.05,0.005), plot = TRUE, se = TRUE)
GA5.Both.L0.se <- summ.results(GA5.Both.clean, par = "L0", true.value = c(2, 2/15) , ylim.ceiling = c(0.55,0.05), plot = TRUE, se = TRUE)
GA5.Both.Sig.se <- summ.results(GA5.Both.clean, par = "Sigma", true.value = c(200, 3000) , ylim.ceiling = c(30,750), plot = TRUE, se = TRUE)

######################
#GA5 with both (max = T)
######################

GA5.BothMax.data.summ <- CH.data.summs(G5.BothMax.Data200)
boxplot(GA5.BothMax.data.summ)

#filter estimates, NAs and infinite, and wildly out values (10 fold)
#one strata at a time
GA5.BothMax.red1 <- find.rogue(G5.BothMax.Results200[,c(1:6)], mag = mag.factor, true = c(0.05, 2, 200))
GA5.BothMax.red2 <- find.rogue(G5.BothMax.Results200[,c(7:12)], mag = mag.factor, true = c(0.05/15, 2/15, 3000))

GA5.BothMax.clean <- merge(GA5.BothMax.red1, GA5.BothMax.red2, by = "row", all.x = TRUE, all.y = TRUE)
summary(GA5.BothMax.clean)

#plot and summarise
#estimates first
GA5.BothMax.D.est <- summ.results(GA5.BothMax.clean, par = "D", true.value = c(0.05,0.05/15) , ylim.ceiling = c(D1.ylim, D2.ylim), plot = TRUE, se = FALSE)
GA5.BothMax.L0.est <- summ.results(GA5.BothMax.clean, par = "L0", true.value = c(2, 2/15) , ylim.ceiling = c(L1.ylim, L2.ylim), plot = TRUE, se = FALSE)
GA5.BothMax.Sig.est <- summ.results(GA5.BothMax.clean, par = "Sigma", true.value = c(200, 3000) , ylim.ceiling = c(Sig1.ylim, Sig2.ylim), plot = TRUE, se = FALSE)

#standard errors
GA5.BothMax.D.se <- summ.results(GA5.BothMax.clean, par = "D", true.value = c(0.05,0.05/15) , ylim.ceiling = c(0.05,0.01), plot = TRUE, se = TRUE)
GA5.BothMax.L0.se <- summ.results(GA5.BothMax.clean, par = "L0", true.value = c(2, 2/15) , ylim.ceiling = c(0.7,0.05), plot = TRUE, se = TRUE)
GA5.BothMax.Sig.se <- summ.results(GA5.BothMax.clean, par = "Sigma", true.value = c(200, 3000) , ylim.ceiling = c(100,800), plot = TRUE, se = TRUE)

########################################################################
#save R workspace

save.image(file = "15FoldScen/SimResultsObjs.RData")

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


#save all results objects into a list for easy loading
#version 3 has 200 reps and includes 100 reps for 1600 m grid spacing

datasumms <- list("Grid.800" = Grid.800.data.summ, "Grid.1600" = Grid.1600.data.summ, 
                  "G4.S1" = GA4.S1.data.summ, "G4.S2" = GA4.S2.data.summ, "G4.Avg" = GA4.Avg.data.summ, "G4.Both" = GA4.Both.data.summ, "G4.BothMax" = GA4.BothMax.data.summ,
                  "G5.S1" = GA5.S1.data.summ, "G5.S2" = GA5.S2.data.summ, "G5.Avg" = GA5.Avg.data.summ, "G5.Both" = GA5.Both.data.summ, "G5.BothMax" = GA5.BothMax.data.summ)
save(datasumms, file = "15FoldScen/DataSumms3.RData")

simresults <- list("Grid800.D" = Grid800.D.est, "Grid800.L0" = Grid800.L0.est, "Grid800.Sig" = Grid800.Sig.est,
                   "Grid1600.D" = Grid1600.D.est, "Grid1600.L0" = Grid1600.L0.est, "Grid1600.Sig" = Grid1600.Sig.est,
                   "GA4.S1.D" = GA4.S1.D.est, "GA4.S1.L0" = GA4.S1.L0.est, "GA4.S1.Sig" = GA4.S1.Sig.est,
                   "GA4.S2.D" = GA4.S2.D.est, "GA4.S2.L0" = GA4.S2.L0.est, "GA4.S2.Sig" = GA4.S2.Sig.est,
                   "GA4.Avg.D" = GA4.Avg.D.est, "GA4.Avg.L0" = GA4.Avg.L0.est, "GA4.Avg.Sig" = GA4.Avg.Sig.est,
                   "GA4.Both.D" = GA4.Both.D.est, "GA4.Both.L0" = GA4.Both.L0.est, "GA4.Both.Sig" = GA4.Both.Sig.est,
                   "GA4.BothMax.D" = GA4.BothMax.D.est, "GA4.BothMax.L0" = GA4.BothMax.L0.est, "GA4.BothMax.Sig" = GA4.BothMax.Sig.est,
                   "GA5.S1.D" = GA5.S1.D.est, "GA5.S1.L0" = GA5.S1.L0.est, "GA5.S1.Sig" = GA5.S1.Sig.est,
                   "GA5.S2.D" = GA5.S2.D.est, "GA5.S2.L0" = GA5.S2.L0.est, "GA5.S2.Sig" = GA5.S2.Sig.est,
                   "GA5.Avg.D" = GA5.Avg.D.est, "GA5.Avg.L0" = GA5.Avg.L0.est, "GA5.Avg.Sig" = GA5.Avg.Sig.est,
                   "GA5.Both.D" = GA5.Both.D.est, "GA5.Both.L0" = GA5.Both.L0.est, "GA5.Both.Sig" = GA5.Both.Sig.est,
                   "GA5.BothMax.D" = GA5.BothMax.D.est, "GA5.BothMax.L0" = GA5.BothMax.L0.est, "GA5.BothMax.Sig" = GA5.BothMax.Sig.est)
save(simresults, file = "15FoldScen/SimResults3.RData")

seresults <- list("Grid800.D" = Grid800.D.se, "Grid800.L0" = Grid800.L0.se, "Grid800.Sig" = Grid800.Sig.se,
                  "Grid1600.D" = Grid1600.D.se, "Grid1600.L0" = Grid1600.L0.se, "Grid1600.Sig" = Grid1600.Sig.se,
                  "GA4.S1.D" = GA4.S1.D.se, "GA4.S1.L0" = GA4.S1.L0.se, "GA4.S1.Sig" = GA4.S1.Sig.se,
                  "GA4.S2.D" = GA4.S2.D.se, "GA4.S2.L0" = GA4.S2.L0.se, "GA4.S2.Sig" = GA4.S2.Sig.se,
                  "GA4.Avg.D" = GA4.Avg.D.se, "GA4.Avg.L0" = GA4.Avg.L0.se, "GA4.Avg.Sig" = GA4.Avg.Sig.se,
                  "GA4.Both.D" = GA4.Both.D.se, "GA4.Both.L0" = GA4.Both.L0.se, "GA4.Both.Sig" = GA4.Both.Sig.se,
                  "GA4.BothMax.D" = GA4.BothMax.D.se, "GA4.BothMax.L0" = GA4.BothMax.L0.se, "GA4.BothMax.Sig" = GA4.BothMax.Sig.se,
                  "GA5.S1.D" = GA5.S1.D.se, "GA5.S1.L0" = GA5.S1.L0.se, "GA5.S1.Sig" = GA5.S1.Sig.se,
                  "GA5.S2.D" = GA5.S2.D.se, "GA5.S2.L0" = GA5.S2.L0.se, "GA5.S2.Sig" = GA5.S2.Sig.se,
                  "GA5.Avg.D" = GA5.Avg.D.se, "GA5.Avg.L0" = GA5.Avg.L0.se, "GA5.Avg.Sig" = GA5.Avg.Sig.se,
                  "GA5.Both.D" = GA5.Both.D.se, "GA5.Both.L0" = GA5.Both.L0.se, "GA5.Both.Sig" = GA5.Both.Sig.se,
                  "GA5.BothMax.D" = GA5.BothMax.D.se, "GA5.BothMax.L0" = GA5.BothMax.L0.se, "GA5.BothMax.Sig" = GA5.BothMax.Sig.se)

save(seresults, file = "15FoldScen/seResults3.RData")

excluded <- list("Grid800.S1" = Grid.800.red1, "Grid800.S2" = Grid.800.red2,
                 "Grid1600.S1" = Grid.1600.red1, "Grid1600.S2" = Grid.1600.red2,
                 "GA4.S1.S1" = GA4.S1.red1, "GA4.S1.S2" = GA4.S1.red2,
                 "GA4.S2.S1" = GA4.S2.red1, "GA4.S2.S2" = GA4.S2.red2,
                 "GA4.Avg.S1" = GA4.Avg.red1, "GA4.Avg.S2" = GA4.Avg.red2,
                 "GA4.Both.S1" = GA4.Both.red1, "GA4.Both.S2" = GA4.Both.red2,
                 "GA4.BothMax.S1" = GA4.BothMax.red1, "GA4.BothMax.S2" = GA4.BothMax.red2,
                 "GA5.S1.S1" = GA5.S1.red1, "GA5.S1.S2" = GA5.S1.red2,
                 "GA5.S2.S1" = GA5.S2.red1, "GA5.S2.S2" = GA5.S2.red2,
                 "GA5.Avg.S1" = GA5.Avg.red1, "GA5.Avg.S2" = GA5.Avg.red2,
                 "GA5.Both.S1" = GA5.Both.red1, "GA5.Both.S2" = GA5.Both.red2,
                 "GA5.BothMax.S1" = GA5.BothMax.red1, "GA5.BothMax.S2" = GA5.BothMax.red2)
save(excluded, file = "15FoldScen/Excluded3.RData")
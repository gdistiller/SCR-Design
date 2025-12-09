#September 2024
#generates both expected and simulated numbers using CLT S values
#excluded both Max = T
#1000 reps

library(secrdesign)
library(dplyr)
library(tidyr)
library(sf)

##########################################################
#Set up
#large area created and objects with and without a buffer are in LargeSCRObjs.RData
##########################################################

load("LargeSCRObjs.RData")

locs <- scr.objs$`SF traps polygon`
mask <- scr.objs$`SCR mask (buff)`

#Use CLT S values with 40 traps
sigmaF = 2000 ; sigmaM = 6500 ; sigma <- c(sigmaF, sigmaM)
DF = 0.00009 ; DM = 0.00006 ; D <- c(DF, DM)
lambda0F = 6 ; lambda0M = 3.25 ; lambda0 <- c(lambda0F,lambda0M) 
nT <- 40 ; nreps = 1000

Scen.F <- make.scenarios (trapsindex = 1, D = DF, lambda0 = lambda0F, sigma = sigmaF, noccasions = 1, detectfn = 14) 
Scen.M <- make.scenarios (trapsindex = 1, D = DM, lambda0 = lambda0M, sigma = sigmaM, noccasions = 1, detectfn = 14) 

#function to extract expected and data summaries
get.data.summary <- function(scenario, traps.list, mask.obj, nreps){
  expected <- scenarioSummary(scenario, trapset = traps.list[[1]], maskset = mask.obj)
  data <- run.scenarios(nreps, scenario, trapset = traps.list[[1]], maskset = mask, fit = FALSE, extract = identity, 
                        det.args = list(savepopn = TRUE), fit.args = list(detectfn = 'HHN'))
  results <- list("Expected" = expected, "Simulated data" = data)
  return(results)
}

###########################
#Grids
###########################
load("LGridDesignsCLTS40.RData")

###################
#Grid with F sigma
###################
grid.Fsig.F <- get.data.summary(Scen.F, L.grid.designs.CLT.S$`2*F Sigma`, mask, nreps)
grid.Fsig.M <- get.data.summary(Scen.M, L.grid.designs.CLT.S$`2*F Sigma`, mask, nreps)

grid.Fsig.F$Expected
grid.Fsig.M$Expected

#F: En = 8 ; Er = 46.3  ; Em = 28.3
#M: En = 14.9 ; Er = 192.2 ; Em = 178.4

summary(grid.Fsig.F$`Simulated data`)
summary(grid.Fsig.M$`Simulated data`)

#F: n = 8 ; r =  45.5 ; nmov = 13.1 
#M: n = 15 ; r = 192.9 ; nmov = 98

###################
#Grid with M sigma
###################
grid.Msig.F <- get.data.summary(Scen.F, L.grid.designs.CLT.S$`2*M Sigma`, mask, nreps)
grid.Msig.M <- get.data.summary(Scen.M, L.grid.designs.CLT.S$`2*M Sigma`, mask, nreps)

grid.Msig.F$Expected
grid.Msig.M$Expected

#F: En = 21.4 ; Er = 32.9 ; Em = 0.01
#M: En = 51.7 ; Er = 155.4 ; Em = 95.9 

summary(grid.Msig.F$`Simulated data`)
summary(grid.Msig.M$`Simulated data`)

#F: n = 21.3 ; r = 32.1  ; nmov = 0.01  
#M: n = 51.7 ; r = 154.9 ; nmov = 60.3  

###################
#Grid with Avg sigma
###################
grid.Avgsig.F <- get.data.summary(Scen.F, L.grid.designs.CLT.S$`2*Avg Sigma`, mask, nreps)
grid.Avgsig.M <- get.data.summary(Scen.M, L.grid.designs.CLT.S$`2*Avg Sigma`, mask, nreps)

grid.Avgsig.F$Expected
grid.Avgsig.M$Expected

#F: En = 20 ; Er = 34.3 ; Em = 1.9 
#M: En = 30.9 ; Er = 176.2 ; Em = 140.8 

summary(grid.Avgsig.F$`Simulated data`)
summary(grid.Avgsig.M$`Simulated data`)

#F: n = 19.6 ; r = 33.6 ; nmov = 1.5 
#M: n = 31 ; r = 176.6 ; nmov = 81.7   

###########################
#GA Designs
###########################

load("Cluster/ProposedDesigns/LargeArea/Results/GADesigns.RData")

###########################
#GA Crit 4
###########################

#######################
#GA crit 4 with F sigma
#######################
Crit4.Fsig.F <- get.data.summary(Scen.F, GA.designs$`F4 Sigma`, mask, nreps)
Crit4.Fsig.M <- get.data.summary(Scen.M, GA.designs$`F4 Sigma`, mask, nreps)

Crit4.Fsig.F$Expected
Crit4.Fsig.M$Expected

#F: En =  ; Er =  ; Em = 
#M: En =  ; Er =  ; Em = 

summary(Crit4.Fsig.F$`Simulated data`)
summary(Crit4.Fsig.M$`Simulated data`)

#F: n =  ; r =  ; nmov =  
#M: n =  ; r =  ; nmov =  

###############################
#GA crit 4 with M sigma
###############################
Crit4.Msig.F <- get.data.summary(Scen.F, GA.designs$`M4 Sigma`, mask, nreps)
Crit4.Msig.M <- get.data.summary(Scen.M, GA.designs$`M4 Sigma`, mask, nreps)

Crit4.Msig.F$Expected
Crit4.Msig.M$Expected

#F: En =  ; Er =  ; Em = 
#M: En =  ; Er =  ; Em = 

summary(Crit4.Msig.F$`Simulated data`)
summary(Crit4.Msig.M$`Simulated data`)

#F: n =  ; r =  ; nmov = 
#M: n =  ; r =   ; nmov = 

###############################
#GA crit 4 with Avg sigma
###############################
Crit4.Avgsig.F <- get.data.summary(Scen.F, GA.designs$`Avg4 Sigma`, mask, nreps)
Crit4.Avgsig.M <- get.data.summary(Scen.M, GA.designs$`Avg4 Sigma`, mask, nreps)

Crit4.Avgsig.F$Expected
Crit4.Avgsig.M$Expected

#F: En = ; Er =  ; Em =  
#M: En =  ; Er =  ; Em =

summary(Crit4.Avgsig.F$`Simulated data`)
summary(Crit4.Avgsig.M$`Simulated data`)

#F: n = ; r =  ; nmov =  
#M: n = ; r =  ; nmov =  

###############################
#GA crit 4 with both
###############################
Crit4.both.F <- get.data.summary(Scen.F, GA.designs$`Both4`, mask, nreps)
Crit4.both.M <- get.data.summary(Scen.M, GA.designs$`Both4`, mask, nreps)

Crit4.both.F$Expected
Crit4.both.M$Expected

#F: En =  ; Er =  ; Em =  
#M: En =  ; Er =  ; Em = 

summary(Crit4.both.F$`Simulated data`)
summary(Crit4.both.M$`Simulated data`)

#F: n =  ; r =  ; nmov =  
#M: n =  ; r =  ; nmov = 

###########################
#GA Crit 5
###########################

#######################
#GA crit 5 with F sigma
#######################
Crit5.Fsig.F <- get.data.summary(Scen.F, GA.designs$`F5 Sigma`, mask, nreps)
Crit5.Fsig.M <- get.data.summary(Scen.M, GA.designs$`F5 Sigma`, mask, nreps)

Crit5.Fsig.F$Expected
Crit5.Fsig.M$Expected

#F: En =  ; Er =  ; Em =
#M: En =  ; Er =  ; Em = 

summary(Crit5.Fsig.F$`Simulated data`)
summary(Crit5.Fsig.M$`Simulated data`)

#F: n =  ; r =  ; nmov = 
#M: n =  ; r =  ; nmov = 

###############################
#GA crit 5 with M sigma
###############################
Crit5.Msig.F <- get.data.summary(Scen.F, GA.designs$`M5 Sigma`, mask, nreps)
Crit5.Msig.M <- get.data.summary(Scen.M, GA.designs$`M5 Sigma`, mask, nreps)

Crit5.Msig.F$Expected
Crit5.Msig.M$Expected

#F: En =  ; Er =  ; Em = 
#M: En =  ; Er =  ; Em = 

summary(Crit5.Msig.F$`Simulated data`)
summary(Crit5.Msig.M$`Simulated data`)

#F: n =  ; r = ; nmov = 
#M: n =  ; r =  ; nmov = 

###############################
#GA crit 5 with Avg sigma
###############################
Crit5.Avgsig.F <- get.data.summary(Scen.F, GA.designs$`Avg5 Sigma`, mask, nreps)
Crit5.Avgsig.M <- get.data.summary(Scen.M, GA.designs$`Avg5 Sigma`, mask, nreps)

Crit5.Avgsig.F$Expected
Crit5.Avgsig.M$Expected

#F: En =  ; Er =  ; Em = 
#M: En =  ; Er =  ; Em = 

summary(Crit5.Avgsig.F$`Simulated data`)
summary(Crit5.Avgsig.M$`Simulated data`)

#F: n =  ; r =  ; nmov =  
#M: n =  ; r =   ; nmov =  

###############################
#GA crit 5 with both
###############################
Crit5.both.F <- get.data.summary(Scen.F, GA.designs$`Both5`, mask, nreps)
Crit5.both.M <- get.data.summary(Scen.M, GA.designs$`Both5`, mask, nreps)

Crit5.both.F$Expected
Crit5.both.M$Expected

#F: En =  ; Er =  ; Em =  
#M: En =  ; Er =  ; Em = 

summary(Crit5.both.F$`Simulated data`)
summary(Crit5.both.M$`Simulated data`)

#F: n =  ; r =  ; nmov = 
#M: n =  ; r =  ; nmov =  

##############################################################
#dataframe for plotting

Design <- rep(c("Grid.FSigma", "Grid.MSigma","Grid.AvgSigma", 
            "GA4.F","GA4.M", "GA4.Avg","GA4.both",
            "GA5.F","GA5.M", "GA5.Avg","GA5.both"),each = 2)
DesignCat <- c(rep("Grid",6),rep("GA4",8),rep("GA5",8))
Strata <- rep(c("F","M"),11)

Individuals <- c(summary(grid.Fsig.F$`Simulated data`)[[2]][[1]][1,2], summary(grid.Fsig.M$`Simulated data`)[[2]][[1]][1,2], summary(grid.Msig.F$`Simulated data`)[[2]][[1]][1,2], summary(grid.Msig.M$`Simulated data`)[[2]][[1]][1,2], summary(grid.Avgsig.F$`Simulated data`)[[2]][[1]][1,2], summary(grid.Avgsig.M$`Simulated data`)[[2]][[1]][1,2],
                 summary(Crit4.Fsig.F$`Simulated data`)[[2]][[1]][1,2], summary(Crit4.Fsig.M$`Simulated data`)[[2]][[1]][1,2], summary(Crit4.Msig.F$`Simulated data`)[[2]][[1]][1,2], summary(Crit4.Msig.M$`Simulated data`)[[2]][[1]][1,2], summary(Crit4.Avgsig.F$`Simulated data`)[[2]][[1]][1,2], summary(Crit4.Avgsig.M$`Simulated data`)[[2]][[1]][1,2], summary(Crit4.both.F$`Simulated data`)[[2]][[1]][1,2], summary(Crit4.both.M$`Simulated data`)[[2]][[1]][1,2],
                 summary(Crit5.Fsig.F$`Simulated data`)[[2]][[1]][1,2], summary(Crit5.Fsig.M$`Simulated data`)[[2]][[1]][1,2], summary(Crit5.Msig.F$`Simulated data`)[[2]][[1]][1,2], summary(Crit5.Msig.M$`Simulated data`)[[2]][[1]][1,2], summary(Crit5.Avgsig.F$`Simulated data`)[[2]][[1]][1,2], summary(Crit5.Avgsig.M$`Simulated data`)[[2]][[1]][1,2], summary(Crit5.both.F$`Simulated data`)[[2]][[1]][1,2], summary(Crit5.both.M$`Simulated data`)[[2]][[1]][1,2])
Recaptures <- c(summary(grid.Fsig.F$`Simulated data`)[[2]][[1]][2,2], summary(grid.Fsig.M$`Simulated data`)[[2]][[1]][2,2], summary(grid.Msig.F$`Simulated data`)[[2]][[1]][2,2], summary(grid.Msig.M$`Simulated data`)[[2]][[1]][2,2], summary(grid.Avgsig.F$`Simulated data`)[[2]][[1]][2,2], summary(grid.Avgsig.M$`Simulated data`)[[2]][[1]][2,2], 
                summary(Crit4.Fsig.F$`Simulated data`)[[2]][[1]][2,2], summary(Crit4.Fsig.M$`Simulated data`)[[2]][[1]][2,2], summary(Crit4.Msig.F$`Simulated data`)[[2]][[1]][2,2], summary(Crit4.Msig.M$`Simulated data`)[[2]][[1]][2,2], summary(Crit4.Avgsig.F$`Simulated data`)[[2]][[1]][2,2], summary(Crit4.Avgsig.M$`Simulated data`)[[2]][[1]][2,2], summary(Crit4.both.F$`Simulated data`)[[2]][[1]][2,2], summary(Crit4.both.M$`Simulated data`)[[2]][[1]][2,2],
                summary(Crit5.Fsig.F$`Simulated data`)[[2]][[1]][2,2], summary(Crit5.Fsig.M$`Simulated data`)[[2]][[1]][2,2], summary(Crit5.Msig.F$`Simulated data`)[[2]][[1]][2,2], summary(Crit5.Msig.M$`Simulated data`)[[2]][[1]][2,2], summary(Crit5.Avgsig.F$`Simulated data`)[[2]][[1]][2,2], summary(Crit5.Avgsig.M$`Simulated data`)[[2]][[1]][2,2], summary(Crit5.both.F$`Simulated data`)[[2]][[1]][2,2], summary(Crit5.both.M$`Simulated data`)[[2]][[1]][2,2])
SpatialRCS <- c(summary(grid.Fsig.F$`Simulated data`)[[2]][[1]][3,2], summary(grid.Fsig.M$`Simulated data`)[[2]][[1]][3,2], summary(grid.Msig.F$`Simulated data`)[[2]][[1]][3,2], summary(grid.Msig.M$`Simulated data`)[[2]][[1]][3,2], summary(grid.Avgsig.F$`Simulated data`)[[2]][[1]][3,2], summary(grid.Avgsig.M$`Simulated data`)[[2]][[1]][3,2], 
                summary(Crit4.Fsig.F$`Simulated data`)[[2]][[1]][3,2], summary(Crit4.Fsig.M$`Simulated data`)[[2]][[1]][3,2], summary(Crit4.Msig.F$`Simulated data`)[[2]][[1]][3,2], summary(Crit4.Msig.M$`Simulated data`)[[2]][[1]][3,2], summary(Crit4.Avgsig.F$`Simulated data`)[[2]][[1]][3,2], summary(Crit4.Avgsig.M$`Simulated data`)[[2]][[1]][3,2], summary(Crit4.both.F$`Simulated data`)[[2]][[1]][3,2], summary(Crit4.both.M$`Simulated data`)[[2]][[1]][3,2],
                summary(Crit5.Fsig.F$`Simulated data`)[[2]][[1]][3,2], summary(Crit5.Fsig.M$`Simulated data`)[[2]][[1]][3,2], summary(Crit5.Msig.F$`Simulated data`)[[2]][[1]][3,2], summary(Crit5.Msig.M$`Simulated data`)[[2]][[1]][3,2], summary(Crit5.Avgsig.F$`Simulated data`)[[2]][[1]][3,2], summary(Crit5.Avgsig.M$`Simulated data`)[[2]][[1]][3,2], summary(Crit5.both.F$`Simulated data`)[[2]][[1]][3,2], summary(Crit5.both.M$`Simulated data`)[[2]][[1]][3,2])

En <- c(grid.Fsig.F$Expected$En, grid.Fsig.M$Expected$En, grid.Msig.F$Expected$En, grid.Msig.M$Expected$En, grid.Avgsig.F$Expected$En, grid.Avgsig.M$Expected$En, 
        Crit4.Fsig.F$Expected$En, Crit4.Fsig.M$Expected$En,  Crit4.Msig.F$Expected$En, Crit4.Msig.M$Expected$En, Crit4.Avgsig.F$Expected$En, Crit4.Avgsig.M$Expected$En, Crit4.both.F$Expected$En, Crit4.both.M$Expected$En, 
        Crit5.Fsig.F$Expected$En, Crit5.Fsig.M$Expected$En,  Crit5.Msig.F$Expected$En, Crit5.Msig.M$Expected$En, Crit5.Avgsig.F$Expected$En, Crit5.Avgsig.M$Expected$En, Crit5.both.F$Expected$En, Crit5.both.M$Expected$En)
Er <- c(grid.Fsig.F$Expected$Er, grid.Fsig.M$Expected$Er, grid.Msig.F$Expected$Er, grid.Msig.M$Expected$Er, grid.Avgsig.F$Expected$Er, grid.Avgsig.M$Expected$Er, 
        Crit4.Fsig.F$Expected$Er, Crit4.Fsig.M$Expected$Er, Crit4.Msig.F$Expected$Er, Crit4.Msig.M$Expected$Er, Crit4.Avgsig.F$Expected$Er, Crit4.Avgsig.M$Expected$Er, Crit4.both.F$Expected$Er, Crit4.both.M$Expected$Er, 
        Crit5.Fsig.F$Expected$Er, Crit5.Fsig.M$Expected$Er,  Crit5.Msig.F$Expected$Er, Crit5.Msig.M$Expected$Er, Crit5.Avgsig.F$Expected$Er, Crit5.Avgsig.M$Expected$Er, Crit5.both.F$Expected$Er, Crit5.both.M$Expected$Er)
Em <- c(grid.Fsig.F$Expected$Em, grid.Fsig.M$Expected$Em, grid.Msig.F$Expected$Em, grid.Msig.M$Expected$Em, grid.Avgsig.F$Expected$Em, grid.Avgsig.M$Expected$Em, 
        Crit4.Fsig.F$Expected$Em, Crit4.Fsig.M$Expected$Em, Crit4.Msig.F$Expected$Em, Crit4.Msig.M$Expected$Em, Crit4.Avgsig.F$Expected$Em, Crit4.Avgsig.M$Expected$Em, Crit4.both.F$Expected$Em, Crit4.both.M$Expected$Em, 
        Crit5.Fsig.F$Expected$Em, Crit5.Fsig.M$Expected$Em, Crit5.Msig.F$Expected$Em, Crit5.Msig.M$Expected$Em, Crit5.Avgsig.F$Expected$Em, Crit5.Avgsig.M$Expected$Em, Crit5.both.F$Expected$Em, Crit5.both.M$Expected$Em)

df <- data.frame("Type" = as.factor(DesignCat), Design = as.factor(Design), Strata = as.factor(Strata), Individuals, Recaptures, SpatialRCS, En, Er, Em)

#################################################

#plot simulated data summary number
pdf("SimulatedNumbers.pdf", height = 12, width = 8)
par(mfrow=c(3,1))
#get axis labels
Designlabs <- levels(df$Design)

#Individuals
plot(df$Individuals ~ as.numeric(df$Design), type='n', xlab = "", ylab = "Number of individuals (simulated)", xaxt = 'n')
axis(side = 1, at = 1:11, labels = Designlabs, las = 2,  cex.axis = 0.7)
points(df$Individuals[df$Strata=="F"] ~ as.numeric(df$Design[df$Strata=="F"]), pch = c(15,16,17)[as.numeric(df$Type[df$Strata=="F"])], col = c("red", "green","blue")[as.numeric(df$Type[df$Strata=="F"])])
points(df$Individuals[df$Strata=="M"] ~ as.numeric(df$Design[df$Strata=="M"]), pch = c(0,1,2)[as.numeric(df$Type[df$Strata=="M"])], col = c("red", "green","blue")[as.numeric(df$Type[df$Strata=="M"])])
legend("left", legend=c("Space filling", "Clustered","Grid", "F", "M"), pch = c(15,16,17,15,0),col = c("red", "green","blue","black","black"), cex = 0.85 )

#Recaps
plot(df$Recaptures ~ as.numeric(df$Design), type='n', xlab = "", ylab = "Number of recaptures (simulated)", xaxt = 'n')
axis(side = 1, at = 1:11, labels = Designlabs, las = 2,  cex.axis = 0.7)
points(df$Recaptures[df$Strata=="F"] ~ as.numeric(df$Design[df$Strata=="F"]), pch = c(15,16,17)[as.numeric(df$Type[df$Strata=="F"])], col = c("red", "green","blue")[as.numeric(df$Type[df$Strata=="F"])])
points(df$Recaptures[df$Strata=="M"] ~ as.numeric(df$Design[df$Strata=="M"]), pch = c(0,1,2)[as.numeric(df$Type[df$Strata=="M"])], col = c("red", "green","blue")[as.numeric(df$Type[df$Strata=="M"])])
legend("topleft", legend=c("Space filling", "Clustered","Grid", "F", "M"), pch = c(15,16,17,15,0),col = c("red", "green","blue","black","black"), cex = 0.85 )

#SRCs
plot(df$SpatialRCS ~ as.numeric(df$Design), type='n', xlab = "", ylab = "Spatial recaptures (simulated)", xaxt = 'n', ylim=c(0,160))
axis(side = 1, at = 1:11, labels = Designlabs, las = 2,  cex.axis = 0.7)
points(df$SpatialRCS[df$Strata=="F"] ~ as.numeric(df$Design[df$Strata=="F"]), pch = c(15,16,17)[as.numeric(df$Type[df$Strata=="F"])], col = c("red", "green","blue")[as.numeric(df$Type[df$Strata=="F"])])
points(df$SpatialRCS[df$Strata=="M"] ~ as.numeric(df$Design[df$Strata=="M"]), pch = c(0,1,2)[as.numeric(df$Type[df$Strata=="M"])], col = c("red", "green","blue")[as.numeric(df$Type[df$Strata=="M"])])
legend("topleft", legend=c("Space filling", "Clustered","Grid", "F", "M"), pch = c(15,16,17,15,0),col = c("red", "green","blue","black","black"), cex = 0.85 )

dev.off()

#plot expected data summary number
pdf("ExpectedNumbers.pdf", height = 12, width = 8)
par(mfrow=c(3,1))

#Individuals
plot(df$En ~ as.numeric(df$Design), type='n', xlab = "", ylab = "E(n)", xaxt = 'n')
axis(side = 1, at = 1:11, labels = Designlabs, las = 2,  cex.axis = 0.7)
points(df$En[df$Strata=="F"] ~ as.numeric(df$Design[df$Strata=="F"]), pch = c(15,16,17)[as.numeric(df$Type[df$Strata=="F"])], col = c("red", "green","blue")[as.numeric(df$Type[df$Strata=="F"])])
points(df$En[df$Strata=="M"] ~ as.numeric(df$Design[df$Strata=="M"]), pch = c(0,1,2)[as.numeric(df$Type[df$Strata=="M"])], col = c("red", "green","blue")[as.numeric(df$Type[df$Strata=="M"])])
legend("left", legend=c("Space filling", "Clustered","Grid", "F", "M"), pch = c(15,16,17,15,0),col = c("red", "green","blue","black","black"), cex = 0.85 )

#Recaps
plot(df$Er ~ as.numeric(df$Design), type='n', xlab = "", ylab = "E(r)", xaxt = 'n')
axis(side = 1, at = 1:11, labels = Designlabs, las = 2,  cex.axis = 0.7)
points(df$Er[df$Strata=="F"] ~ as.numeric(df$Design[df$Strata=="F"]), pch = c(15,16,17)[as.numeric(df$Type[df$Strata=="F"])], col = c("red", "green","blue")[as.numeric(df$Type[df$Strata=="F"])])
points(df$Er[df$Strata=="M"] ~ as.numeric(df$Design[df$Strata=="M"]), pch = c(0,1,2)[as.numeric(df$Type[df$Strata=="M"])], col = c("red", "green","blue")[as.numeric(df$Type[df$Strata=="M"])])
legend("topleft", legend=c("Space filling", "Clustered","Grid", "F", "M"), pch = c(15,16,17,15,0),col = c("red", "green","blue","black","black"), cex = 0.85 )

#SRCs
plot(df$Em ~ as.numeric(df$Design), type='n', xlab = "", ylab = "E(m)", xaxt = 'n', , ylim=c(0,160))
axis(side = 1, at = 1:11, labels = Designlabs, las = 2,  cex.axis = 0.7)
points(df$Em[df$Strata=="F"] ~ as.numeric(df$Design[df$Strata=="F"]), pch = c(15,16,17)[as.numeric(df$Type[df$Strata=="F"])], col = c("red", "green","blue")[as.numeric(df$Type[df$Strata=="F"])])
points(df$Em[df$Strata=="M"] ~ as.numeric(df$Design[df$Strata=="M"]), pch = c(0,1,2)[as.numeric(df$Type[df$Strata=="M"])], col = c("red", "green","blue")[as.numeric(df$Type[df$Strata=="M"])])
legend("topleft", legend=c("Space filling", "Clustered","Grid", "F", "M"), pch = c(15,16,17,15,0),col = c("red", "green","blue","black","black"), cex = 0.85 )

dev.off()

#plot together
pdf("DataSumms.pdf", height = 12, width = 8, pointsize = 11)
par(mfrow=c(3,2))

############################################################################
#for a smaller area

#create mask and potential traplocs
mask.r <- rast(ncol = 60, nrow = 60, xmin = 0, xmax = 75000, ymin = 0, ymax = 75000)  #ncol/row determine the resolution
traplocs.r <- rast(ncol = 30, nrow = 30, xmin = 20000, xmax = 55000, ymin = 20000, ymax = 55000) #res of 1.5K

values(mask.r) <- 1:ncell(mask.r) #numbers each cell, seems needed
values(traplocs.r) <- 1:ncell(traplocs.r) 

#create polyon
msk.p1 <- rbind(c(0, 0), c(75000, 0), c(75000, 75000), c(0, 75000), c(0,0))
mask.v <- vect(msk.p1, type = "polygons")

traplocs.p1 <- rbind(c(20000, 20000), c(55000, 20000), c(55000, 55000), c(20000, 55000), c(20000,20000))
traplocs.v <- vect(traplocs.p1, type = "polygons")

#create mask and traps obj
#get coordinates
msk.coords <- crds(mask.r, df = T)
traploc.coords <- crds(traplocs.r, df=T)

mask <- read.mask(data = msk.coords, spacing = 500)
trap.locs <- read.traps(data = traploc.coords, detector = "count", spacing = 1500)

# 40 traps
D4 = 0.008 ; L0 = 3 ; sig = 3000

nT <- 40

######
#GA
######
GA.F4 <- Proposed.traps(poly = NULL, alltraps = trap.locs, mask.buff = mask, D = D4, sigma = sig, 
                        lambda0 = L0, sigma.buff = sig, grid.spacing = NULL, 
                        criterion = 4, n.reps = 1, grid = FALSE, nT = nT, GA.ngen = 50)
plot(mask.v)
plot(GA.F4[[1]], add = T)

#Enrm numbers
scenF <- make.scenarios (trapsindex = 1, detectfn = 'HHN', D = D4, lambda0 = L0, 
                         sigma = sig, noccasions = 1)

scenarioSummary(scenF, traps = GA.F4[[1]], mask = mask)




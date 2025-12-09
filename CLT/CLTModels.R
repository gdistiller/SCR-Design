#Sept 2024, fitting models to CLT Boland Survey
#redoing it after getting revised Mask from Anita
#setup capture history objs

library(secr)
library(secrdesign)
load("CLT/SepAreas.RData")

#using CLTdat that includes the 4 captures at traps that were corrected to be R5 (L2-1.Ttol.R6 / R7 / R8)
#also includes the captures with NA sex
#there were spaces in the trap names, replaced with fullstops and named BolandCHs.csv
#CLTTrapFile txt files have had spaces removed too, two stations needed upper case (L2-1 and L22-2)

CH <- read.capthist("CLT/BolandCHs.csv", "CLT/CLTTrapFile.txt", detector = "count", fmt = "trapID", covnames = "Sex")

summary(CH)
par(mar = c(1,1,3,1)) # reduce margins
plot (CH, tracks = TRUE)

m <- unlist(moves(CH))
par(mar = c(3.2,4,1,1), mgp = c(2.1,0.6,0)) # reduce margins
hist(m, xlab = "Movement m", main = "")
RPSV(CH, CC = TRUE)

# as 2 seperate sessions, 1 is the Northern section and 2 the Southern
NorthernTraps <- read.traps("CLT/CLTTrapFile1.txt", detector = "count")
SouthernTraps <- read.traps("CLT/CLTTrapFile2.txt", detector = "count")

caps.session1 <- read.csv("CLT/CLTNorthernSession.csv")
caps.session2 <- read.csv("CLT/CLTSouthernSession.csv") #had to renumber occs to start at 1

CH.Northern <- make.capthist(caps.session1, NorthernTraps)
CH.Southern <- make.capthist(caps.session2, SouthernTraps)

#MS object (with occs)
trapsession.list <- list(Northern = traps(CH.Northern), Southern = traps(CH.Southern))

#read in data split into two sessions
allcaps <- read.csv("CLT/BolandCHs.csv")
allcaps[allcaps$Occasion<131,1] <- "Northern"
allcaps[allcaps$Occasion>130,1] <- "Southern"
#renumber occasions in Southern session
allcaps[allcaps$X.Session=="Southern",3] <- allcaps[allcaps$X.Session=="Southern",3] - 130  
MS <- make.capthist(allcaps, traps = trapsession.list)

summary(MS)

####################################################################################
#reducing to 1 occasions
#alltogether
CH.red <- reduce(CH, by = 'all')
MS.red <- reduce(MS, by = "all")

#after using reduce no need to manually sort out usage any longer
# sum usage info
usage.occ1 <- usage(session1)
usage.occ2 <- usage(session2) 

usage.sum1 <- apply(usage.occ1,1,sum)
usage.sum2 <- apply(usage.occ2,1,sum)

usage(session1) <- NULL ; usage(session1) <- matrix(usage.sum1,ncol = 1)
usage(session2) <- NULL ; usage(session2) <- matrix(usage.sum2,ncol = 1)

save(CH, CH.red, CH.Northern ,CH.Southern, MS, MS.red, file = "CLTCapthists.RData")

####################################################################################
#load results from cluster
#1st for multi occs
load("CLT/Cluster/CLTResults2.RData")
load("CLT/Cluster/CLTResultsRed.RData")

#all data as 1 session
fit.null.both
fit.sex.both

#N fits
fit.null.N
fit.sex.N

#S fits
fit.null.S
fit.sex.S

#using MS obj
#note need to specify pmix ~ session to est pmix per area (b)
fit.MS.null
fit.MS.sexa
fit.MS.sexb

###################################################################
#refitting to one occasion
#both sessions together
fit.null.both.red
fit.sex.both.red

#seperately using MS obj
#note need to specify pmix ~ session to est pmix per area (b)
fit.MS.null.red
fit.MS.sexa.red
fit.MS.sexb.red

AIC(fit.MS.null.red,fit.MS.sexa.red,fit.MS.sexb.red)

############################################
#get Enrm numbers for CLT areas
load("CLT/CLTCapthists.RData")
plot(CLT_Both_msk)
plot(traps(CH), add = T)

#together
#data summaries, note that grouped scenarios are excluded from scenarioSummary so will do one for each sex
CLTscen.F <- make.scenarios (trapsindex = 1, D = 0.000084, lambda0 = 4.16, sigma = 2363, noccasions = 1, detectfn = 14) 
CLTscen.M <- make.scenarios (trapsindex = 1, D = 0.000081, lambda0 = 3.12, sigma = 4854, noccasions = 1, detectfn = 14) 

CLTTraps <- traps(CH)
usage(CLTTraps) <- NULL

scenarioSummary(CLTscen.F , trapset = CLTTraps, maskset = CLT_Both_msk)
#En = 19.3 ; Er = 82
scenarioSummary(CLTscen.M , trapset = CLTTraps, maskset = CLT_Both_msk)
#En = 22.8 ; Er = 249.8

dat <- read.csv("CLT/CLTdat.csv", sep = ";")
F <- as.matrix(table(dat$ID[dat$Sex=="F"], dat$Detector[dat$Sex=="F"]))
M <- as.matrix(table(dat$ID[dat$Sex=="M"], dat$Detector[dat$Sex=="M"]))

sum(F)
sum(M)

#North
#data summaries, note that grouped scenarios are excluded from scenarioSummary so will do one for each sex
CLTscenN.F <- make.scenarios (trapsindex = 1, D = 0.00008, lambda0 = 3.9, sigma = 2578, noccasions = 1, detectfn = 14) 
CLTscenN.M <- make.scenarios (trapsindex = 1, D = 0.0001, lambda0 = 4.225, sigma = 3727, noccasions = 1, detectfn = 14) 

CLTTrapsN <- traps(CH.Northern)
usage(CLTTrapsN) <- NULL

scenarioSummary(CLTscenN.F , trapset = CLTTrapsN, maskset = CLT_Both_msk)
#En = 9.5 ; Er = 43.4
scenarioSummary(CLTscenN.M , trapset = CLTTrapsN, maskset = CLT_Both_msk)
#En = 14.4 ; Er = 129.4
#actual 9 * F, 13 * M, 147 detections so quite close

#South
#data summaries, note that grouped scenarios are excluded from scenarioSummary so will do one for each sex
CLTscenS.F <- make.scenarios (trapsindex = 1, D = 0.00009, lambda0 = 6, sigma = 1962, noccasions = 1, detectfn = 14) 
CLTscenS.M <- make.scenarios (trapsindex = 1, D = 0.000056, lambda0 = 3.25, sigma = 6483, noccasions = 1, detectfn = 14) 

CLTTrapsS <- traps(CH.Southern)
usage(CLTTrapsS) <- NULL

scenarioSummary(CLTscenS.F , trapset = CLTTrapsS, maskset = CLT_Both_msk)
#En = 11.2 ; Er = 45.1
scenarioSummary(CLTscenS.M , trapset = CLTTrapsS, maskset = CLT_Both_msk)
#En = 11.7 ; Er = 142.8 

#actual 10 * F, 11 * M, 165 detections so quite close
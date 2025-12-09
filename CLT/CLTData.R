#Sept 2023, data from CLT Boland Survey
#note that the full data have 312 detections including 4 with NA for sex for 3 individuals
dat <- read.csv("CLT/CLT Boland cam survey 2010 - 2011 Capture history files.csv")

F <- as.matrix(table(dat$ID[dat$Sex=="F"], dat$Detector[dat$Sex=="F"]))
M <- as.matrix(table(dat$ID[dat$Sex=="M"], dat$Detector[dat$Sex=="M"]))

table(dat$ID[dat$Sex=="F"], dat$Detector[dat$Sex=="F"])
table(dat$ID[dat$Sex=="M"], dat$Detector[dat$Sex=="M"])

n.caps.F <- apply(F,1,sum)
detectors.F <- apply(F,1,function(x) which(x>0))
n.detectors.F <- lapply(detectors.F,function(x) length(unique(x)))
n.caps.M <- apply(M,1,sum)
detectors.M <- apply(M,1,function(x) which(x>0))
n.detectors.M <- lapply(detectors.M,function(x) length(unique(x)))

mean(as.numeric(n.detectors.F[1:13]))
mean(as.numeric(n.detectors.M[1:21]))

#F: 3 caught at 1, 6 at 2, 3 at 3, 1 at 4 (13 at 78 caps)
#M: 6 caught at 1, 4 at 2, 0 at 3, 2 at 4, 3 at 5, 1 at 6, 1 at 7, 3 at 9
# 1 at 12 (21 at 222)

#######################################################
#setup capture history objs
library(secr)
load("CLTAreas.RData")

#using CLTdat that includes the 4 captures at traps that were corrected to be R5 (L2-1.Ttol.R6 / R7 / R8)
#also includes the captures with NA sex
#there were spaces in the trap names, replaced with fullstops and named BolandCHs.csv
#CLTTrapFile txt files have had spaces removed too, two stations needed upper case (L2-1 and L22-2)

CH <- read.capthist("BolandCHs.csv", "CLTTrapFile.txt", detector = "count", fmt = "trapID", covnames = "Sex")

summary(CH)
par(mar = c(1,1,3,1)) # reduce margins
plot (CH, tracks = TRUE)

m <- unlist(moves(CH))
par(mar = c(3.2,4,1,1), mgp = c(2.1,0.6,0)) # reduce margins
hist(m, xlab = "Movement m", main = "")
RPSV(CH, CC = TRUE)

# as 2 seperate sessions, 1 is the Northern section and 2 the Southern
NorthernTraps <- read.traps("CLTTrapFile1.txt", detector = "count")
SouthernTraps <- read.traps("CLTTrapFile2.txt", detector = "count")

caps.session1 <- read.csv("CLTNorthernSession.csv")
caps.session2 <- read.csv("CLTSouthernSession.csv") #had to renumber occs to start at 1

CH.Northern <- make.capthist(caps.session1, NorthernTraps)
CH.Southern <- make.capthist(caps.session2, SouthernTraps)

#MS object (with occs)
trapsession.list <- list(Northern = traps(CH.Northern), Southern = traps(CH.Southern))

#read in data split into two sessions
allcaps <- read.csv("BolandCHs.csv")
allcaps[allcaps$Occasion<131,1] <- "Northern"
allcaps[allcaps$Occasion>130,1] <- "Southern"
#renumber occasions in Southern session
allcaps[allcaps$X.Session=="Southern",3] <- allcaps[allcaps$X.Session=="Southern",3] - 130  
MS <- make.capthist(allcaps, traps = trapsession.list)

summary(MS)

#fit models on cluster

#all together
fit.all.null <- secr.fit(CH, mask = CLTsimp.nobuff, list(D ~ 1, lambda0 ~ 1, sigma ~ 1), detectfn = 'HHN')
fit.all.sex <- secr.fit(CH, mask = CLTsimp.nobuff, hcov = "Sex", list(D ~ 1, lambda0 ~ h2, sigma ~ h2), detectfn = 'HHN')

#by section
fit.N.null <- secr.fit(CH.Northern, mask = CLTsimp.nobuff, list(D ~ 1, lambda0 ~ 1, sigma ~ 1), detectfn = 'HHN')
fit.N.sex <- secr.fit(CH.Northern, mask = CLTsimp.nobuff, hcov = "Sex", list(D ~ 1, lambda0 ~ h2, sigma ~ h2),  detectfn = 'HHN')

fit.S.null <- secr.fit(CH.Southern, mask = CLTsimp.nobuff, list(D ~ 1, lambda0 ~ 1, sigma ~ 1), detectfn = 'HHN')
fit.S.sex <- secr.fit(CH.Southern, mask = CLTsimp.nobuff, hcov = "Sex", list(D ~ 1, lambda0 ~ h2, sigma ~ h2), detectfn = 'HHN')

#using MS obj
fit.MS.NS.null <- secr.fit(MS, mask = CLTsimp.nobuff, hcov = "Sex", model = list(D ~ session, lambda0 ~ session, sigma ~ session, pmix ~ session), detectfn = 'HHN', trace = FALSE, ncores = 1)
fit.MS.NS.sexa <- secr.fit(MS, mask = CLTsimp.nobuff, hcov = "Sex", model = list(D ~ session, lambda0 ~ session*h2, sigma ~ session*h2), detectfn = 'HHN', trace = FALSE, ncores = 1)
fit.MS.NS.sexb <- secr.fit(MS, mask = CLTsimp.nobuff, hcov = "Sex", model = list(D ~ session, lambda0 ~ session*h2, sigma ~ session*h2, pmix ~ session*h2), detectfn = 'HHN', trace = FALSE, ncores = 1)

####################################################################################
#refit after reducing to 1 occasions
#alltogether
CH.red <- reduce(CH, by = 'all')
MS.red <- reduce(MS, by = "all")

load("Cluster/CLTResultsRed.RData")

fit.all.null.1oc <- secr.fit(CH.red, mask = CLTsimp.nobuff, list(D ~ 1, lambda0 ~ 1, sigma ~ 1), detectfn = 'HHN')
fit.all.sex.1oc <- secr.fit(CH.red, mask = CLTsimp.nobuff, hcov = "Sex", list(D ~ 1, lambda0 ~ h2, sigma ~ h2), detectfn = 'HHN')

fit.all.null.1oc
fit.all.sex.1oc

#null: D = 0.00016 ; L0 = 0.02 (2.6) ; Sig = 4326
#sex: Df = 0.00007 ; Dm = 0.00009 ; L0f = 0.033 ; Lom = 0.024 ; Sigf = 23786 ; sigm = 4735

#using MS obj
fit.MS.NS.null <- secr.fit(MS.red, mask = CLTsimp.nobuff, hcov = "Sex", model = list(D ~ session, lambda0 ~ session, sigma ~ session), detectfn = 'HHN', trace = FALSE, ncores = 1)
fit.MS.NS.sexa <- secr.fit(MS.red, mask = CLTsimp.nobuff, hcov = "Sex", model = list(D ~ session, lambda0 ~ session*h2, sigma ~ session*h2), detectfn = 'HHN', trace = FALSE, ncores = 1)
fit.MS.NS.sexb <- secr.fit(MS.red, mask = CLTsimp.nobuff, hcov = "Sex", model = list(D ~ session, lambda0 ~ session*h2, sigma ~ session*h2, pmix ~ session), detectfn = 'HHN', trace = FALSE, ncores = 1)

#after using reduce no need to manually sort out usage any longer
# sum usage info
usage.occ1 <- usage(session1)
usage.occ2 <- usage(session2) 

usage.sum1 <- apply(usage.occ1,1,sum)
usage.sum2 <- apply(usage.occ2,1,sum)

usage(session1) <- NULL ; usage(session1) <- matrix(usage.sum1,ncol = 1)
usage(session2) <- NULL ; usage(session2) <- matrix(usage.sum2,ncol = 1)


save(CH, CH.red, CH.Northern ,CH.Southern, MS, MS.red, file = "CLTCapthists.RData")

################################################

#load results from cluster
#1st for multi occs
load("Cluster/CLTResults.RData")

#both areas together
fit.all.null
fit.all.sex

#by area
fit.N.null
fit.N.sex

fit.S.null
fit.S.sex

#using MS object
fit.MS.NS.null
fit.MS.NS.sexa
fit.MS.NS.sexb

AIC(fit.MS.NS.null, fit.MS.NS.sexa, fit.MS.NS.sexb)

#so a single pmix estimate seems best
fit.MS.NS.sexa

load("Cluster/CLTResultsRed.RData")

#both areas together
fit.all.null.1occ
fit.all.sex.1occ

#using MS object
fit.MS.NS.null.1occ
fit.MS.NS.sexa.1occ
fit.MS.NS.sexb.1occ

AIC(fit.MS.NS.null.1occ, fit.MS.NS.sexa.1occ, fit.MS.NS.sexb.1occ)
fit.MS.NS.sexa.1occ

#and again using the complex mask
load("CLT/Cluster/CLTResultsRed2.RData")

fit.all.null
fit.all.sex

fit.MS.NS.null
fit.MS.NS.sexa
fit.MS.NS.sexb
fit.MS.NS.sexc
fit.MS.NS.sexd

AIC(fit.MS.NS.null, fit.MS.NS.sexa, fit.MS.NS.sexb, fit.MS.NS.sexc, fit.MS.NS.sexd)
#model with pmix ~ session + h2 best

############################################
#get Enrm numbers for CLT
plot(CLTsimp.nobuff)
plot(traps(CH), add = T)

plot(CLTsimp.buff)
plot(traps(CH), add = T)

#data summaries, note that grouped scenarios are excluded from scenarioSummary so will do one for each sex
CLTscen.F <- make.scenarios (trapsindex = 1, D = 0.00009, lambda0 = 3.8, sigma = 2388, noccasions = 1, detectfn = 14) 
CLTscen.M <- make.scenarios (trapsindex = 1, D = 0.000087, lambda0 = 3.12, sigma = 4785, noccasions = 1, detectfn = 14) 

CLTTraps <- traps(CH)
usage(CLTTraps) <- NULL

scenarioSummary(CLTscen.F , trapset = CLTTraps, maskset = CLTsimp.nobuff)
#En = 20 ; Er = 78
scenarioSummary(CLTscen.M , trapset = CLTTraps, maskset = CLTsimp.nobuff)
#En = 22 ; Er = 253
#so these numbers are fairly close to the actual data

#now look by area
#first Northern area
dat.N <- read.csv("CLT/CLTdatSession1.csv")
F <- as.matrix(table(dat.N$ID[dat.N$Sex=="F"], dat.N$Detector[dat.N$Sex=="F"]))
M <- as.matrix(table(dat.N$ID[dat.N$Sex=="M"], dat.N$Detector[dat.N$Sex=="M"]))

sum(F)
sum(M)

dat.S <- read.csv("CLT/CLTdatSession2.csv")
F <- as.matrix(table(dat.S$ID[dat.S$Sex=="F"], dat.S$Detector[dat.S$Sex=="F"]))
M <- as.matrix(table(dat.S$ID[dat.S$Sex=="M"], dat.S$Detector[dat.S$Sex=="M"]))
sum(F)
sum(M)

#data summaries
#for North
CLTscenN.F <- make.scenarios (trapsindex = 1, D = 0.00011, lambda0 = 3.1, sigma = 2600, noccasions = 1, detectfn = 14) 
CLTscenN.M <- make.scenarios (trapsindex = 1, D = 0.0001, lambda0 = 4.2, sigma = 3700, noccasions = 1, detectfn = 14) 

CLTNTraps <- traps(CH.Northern)
usage(CLTNTraps) <- NULL

scenarioSummary(CLTscenN.F , trapset = CLTNTraps, maskset = CLTsimp.nobuff)
#Females: En = 11 ; Er = 43
scenarioSummary(CLTscenN.M , trapset = CLTNTraps, maskset = CLTsimp.nobuff)
#Males: En = 12 ; Er = 116

#for South
CLTscenS.F <- make.scenarios (trapsindex = 1, D = 0.00008, lambda0 = 6.2, sigma = 1950, noccasions = 1, detectfn = 14) 
CLTscenS.M <- make.scenarios (trapsindex = 1, D = 0.00007, lambda0 = 3.3, sigma = 6500, noccasions = 1, detectfn = 14) 

CLTSTraps <- traps(CH.Southern)
usage(CLTSTraps) <- NULL

scenarioSummary(CLTscenS.F , trapset = CLTSTraps, maskset = CLTsimp.nobuff)
#Females: En = 10 ; Er = 42
scenarioSummary(CLTscenS.M , trapset = CLTSTraps, maskset = CLTsimp.nobuff)
#Males: En = 13 ; Er = 195

#so expected numbers are fairly similar to observed numbers, main diff is less recaps observed than expected in Sout

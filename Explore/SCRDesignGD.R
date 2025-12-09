#SCR design, Feb 2023
#this script firstly explores how to use secrdesign
#it then plays with generating different designs
#sims were then run on the cluster, and the results are explored here too

#Useful packages: elevatr (generates elevation from (x,y)), 
#secrdesign (for lots of things, including optimized grid), 
#"terrain" function from raster package

tri <- terrain(elev, opt = "TRI", filename = "TRI.tiff")

#takes an elevation raster and makes a terrain ruggedness raster

#also 'extract' function, extract value of a raster at a point

#get proposed design from GAoptim for Enrm, plus other designs
#then can use std secrdesign tools to run the simulations

#hint: first just generate raw data summaries and check 

############################
#exploring secrdesign
############################

#code to test expected numbers

library(secr) ; library(secrdesign) ; library(kofnGA)

# generate scenario dataframe
scen <- make.scenarios (trapsindex = 1:2, detectfn = 'HHN', D = c(20),
                        lambda0 = c(0.05,0.2), sigma = 1:2, noccasions = 1)

# set sigma = k/sqrt(D) for k = 50, 100
scen$sigma <- c(50,100)[scen$sigma] / sqrt(scen$D)

# detector layouts
traplist <- list(
#  make.grid(6,6, detector = 'multi', spacing = 20),
  make.grid(6,6, detector = 'proximity', spacing = 20),
#  make.grid(6,6, detector = 'count', spacing = 20),
#  make.grid(10,10, detector = 'multi', spacing = 20),
  make.grid(10,10, detector = 'proximity', spacing = 20)
#  make.grid(10,10, detector = 'count', spacing = 20)
)

# deterministic summary: expected counts
nrm <- scenarioSummary(scen, traplist)

# stochastic simulations
sumnrm <- function(CH) {
  c(
    n = nrow(CH),
    r = sum(CH) - nrow(CH),
    # moves(CH, names = TRUE) to dodge bug in secr < 4.5.2
    m = sum(unlist(secr::moves(CH, names = TRUE))>0, na.rm = TRUE),
    n2 = sum( apply(apply(CH, c(1,3), sum)>0, 1, sum) >1 )
  )
}

sim <- run.scenarios(scen, trapset = traplist, nrepl = 1000,
                     extractfn = sumnrm, fit = FALSE)

get.est <- function(df) {
 est <- df[4,2]
}

vec <- lapply(sims2$output[[1]],get.est)
mean(as.numeric(vec))

# collate deterministic and stochastic results
out <- data.frame(nrm, simout)
plotc <- function (v = 'n') {
  dtype <- c('multi','proximity','count')
  detector <- rep(dtype,2)[out$trapsindex]
  out2 <- matrix(nrow=3, ncol = 2, dimnames = list(dtype, c('mean','sd')))
  for (d in 1:3) {
    OK <- (detector == dtype[d])
    RB <- (out[OK,v] - out[OK, paste0('E',v)]) / out[OK, paste0('E', v)]
    # use log scales to spread values
    plot(out[OK, paste0('E',v)], out[OK,v], log='xy',
         xlab = 'Expected', ylab = 'simulated')
    abline(0,1) # y = x line
    # return mean and sd of estimated relative bias
    out2[d,] <- round(c(mean = mean(RB), sd = sd(RB)),5)
    mtext (side=3, line=0.2, paste(v, " ", dtype[d]), cex = 0.9)
  }
  out2
}

par(mfrow=c(4,3), mgp=c(2.3,0.6,0), mar=c(4,4,2,1), pty='s')
cat("Relative discrepancy between expected and simulated counts\n")
cat("Number of individuals\n")
plotc('n')
cat("Number of recaptures\n")
plotc('r')
cat("Number of movements\n")
plotc('m')
cat("Number of individuals at two or more detectors\n")
plotc('n2')

####################################################################
##scenarios with a regular box as mask, count detectors, 2 levels of D, 2 levels of lambda0
####################################################################

# generate scenario dataframe
scen <- make.scenarios (trapsindex = 1, detectfn = 'HHN', D = c(0.001),
                        lambda0 = c(0.5), sigma = c(3000), noccasions = 1)

# create trap array
traps <- make.grid(10,10, detector = 'count', spacing = 4000)

# deterministic summary: expected counts ; note that this function excludes grouped scenarios
nrm <- scenarioSummary(scen, traps)

# stochastic simulations
raw <- run.scenarios(scen, trapset = traps, nrepl = 50,
                     extractfn = identity, fit = FALSE)
summary(raw)
raw$output[[1]]
str(raw$output[[1]])  #output from 1st scenario, CHs for each of the 5 reps

#fit models as a 2nd step
sims1 <- fit.models(raw, fit.args = list(detectfn = 'HHN'), fit = TRUE, ncores = 4, byscenario = FALSE)
summary(sims1)

D.ests <- select.stats(sims1, parameter = "D", statistics = c("estimate","SE.estimate","RB", "RSE", "COV","lcl","ucl"))
L0.ests <- select.stats(sims1, parameter = "lambda0", statistics = c("estimate","SE.estimate","RB", "RSE", "COV","lcl","ucl"))
Sig.ests <- select.stats(sims1, parameter = "sigma", statistics = c("estimate","SE.estimate","RB", "RSE", "COV","lcl","ucl"))

summary(D.ests, dec = 6, fields = c("n", "mean", "se"), alpha = 0.05, type = "list")
summary(L0.ests, dec = 3, fields = c("n", "mean", "se"), alpha = 0.05, type = "list")
summary(Sig.ests, dec = 2, fields = c("n", "mean", "se"), alpha = 0.05, type = "list")

##fit model manually, to 1st CH from 7th scen
CH1 <- (raw$output[[7]][[1]])
mod1 <- secr.fit(CH1, detectfn = 14, buffer = 4000)
summary(mod1)
exampleOutput <- predict(mod1)

#sim data and fit model
sims2 <- run.scenarios(scen, trapset = traps, nrepl = 5, fit = TRUE)

#Note that the output to run.scenarios includes the result of extractfn for each rep
#when fit = T, the default output from each rep is the result of applying predict 
#i.e. a data frame of 'real' estimates with se's

sims3 <- run.scenarios(scen, trapset = traps, nrepl = 5, fit = TRUE, extractfn = predict)
find.param(sims3)
Dests <- select.stats(sims3, parameter = "sigma", statistics = c("estimate","SE.estimate","RB", "RSE", "COV","lcl","ucl"))
#note that validate() can be used to remove rouge values
apply(as.matrix(Dests$output[[1]]),2,mean)
#summary can be used on the 'selectedstatistics' object
summary(Dests, dec = 2, fields = c("n", "mean", "se"), alpha = 0.05, type = "list")

#plot estimates from 'selectedstatistics' object
par(mfrow = c(2,2))
plot(Dests, type = "hist", statistic = "estimate")
plot(Dests, type = "CI")
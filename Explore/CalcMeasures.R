
Dests <- c(res$output[[1]][[1]][1,2], res$output[[1]][[2]][1,2])
se.ests <- c(res$output[[1]][[1]][1,3], res$output[[1]][[2]][1,3])

#estimate
mean(Dests)
sd(Dests) / sqrt(2)

#SE
mean(se.ests)
sd(se.ests) / sqrt(2)

#RB
RBS <- (Dests -0.0008) / 0.0008
mean(RBS)
sd(RBS) / sqrt(2)

#RSE
RSEs <- se.ests / Dests
mean(RSEs)
sd(RSEs) / sqrt(2)



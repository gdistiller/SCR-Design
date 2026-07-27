#July 26, processing results vs varying K factor
library(dplyr)

#load results from K = 2 to 12
#K = 2 to 6
load("Cluster/Sims/GridResults40Ka.RData")

#K = 7 - 9
load("Cluster/Sims/GridResults40Kb.RData")

#K = 10-12
load("Cluster/Sims/GridResults40Kc.RData")

#combine results
KResults <- list(
  K2 = Grid.K2.results,
  K3 = Grid.K3.results,
  K4 = Grid.K4.results,
  K5 = Grid.K5.results,
  K6 = Grid.K6.results,
  K7 = Grid.K7.results,
  K8 = Grid.K8.results,
  K9 = Grid.K9.results,
  K10 = Grid.K10.results,
  K11 = Grid.K11.results,
  K12 = Grid.K12.results,
  K13 = Grid.K13.results,
  K14 = Grid.K14.results
)

K.results <- do.call(
  rbind,
  Map(function(x, nm) {
    df <- as.data.frame(x)
    df$K <- as.numeric(sub("K", "", nm))
    df
  },
  KResults,
  names(KResults))
)

CV.df <- K.results %>%
  group_by(K) %>%
  summarise(
    CVG1 = sd(D1) / 0.05,
    CVG2 = sd(D2) / (0.05 / first(K))
  )

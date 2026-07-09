library(secrdesign)
library(dplyr)
library(tidyr)
library(ggplot2)

source('optimalSpacing2G/optimalSpacing2G.R')

## User settings
T <- 40                  # number of traps
noccasions <- 1           # change if needed
H <- seq(100, 3000, 50)   # candidate spacings, metres

## Group-specific SCR parameters
sigma <- c(200, 3000)
lambda0 <- c(2, 2 / 15)
D <- c(1 / 20, 1 / 300)   # animals / hectare

## Fixed-T near-square rectangular trap grid
nxtraps <- floor(sqrt(T))
while (T %% nxtraps != 0) nxtraps <- nxtraps - 1
nytraps <- T / nxtraps

traps <- make.grid(
  nx = nxtraps,
  ny = nytraps,
  spacing = 3000,
  detector = "count"
)

mask <- make.mask(traps, buffer = 9000, spacing = 200)

## optimalSpacing per group

os.G1 = optimalSpacing(D = D[1], 
                       traps = traps, 
                       detectpar = list(lambda0 = lambda0[1], sigma = sigma[1]), 
                       noccasions = 1, 
                       plt = FALSE)

os.G2 = optimalSpacing(D = D[2], 
                       traps = traps, 
                       detectpar = list(lambda0 = lambda0[2], sigma = sigma[2]), 
                       noccasions = 1, 
                       plt = FALSE)

os.G1$rotRSE$optimum.spacing
os.G2$rotRSE$optimum.spacing

# max(min(n1 + n2, r1 + r2))
os2.sum_min = optimalSpacing2G(D1 = D[1], D2 = D[2],
                               traps0 = traps,
                               detectpar1 = list(lambda0 = lambda0[1], sigma = sigma[1]),
                               detectpar2 = list(lambda0 = lambda0[2], sigma = sigma[2]),
                               noccasions = 1,
                               criterion = c("sum_min"), 
                               spacing_m = seq(100,2000,50))

# max(min(n1, n2, r1, r2))
os2.all_min = optimalSpacing2G(D1 = D[1], D2 = D[2],
                               traps0 = traps,
                               detectpar1 = list(lambda0 = lambda0[1], sigma = sigma[1]),
                               detectpar2 = list(lambda0 = lambda0[2], sigma = sigma[2]),
                               noccasions = 1,
                               criterion = c("all_min"), 
                               spacing_m = seq(100,2000,50))

# min(mean(CV1, CV2))
os2.min_mean_CV = optimalSpacing2G(D1 = D[1], D2 = D[2],
                                   traps0 = traps,
                                   detectpar1 = list(lambda0 = lambda0[1], sigma = sigma[1]),
                                   detectpar2 = list(lambda0 = lambda0[2], sigma = sigma[2]),
                                   noccasions = 1,
                                   criterion = c("min_mean_CV"),
                                   spacing_m = seq(100,2000,50))

# min(max(CV1, CV2))
os2.min_max_CV = optimalSpacing2G(D1 = D[1], D2 = D[2],
                                  traps0 = traps,
                                  detectpar1 = list(lambda0 = lambda0[1], sigma = sigma[1]),
                                  detectpar2 = list(lambda0 = lambda0[2], sigma = sigma[2]),
                                  noccasions = 1,
                                  criterion = c("min_max_CV"),
                                  spacing_m = seq(100,2000,50))

# max((En2_1 + En2_2))
os2.max_mean_En2 = optimalSpacing2G(D1 = D[1], D2 = D[2],
                                    traps0 = traps,
                                    detectpar1 = list(lambda0 = lambda0[1], sigma = sigma[1]),
                                    detectpar2 = list(lambda0 = lambda0[2], sigma = sigma[2]),
                                    noccasions = 1,
                                    criterion = c("max_mean_En2"),
                                    spacing_m = seq(100,2000,50))

# max(min(En2_1, En2_2))
os2.max_min_En2 = optimalSpacing2G(D1 = D[1], D2 = D[2],
                                   traps0 = traps,
                                   detectpar1 = list(lambda0 = lambda0[1], sigma = sigma[1]),
                                   detectpar2 = list(lambda0 = lambda0[2], sigma = sigma[2]),
                                   noccasions = 1,
                                   criterion = c("max_min_En2"), 
                                   spacing_m = seq(100,2000,50))

# compare
os2.sum_min$optimum.spacing
os2.all_min$optimum.spacing
os2.min_mean_CV$optimum.spacing
os2.min_max_CV$optimum.spacing
os2.max_mean_En2$optimum.spacing
os2.max_min_En2$optimum.spacing

# hold G1 fixed (D[1], lambda0[1], sigma[1]) and vary G2 as factor of G1.
# for each value, find optimal spacing using criterion = "min_mean_CV"

os.by_k = list()
for(k in 1:15){
  
  st = k * sigma[1]
  lt = lambda0[1] / k
  Dt = D[1] / k 
  
  os.by_k[[k]] = optimalSpacing2G(D1 = D[1], D2 = Dt,
                                  traps0 = traps,
                                  detectpar1 = list(lambda0 = lambda0[1], sigma = sigma[1]),
                                  detectpar2 = list(lambda0 = lt, sigma = st),
                                  noccasions = 1,
                                  criterion = c("min_mean_CV"),
                                  spacing_m = seq(0,1500,50))
  
}

os.by_k_all <- sapply(os.by_k, function(z) z$optimum.spacing)
os.by_k_all

os.t = os.by_k[[15]]$values

# plot how CV1 and CV2 change with spacing and k
os_df <- do.call(rbind, lapply(os.by_k, `[[`, "values"))
os_long <- os_df %>%
  mutate(k = rep(1:15, each = 31)) |> 
  pivot_longer(
    cols = c(CV1, CV2),
    names_to = "CV",
    values_to = "value"
  )

ggplot(os_long, aes(x = spacing, y = value, colour = CV)) +
  geom_line() +
  facet_wrap(~ k) +
  labs(
    x = "Spacing",
    y = "CV",
    colour = "Criterion"
  ) +
  theme_bw()

# manually check results for k = 10, spacing = 50 (note: works out)
kk = 10; ss = 50
xx = os.by_k[[kk]]$values
xx[xx$spacing == ss, ]

# verify
traps2 <- make.grid(
  nx = nxtraps, 
  ny = nytraps,
  spacing = ss,
  detector = "count"
)

# group 1
Enrm(D = D[1] , 
     traps = traps2, 
     mask = mask,
     detectpar = list(lambda0 = lambda0[1], sigma = sigma[1]), 
     noccasions = 1)

# group 2
Enrm(D = D[1] / kk, 
     traps = traps2, 
     mask = mask,
     detectpar = list(lambda0 = lambda0[1] / kk, sigma = sigma[1] * kk), 
     noccasions = 1)

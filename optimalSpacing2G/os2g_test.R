library(secrdesign)
source('optimalSpacing2G/optimalSpacing2G.R')
grid <- make.grid(7,7,50)
os1 = optimalSpacing(D = 5, traps = grid, detectpar = list(lambda0 = 0.2, sigma = 40), 
               noccasions = 5, plt = FALSE)
os1$rotRSE$optimum.spacing

os2 = optimalSpacing2G(D1 = 5, D2 = 5,
                 traps0 = grid,
                 detectpar1 = list(lambda0 = 0.2, sigma = 40),
                 detectpar2 = list(lambda0 = 0.2, sigma = 40),
                 noccasions = 5,
                 criterion = c("sum_min"), # all_min
                 spacing_m = seq(0,120,1))

os2$optimum.spacing
xx=os2$values
xx[1,]
xx[xx$spacing == round(os2$optimum.spacing), ]

os2 = optimalSpacing2G(D1 = 5, D2 = 5,
                       traps0 = grid,
                       detectpar1 = list(lambda0 = 0.2, sigma = 80),
                       detectpar2 = list(lambda0 = 0.2, sigma = 40),
                       noccasions = 5,
                       criterion = c("sum_min"), # all_min
                       spacing_m = seq(0,120,1))

os2$optimum.spacing
xx=os2$values
xx[xx$spacing == round(os2$optimum.spacing), ]

os2_all = optimalSpacing2G(D1 = 5, D2 = 5,
                           traps0 = grid,
                           detectpar1 = list(lambda0 = 0.2, sigma = 40),
                           detectpar2 = list(lambda0 = 0.2, sigma = 40),
                           noccasions = 5,
                           criterion = "all_min",
                           spacing_m = seq(40,120,1))

os2_mean = optimalSpacing2G(D1 = 5, D2 = 5,
                            traps0 = grid,
                            detectpar1 = list(lambda0 = 0.2, sigma = 40),
                            detectpar2 = list(lambda0 = 0.2, sigma = 40),
                            noccasions = 5,
                            criterion = "min_mean",
                            spacing_m = seq(40,120,1))

os2_max = optimalSpacing2G(D1 = 5, D2 = 5,
                           traps0 = grid,
                           detectpar1 = list(lambda0 = 0.2, sigma = 40),
                           detectpar2 = list(lambda0 = 0.2, sigma = 40),
                           noccasions = 5,
                           criterion = "min_max",
                           spacing_m = seq(40,120,1))

stopifnot(isTRUE(all.equal(os2_mean$values$CV1,
                           1 / sqrt(pmin(os2_mean$values$En1, os2_mean$values$Er1)))))
stopifnot(isTRUE(all.equal(os2_mean$values$CV2,
                           1 / sqrt(pmin(os2_mean$values$En2, os2_mean$values$Er2)))))
stopifnot(isTRUE(all.equal(os2_mean$values$crit,
                           -rowMeans(os2_mean$values[, c("CV1", "CV2")]))))
stopifnot(isTRUE(all.equal(os2_max$values$crit,
                           -pmax(os2_max$values$CV1, os2_max$values$CV2))))
stopifnot(abs(os2_all$optimum.spacing - os2_mean$optimum.spacing) < 1)
stopifnot(abs(os2_mean$optimum.spacing - os2_max$optimum.spacing) < 1)

os2_en2_sum = optimalSpacing2G(D1 = 5, D2 = 5,
                               traps0 = grid,
                               detectpar1 = list(lambda0 = 0.2, sigma = 40),
                               detectpar2 = list(lambda0 = 0.2, sigma = 40),
                               noccasions = 5,
                               criterion = "max_mean_En2",
                               spacing_m = seq(40,120,1))

os2_en2_min = optimalSpacing2G(D1 = 5, D2 = 5,
                               traps0 = grid,
                               detectpar1 = list(lambda0 = 0.2, sigma = 40),
                               detectpar2 = list(lambda0 = 0.2, sigma = 40),
                               noccasions = 5,
                               criterion = "max_min_En2",
                               spacing_m = seq(40,120,1))

stopifnot(isTRUE(all.equal(os2_en2_sum$values$crit,
                           os2_en2_sum$values$En2plus1 + os2_en2_sum$values$En2plus2)))
stopifnot(isTRUE(all.equal(os2_en2_min$values$crit,
                           pmin(os2_en2_min$values$En2plus1, os2_en2_min$values$En2plus2))))

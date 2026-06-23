library(secrdesign)
source('optimalSpacing2G/optimalSpacing2G.R')
grid <- make.grid(7,7,1000)
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

os3 = optimalSpacing2G(D1 = 5, D2 = 5/15,
                       traps0 = grid,
                       detectpar1 = list(lambda0 = 0.2, sigma = 40),
                       detectpar2 = list(lambda0 = 0.2/15, sigma = 40*15),
                       noccasions = 5,
                       criterion = c("sum_min"), # all_min
                       spacing_m = seq(100,4000,100))

os3$optimum.spacing

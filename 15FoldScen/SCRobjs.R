#recreating mask as I can't find the code
#to be safe I am going to reuse the existing objs and replace the cluster OS masks

load("15FoldScen/Cluster/ProposedDesigns/SCRObjs.RData")

#store existing slots
Fullmask <- SCR.objs$Full
TwoSig <- SCR.objs$`2 sig`

#now recreate the cluster OS mask for 40 and 120 traps

ClustOS40 <- create.extent(sigma = 3000, buff.factor = 3.12, res = 200)
mask1 <- ClustOS40 [[1]]
clust.locs1 <- ClustOS40 [[2]]
clustlocs.sf1 <- ClustOS40 [[3]]

ClustOS40 <- list("SCR Mask" = mask1, "SCR trap locs" = clust.locs1, "SF traps polygon" = clustlocs.sf1)

ClustOS120 <- create.extent(sigma = 3000, buff.factor = 3.14, res = 200)
mask2 <- ClustOS120[[1]]
clust.locs2 <- ClustOS120[[2]]
clustlocs.sf2 <- ClustOS120[[3]]

ClustOS120 <- list("SCR Mask" = mask2, "SCR trap locs" = clust.locs2, "SF traps polygon" = clustlocs.sf2)


SCR.objs <- list("Full" = Fullmask, "2 sig" = TwoSig, "OS (40)" = ClustOS40, "OS (120)" = ClustOS120)

save(SCR.objs, file = "15FoldScen/Cluster/ProposedDesigns/SCRObjs.RData")
save(SCR.objs, file = "15FoldScen/Cluster/Sims/SCRObjs.RData")

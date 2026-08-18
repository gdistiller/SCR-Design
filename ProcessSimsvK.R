#July 26, processing results vs varying K factor
library(dplyr)
library(tidyr)
library(ggplot2)

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

#get vectors of true values for the different Ks
#spacings used
sigma1 = 200 ; L01 = 2; D1 = 0.05
L0.vec <- NULL
D.vec <- NULL 
Sig.vec <- NULL

for (k in 2:14){
  Sig.vec[k-1] <-  k * sigma1
  L0.vec[k-1] <- L01 / k
  D.vec[k-1] <- D1 / k
}

Sig.vec
L0.vec
D.vec

#filter using find.rouge() based on the same mag factor, and construct_df 
mag.factor <- 10

K.values <- sort(unique(K.results$K))

g1 <- find.rogue(
  K.results[, 1:6],
  mag = mag.factor,
  true = c(D1, L01, sigma1)
)

g1$K <- K.results[g1$row, "K"]

g2 <- do.call(
  rbind,
  lapply(seq_along(K.values), function(i) {
    
    k <- K.values[i]
    
    dat.k <- K.results[K.results$K == k, ]
    
    out <- find.rogue(
      dat.k[, 7:12],
      mag = mag.factor,
      true = c(D.vec[i], L0.vec[i], Sig.vec[i])
    )
    
    out$K <- k
    out
  })
)

K.clean <- merge(
  g1,
  g2,
  by = c("row", "K"),
  all = TRUE
)

K.clean <- K.clean[match(rownames(K.results), K.clean$row),]


CV.df <- K.clean %>%
  group_by(K) %>%
  summarise(
    n1 = sum(!is.na(D.x)),
    n2 = sum(!is.na(D.y)),
    CVG1 = sd(D.x, na.rm = T) / 0.05,
    CVG2 = sd(D.y, na.rm = T) / (0.05 / first(K))
  ) %>%
  mutate(
    SE.CVG1 = CVG1 / sqrt(2 * (n1 - 1)),
    SE.CVG2 = CVG2 / sqrt(2 * (n2 - 1))
  )

#calculating for K = 15 from main results file
load("15FoldScen/AllResults4.RData")
results.df <- Allresults_4[[1]]
RSE.df <- results.df %>%
  filter(design == "Grid 800m",
         param == "D") %>%
  group_by(stratum) %>%
  summarise(
    Truth = first(truth),
    SD = sd(estimate, na.rm = T),
    RSE = SD / Truth,
    .groups = "drop"
  )
RSE.df

0.346 / sqrt(2*487-1) ; 0.586 / sqrt(2*462-1)

#Add results for K = 15
CV.df <- rbind(CV.df, c(15, 487, 462,0.346,0.586, 0.011, 0.0193))

#reshape to long format
CV.long <- CV.df %>%
  select(K, CVG1, CVG2, SE.CVG1, SE.CVG2) %>%
  pivot_longer(
    cols = c(CVG1, CVG2),
    names_to = "Group",
    values_to = "CV"
  ) %>%
  mutate(
    SE = ifelse(Group == "CVG1", SE.CVG1, SE.CVG2),
    Group = factor(Group,
                   levels = c("CVG1", "CVG2"),
                   labels = c("Group 1", "Group 2"))
  )

#Plot, with and without SEs
p1.RSE.K <- ggplot(CV.long,
       aes(x = K,
           y = CV,
           colour = Group)) +
  
  geom_line(linewidth = 1.2) +
  
  geom_point(size = 3) +
  
  geom_errorbar(aes(ymin = CV - 1.96 * SE,
                    ymax = CV + 1.96 * SE),
                width = 0.15,
                linewidth = 0.5) +
  
  scale_colour_manual(
    values = c(
      "Group 1" = "#0072B2",  # blue
      "Group 2" = "#D55E00"   # orange
    )
  ) +
  
  scale_x_continuous(breaks = unique(CV.long$K)) +
  
  labs(x = TeX("Between-group difference factor ($R$)"),
               y = expression(CV(hat(D)^{(g)}))) +
  
  theme_classic(base_size = 13) +
  
  theme(
    axis.title = element_text(face = "bold"),
    legend.position = c(0.15, 0.85),
    legend.background = element_blank(),
    legend.title = element_blank()
  )

#and without the SE bar
p2.RSE.K <- ggplot(CV.long,
       aes(x = K,
           y = CV,
           colour = Group)) +
  
  geom_line(linewidth = 1.2) +
  geom_point(size = 3) +
  
  scale_colour_manual(
    values = c(
      "Group 1" = "#0072B2",
      "Group 2" = "#D55E00"
    )
  ) +
  
  scale_x_continuous(breaks = unique(CV.long$K)) +
  
  labs(
    x = TeX("Between-group difference factor ($k$)"),
  ) +
  
  ylab(TeX(r'($CV(\hat{D}^{(g)})$)')) +
  
  theme_classic(base_size = 13) +
  
  theme(
    axis.title = element_text(face = "bold"),
    legend.position = "top"
  )

ggsave("15FoldScen/figures/CV_vs_Ka.pdf", p1.RSE.K, width = 6, height = 4)
ggsave("15FoldScen/figures/CV_vs_Kb.pdf", p2.RSE.K, width = 6, height = 4)

#boxplots of the estimates, 1st make long format
D.plot <- K.clean %>%
  select(K, D.x, D.y) %>%
  pivot_longer(
    cols = c(D.x, D.y),
    names_to = "Group",
    values_to = "Estimate"
  ) %>%
  mutate(
    Group = recode(Group,
                   "D.x" = "Group 1",
                   "D.y" = "Group 2")
  )

ggplot(D.plot,
       aes(x = factor(K),
           y = Estimate,
           fill = Group)) +
  geom_boxplot(position = position_dodge(width = 0.8)) +
  labs(
    x = "K",
    y = "D estimate"
  ) +
  theme_bw()



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

CV.df <- K.results %>%
  group_by(K) %>%
  summarise(
    n = n(),
    CVG1 = sd(D1) / 0.05,
    CVG2 = sd(D2) / (0.05 / first(K))
  ) %>%
  mutate(
    SE.CVG1 = CVG1 / sqrt(2 * (n - 1)),
    SE.CVG2 = CVG2 / sqrt(2 * (n - 1))
  )

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

#PLot, with and without SEs
p.RSE.K <- ggplot(CV.long,
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
  
  labs(
    x = "Difference factor (K)",
    y = "Relative standard error (RSE) of density",
    colour = "Group"
  ) +
  
  theme_classic(base_size = 13) +
  
  theme(
    axis.title = element_text(face = "bold"),
    legend.position = c(0.15, 0.85),
    legend.background = element_blank(),
    legend.title = element_blank()
  )

#and witout the SE bar
p.RSE.K <- ggplot(CV.long,
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
    x = expression(K),
    y = "Relative standard error (RSE) of density",
    colour = NULL
  ) +
  
  theme_classic(base_size = 13) +
  
  theme(
    axis.title = element_text(face = "bold"),
    legend.position = "top"
  )

ggsave("15FoldScen/figures/CV_vs_K.pdf", p.RSE.K, width = 6, height = 4)

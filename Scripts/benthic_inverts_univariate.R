library(tidyverse) # for figures, pipes etc.
library(vegan) # for NMDS
library(ggpubr) # for arranging figures
library(agricolae)
library(Hmisc)
library(car)
library(viridis) # colours


# Load Data ---------------------------------------------------------------

benthic <- read.csv("Data/Benthic_vegetation_QCC.csv") # occurrences = 1 removed
str(benthic)

benthic.data <- benthic %>% select(Oligochaeta:Hydroptilidae)
benthic.env <- benthic %>% select(Site.ID:Collection.date)


# Univariate Data ---------------------------------------------------------
richness <- rowSums(benthic.data > 0) # species richness
rich <- specnumber(benthic.data) # species richness in vegan
abundance <- rowSums(benthic.data) # abundance
H <- diversity(benthic.data) # Shannon Weiner
D1 <- diversity(benthic.data, index = "simpson") #default is base log, but can change it
J <- H/log(specnumber(benthic.data))


# make new little data frame
benthic.uni<- benthic.env
benthic.uni$rich <- richness
benthic.uni$abundance <- abundance
benthic.uni$H <- H
benthic.uni$D1 <- D1
benthic.uni$J <- J

write.csv(benthic.uni, "Data/benthic_invertebrates_univariate.csv")

colnames(benthic.uni)


# Univariate Analyses -----------------------------------------------------

benthic.uni <- read.csv("Data/benthic_invertebrates_univariate.csv")


abundance.lm <- lm(abundance ~ Habitat, data = benthic.uni)
Anova(abundance.lm, type = 3)

abundance.hsd <- HSD.test(abundance.lm, "Habitat")

#Anova Table (Type III tests)

#Response: abundance
#               Sum Sq Df F value   Pr(>F)   
#(Intercept)  2726112  1  0.8434 0.368396   
#Habitat     47207700  2  7.3022 0.003696 **
#Residuals   71113597 22  

#          abundance       std r Min  Max  Q25    Q50     Q75
#Invaded      583.75  622.2551 8  73 1767  233  320.5  661.25
#Treated   3375.00 2920.3388 9  42 9851 2194 2472.0 3180.00
#Uninvaded    445.25  158.6368 8 234  781  381  406.5  476.00


#           abundance groups
#Treated   3375.00      a
#Invaded      583.75      b
#Uninvaded    445.25      b


rich.lm <- lm(rich ~ Habitat, data = benthic.uni)
Anova(rich.lm, type = 3)

#Anova Table (Type III tests)

#Response: rich
#             Sum Sq Df  F value   Pr(>F)    
#(Intercept) 1800.00  1 155.3703 1.91e-11 ***
#Habitat       74.57  2   3.2181  0.05944 .  
#Residuals    254.87 22 


benthic.uni %>% 
  group_by(Habitat) %>% 
  mutate(mean.s = mean(rich),
         median.s = median(rich),
         N = length(abundance),
         sd.s = sd(rich),
         sterr.s = (sd.s/(sqrt(N))),
         mean.ab = mean(abundance),
         median.ab = median(abundance),
         sd.ab = sd(abundance),
         sterr.ab = (sd.ab/sqrt(N))) -> benthic.uni


# Univariate Figures ------------------------------------------------------

ggplot(data = benthic.uni, 
       aes(y = abundance, x = Habitat)) +
  geom_point()

ggplot(data = benthic.uni, 
       aes(y = rich, x = Habitat)) +
  geom_point()


# Density plots 

abundance <- ggplot(data = benthic.uni, 
                    aes(x = abundance, group = Habitat, fill = Habitat)) +
  geom_density(adjust = 1.5, alpha = 0.6) +
  theme_classic() +
  xlab("Abundance") +
  ylab("Density") +
  theme(legend.position = "none") +
  geom_vline(aes(xintercept = mean.ab, colour = Habitat),
             size = 1,
             show.legend = FALSE) +
  scale_fill_viridis(discrete = TRUE) +
  scale_colour_viridis(discrete = TRUE)


richness <- ggplot(data = benthic.uni, 
                   aes(x = rich, group = Habitat, fill = Habitat)) +
  geom_density(adjust = 1.5, alpha = 0.6) +
  xlim(0, 42) +
  theme_classic() +
  xlab("Species Richness") +
  ylab("Density") +
  theme(legend.position = c(0.8, 0.8)) +
  geom_vline(aes(xintercept = mean.s, colour = Habitat),
             size = 1,
             show.legend = FALSE) +
  scale_fill_viridis(discrete = TRUE) +
  scale_colour_viridis(discrete = TRUE)



ggarrange(abundance, richness,
          labels = "AUTO",
          hjust = c(-5.5),
          vjust = 2)

# Violin Plots

(violin.s <- ggplot(data = benthic.uni, 
                    aes(x = Habitat, y = rich, group = Habitat, colour = Habitat)) +
    geom_violin(trim = FALSE, lwd = 0.75, colour = "black") +
    geom_point(size = 3) +
    theme_classic(base_size = 14) +
    xlab(" ") +
    ylab("Species Richness") +
    theme(legend.position = "none",
          axis.text = element_text(size = 14)) +
    ylim(0, 30) +
    scale_colour_manual(values = c("#969696","#35978f", "#2166ac")) +
    stat_summary(     
      aes(shape = Habitat),
      colour = "black",
      fun.data = "mean_se", fun.args = list(mult = 1), 
      geom = "pointrange", size = 1,
      position = position_dodge(0.8)) +
    theme(panel.border = element_rect(fill = NA)))



(violin.ab <- ggplot(data = benthic.uni, 
                     aes(x = Habitat, y = abundance, group = Habitat, colour = Habitat)) +
    geom_violin(trim = FALSE, lwd = 0.75, colour = "black") +
    geom_point(size = 3) +
    theme_classic(base_size = 14) +
    xlab(" ") +
    ylab("Abundance") +
    theme(legend.position = "none",
          axis.text = element_text(size = 14)) +
    scale_colour_manual(values = c("#969696","#35978f", "#2166ac")) +
    annotate("text", x = 1:3, y = c(4600, 12000, 4600),
             label = c("a", "b", "a"),
             size = 4.5) +
    stat_summary(     
      aes(shape = Habitat),
      colour = "black",
      fun.data = "mean_se", fun.args = list(mult = 1), 
      geom = "pointrange", size = 1,
      position = position_dodge(0.8)) +
    theme(panel.border = element_rect(fill = NA)))

(violin.benthic <- ggarrange(violin.ab, violin.s,
                             labels = "AUTO",
                             hjust = c(-8, -5.5),
                             vjust = 2))

ggsave("Figures/benthic_violin.TIFF", violin.benthic,
       dpi = 300,
       height = 6.1,
       width = 11.9,
       units = "in")

# Ridge Plots

library(ggridges)

ggplot(benthic.uni, aes(x = abundance, y = Habitat, fill = Habitat)) +
  geom_density_ridges() +
  theme_ridges() + 
  theme(legend.position = "none") +
  scale_fill_viridis(discrete = TRUE)


ggplot(benthic.uni, aes(x = rich, y = Habitat, fill = Habitat)) +
  geom_density_ridges() +
  theme_ridges() + 
  theme(legend.position = "none") +
  scale_fill_viridis(discrete = TRUE)



### Additional Diversity Indices ####

(H.vio <- ggplot(data = benthic.uni, 
                 aes(x = Habitat, y = H, group = Habitat, fill = Habitat)) +
   geom_violin(trim = FALSE) +
   theme_classic(base_size = 14) +
   xlab(" ") +
   ylab("Shannon Weiner (H)") +
   theme(legend.position = "none",
         axis.text = element_text(size = 14)) +
   scale_fill_viridis(discrete = TRUE))

(D.vio <- ggplot(data = benthic.uni, 
                 aes(x = Habitat, y = D1, group = Habitat, fill = Habitat)) +
    geom_violin(trim = FALSE) +
    theme_classic(base_size = 14) +
    xlab(" ") +
    ylab("Simpson's Diversity") +
    theme(legend.position = "none",
          axis.text = element_text(size = 14)) +
    scale_fill_viridis(discrete = TRUE))

(J.vio <- ggplot(data = benthic.uni, 
                 aes(x = Habitat, y = J, group = Habitat, fill = Habitat)) +
    geom_violin(trim = FALSE) +
    theme_classic(base_size = 14) +
    xlab(" ") +
    ylab("Pielou's J") +
    theme(legend.position = "none",
          axis.text = element_text(size = 14)) +
    scale_fill_viridis(discrete = TRUE))


ggarrange(H.vio, D.vio, J.vio,
          nrow = 3,
          labels = "AUTO",
          hjust = -6)


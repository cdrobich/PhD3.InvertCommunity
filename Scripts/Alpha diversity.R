
library(tidyverse) # for figures, pipes etc.
library(vegan) # for NMDS
library(ggpubr) # for arranging figures
library(agricolae)
library(Hmisc)
library(car)

data <- read.csv("Data/Inverts_sum_richness.csv")

data$Year <- as.factor(data$Year)
data$Treatment <- as.factor(data$Treatment)
str(data)
unique(data$Treatment)

inverts <- data %>% 
  mutate(Treatment = fct_recode(Treatment,
                               "Uninvaded" = "Open Water",
                               "Uninvaded" = "Open water",
                               "Uninvaded" = "Meadow",
                               "Uninvaded" = "Meadow ",
                               "Uninvaded" = "Typha",
                               "Herbicide-treated" = "Restored"))
unique(inverts$Treatment)


sample.size <- inverts %>% 
  group_by(Treatment) %>% 
  summarise(sample = sum(N))

#Treatment           sample
#1 Invaded               88
#2 Herbicide-treated     83
#3 Uninvaded             88


abundance.lm <- lm(Sum ~ Treatment, data = inverts)
Anova(abundance.lm, type = "2")
tukey1 <- HSD.test(abundance.lm, "Treatment")

#Response: Sum
#Sum Sq Df F value   Pr(>F)    
#Treatment 20656332  2  14.154 1.29e-05 ***
#Residuals 37214525 51 

#Sum       std  r Min  Max    Q25    Q50     Q75
#Herbicide-treated 1756.5000 1343.4022 18 149 5055 542.00 1703.5 2622.25
#Invaded            383.3889  204.4773 18 101  757 268.75  311.0  554.25
#Uninvaded          515.6111  585.2771 18  51 2257 174.25  256.0  568.25

#Sum groups
#Herbicide-treated 1756.5000      a
#Uninvaded          515.6111      b
#Invaded            383.3889      b

rich.lm <- lm(Richness ~ Treatment, data = inverts)
Anova(rich.lm, type = "2")
tukey2 <- HSD.test(rich.lm, "Treatment")

#Response: Richness
#Sum Sq Df F value    Pr(>F)    
#Treatment  942.81  2  10.503 0.0001514 ***
#Residuals 2289.11 51

#Richness      std  r Min Max   Q25  Q50   Q75
#Herbicide-treated 15.33333 4.458963 18   8  23 11.25 15.5 18.75
#Invaded           25.33333 6.782330 18  10  36 20.25 25.0 29.75
#Uninvaded         22.22222 8.292843 18  10  37 17.00 22.5 27.75

#Richness groups
#Invaded           25.33333      a
#Uninvaded         22.22222      a
#Herbicide-treated 15.33333      b


inverts <- inverts %>% 
  mutate(Treatment = fct_relevel(Treatment,
                                 "Invaded",
                                 "Herbicide-treated",
                                 "Uninvaded"))

Abundance <- ggplot(inverts, aes(x = Treatment, y = Sum)) +
  geom_jitter(
    aes(shape = Treatment, color = Treatment),
    position = position_jitterdodge(jitter.width = 0.5, dodge.width = 0.9),
    size = 5) +
  theme_classic(base_size = 20) +
  stat_summary(
    aes(shape = Treatment),
    fun.data = "mean_se", fun.args = list(mult = 1),
    geom = "pointrange", size = 1,
    position = position_dodge(0.6)) +
  labs(x = " ",
       y = "Abundance") +
  scale_color_manual(values = c("#9970ab", "#1b7837", "#2166ac")) +
  theme(panel.border = element_rect(fill = NA)) +
  scale_y_continuous(breaks = seq(0, 6000, by = 500)) +
  theme(legend.position = "blank")
  
Abundance


Rich <- ggplot(inverts, aes(x = Treatment, y = Richness)) +
  geom_jitter(
    aes(shape = Treatment, color = Treatment),
    position = position_jitterdodge(jitter.width = 0.5, dodge.width = 0.9),
    size = 5) +
  theme_classic(base_size = 20) +
  stat_summary(
    aes(shape = Treatment),
    fun.data = "mean_se", fun.args = list(mult = 1),
    geom = "pointrange", size = 1,
    position = position_dodge(0.6)) +
  labs(x = " ",
       y = "Richness (Family)") +
  ylim(0, 50) +
  scale_color_manual(values = c("#9970ab", "#1b7837", "#2166ac")) +
  theme(panel.border = element_rect(fill = NA)) +
  theme(legend.position = "blank")

Rich
 
Panel <- ggarrange(Abundance, Rich,
                   labels = c("A","B"),
                   hjust = c(-9, -7),
                   vjust = 2.5)

Panel

ggsave("Figures/Inverts_AbSpanel.jpeg")















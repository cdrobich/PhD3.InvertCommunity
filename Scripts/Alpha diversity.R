
library(tidyverse) # for figures, pipes etc.
library(vegan) # for NMDS
library(ggpubr) # for arranging figures
library(agricolae)
library(Hmisc)

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















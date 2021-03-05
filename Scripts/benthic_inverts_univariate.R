library(tidyverse) # for figures, pipes etc.
library(vegan) # for NMDS
library(ggpubr) # for arranging figures
library(agricolae)
library(Hmisc)
library(car)
library(viridis) # colours

library(emmeans)
library(rstatix)

library(performance)
library(see)

# Load Data ---------------------------------------------------------------

benthic <- read.csv("Data/Aquatic/aquatic_inverts_rares2.csv") # occurrences <= 2 removed
str(benthic)

benthic.data <- benthic %>% select(Oligochaeta:Leptoceridae)
benthic.env <- benthic %>% select(Site.ID:Collection.date)


# Univariate Data ---------------------------------------------------------
richness <- rowSums(benthic.data > 0) # species richness
rich <- specnumber(benthic.data) # species richness in vegan
abundance <- rowSums(benthic.data) # abundance
H <- diversity(benthic.data) # Shannon Weiner
D1 <- diversity(benthic.data, index = "simpson") # 1 - D
J <- H/log(specnumber(benthic.data))


# make new little data frame
benthic.uni<- benthic.env
benthic.uni$rich <- richness
benthic.uni$abundance <- abundance
benthic.uni$H <- H
benthic.uni$D1 <- D1
benthic.uni$J <- J

benthic.uni1 <- benthic.uni

write.csv(benthic.uni, "Data/Aquatic/benthic_invertebrates_univariate.csv")

colnames(benthic.uni)

install.packages("plotrix")
library(plotrix)

summary <- benthic.uni %>% group_by(Habitat) %>% 
  summarise(across(
    .cols = where(is.numeric),
    .fns = list(Mean = mean, SD = sd, SE = std.error), na.rm = TRUE,
    .names = "{col}_{fn}"
  ))

sum <- summary %>% t %>% as.data.frame

#Habitat           Invaded    Treated  Uninvaded
#Depth_Mean       39.92500   32.05556   22.93750
#Depth_SD         17.86719   11.11925   10.80317
#Depth_SE         6.317005   3.706418   3.819496
#rich_Mean        15.00000   15.66667   18.50000
#rich_SD          1.309307   3.201562   4.440077
#rich_SE          0.462910   1.067187   1.569804
#abundance_Mean    583.875   3374.667    444.625
#abundance_SD     622.1540  2920.5174   158.4297
#abundance_SE    219.96464  973.50579   56.01337
#H_Mean           1.426528   1.766169   1.909067
#H_SD            0.3773256  0.2415291  0.1909684
#H_SE           0.13340473 0.08050969 0.06751754
#D1_Mean         0.6058429  0.7555122  0.7799573
#D1_SD          0.13679184 0.09687545 0.05650638
#D1_SE          0.04836322 0.03229182 0.01997802
#J_Mean          0.5294864  0.6506698  0.6606395
#J_SD           0.15016095 0.10256570 0.05272685
#J_SE           0.05308991 0.03418857 0.01864176

# Histograms --------------------------------------------------------------

benthic.uni <- read.csv("Data/Aquatic/benthic_invertebrates_univariate.csv")

# abundance histogram
ggplot(benthic.uni, aes(x = abundance)) + 
  geom_histogram(binwidth = 100,
                 color="black", fill="white")

#richness histogram
ggplot(benthic.uni, aes(x = rich)) + 
  geom_histogram(binwidth = 1.5,
                 color="black", fill="white")

# shannon # looks ok, bit left skewed

ggplot(benthic.uni, aes(x = H)) + 
  geom_histogram(binwidth = 0.25,
                 color="black", fill="white")

#simp # bit left skewed

ggplot(benthic.uni, aes(x = D1)) + 
  geom_histogram(binwidth = 0.1,
                 color="black", fill="white")


#pie  # looks nice
ggplot(benthic.uni, aes(x = J)) + 
  geom_histogram(binwidth = 0.1,
                 color="black", fill="white")


## Transformations
benthic.uni$logAb <- log(benthic.uni$abundance)
benthic.uni$sqH <- sqrt(benthic.uni$H)
benthic.uni$sqD1 <- log(benthic.uni$D1 + 1)
benthic.uni$sqJ <- sqrt(benthic.uni$J)
benthic.uni$sqrich <- sqrt(benthic.uni$rich)


# abundance histogram
ggplot(benthic.uni, aes(x = logAb)) + 
  geom_histogram(binwidth = 1,
                 color="black", fill="white")

ggplot(benthic.uni, aes(x = sqH)) + 
  geom_histogram(binwidth = 0.1,
                 color="black", fill="white")

ggplot(benthic.uni, aes(x = sqD1)) + 
  geom_histogram(binwidth = 0.1,
                 color="black", fill="white")

ggplot(benthic.uni, aes(x = sqJ)) + 
  geom_histogram(binwidth = 0.1,
                 color="black", fill="white")

ggplot(benthic.uni, aes(x = sqrich)) + 
  geom_histogram(binwidth = 0.2,
                 color="black", fill="white")

# Univariate Analyses -----------------------------------------------------

# ANCOVA with habitat and water depth, richness prob poisson distribution

richness <- glm(rich ~ Habitat * Depth, data = benthic.uni, 
                x = TRUE, family = poisson)

summary(richness)
anova(richness)

#Analysis of Deviance Table

#Model: poisson, link: log

#Response: rich

#Terms added sequentially (first to last)


#              Df Deviance Resid. Df Resid. Dev
#NULL                             24    17.2028
#Habitat        2   3.3466        22    13.8562
#Depth          1   0.1013        21    13.7549
#Habitat:Depth  2   3.8007        19     9.9542

plot(richness)

check_model(richness)




richness2 <- lm(sqrich ~ Habitat * Depth, data = benthic.uni)
summary(richness2)
anova(richness2)


#Anova Table (Type II tests)

#Response: sqrich
#               Sum Sq Df F value  Pr(>F)  
#Habitat       0.67417  2  2.5380 0.10547  
#Depth         0.00843  1  0.0635 0.80378  
#Habitat:Depth 0.99305  2  3.7384 0.04275 *
#Residuals     2.52351 19  

plot(richness2)

# this looks way better than GLM

check_model(richness2)
check_heteroscedasticity(richness2) #p = 0.675
check_normality(richness2) # p = 0.704

# Pairwise comparisons
pwc.rich <- benthic.uni %>% 
  emmeans_test(
    rich ~ Habitat, covariate = Depth,
    p.adjust.method = "bonferroni"
  )

#   term          .y.   group1  group2       df statistic      p p.adj p.adj.signif
#1 Depth*Habitat rich  Invaded Treated      21    -0.496 0.625  1     ns          
#2 Depth*Habitat rich  Invaded Uninvaded    21    -2.04  0.0536 0.161 ns          
#3 Depth*Habitat rich  Treated Uninvaded    21    -1.80  0.0867 0.260 ns 

get_emmeans(pwc.rich)

#   Depth Habitat   emmean    se    df conf.low conf.high method     
#1  31.7 Invaded     14.8  1.25    21     12.2      17.4 Emmeans test
#2  31.7 Treated     15.7  1.10    21     13.4      18.0 Emmeans test
#3  31.7 Uninvaded   18.7  1.26    21     16.1      21.3 Emmeans test

## Shannon Weiner

shannon.lm <- lm(H ~ Habitat * Depth, data = benthic.uni)
summary(shannon.lm)
anova(shannon.lm)

#Analysis of Variance Table

#Response: H
#               Df  Sum Sq Mean Sq F value   Pr(>F)   
#Habitat        2 0.98712 0.49356  6.0367 0.009342 **
#Depth          1 0.00074 0.00074  0.0091 0.924956   
#Habitat:Depth  2 0.16442 0.08221  1.0055 0.384504   
#Residuals     19 1.55343 0.08176  

shannon <- HSD.test(shannon.lm, "Habitat", group = TRUE)
plot(shannon)

# H groups
#Uninvaded 1.909067      a
#Treated   1.766169     ab
#Invaded   1.426528      b


# Residuals
plot(shannon.lm)




# Simpsons

simp.lm <- lm(D1 ~ Habitat * Depth, data = benthic.uni)
summary(simp.lm)
anova(simp.lm)

#Analysis of Variance Table

#Response: D1
#              Df   Sum Sq  Mean Sq F value   Pr(>F)   
#Habitat        2 0.143844 0.071922  6.2542 0.008187 **
#Depth          1 0.000028 0.000028  0.0024 0.961364   
#Habitat:Depth  2 0.009890 0.004945  0.4300 0.656688   
#Residuals     19 0.218496 0.011500  

simp <- HSD.test(simp.lm, "Habitat", group = TRUE)
plot(simp)

# D1 groups
#Uninvaded 0.7799573      a
#Treated   0.7555122      a
#Invaded   0.6058429      b

# Residuals
plot(simp.lm)

# Pielous J

pie.lm <- lm(J ~ Habitat * Depth, data = benthic.uni)
summary(pie.lm)
anova(pie.lm)

#Analysis of Variance Table

#Response: J
#              Df   Sum Sq  Mean Sq F value  Pr(>F)  
#Habitat        2 0.086615 0.043308  3.3800 0.05549 .
#Depth          1 0.000067 0.000067  0.0052 0.94331  
#Habitat:Depth  2 0.017942 0.008971  0.7002 0.50887  
#Residuals     19 0.243448 0.012813 

plot(pie.lm)


pie.lm.t <- lm(sqJ ~ Habitat * Depth, data = benthic.uni)
summary(pie.lm.t)
anova(pie.lm.t)

plot(pie.lm.t)


# Univariate Figures ------------------------------------------------------

# ANCOVA plots ------------------------------------------------------------
unique(benthic.uni$Habitat)

rich <- ggplot(benthic.uni, aes(x = Depth, y = rich, 
                           group = Habitat)) +
  geom_point(data = benthic.uni,
             aes(fill = Habitat, shape = Habitat),
             size = 5,
             stroke = 1.5) +
  geom_smooth(method = lm,
              stat = "smooth",
              level = 0.95,
              colour = "black") +
  theme_classic() +
  labs(x = "Water Depth (cm) ",
       y = "Taxonomic Richness") +
  scale_fill_viridis(discrete = TRUE) +
  theme(panel.border = element_rect(fill = NA)) +
  scale_shape_manual(values = c(21, 24, 22)) +
  scale_x_continuous(breaks = c(0, 10, 20, 30, 40, 50, 60, 70, 80)) 
  #theme(legend.position = "blank")


div <- ggplot(benthic.uni, aes(x = Depth, y = H, 
                               group = Habitat)) +
  geom_point(data = benthic.uni,
             aes(fill = Habitat, shape = Habitat),
             size = 5,
             stroke = 1.5) +
  geom_smooth(method = lm,
              stat = "smooth",
              level = 0.95,
              colour = "black") +
  theme_classic() +
  labs(x = "Water Depth (cm) ",
       y = "Shannon-Weiner (H')") +
  scale_fill_viridis(discrete = TRUE) +
  theme(panel.border = element_rect(fill = NA)) +
  scale_shape_manual(values = c(21, 24, 22)) +
  scale_x_continuous(breaks = c(0, 10, 20, 30, 40, 50, 60, 70, 80)) +
  ylim(0, 3)

domin <- ggplot(benthic.uni, aes(x = Depth, y = D1, 
                                 group = Habitat)) +
  geom_point(data = benthic.uni,
             aes(fill = Habitat, shape = Habitat),
             size = 5,
             stroke = 1.5) +
  geom_smooth(method = lm,
              stat = "smooth",
              level = 0.95,
              colour = "black") +
  theme_classic() +
  labs(x = "Water Depth (cm) ",
       y = "Simpson's Diversity (1-D)") +
  scale_fill_viridis(discrete = TRUE) +
  theme(panel.border = element_rect(fill = NA)) +
  scale_shape_manual(values = c(21, 24, 22)) +
  scale_x_continuous(breaks = c(0, 10, 20, 30, 40, 50, 60, 70, 80)) +
  ylim(0,1)


even <- ggplot(benthic.uni, aes(x = Depth, y = J, 
                                group = Habitat)) +
  geom_point(data = benthic.uni,
             aes(fill = Habitat, shape = Habitat),
             size = 5,
             stroke = 1.5) +
  geom_smooth(method = lm,
              stat = "smooth",
              level = 0.95,
              colour = "black") +
  theme_classic() +
  labs(x = "Water Depth (cm) ",
       y = "Pielou's Evenness (J)") +
  scale_fill_viridis(discrete = TRUE) +
  theme(panel.border = element_rect(fill = NA)) +
  scale_shape_manual(values = c(21, 24, 22)) +
  scale_x_continuous(breaks = c(0, 10, 20, 30, 40, 50, 60, 70, 80)) +
  ylim(0,1)


uni.panel <- ggarrange(rich, div,
          domin, even,
          common.legend = TRUE,
          legend = "right",
          labels = "AUTO")

ggsave("Figures/aquatic_invert_ANCOVApanels.jpeg", uni.panel)

# ANOVA plots -------------------------------------------------------------

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
  scale_colour_viridis(discrete = TRUE) +
  xlim(0, 11000)


richness <- ggplot(data = benthic.uni, 
                   aes(x = rich, group = Habitat, fill = Habitat)) +
  geom_density(adjust = 1.5, alpha = 0.6) +
  xlim(0, 42) +
  theme_classic() +
  xlab("Taxa Richness") +
  ylab("Density") +
  theme(legend.position = c(0.8, 0.8)) +
  geom_vline(aes(xintercept = mean.s, colour = Habitat),
             size = 1,
             show.legend = FALSE) +
  scale_fill_viridis(discrete = TRUE) +
  scale_colour_viridis(discrete = TRUE)

ggsave("Figures/aquatic_invert_richness_density.jpeg", richness)

ggarrange(abundance, richness,
          labels = "AUTO",
          hjust = c(-5.5),
          vjust = 2)

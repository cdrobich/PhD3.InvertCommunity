
# Packages ----------------------------------------------------------------

library(tidyverse) # for figures, pipes etc.
library(vegan) # for NMDS

library(ggpubr) # for arranging figures


library(agricolae) # Tukeys post hoc
library(Hmisc)
library(car)

library(lme4) # GLMM
library(MuMIn)
library(AICcmodavg)

# Data Import -------------------------------------------------------------

invert <- read.csv("Data/emerging_invertebrates.csv")
str(invert)

colnames(invert)

invert$Year <- as.factor(invert$Year)
invert <- invert %>% unite("TrtYr", Treatment,Year, remove = FALSE)

invert.data <- invert %>% select(Araneae:Crambidae)
invert.env <- invert %>% select(Site:N)

unique(invert$Treatment)

# Diversity Metrics -------------------------------------------------------

richness <- rowSums(invert.data > 0) # species richness
abundance <- rowSums(invert.data) # abundance
H <- diversity(invert.data) # Shannon Weiner
D1 <- diversity(invert.data, index = "simpson") #default is base log, but can change it
J <- H/log(specnumber(invert.data)) #Pielous


# make new little data frame
invert.uni <- invert.env
invert.uni$rich <- richness
invert.uni$abundance <- abundance
invert.uni$H <- H
invert.uni$D1 <- D1
invert.uni$J <- J

invert.univariate <- invert.uni

write.csv(invert.univariate, "Data/emerging_invertebrate_univariate.csv")


# Histograms ---------------------------------------------------------------
invert <- read.csv("Data/emerging_invertebrate_univariate.csv")
invert$Year <- as.factor(invert$Year)
invert$logAb <- log(invert$abundance)

invert %>% group_by(Treatment, Year) %>% 
  summarise(ab.mean = mean(abundance),
            ab.N = length(abundance),
            std.ab = sd(abundance),
            str.ab = (std.ab)/(sqrt(ab.N)),
            S.mean = mean(rich),
            S.N = length(rich),
            std.S = sd(rich),
            str.S = (std.S)/(sqrt(S.N)))

#Treatment Year  ab.mean  ab.N std.ab str.ab S.mean   S.N std.S str.S
#1 Invaded   2017     317.     9   196.   65.2   22.3     9  5.57  1.86
#2 Invaded   2018     449      9   202.   67.2   27.1     9  6.43  2.14
#3 Treated   2017     932.     9   860.  287.    13.2     9  4.12  1.37
#4 Treated   2018    2580.     9  1252.  417.    17.1     9  4.14  1.38
#5 Uninvaded 2017     268.     9   202.   67.2   16.8     9  5.09  1.70
#6 Uninvaded 2018     761.     9   741.  247.    26.4     9  7.67  2.56

invert %>% group_by(Treatment) %>% 
  summarise(ab.mean = mean(abundance),
            ab.N = length(abundance),
            std.ab = sd(abundance),
            str.ab = (std.ab)/(sqrt(ab.N)),
            S.mean = mean(rich),
            S.N = length(rich),
            std.S = sd(rich),
            str.S = (std.S)/(sqrt(S.N)))

# Treatment ab.mean  ab.N std.ab str.ab S.mean   S.N std.S str.S
#1 Invaded      383.    18   204.   48.2   24.7    18  6.33  1.49
#2 Treated     1756.    18  1343.  317.    15.2    18  4.48  1.05
#3 Uninvaded    515.    18   585.  138.    21.6    18  8.04  1.89



# abundance histogram
ggplot(invert, aes(x = abundance)) + 
  geom_histogram(binwidth = 100,
                 color="black", fill="white")


#richness histogram
ggplot(invert, aes(x = rich)) + 
  geom_histogram(binwidth = 1,
                 color="black", fill="white")




# abundance histogram
ggplot(invert, aes(x = logAb)) + 
  geom_histogram(binwidth = 1,
                 color="black", fill="white")

# ANOVAs ------------------------------------------------------------------


# abundance two-way ANOVA

abundance.lm <- lm(abundance ~ Treatment * Year, data = invert)
Anova(abundance.lm, type = 3)

#Response: abundance
#                  Sum Sq Df F value   Pr(>F)   
#(Intercept)      901867  1  1.8186 0.183813   
#Treatment       2466358  2  2.4866 0.093849 . 
#Year              78937  1  0.1592 0.691690   
#Treatment:Year  5642917  2  5.6893 0.006063 **
#Residuals      23804326 48



# log Abundance

labundance.lm <- lm(logAb ~ Treatment * Year, data = invert)
Anova(labundance.lm, type = 3)

#Response: logAb
#                Sum Sq Df  F value    Pr(>F)    
#(Intercept)    282.780  1 515.3912 < 2.2e-16 ***
#Treatment        6.784  2   6.1819  0.004085 ** 
#Year             0.799  1   1.4556  0.233543    
#Treatment:Year   1.323  2   1.2056  0.308423    
#Residuals       26.336 48  

labundance.hsd <- HSD.test(labundance.lm, "Treatment")

#logAb groups
#Treated   7.108355      a
#Invaded   5.815990      b
#Uninvaded 5.775302      b

# richness two-way ANOVA

rich.lm <- lm(rich ~ Treatment * Year, data = invert)
Anova(rich.lm, type = 3)

#Response: rich
#Sum Sq Df  F value    Pr(>F)    
#(Intercept)    4489.0  1 140.9132 6.907e-16 ***
#Treatment       379.6  2   5.9573  0.004887 ** 
#Year            102.7  1   3.2245  0.078839 .  
#Treatment:Year   87.1  2   1.3672  0.264552    
#Residuals      1529.1 48 

rich.hsd <- HSD.test(rich.lm, "Treatment")

#rich groups
#Invaded   24.72222      a
#Uninvaded 21.61111      a
#Treated   15.16667      b



# General Linear Mixed Models ---------------------------------------------


# Abundance GLMM ----------------------------------------------------------

# Abundance, with year and number of collections as random

abundance.glmm <- lmer(abundance ~ Treatment + (1|Year) + (1|N),
                    data = invert, REML = FALSE)

summary(abundance.glmm)

#Linear mixed model fit by maximum likelihood  ['lmerMod']
#Formula: abundance ~ Treatment + (1 | Year) + (1 | N)
#Data: invert

#AIC      BIC   logLik deviance df.resid 
#869.0    880.8   -428.5    857.0       47 

#Scaled residuals: 
#  Min      1Q  Median      3Q     Max 
#-2.1973 -0.6076  0.0231  0.2638  3.9185 

#Random effects:
#Groups   Name        Variance Std.Dev.
#N        (Intercept)      0     0.0   
#Year     (Intercept) 128827   358.9   
#Residual             573135   757.1   
#Number of obs: 53, groups:  N, 5; Year, 2

#Fixed effects:
#                    Estimate Std. Error t value
#(Intercept)           406.4      313.3   1.297
#TreatmentTreated     1349.9      256.1   5.271
#TreatmentUninvaded    108.5      256.1   0.424

#Correlation of Fixed Effects:
#             (Intr) TrtmnT
#TretmntTrtd -0.421       
#TrtmntUnnvd -0.421  0.515

r.squaredGLMM(abundance.glmm)

#      R2m       R2c
#0.3540902 0.4726298

# examine the residuals 

plot(abundance.glmm)
qqnorm(resid(abundance.glmm))
qqline(resid(abundance.glmm))


# null model 

ab.null <- lmer(abundance ~ (1|Year) + (1|N), 
                data = invert, REML = FALSE) 

summary(ab.null)

#Linear mixed model fit by maximum likelihood  ['lmerMod']
#Formula: abundance ~ (1 | Year) + (1 | N)
#Data: invert

#AIC      BIC   logLik deviance df.resid 
#891.6    899.5   -441.8    883.6       49 

#Scaled residuals: 

#  Min      1Q  Median      3Q     Max 
#-1.0976 -0.5546 -0.3346  0.2112  3.9187 

#Random effects:

#Groups   Name        Variance Std.Dev.
#N        (Intercept)      0     0.0   
#Year     (Intercept) 121398   348.4   
#Residual             964325   982.0   
#Number of obs: 53, groups:  N, 5; Year, 2

#Fixed effects:
#             Estimate Std. Error t value
#(Intercept)    901.2      280.9   3.208

r.squaredGLMM(ab.null)

#R2m       R2c
# 0   0.1118127

## Comparison of abundance model and null model

models <- list(abundance.glmm, ab.null)
names <- c(1:2)

aictab(cand.set = models, modnames = names,
       sort = TRUE, c.hat = 1, second.ord = TRUE,
       nobs = NULL)

#Model selection based on AICc:
  
# K   AICc Delta_AICc AICcWt Cum.Wt      LL
# 6 870.83       0.00      1      1 -428.50
# 4 892.47      21.64      0      1 -441.82

# or do it this way and get info re: test statistic

anova(abundance.glmm, ab.null)

#               npar    AIC    BIC  logLik deviance  Chisq Df Pr(>Chisq)    
#ab.null           4 891.64 899.52 -441.82   883.64                         
#abundance.glmm    6 869.01 880.83 -428.50   857.01 26.631  2  1.649e-06 ***


# log Abundance GLMMM


labundance.glmm <- lmer(logAb ~ Treatment + (1|Year) + (1|N),
                       data = invert, REML = FALSE)

summary(labundance.glmm)



#Linear mixed model fit by maximum likelihood  ['lmerMod']
#Formula: logAb ~ Treatment + (1 | Year) + (1 | N)
#Data: invert

#AIC      BIC   logLik deviance df.resid 
#133.8    145.6    -60.9    121.8       47 

#Scaled residuals: 
#  Min      1Q  Median      3Q     Max 
#-2.3632 -0.6599  0.1089  0.5649  2.1454 

#Random effects:
#  Groups   Name        Variance Std.Dev.
#N        (Intercept) 0.0000   0.0000  
#Year     (Intercept) 0.1563   0.3953  
#Residual             0.5374   0.7331  
#Number of obs: 53, groups:  N, 5; Year, 2

#Fixed effects:
#  Estimate Std. Error t value
#(Intercept)         5.84197    0.33134  17.631
#TreatmentTreated    1.26638    0.24800   5.106
#TreatmentUninvaded -0.06667    0.24800  -0.269

#Correlation of Fixed Effects:
#  (Intr) TrtmnT
#TretmntTrtd -0.385       
#TrtmntUnnvd -0.385  0.515

r.squaredGLMM(labundance.glmm)

#R2m       R2c
#0.3583682 0.5029123

plot(labundance.glmm)
qqnorm(resid(labundance.glmm))
qqline(resid(labundance.glmm))



# Richness GLMM -----------------------------------------------------------

rich.glmm <- lmer(rich ~ Treatment + (1|Year) + (1|N),
                     data = invert, REML = FALSE)

summary(rich.glmm)

#AIC      BIC   logLik deviance df.resid 
#340.8    352.6   -164.4    328.8       47 

#Scaled residuals: 
#  Min      1Q  Median      3Q     Max 
#-2.4030 -0.4817  0.0774  0.6760  2.2532 

#Random effects:
#  Groups   Name        Variance Std.Dev.
#N        (Intercept) 12.97    3.602   
#Year     (Intercept)  0.00    0.000   
#Residual             25.12    5.012   
#Number of obs: 53, groups:  N, 5; Year, 2

#Fixed effects:
#                   Estimate Std. Error t value
#(Intercept)          23.885      2.209  10.812
#TreatmentTreated     -9.432      1.711  -5.513
#TreatmentUninvaded   -3.630      1.696  -2.141

#Correlation of Fixed Effects:
#(Intr) TrtmnT
#TretmntTrtd -0.440       
#TrtmntUnnvd -0.393  0.508

r.squaredGLMM(rich.glmm)

#      R2m       R2c
# 0.2865247 0.5295009

# examine the residuals 

plot(rich.glmm)
qqnorm(resid(rich.glmm))
qqline(resid(rich.glmm))

### null abundance model 

rich.null <- lmer(rich ~ (1|Year) + (1|N), 
                data = invert, REML = FALSE) 

summary(rich.null)

#Formula: rich ~ (1 | Year) + (1 | N)
#Data: invert

#AIC      BIC   logLik deviance df.resid 
#361.0    368.9   -176.5    353.0       49 

#Scaled residuals: 
#  Min       1Q   Median       3Q      Max 
#-2.03558 -0.82443 -0.03973  0.78934  1.88792 

#Random effects:
#Groups   Name        Variance Std.Dev.
#N        (Intercept) 14.48    3.806   
#Year     (Intercept)  0.00    0.000   
#Residual             40.60    6.372   
#Number of obs: 53, groups:  N, 5; Year, 2

#Fixed effects:
#Estimate Std. Error t value
#(Intercept)   19.041      2.146   8.874

r.squaredGLMM(rich.null)

#R2m       R2c
#0     0.2629343


## Comparison 

models <- list(rich.glmm, rich.null)
names <- c(1:2)

aictab(cand.set = models, modnames = names,
       sort = TRUE, c.hat = 1, second.ord = TRUE,
       nobs = NULL)

#K   AICc Delta_AICc AICcWt Cum.Wt      LL
#6  342.63       0.00      1      1 -164.40
#4  361.86      19.23      0      1 -176.51
# or do it this way and get info re: test statistic

anova(rich.glmm, rich.null)

#npar    AIC    BIC  logLik deviance  Chisq Df Pr(>Chisq)    
#rich.null    4 361.02 368.91 -176.51   353.02                         
#rich.glmm    6 340.80 352.62 -164.40   328.80 24.221  2  5.502e-06 ***


# Univariate Figures ------------------------------------------------------

colnames(invert)

# Abundance

abundance <- ggplot(data = invert, 
                    aes(x = Treatment, 
                        y = abundance,
                        shape = Year, colour = Treatment),
                    size = 4) +
  geom_violin(trim = FALSE, 
              lwd = 0.75, 
              position = position_dodge(0.5),
              colour = "black") +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.1, 
                                              dodge.width = 0.5), size = 3) +
  stat_summary(aes(shape = Year), 
               colour = "black", 
               fun.data = "mean_se", 
               fun.args = list(mult = 1), 
               geom = "pointrange", 
               size = 1, position = position_dodge(0.5)) +
  theme_classic() +
  labs(x = " ",
       y = expression(paste("Invertebrate Abundance"))) +
  theme(panel.border = element_rect(fill = NA)) +
  theme(text = element_text(size = 14),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14)) +
  scale_color_manual(values = c("#969696","#35978f", "#2166ac")) +
  scale_shape_manual(values = c(15, 16)) +
  guides(color = "none") 


richness <- ggplot(data = invert.univariate, 
                   aes(x = Treatment, 
                       y = rich, 
                       shape = Year, colour = Treatment), size = 4) +
  geom_violin(trim = FALSE,
              lwd = 0.75, 
              position = position_dodge(0.8),
              colour = "black") +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.1, 
                                              dodge.width = 0.8), size = 3) +
  stat_summary(aes(shape = Year), 
               colour = "black", 
               fun.data = "mean_se", 
               fun.args = list(mult = 1),
               geom = "pointrange", size = 1,
               position = position_dodge(0.8)) +
  theme_classic() +
  labs(x = " ",
       y = expression(paste("Taxonomic Richness"))) +
  theme(panel.border = element_rect(fill = NA)) +
  theme(text = element_text(size = 14),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14)) +
  scale_colour_manual(values = c("#969696","#35978f", "#2166ac")) +
  scale_shape_manual(values = c(15, 16)) +
  guides(color = "none") +
  theme(legend.position = "right") 

richness +
  annotate("text", x = 1:3, y = c(4600, 12000, 4600),
           label = c("a", "b", "a"),
           size = 4.5)


uni.panel <- ggarrange(abundance, richness,
                       common.legend = TRUE,
                       legend = "bottom",
                       labels = "AUTO",
                       hjust = c(-7,-6),
                       vjust = 2)


ggsave("Figures/emerging_invert_unipanel.jpeg", uni.panel,
       height = 5.49,
       width = 9.58,
       units = "in")


### without year

ab.hab <- ggplot(data = invert, 
                 aes(x = Treatment, 
                     y = abundance,
                     colour = Treatment),
                 size = 4) +
  geom_violin(trim = FALSE, lwd = 1) +
  geom_point() +
  stat_summary(aes(shape = Treatment), 
               colour = "black", 
               fun.data = "mean_se", 
               fun.args = list(mult = 1), 
               geom = "pointrange", 
               size = 1) +
  theme_classic() +
  labs(x = " ",
       y = expression(paste("Invertebrate Abundance"))) +
  theme(panel.border = element_rect(fill = NA)) +
  theme(text = element_text(size = 14),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14)) +
  scale_colour_manual(values = c("#969696","#35978f", "#2166ac")) +
  scale_shape_manual(values = c(17, 18, 15)) +
  guides(color = "none") +
  annotate("text", x = 1:3, y = c(4600, 8000, 4600),
           label = c("a", "b", "a"),
           size = 4.5)

ab.hab

rich.hab <- ggplot(data = invert, 
                 aes(x = Treatment, 
                     y = rich,
                     colour = Treatment),
                 size = 4) +
  geom_violin(trim = FALSE, lwd = 1) +
  geom_point() +
  stat_summary(aes(shape = Treatment), 
               colour = "black", 
               fun.data = "mean_se", 
               fun.args = list(mult = 1), 
               geom = "pointrange", 
               size = 1) +
  theme_classic() +
  labs(x = " ",
       y = expression(paste("Invertebrate Abundance"))) +
  theme(panel.border = element_rect(fill = NA)) +
  theme(text = element_text(size = 14),
        axis.text.x = element_text(size = 14),
        axis.text.y = element_text(size = 14)) +
  scale_colour_manual(values = c("#969696","#35978f", "#2166ac")) +
  scale_shape_manual(values = c(17, 18, 15)) +
  guides(color = "none") +
  annotate("text", x = 1:3, y = c(55, 35, 55),
           label = c("a", "b", "a"),
           size = 4.5)

rich.hab



uni.hab.panel <- ggarrange(ab.hab, rich.hab,
                       common.legend = TRUE,
                       legend = "bottom",
                       labels = "AUTO",
                       hjust = c(-7,-6),
                       vjust = 2)


ggsave("Figures/emerging_invert_unipanel.jpeg",
       uni.hab.panel,
       height = 5.49,
       width = 9.58,
       units = "in")

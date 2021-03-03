
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
library(viridis)

library(effects)
library(lmerTest)

library(sjPlot)
# Data Import -------------------------------------------------------------

invert <- read.csv("Data/Emerging/emerging_invertebrates.csv")
str(invert)
unique(invert$Treatment)

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
D1 <- diversity(invert.data, index = "simpson") # 1 - D
J <- H/log(specnumber(invert.data)) #Pielous


# make new little data frame
invert.uni <- invert.env
invert.uni$rich <- richness
invert.uni$abundance <- abundance
invert.uni$H <- H
invert.uni$D1 <- D1
invert.uni$J <- J

invert.univariate <- invert.uni

write.csv(invert.univariate, "Data/Emerging/emerging_invertebrate_univariate.csv")


# Histograms ---------------------------------------------------------------
invert <- read.csv("Data/Emerging/emerging_invertebrate_univariate.csv")

invert$Year <- as.factor(invert$Year)
invert$logAb <- log(invert$abundance)

install.packages("plotrix")
library(plotrix)

invert.hab <- invert %>% group_by(Treatment) %>% 
  summarise(across(
    .cols = where(is.numeric),
    .fns = list(Mean = mean, SD = sd, SE = std.error), na.rm = TRUE,
    .names = "{col}_{fn}"
  ))

invert.hab %>% t %>% as.data.frame

#Treatment         Invaded    Treated  Uninvaded
#Depth_Mean       43.06111   41.36111   33.15556
#Depth_SD         18.53579   15.33282   18.69666
#Depth_SE         4.368927   3.613981   4.406844

#N_Mean           4.823529   4.611111   4.888889
#N_SD             1.074436   1.195033   1.078610
#N_SE            0.2605889  0.2816720  0.2542307

#rich_Mean        24.72222   15.16667   21.61111
#rich_SD          6.332043   4.475423   8.037697
#rich_SE          1.492477   1.054867   1.894503

#abundance_Mean   382.7778  1756.3333   514.8889
#abundance_SD     204.4276  1343.4526   584.5660
#abundance_SE     48.18405  316.65482  137.78353

#H_Mean          1.6472675  0.3933663  1.3781480
#H_SD            0.4860376  0.2538906  0.7731806
#H_SE           0.11456015 0.05984259 0.18224040

#D1_Mean         0.6171321  0.1552774  0.5095815
#D1_SD           0.1673317  0.1329803  0.2809856
#D1_SE          0.03944046 0.03134375 0.06622894

#J_Mean          0.5181020  0.1478616  0.4479636
#J_SD           0.13857407 0.09493142 0.23375547
#J_SE           0.03266222 0.02237555 0.05509669




invert.habyr <- invert %>% group_by(Treatment, Year) %>% 
  summarise(across(
    .cols = where(is.numeric),
    .fns = list(Mean = mean, SD = sd, SE = std.error), na.rm = TRUE,
    .names = "{col}_{fn}"
  ))

invert.habyr %>% t %>% as.data.frame


#Treatment         Invaded    Invaded    Treated    Treated  Uninvaded  Uninvaded
#Year                 2017       2018       2017       2018       2017       2018

#Depth_Mean       42.07778   44.04444   50.66667   32.05556   37.58889   28.72222
#Depth_SD         17.20064   20.78606   13.45753   11.11925   17.18476   20.08201
#Depth_SE         5.733546   6.928687   4.485842   3.706418   5.728255   6.694002

#N_Mean           3.888889   5.875000   3.666667   5.555556   3.888889   5.888889
#N_SD            0.3333333  0.3535534  0.7071068  0.7264832  0.3333333  0.3333333
#N_SE            0.1111111  0.1250000  0.2357023  0.2421611  0.1111111  0.1111111

#rich_Mean        22.33333   27.11111   13.22222   17.11111   16.77778   26.44444
#rich_SD          5.567764   6.431260   4.116363   4.136558   5.093569   7.666667
#rich_SE          1.855921   2.143753   1.372121   1.378853   1.697856   2.555556

#abundance_Mean   316.5556   449.0000   932.2222  2580.4444   268.3333   761.4444
#abundance_SD     195.6624   201.6290   859.8657  1252.1426   201.6191   740.7581
#abundance_SE     65.22080   67.20966  286.62188  417.38086   67.20636  246.91936

#H_Mean          1.6687518  1.6257832  0.4041136  0.3826190  1.2662908  1.4900052
#H_SD            0.4602864  0.5376727  0.2232105  0.2947806  0.7483638  0.8259188
#H_SE           0.15342881 0.17922423 0.07440348 0.09826022 0.24945459 0.27530628

#D1_Mean         0.6450358  0.5892285  0.1465689  0.1639859  0.4888948  0.5302681
#D1_SD          0.13468582 0.19901665 0.09529313 0.16830497 0.29554762 0.28189305
#D1_SE          0.04489527 0.06633888 0.03176438 0.05610166 0.09851587 0.09396435

#J_Mean          0.5421213  0.4940827  0.1599113  0.1358119  0.4477526  0.4481745
#J_SD           0.12535814 0.15425028 0.08821076 0.10508386 0.25459609 0.22648218
#J_SE           0.04178605 0.05141676 0.02940359 0.03502795 0.08486536 0.07549406

#richness histogram
ggplot(invert, aes(x = rich)) + 
  geom_histogram(binwidth = 1,
                 color="black", fill="white")


#Shannon histogram
ggplot(invert, aes(x = H)) + 
  geom_histogram(binwidth = .1,
                 color="black", fill="white")

#Simpsons
ggplot(invert, aes(x = D1)) + 
  geom_histogram(binwidth = .1,
                 color="black", fill="white")

#Pielou
ggplot(invert, aes(x = J)) + 
  geom_histogram(binwidth = .1,
                 color="black", fill="white")

# General Linear Mixed Models (glmer) ---------------------------------------------

invert$depthst <- scale(invert$Depth,
                           center = TRUE,
                           scale = TRUE)


# Richness GLMM -----------------------------------------------------------

rich.glmm <- glmer(rich ~ Treatment * depthst + (1|Year) + (1|N),
                     data = invert, family = poisson)

summary(rich.glmm)


anova(rich.glmm)

#Response: rich
#                     Chisq Df Pr(>Chisq)    
#(Intercept)       969.5694  1  < 2.2e-16 ***
#Treatment          39.6564  2  2.448e-09 ***
#depthst             3.4350  1    0.06383 .  
#Treatment:depthst   1.5353  2    0.46410  



anova(rich.glmm)

#  Analysis of Variance Table
#                npar Sum Sq Mean Sq F value
#Treatment          2 37.953 18.9765 18.9765
#Depth              1  3.867  3.8673  3.8673
#Treatment:Depth    2  1.570  0.7849  0.7849

r.squaredGLMM(rich.glmm)

#            R2m       R2c
#delta     0.3594343 0.5813962
#lognormal 0.3633659 0.5877556
#trigamma  0.3553815 0.5748407

# examine the residuals 

plot(rich.glmm)
qqnorm(resid(rich.glmm))
qqline(resid(rich.glmm))

tab_model(rich.glmm,
          string.ci = "95% CI",
          string.pred = "Coefficients",
          string.p = "P value",
          dv.labels = c("Taxa Richness"))

plot(allEffects(rich.glmm))


# Linear models -----------------------------------------------------------

## Output table 
# https://link.springer.com/article/10.3758/s13428-016-0809-y

rich.glmm <- glmer(rich ~ Treatment * depthst + (1|Year) + (1|N),
                   data = invert, family = poisson)

summary(rich.glmm)







# Linear model for this one (non-integer)

shannon.lmm <- lmer(H ~ Treatment * depthst + (1|Year) + (1|N),
                   data = invert, REML = TRUE)

summary(shannon.lmm)

anova(shannon.lmm)

#Type III Analysis of Variance Table with Satterthwaite's method
#                   Sum Sq Mean Sq NumDF DenDF F value    Pr(>F)    
#Treatment         14.9768  7.4884     2    47 23.6600 7.782e-08 ***
#depthst            0.0190  0.0190     1    47  0.0599    0.8077    
#Treatment:depthst  0.2158  0.1079     2    47  0.3409    0.7129 

library(lmerTest)
anova(shannon.lmm)

#Type III Analysis of Variance Table with Satterthwaite's method
#                   Sum Sq Mean Sq NumDF DenDF F value    Pr(>F)    
#Treatment         14.9768  7.4884     2    47 23.6600 7.782e-08 ***
#depthst            0.0190  0.0190     1    47  0.0599    0.8077    
#Treatment:depthst  0.2158  0.1079     2    47  0.3409    0.7129  


r.squaredGLMM(shannon.lmm)

#            R2m       R2c
#      0.4821051 0.4821051

# examine the residuals 

plot(shannon.lmm)
qqnorm(resid(shannon.lmm))
qqline(resid(shannon.lmm))

plot(allEffects(shannon.lmm))


### Linear model for Diversity


simp.lmm <- lmer(D1 ~ Treatment * depthst + (1|Year) + (1|N),
                    data = invert, REML = TRUE)

summary(simp.lmm)


anova(simp.lmm)
#Type III Analysis of Variance Table with Satterthwaite's method
#                   Sum Sq Mean Sq NumDF DenDF F value    Pr(>F)    
#Treatment         1.99192 0.99596     2    47 22.5335 1.374e-07 ***
#depthst           0.00005 0.00005     1    47  0.0012    0.9730    
#Treatment:depthst 0.01207 0.00604     2    47  0.1365    0.8727

r.squaredGLMM(simp.lmm)

#           R2m       R2c
#    0.4671788 0.4671788

plot(simp.lmm)
qqnorm(resid(simp.lmm))
qqline(resid(simp.lmm))

plot(allEffects(simp.lmm))



## Linear model for Pielous 

pie.lmm <- lmer(J ~ Treatment * depthst + (1|Year) + (1|N),
                 data = invert, REML = TRUE)

summary(pie.lmm)
anova(pie.lmm)

#Type III Analysis of Variance Table with Satterthwaite's method
#                   Sum Sq Mean Sq NumDF DenDF F value    Pr(>F)    
#Treatment         1.33479 0.66739     2    47 22.9776 1.096e-07 ***
#depthst           0.00834 0.00834     1    47  0.2871    0.5946    
#Treatment:depthst 0.01340 0.00670     2    47  0.2306    0.7949    


r.squaredGLMM(pie.lmm)

#           R2m       R2c
#    0.4710541 0.4710541

plot(pie.lmm)
qqnorm(resid(pie.lmm))
qqline(resid(pie.lmm))

plot(allEffects(pie.lmm))


tab_model(shannon.lmm, simp.lmm, pie.lmm,
          string.ci = "95% CI",
          string.pred = "Coefficients",
          string.p = "P value",
          dv.labels = c("Shannon-Weiner (H')",
                        "Simpson's Diversity (1-D)",
                        "Pielou's Evenness (J)"),
          file = "Data/Emerging/lmm_output.html")













# Null models -------------------------------------------------------------
 

rich.null <- glmer(rich ~ (1|Year) + (1|N),
                   family = poisson, data = invert) 

summary(rich.null)

#Generalized linear mixed model fit by maximum likelihood (Laplace
#Approximation) [glmerMod]
#Family: poisson  ( log )
#Formula: rich ~ (1 | Year) + (1 | N)
#Data: invert

#AIC      BIC   logLik deviance df.resid 
#372.9    378.8   -183.4    366.9       50 

#Scaled residuals: 
#  Min       1Q   Median       3Q      Max 
#-2.65232 -1.13678  0.08907  1.05232  2.72149 

#Random effects:
#Groups Name        Variance  Std.Dev. 
#N      (Intercept) 3.803e-02 1.950e-01
#Year   (Intercept) 2.037e-09 4.513e-05
#Number of obs: 53, groups:  N, 5; Year, 2

#Fixed effects:
#  Estimate Std. Error z value Pr(>|z|)    
#(Intercept)   2.9085     0.1038   28.03   <2e-16 ***

r.squaredGLMM(rich.null)

#          R2m       R2c
#delta       0 0.4153653
#lognormal   0 0.4217383
#trigamma    0 0.4088534

## Comparison 

models <- list(rich.glmm, rich.null)
names <- c(1:2)

aictab(cand.set = models, modnames = names,
       sort = TRUE, c.hat = 1, second.ord = TRUE,
       nobs = NULL)

#Model selection based on AICc:
  
#  K   AICc Delta_AICc AICcWt Cum.Wt      LL
#1 8 341.62       0.00      1      1 -161.17
#2 3 373.39      31.76      0      1 -183.45

anova(rich.glmm, rich.null)

#Data: invert
#Models:
#  rich.null: rich ~ (1 | Year) + (1 | N)
#rich.glmm: rich ~ Treatment * Depth + (1 | Year) + (1 | N)

#          npar    AIC    BIC  logLik deviance  Chisq Df Pr(>Chisq)    
#rich.null    3 372.90 378.81 -183.45   366.90                         
#rich.glmm    8 338.35 354.11 -161.17   322.35 44.547  5  1.793e-08 ***



# Univariate Figures ------------------------------------------------------

rich <- ggplot(invert, aes(x = Depth, y = rich, 
                                group = Treatment)) +
  geom_point(data = invert,
             aes(fill = Treatment, shape = Treatment),
             size = 5,
             stroke = 1.5) +
  geom_smooth(method = lm,
              stat = "smooth",
              level = 0.95,
              colour = "black") +
  theme_classic(14) +
  labs(x = "Water Depth (cm) ",
       y = "Taxonomic Richness") +
  scale_fill_viridis(discrete = TRUE) +
  scale_shape_manual(values = c(21, 24, 22)) +
  theme(panel.border = element_rect(fill = NA)) +
  scale_x_continuous(breaks = c(0, 10, 20, 30, 40, 50, 60, 70, 80)) +
  ylim(0, 50)




shannon <- ggplot(invert, aes(x = Depth, y = H, 
                              group = Treatment)) +
  geom_point(data = invert,
             aes(fill = Treatment, shape = Treatment),
             size = 5,
             stroke = 1.5) +
  geom_smooth(method = lm,
              stat = "smooth",
              level = 0.95,
              colour = "black") +
  theme_classic(14) +
  labs(x = "Water Depth (cm) ",
       y = "Shannon-Weiner (H')") +
  scale_fill_viridis(discrete = TRUE) +
  scale_shape_manual(values = c(21, 24, 22)) +
  theme(panel.border = element_rect(fill = NA)) +
  scale_x_continuous(breaks = c(0, 10, 20, 30, 40, 50, 60, 70, 80)) +
  ylim(0, 3)



simpson <- ggplot(invert, aes(x = Depth, y = D1, 
                              group = Treatment)) +
  geom_point(data = invert,
             aes(fill = Treatment, shape = Treatment),
             size = 5,
             stroke = 1.5) +
  geom_smooth(method = lm,
              stat = "smooth",
              level = 0.95,
              colour = "black") +
  theme_classic(14) +
  labs(x = "Water Depth (cm) ",
       y = "Simpson's Diversity (1-D)") +
  scale_fill_viridis(discrete = TRUE) +
  scale_shape_manual(values = c(21, 24, 22)) +
  theme(panel.border = element_rect(fill = NA)) +
  scale_x_continuous(breaks = c(0, 10, 20, 30, 40, 50, 60, 70, 80)) +
  ylim(0, 1)



pielou <- ggplot(invert, aes(x = Depth, y = J, 
                             group = Treatment)) +
  geom_point(data = invert,
             aes(fill = Treatment, shape = Treatment),
             size = 5,
             stroke = 1.5) +
  geom_smooth(method = lm,
              stat = "smooth",
              level = 0.95,
              colour = "black") +
  theme_classic(14) +
  labs(x = "Water Depth (cm) ",
       y = "Pielou's Evenness (J)") +
  scale_fill_viridis(discrete = TRUE) +
  scale_shape_manual(values = c(21, 24, 22)) +
  theme(panel.border = element_rect(fill = NA)) +
  scale_x_continuous(breaks = c(0, 10, 20, 30, 40, 50, 60, 70, 80)) +
  ylim(0, 1)


panel <- ggarrange(rich, shannon,
          simpson, pielou,
          common.legend = TRUE,
          legend = "right",
          align = "hv",
          labels = "AUTO")

ggsave("Figures/invert_hab_depth.jpeg", panel)

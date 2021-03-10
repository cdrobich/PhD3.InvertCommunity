
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

library(performance)
library(see)

citation("performance")

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
invert$Factor <- as.factor(invert$Factor)
invert$logAb <- log(invert$abundance + 1)
invert$sqRich <- sqrt(invert$rich)

test.year <- lm(rich ~ Treatment * Year, data = invert)
anova(test.year)

#Response: rich
#                Df  Sum Sq Mean Sq F value    Pr(>F)    
#Treatment       2  855.11  427.56 13.4213 2.345e-05 ***
#Year            1  504.17  504.17 15.8262 0.0002336 ***
#Treatment:Year  2   87.11   43.56  1.3672 0.2645521    
#Residuals      48 1529.11   31.86 
                

year.hsd <- HSD.test(test.year, "Year")

#rich groups
#2018 23.55556      a
#2017 17.44444      b

trt.hsd <- HSD.test(test.year, "Treatment")

#rich groups
#Invaded   24.72222      a
#Uninvaded 21.61111      a
#Treated   15.16667      b

invert.2017 <- invert %>% filter(Year == "2017")

ab.water <- lm(logAb ~ Treatment * Depth, data = invert.2017)
Anova(ab.water)

#Response: logAb
#                 Sum Sq Df F value   Pr(>F)   
#Treatment        7.9593  2  7.2127 0.004126 **
#Depth            1.2544  1  2.2735 0.146498   
#Treatment:Depth  0.3602  2  0.3264 0.725086   
#Residuals       11.5869 21

s.water <- lm(sqRich ~ Treatment * Depth, data = invert.2017)
Anova(s.water)

#Response: sqRich
#                 Sum Sq Df F value   Pr(>F)   
#Treatment       4.3634  2  6.3206 0.007098 **
#Depth           1.3329  1  3.8615 0.062773 . 
#Treatment:Depth 0.6673  2  0.9666 0.396662   
#Residuals       7.2488 21 

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


invert.2017 %>% group_by(TrtYr) %>% 
  summarise(depth.min = min(Depth),
            depth.max = max(Depth),
            range = depth.max - depth.min)

#TrtYr          depth.min depth.max range
#1 Invaded_2017        22.3        75  52.7
#2 Treated_2017        32.2        75  42.8
#3 Uninvaded_2017      17          75  58  


water.depth <- ggplot(invert.2017, aes(y = Depth, x = Treatment)) +
  geom_violin(aes(fill = Treatment),
              trim = FALSE,
              size = 1,
              alpha = 0.8) +
  scale_fill_viridis(discrete = TRUE) +
  theme_classic(18) +
  labs(x = " ",
       y = "Water Depth (cm)") +
  theme(panel.border = element_rect(fill = NA)) +
  theme(legend.position = "none") +
  theme(axis.text = element_text(size=18))

ggsave("Figures/water_depth_2017.jpeg",
       water.depth)



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
                 color="black", fill = "white")

#Pielou
ggplot(invert, aes(x = J)) + 
  geom_histogram(binwidth = .1,
                 color="black", fill="white")

# General Linear Mixed Models (glmer) ---------------------------------------------

invert$depthst <- scale(invert$Depth,
                           center = TRUE,
                           scale = TRUE)


# Richness GLMM -----------------------------------------------------------

# Factor; a = Uninvaded, b = invaded, c = treated

rich.glmm <- glmer(rich ~ Factor  + (1|Year) + (1|N),
                     data = invert, family = poisson(link = "log"))

summary(rich.glmm)

#Family: poisson  ( log )
#Formula: rich ~ Factor + (1 | Year) + (1 | N)
#Data: invert
#
#AIC      BIC   logLik deviance df.resid 
#337.9    347.7   -163.9    327.9       48 
#
#Scaled residuals: 
#  Min       1Q   Median       3Q      Max 
#-2.46402 -0.61819  0.01305  0.77866  2.17309 
#
#Random effects:
#  Groups Name        Variance  Std.Dev. 
#N      (Intercept) 3.279e-02 1.811e-01
#Year   (Intercept) 1.981e-10 1.407e-05
#Number of obs: 53, groups:  N, 5; Year, 2
#
#Fixed effects:
#              Estimate Std. Error z value Pr(>|z|)    
#  (Intercept)  2.98531    0.10682  27.946  < 2e-16 ***
#  Factorb      0.15751    0.07016   2.245   0.0248 *  
#  Factorc     -0.32178    0.07984  -4.031 5.57e-05 ***
#  ---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#Correlation of Fixed Effects:
#  (Intr) Factrb
#Factorb -0.347       
#Factorc -0.361  0.461


anova(rich.glmm)

#npar Sum Sq Mean Sq F value
#Factor    2 37.533  18.767  18.767

?r2_nakagawa
r2_nakagawa(rich.glmm)

# Conditional R2: 0.580
# Marginal R2: 0.320

# examine the residuals 

plot(rich.glmm)
qqnorm(resid(rich.glmm))
qqline(resid(rich.glmm))

# Default if ""Incidence Rate Ratios"" for Poisson (estimates)
tab_model(rich.glmm,
          string.ci = "95% CI",
          string.pred = "Coefficients",
          string.p = "P value",
          dv.labels = c("Taxa Richness"))

          file = "Data/Emerging/glmer_output.html")

plot(allEffects(rich.glmm))


plot_model(
  rich.glmm, 
  type = "pred", 
  terms = c("Factor"), 
  colors = "bw",
  ci.lvl = NA
)

check_overdispersion(rich.glmm) # no overdispersion detected
check_singularity(rich.glmm) # FALSE
check_model(rich.glmm)



rich.null <- glmer(rich ~ 1 + (1|Year) + (1|N),
                   data = invert, family = poisson(link = "log"))

summary(rich.null)

anova(rich.glmm,rich.null)

compare_performance(rich.glmm,rich.null)

plot(compare_performance(rich.glmm,rich.null))

### glmmTMB

library(glmmTMB)

rich.glm <- glmmTMB(rich ~ Treatment * depthst + (1|Year) + (1|N),
                    data = invert, family = poisson,
                    REML = TRUE)
summary(rich.glm)


# Abundance GLMM ----------------------------------------------------------
 # possibly better to use logAb w/ gaussian because of overdispersion of year

# Factor; a = Uninvaded, b = invaded, c = treated

abun.glmm <- glmer(abundance ~ Factor + (1|Year) + (1|N),
                   data = invert, family = Poisson)

summary(abun.glmm)

#Generalized linear mixed model fit by maximum likelihood (Laplace Approximation) ['glmerMod']
#Family: poisson  ( log )
#Formula: abundance ~ Factor + (1 | Year) + (1 | N)
#Data: invert
#
#AIC      BIC   logLik deviance df.resid 
#18068.4  18078.2  -9029.2  18058.4       48 
#
#Scaled residuals: 
#  Min      1Q  Median      3Q     Max 
#-32.248 -12.097  -2.631   5.976  61.875 
#
#Random effects:
#  Groups Name        Variance Std.Dev.
#N      (Intercept) 0.09364  0.3060  
#Year   (Intercept) 0.13718  0.3704  
#Number of obs: 53, groups:  N, 5; Year, 2
#
#Fixed effects:
#  Estimate Std. Error z value Pr(>|z|)    
#(Intercept)  6.04145    0.29582   20.42   <2e-16 ***
#  Factorb     -0.25437    0.01613  -15.78   <2e-16 ***
#  Factorc      1.26682    0.01204  105.23   <2e-16 ***


r2_nakagawa(abun.glmm)

# Conditional R2: 0.998
# Marginal R2: 0.661


# Default if ""Incidence Rate Ratios"" for Poisson (estimates)
tab_model(abun.glmm, rich.glmm,
          string.ci = "95% CI",
          string.pred = "Coefficients",
          string.p = "P value",
          dv.labels = c("Abundance", 
                        "Taxa Richness"),
          file = "Data/Emerging/glmer_output.html")

plot(allEffects(abun.glmm))


plot_model(
  abun.glmm, 
  type = "pred", 
  terms = c("Factor"), 
  colors = "bw",
  ci.lvl = NA
)


check_overdispersion(abun.glmm) # overdispersion detected
check_singularity(abun.glmm) # FALSE
check_model(abun.glmm)


# Linear models -----------------------------------------------------------

## Linear model for Abundance

ab.lmm <- lmer(logAb ~ Factor + (1|Year) + (1|N),
                data = invert, REML = TRUE)

summary(ab.lmm)
anova(ab.lmm)

#Type III Analysis of Variance Table with Satterthwaite's method
#       Sum Sq Mean Sq NumDF  DenDF F value    Pr(>F)    
#Factor 20.043  10.021     2 49.002  18.015 1.365e-06 ***


r2_nakagawa(ab.lmm)

#Conditional R2: NA
#Marginal R2: 0.409

plot(ab.lmm)
qqnorm(resid(ab.lmm))
qqline(resid(ab.lmm))

plot(allEffects(ab.lmm))

check_model(ab.lmm)
check_singularity(ab.lmm) # TRUE

# Just ANOVA

ab.lm <- lm(logAb ~ Factor, data = invert)
Anova(ab.lm)

#Response: logAb
#         Sum Sq Df F value    Pr(>F)    
#Factor    20.608  2  14.373 1.121e-05 ***
#Residuals 36.564 51

ab.hsd <- HSD.test(ab.lm, "Factor")

#logAb groups
#treated 7.109676      a
#invaded  5.819395      b
#uninvaded 5.779883      b

check_model(ab.lm)
check_normality(ab.lm) # good
check_heteroscedasticity(ab.lm) # good


## Linear model for Pielous 

pie.lmm <- lmer(J ~ Factor + (1|Year) + (1|N),
                 data = invert, REML = TRUE)

summary(pie.lmm)
anova(pie.lmm)

#Type III Analysis of Variance Table with Satterthwaite's method
#       Sum Sq Mean Sq NumDF DenDF F value    Pr(>F)    
#Factor  1.326 0.66298     2    50  23.948 5.072e-08 ***  


r2_nakagawa(pie.lmm)

#Conditional R2: NA
#Marginal R2: 0.479

plot(pie.lmm)
qqnorm(resid(pie.lmm))
qqline(resid(pie.lmm))

plot(allEffects(pie.lmm))

check_model(pie.lmm)
check_singularity(pie.lmm) #TRUE

## ANCOVA Pielous

pie.lm <- lm(J ~ Treatment,
              data = invert)

anova(pie.lm)

#Response: J
#          Df Sum Sq Mean Sq F value    Pr(>F)    
#Treatment  2 1.3923 0.69618  25.207 2.441e-08 ***
#Residuals 51 1.4086 0.02762            

pie.hsd <- HSD.test(pie.lm, "Treatment")

#J groups
#Invaded   0.5181020      a
#Uninvaded 0.4479636      a
#Treated   0.1478616      b

check_model(pie.lm)




## all LMM results

tab_model(rich.glmm, pie.lmm,
          string.ci = "95% CI",
          string.pred = "Coefficients",
          string.p = "P value",
          dv.labels = c("Taxa Richness",
                        "Abundance",
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

riche <- ggplot(invert, aes(x = Depth, y = rich, 
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


panel <- ggarrange(riche, shannon,
          simpson, pielou,
          legend = "none",
          ncol = 4,
          labels = c("E","F","G","H"),
          align = "hv")

ggsave("Figures/invert_hab_depth.jpeg", panel)

panel2 <- ggarrange(riche, simpson, 
                   legend = "none",
                   ncol = 2,
                   labels = c("C","D"),
                   align = "hv")



# Boxplots ----------------------------------------------------------------

abund.box <- ggplot(invert, aes(x = Treatment, y = abundance)) +
  stat_boxplot(geom = "errorbar") +
  geom_boxplot(size = 1) +
  geom_jitter(data = invert,
              aes(fill = Treatment, shape = Treatment),
              size = 5,
              stroke = 1.5) +
  theme_classic(14) +
  labs(x = " ",
       y = "Abundance") +
  scale_fill_viridis(discrete = TRUE) +
  scale_shape_manual(values = c(21, 24, 22)) +
  theme(panel.border = element_rect(fill = NA)) +
  theme(legend.position = "none")


pielou.box <- ggplot(invert, aes(x = Treatment, y = J)) +
  stat_boxplot(geom = "errorbar") +
  geom_boxplot(size = 1) +
  geom_jitter(data = invert,
              aes(fill = Treatment, shape = Treatment),
              size = 5,
              stroke = 1.5) +
  theme_classic(14) +
  labs(x = " ",
       y = "Pielou's Evenness (J)") +
  scale_fill_viridis(discrete = TRUE) +
  scale_shape_manual(values = c(21, 24, 22)) +
  theme(panel.border = element_rect(fill = NA)) +
  theme(legend.position = "none") +
  ylim(0, 1)

riche.box <- ggplot(invert, aes(x = Treatment, y = rich)) +
  stat_boxplot(geom = "errorbar") +
  geom_boxplot(size = 1) +
  geom_jitter(data = invert,
              aes(fill = Treatment, shape = Treatment),
              size = 5,
              stroke = 1.5) +
  theme_classic(14) +
  labs(x = " ",
       y = "Taxonomic Richness") +
  scale_fill_viridis(discrete = TRUE) +
  scale_shape_manual(values = c(21, 24, 22)) +
  theme(panel.border = element_rect(fill = NA)) +
  ylim(0, 50) +
  theme(legend.position = "none")

emerg.boxplots <- ggarrange(abund.box, riche.box, pielou.box,
          nrow = 1,
          labels = c("D","E","F"))


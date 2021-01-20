library(tidyverse) # for figures, pipes etc.
library(vegan) # for NMDS
library(ggpubr) # for arranging figures
library(agricolae)
library(Hmisc)
library(car)
library(viridis)
library(lme4)
library(MuMIn)
library(AICcmodavg)
library(ggrepel)



# Data Import -------------------------------------------------------------

invert <- read.csv("Data/emerging_invertebrates.csv")
str(invert)

colnames(invert)

invert$Year <- as.factor(invert$Year)
invert <- invert %>% unite("TrtYr", Treatment,Year, remove = FALSE)

invert.data <- invert %>% select(Araneae:Crambidae)
invert.env <- invert %>% select(ID:N)

unique(invert$Habitat)

# Diversity Metrics -------------------------------------------------------

richness <- rowSums(invert.data > 0) # species richness
abundance <- rowSums(invert.data) # abundance
H <- diversity(invert.data) # Shannon Weiner
D1 <- diversity(invert.data, index = "simpson") #default is base log, but can change it
J <- H/log(specnumber(invert.data))


# make new little data frame
invert.uni<- invert.env
invert.uni$rich <- richness
invert.uni$abundance <- abundance
invert.uni$H <- H
invert.uni$D1 <- D1
invert.uni$J <- J

invert.univariate <- invert.uni

write.csv(invert.univariate, "Data/emerging_invertebrate_univariate.csv")


colnames(invert.uni)


# Univariate Analyses -----------------------------------------------------

invert.uni <- read.csv("Data/emerging_invertebrate_univariate.csv")
invert.uni$Year <- as.factor(invert.uni$Year)
invert.uni$Treatment <- as.factor(invert.uni$Treatment)


abundance.lm <- lm(abundance ~ Treatment * Year, data = invert.uni)
Anova(abundance.lm, type = 3)

#Response: abundance
#                Sum Sq Df F value   Pr(>F)   
#(Intercept)    901867  1  1.8186 0.183813   
#Habitat       2466358  2  2.4866 0.093849 . 
#Year            78937  1  0.1592 0.691690   
#Habitat:Year  5642917  2  5.6893 0.006063 **
#Residuals    23804326 48   


## GLMM abundance #

ggplot(invert.data, aes(x = abundance)) + 
  geom_histogram(binwidth = 100,
                 color="black", fill="white")


invert.glmm <- lmer(abundance ~ Treatment + (1|Year) + (1|N), data = invert.uni, REML = FALSE)
summary(invert.glmm)

#Linear mixed model fit by maximum likelihood  ['lmerMod']
#Formula: abundance ~ Treatment + (1 | Year) + (1 | N)
#Data: invert.uni

#AIC      BIC   logLik deviance df.resid 
#884.4    896.4   -436.2    872.4       48 

#Scaled residuals: 
#  Min      1Q  Median      3Q     Max 
#-2.1987 -0.5641 -0.0023  0.2849  3.9533 

#Random effects:
#Groups   Name        Variance  Std.Dev. 
#N        (Intercept) 8.354e-05   0.00914
#Year     (Intercept) 1.226e+05 350.20975
#Residual             5.664e+05 752.60881
#Number of obs: 54, groups:  N, 5; Year, 2

#Fixed effects:
#                    Estimate Std. Error t value
#(Intercept)           382.4      304.6   1.255
#TreatmentTreated     1373.7      250.9   5.476
#TreatmentUninvaded    132.4      250.9   0.528

#Correlation of Fixed Effects:
#             (Intr) TrtmnT
#TretmntTrtd -0.412       
#TrtmntUnnvd -0.412  0.500

#optimizer (nloptwrap) convergence code: 0 (OK)
#boundary (singular) fit: see ?isSingular


VarCorr(invert.glmm)

#Groups   Name        Std.Dev.  
#N        (Intercept) 9.1398e-03
#Year     (Intercept) 3.5021e+02
#Residual             7.5261e+02

r.squaredGLMM(invert.glmm) # marginal and conditional r2 for model 1

#       R2m       R2c
#[1,] 0.3134326 0.5267209

## coefficients

coef(invert.glmm)

#     (Intercept) HabitatRestored HabitatUninvaded
#2017    32.59426        1373.556         132.1111
#2018   732.96130        1373.556         132.1111

# examine the residuals 
plot(invert.glmm)
qqnorm(resid(invert.glmm))
qqline(resid(invert.glmm))


# null model 

ab.null <- lmer(abundance ~ (1|Year) + (1|N), data = invert.uni, REML = FALSE) #boundary(singular) = variances of one or more linear combination of effects are 0
summary(ab.null)

#Linear mixed model fit by maximum likelihood  ['lmerMod']
#Formula: abundance ~ (1 | Year) + (1 | N)
#Data: invert.uni

#AIC      BIC   logLik deviance df.resid 
#908.1    916.0   -450.0    900.1       50 

#Scaled residuals: 
#  Min      1Q  Median      3Q     Max 
#-1.0595 -0.5600 -0.3445  0.1872  3.9578 

#Random effects:
#Groups   Name        Variance  Std.Dev. 
#N        (Intercept) 2.753e-09 5.247e-05
#Year     (Intercept) 1.079e+05 3.285e+02
#Residual             9.639e+05 9.818e+02
#Number of obs: 54, groups:  N, 5; Year, 2

#              Fixed effects:
#              Estimate Std. Error t value
#(Intercept)    884.5      268.0     3.3

#optimizer (nloptwrap) convergence code: 0 (OK)
#boundary (singular) fit: see ?isSingular

confint(ab.null)
r.squaredGLMM(ab.null)

#  R2m       R2c
# 0      0.1006882


models <- list(invert.glmm, ab.null)
names <- c(1:2)

#Model selection based on AICc:
  
#  K   AICc Delta_AICc AICcWt Cum.Wt      LL
#1 6 886.22       0.00      1      1 -436.22
#2 4 908.90      22.68      0      1 -450.04


# or do it this way and get info re: test statistic

anova(invert.glmm, ab.null)

#Data: invert.uni
#Models:
#  ab.null: abundance ~ (1 | Year) + (1 | N)
#invert.glmm: abundance ~ Treatment + (1 | Year) + (1 | N)

#             npar    AIC    BIC  logLik deviance  Chisq Df Pr(>Chisq)    
#ab.null        4 908.08 916.04 -450.04   900.08                         
#invert.glmm    6 884.44 896.37 -436.22   872.44 27.649  2  9.913e-07 ***


## Richness 
# ANOVA

rich.lm <- lm(rich ~ Treatment * Year, data = invert.uni)
Anova(rich.lm, type = 2)

#Response: rich Type III SS

#Sum Sq Df  F value    Pr(>F)    
#(Intercept)  4489.0  1 140.9132 6.907e-16 ***
#Habitat       379.6  2   5.9573  0.004887 ** 
#Year          102.7  1   3.2245  0.078839 .  
#Habitat:Year   87.1  2   1.3672  0.264552    
#Residuals    1529.1 48  

#Anova Table (Type II tests)

#Response: rich
#                   Sum Sq Df F value    Pr(>F)    
#Treatment       840.70  2  13.034  3.01e-05 ***
#Year            498.07  1  15.444 0.0002719 ***
#Treatment:Year   93.59  2   1.451 0.2444183    
#Residuals      1548.00 48 

# GLMM 

ggplot(invert.uni, aes(x = rich)) + 
  geom_histogram(binwidth = 1,
                 color="black", fill="white")


rich.in.glmm <- lmer(rich ~ Treatment + (1|Year) + (1|N),
                     data = invert.uni, REML = FALSE)
summary(rich.in.glmm)

#Linear mixed model fit by maximum likelihood  ['lmerMod']
#Formula: rich ~ Treatment + (1 | Year) + (1 | N)
#Data: invert.uni

#AIC      BIC   logLik deviance df.resid 
#349.0    361.0   -168.5    337.0       48 

#Scaled residuals: 
#  Min      1Q  Median      3Q     Max 
#-2.2232 -0.6037  0.1341  0.6268  2.2744 

#Random effects:
#Groups   Name        Variance Std.Dev.
#N        (Intercept) 12.35    3.514   
#Year     (Intercept)  0.00    0.000   
#Residual             26.28    5.127   
#Number of obs: 54, groups:  N, 5; Year, 2

#Fixed effects:
#                   Estimate Std. Error t value
#(Intercept)          23.029      2.181  10.557
#TreatmentTreated     -8.814      1.727  -5.104
#TreatmentUninvaded   -2.944      1.709  -1.723

#Correlation of Fixed Effects:
#  (Intr) TrtmnT
#TretmntTrtd -0.440       
#TrtmntUnnvd -0.392  0.495
#optimizer (nloptwrap) convergence code: 0 (OK)
#boundary (singular) fit: see ?isSingular

VarCorr(rich.in.glmm)

#Groups   Name        Std.Dev.
#N        (Intercept) 3.5143  
#Year     (Intercept) 0.0000  
#Residual             5.1266  

r.squaredGLMM(rich.in.glmm) # marginal and conditional r2 for model 1

#       R2m       R2c
#[1,] 0.2614732 0.4975748

## coefficients
coef(rich.in.glmm)

#(Intercept) HabitatRestored HabitatUninvaded
#2017    22.04341       -9.555556        -3.111111
#2018    27.40104       -9.555556        -3.111111


confint(rich.in.glmm)

#                     2.5 %     97.5 %
#.sig01               1.595810  8.6648467
#.sig02               0.000000  6.7664837
#.sigma               4.265627  6.3189536
#(Intercept)         17.983681 27.8249298
#TreatmentTreated   -12.271214 -5.3680080
#TreatmentUninvaded  -6.358896  0.4700068

# examine the residuals 
plot(rich.in.glmm)
qqnorm(resid(rich.in.glmm))
qqline(resid(rich.in.glmm))

# null model 

rich.null <- lmer(rich ~ (1|Year) + (1|N), data = invert.uni, REML = FALSE) #boundary(singular) = variances of one or more linear combination of effects are 0
summary(rich.null)

#Linear mixed model fit by maximum likelihood  ['lmerMod']
#Formula: rich ~ (1 | Year) + (1 | N)
#Data: invert.uni

#AIC      BIC   logLik deviance df.resid 
#366.8    374.7   -179.4    358.8       50 

#Scaled residuals: 
#  Min       1Q   Median       3Q      Max 
#-1.98998 -0.83028 -0.00916  0.70006  1.96595 

#Random effects:
#  Groups   Name        Variance Std.Dev.
#N        (Intercept) 14.46    3.802   
#Year     (Intercept)  0.00    0.000   
#Residual             39.94    6.320   
#Number of obs: 54, groups:  N, 5; Year, 2

#Fixed effects:
#            Estimate Std. Error t value
#(Intercept)   18.665      2.138   8.731

#optimizer (nloptwrap) convergence code: 0 (OK)
#boundary (singular) fit: see ?isSingular

confint(rich.null)
r.squaredGLMM(rich.null)

#  R2m       R2c
#  0    0.2657744

models <- list(rich.in.glmm, rich.null)
names <- c(1:2)

aictab(cand.set = models, modnames = names,
       sort = TRUE, c.hat = 1, second.ord = TRUE,
       nobs = NULL)

#Model selection based on AICc:
  
#  K   AICc Delta_AICc AICcWt Cum.Wt      LL
#1 6 350.81       0.00      1      1 -168.51
#2 4 367.58      16.77      0      1 -179.38


# or do it this way and get info re: test statistic

anova(rich.in.glmm, rich.null)

#Data: invert.uni
#Models:
#rich.null: rich ~ (1 | Year) + (1 | N)
#rich.in.glmm: rich ~ Treatment + (1 | Year) + (1 | N)

#             npar    AIC    BIC  logLik deviance  Chisq Df Pr(>Chisq)    
#rich.null       4 366.76 374.72 -179.38   358.76                         
#rich.in.glmm    6 349.02 360.96 -168.51   337.02 21.741  2  1.901e-05 ***
 

# Univariate Figures ------------------------------------------------------

# Abundance

abundance <- ggplot(data = invert.uni, aes(x = Treatment, y = abundance, 
                              shape = Year, colour = Treatment),
       size = 4) +
  geom_violin(trim = FALSE, lwd = 1, position = position_dodge(0.5),
              colour = "black") +
  geom_jitter(position = position_jitterdodge(
    jitter.width = 0.1, dodge.width = 0.5), size = 3) +
  stat_summary(     
    aes(shape = Year),
    colour = "black",
    fun.data = "mean_se", fun.args = list(mult = 1), 
    geom = "pointrange", size = 1,
    position = position_dodge(0.5)) +
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

abundance


richness <- ggplot(data = invert.uni, aes(x = Treatment, y = rich, 
                              shape = Year, colour = Treatment),
       size = 4) +
  geom_violin(trim = FALSE, lwd = 1, position = position_dodge(0.8),
              colour = "black") +
  geom_jitter(position = position_jitterdodge(
    jitter.width = 0.1, dodge.width = 0.8), size = 3) +
  stat_summary(     
    aes(shape = Year),
    colour = "black",
    fun.data = "mean_se", fun.args = list(mult = 1), 
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
  ylim(0, 40) +
  guides(color = "none") +
  theme(legend.position = "right") 

richness

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


# Multivariate Analyses ---------------------------------------------------

invert.rel <- decostand(invert.data, "max", 2, na.rm = NULL) # divide by column max

write.csv(invert.rel, "Data/emerging_invert_relativized.csv")

##### perMANVOA ####

(per.inv <- adonis2(invert.rel ~ Treatment * Year,
                    data = invert,
                    permutations = 999, method = "bray"))

#Permutation test for adonis under reduced model
#Terms added sequentially (first to last)
#Permutation: free
#Number of permutations: 999

#                Df SumOfSqs    R2     F     Pr(>F)    
#Treatment       2   2.4582 0.12278 3.8561  0.001 ***
#Year            1   1.2360 0.06174 3.8779  0.001 ***
#Treatment:Year  2   1.0278 0.05133 1.6123  0.007 ** 
#Residual       48  15.2994 0.76415                  
#Total          53  20.0215 1.00000 

bugs.b <- vegdist(invert.rel, method = "bray")

groups <- factor(invert$TrtYr)

(dispersion <- betadisper(bugs.b, groups)) 

#Average distance to median:
#Invaded_2017   Invaded_2018   Treated_2017   Treated_2018 Uninvaded_2017 Uninvaded_2018 
#0.5348         0.4988         0.5285         0.4853         0.5607         0.5435

plot(dispersion)
boxplot(dispersion)

library(devtools)
install_github("GuillemSalazar/EcolUtils")
library(EcolUtils)

(adonis.pair(bugs.b, groups, nper = 1000, corr.method = "bonferroni"))

#                         combination SumsOfSqs   MeanSqs  F.Model         R2     P.value P.value.corrected
#1      Invaded_2017 <-> Invaded_2018 0.8250993 0.8250993 2.653724 0.14226245 0.000999001        0.01498501
#2      Invaded_2017 <-> Treated_2017 1.0810539 1.0810539 3.296064 0.17081534 0.000999001        0.01498501
#3      Invaded_2017 <-> Treated_2018 1.3744035 1.3744035 4.582915 0.22265626 0.000999001        0.01498501
#4    Invaded_2017 <-> Uninvaded_2017 0.5422269 0.5422269 1.570793 0.08939794 0.051948052        0.77922078
#5    Invaded_2017 <-> Uninvaded_2018 0.7485096 0.7485096 2.236271 0.12262766 0.000999001        0.01498501
#6      Invaded_2018 <-> Treated_2017 1.2823550 1.2823550 4.168120 0.20666875 0.000999001        0.01498501
#7      Invaded_2018 <-> Treated_2018 1.2592004 1.2592004 4.504033 0.21966569 0.000999001        0.01498501
#8    Invaded_2018 <-> Uninvaded_2017 1.0150108 1.0150108 3.124382 0.16337165 0.000999001        0.01498501
#9    Invaded_2018 <-> Uninvaded_2018 0.4905781 0.4905781 1.560423 0.08886024 0.022977023        0.34465534
#10     Treated_2017 <-> Treated_2018 0.6968943 0.6968943 2.349342 0.12803411 0.000999001        0.01498501
#11   Treated_2017 <-> Uninvaded_2017 0.8748130 0.8748130 2.558458 0.13785941 0.000999001        0.01498501
#12   Treated_2017 <-> Uninvaded_2018 1.0702554 1.0702554 3.229011 0.16792394 0.002997003        0.04495504
#13   Treated_2018 <-> Uninvaded_2017 1.1827088 1.1827088 3.768464 0.19063010 0.000999001        0.01498501
#14   Treated_2018 <-> Uninvaded_2018 0.9810781 0.9810781 3.233999 0.16813971 0.000999001        0.01498501
#15 Uninvaded_2017 <-> Uninvaded_2018 0.7418541 0.7418541 2.127731 0.11737438 0.001998002        0.02997003


############ NMDS ############

library(ggrepel)

## stress plot/ see how many axes make sense

k_vec <- 1:10 #dimensions 1 - 10
stress <- numeric(length(k_vec)) # stress of each model put here
dune_dij <- metaMDSdist(invert.rel)

set.seed(25)

for(i in seq_along(k_vec)) {
  sol <- metaMDSiter(dune_dij, k = i, # you can tell I lifted this from a tutorial on the dune package
                     trace = FALSE)
  stress[i] <- sol$stress
}

plot(stress) # 3D makes sense

#### NMDS analysis ####
set.seed(120) 

nms.invert <- metaMDS(invert.rel, distance = "bray", # species data, bray-curtis dissimilarity
                      autotransform = FALSE,  # NMDS will do autotransformations for you
                      k = 3, trymax = 1000)   # k = number of axes
nms.invert


#global Multidimensional Scaling using monoMDS

#Data:     invert.rel 
#Distance: bray 

#Dimensions: 3 
#Stress:     0.2017306 
#Stress type 1, weak ties
#Two convergent solutions found after 22 tries
#Scaling: centring, PC rotation, halfchange scaling 
#Species: expanded scores based on ‘invert.rel’ 


# look at points and the stress real quick
layout(matrix(1:2, ncol = 2))
plot(nms.invert, main = "Invertebrate NMDS plot"); stressplot(nms.invert, main = "Shepard plot")
layout(1)

ordiplot(nms.invert, type = "n")
orditorp(nms.invert, display = "species")
orditorp(nms.invert, display = "sites")


ordiplot(nms.invert, type = "n", choices = c(1,3))
orditorp(nms.invert, display = "species", choices = c(1,3))
orditorp(nms.invert, display = "sites", choices = c(1,3))


# how many iterations of the NMDS
nms.invert$iters # 136

# Goodness of fit
(g <- goodness(nms.invert)) # smaller the number the better the fit
sum(g^2)
nms.invert$stress^2  # 0.04069524

1-nms.invert$stress^2 # 0.9593048 #analogous to square correlation coefficient


## extract the scores for plotting 
scr <- as.data.frame(scores(nms.invert, display = "sites")) # extract NMDS scores
colnames(scr)
 
# adding categorical info to scores
invert.env$NMDS1 <- scr$NMDS1
invert.env$NMDS2 <- scr$NMDS2
invert.env$NMDS3 <- scr$NMDS3

scores <- invert.env %>% select(ID:Year,NMDS1:NMDS3)

write.csv(scores,"Data/NMDS_emerging_inverts_NMDSscores.csv") # save this as a csv

### Vectors correlated with Axis 1 & 2 ####

alltaxa12 <- envfit(nms.invert, invert.rel,
                  choices = c(1,2)) #produces a list with r2, p value, and NMDS coordinates

all.taxa.df <- data.frame((alltaxa12$vectors)$arrows,
                          (alltaxa12$vectors)$r,
                          (alltaxa12$vectors)$pvals) #take list and make into data frame

write.csv(all.taxa.df, "Data/NMDS_emerg_vectors_axis12.csv") # save vector scores as csv


#### Vectors correlated with axis 1 & 3 ###

alltaxa.13 <- envfit(nms.invert, invert.rel, 
                     permutations = 999, choices = c(1,3)) 


all.taxa13.df <- data.frame((alltaxa.13$vectors)$arrows,
                            (alltaxa.13$vectors)$r,
                            (alltaxa.13$vectors)$pvals)

write.csv(all.taxa13.df, "Data/NMDS_emerging_vectors_axis13.csv")



## Picking the vectors we want for the figure based on fit (r > 0.2)

nmds.axis12 <- read.csv("Data/NMDS_emerg_vectors_axis12.csv")
nmds.axis13 <- read.csv("Data/NMDS_emerging_vectors_axis13.csv")

colnames(nmds.axis13)

# axis 1, 2
corr.sp <- nmds.axis12 %>% filter(X.alltaxa12.vectors..r > 0.2) 
target12 <- corr.sp$X # string of the Family names

axis12.vectors <- invert.rel %>% select(all_of(target12)) # make a matrix of just those

# axis 1, 3

corr.sp13 <- nmds.axis13 %>% filter(X.alltaxa.13.vectors..r > 0.2) 
target13 <- corr.sp13$X # string of the Family names

axis13.vectors <- invert.rel %>% select(all_of(target13)) # make a matrix of just those


# fit them to the nms
# axis 1, 2
(nmds.vectors.12 <- envfit(nms.invert$points, axis12.vectors,
                         permutations = 999, choices = c(1,2)))                        


corr.vectors.12 <- as.data.frame(nmds.vectors.12$vectors$arrows*sqrt(nmds.vectors.12$vectors$r)) #scaling vectors
corr.vectors.12$species <- rownames(corr.vectors.12) # add Family as a column

write.csv(corr.vectors.12, "Data/NMDS_emerging_correlatedvectors_axis12.csv")

# axis 1, 3
(nmds.vectors.13 <- envfit(nms.invert$points, axis13.vectors,
                           permutations = 999, choices = c(1,3)))                        


corr.vectors.13 <- as.data.frame(nmds.vectors.13$vectors$arrows*sqrt(nmds.vectors.13$vectors$r)) #scaling vectors
corr.vectors.13$species <- rownames(corr.vectors.13) # add Family as a column

write.csv(corr.vectors.13, "Data/NMDS_emerging_correlatedvectors_axis13.csv")

## plotting with base R ##
# I don't like base R but it has some nice things built in for ordinations

# axis 1, 2 

col_vec <- c("#9970ab", "#1b7837", "#2166ac")

ordiplot(nms.invert, choices = c(1,2), 
         type = "points",
         display = "sites")


ordihull(nms.invert, groups = scores$Treatment, # ellipse hull
         col = col_vec)

plot(alltaxa12, p.max = 0.001, col = "black")


ordiplot(nms.invert, choices = c(1,3), 
         type = "points",
         display = "sites")


ordihull(nms.invert, groups = scores$Treatment, # ellipse hull
         col = col_vec)

plot(alltaxa13, p.max = 0.013, col = "black")

#### ggPlot Figures ####

# NMDS Figures ------------------------------------------------------------
# load the data so you don't have to redo it every time

nmds.scores <- read.csv("Data/NMDS_emerging_inverts_NMDSscores.csv") #points
nmds.scores$Year <- as.factor(nmds.scores$Year)
nmds.scores <- nmds.scores %>% unite("HabYr", Habitat,Year, remove = FALSE)

#vectors r > 0.2
vectors.12 <- read.csv("Data/NMDS_emerging_correlatedvectors_axis12.csv") 
vectors.13 <- read.csv("Data/NMDS_emerging_correlatedvectors_axis13.csv")

## NMDS Axis 1, 2 

invert.12 <- ggplot(data = nmds.scores,
                    aes(x = NMDS1, y = NMDS2)) +
  geom_point(data = nmds.scores, 
             aes(x = NMDS1, y = NMDS2, colour = Treatment, shape = HabYr),
             size = 4, stroke = 1.5) + # sites as points
  stat_ellipse(data = nmds.scores, aes(x = NMDS1,y = NMDS2,
                                  linetype = Year, colour = Treatment), 
               size = 1, level = 0.9) + 
  geom_segment(data = vectors.12, aes(x = 0, xend = MDS1, y = 0, yend = MDS2), # adding in the vectors, c
               arrow = arrow(length = unit(0.5, "cm")), colour = "black") + # can add in geom_label or geom_text for labels
  theme_minimal() + # no background
  theme(panel.border = element_rect(fill = NA)) + # full square around figure
  xlab("NMDS 1") +
  ylab("NMDS 2") +
  #ylim(-1, 1.5) +
  #xlim(-1.45, 1) +
  #theme(legend.position = "none") +
  geom_text_repel(data = vectors.12, 
                  aes(x = MDS1, y = MDS2, label = X),
                  color="black",
                  size = 6) +
  scale_color_manual(values = c("#969696","#35978f", "#2166ac")) +
  scale_shape_manual(values = c(15,0, 16, 1, 17, 2, 18, 5))

invert.12

## NMDS Axis 1, 3

invert.13 <- ggplot(data = nmds.scores,
                    aes(x = NMDS1, y = NMDS3)) +
  geom_point(data = nmds.scores, 
             aes(x = NMDS1, y = NMDS3, 
                 colour = Treatment, shape = HabYr), 
             size = 4, stroke = 1.5) +
  stat_ellipse(data = nmds.scores, 
               aes(x = NMDS1,y = NMDS3,linetype = Year, 
                   colour = Treatment), 
               size = 1, level = 0.9) +
  geom_segment(data = vectors.13, 
               aes(x = 0, xend = MDS1, y = 0, yend = MDS3),
               arrow = arrow(length = unit(0.5, "cm")), colour = "black") +
  theme_minimal() +
  theme(panel.border = element_rect(fill = NA)) +
  geom_text_repel(data = vectors.13, 
                  aes(x = MDS1, y = MDS3, label = X),
                  color="black",
                  size = 6) +
  scale_color_manual(values = c("#969696","#35978f", "#2166ac")) +
  scale_shape_manual(values = c(15,0, 16, 1, 17, 2, 18, 5))

invert.13


(NMS.emerging.panel <- ggarrange(invert.12, invert.13,
                                 common.legend = TRUE,
                                 legend = "bottom"))
  

ggsave("Figures/Emerging_NMDS_panel.TIFF", NMS.emerging.panel,
       dpi = 300)




# NMDS Clusters -----------------------------------------------------------

# added the grouping variable derived from the ISA and cluster analysis

inverts <- read.csv("Data/emerging_invert_relativized_10 groups.csv")

str(inverts)
dim(inverts)
colnames(inverts)

invert.sp <- inverts %>% select(Araneae:Crambidae)# just the Families
env.group <- inverts %>% select(Site:Group4) # Categorical variables

env.group$Group4 <- as.factor(env.group$Group4) # adding in the cluster groups
env.group$Group5 <- as.factor(env.group$Group5)


scores$ClusterGroups4 <- env.group$Group4
scores$ClusterGroups5 <- env.group$Group5

colnames(scores)

col_order <- c("ID", "Treatment", "Habitat", "Year", "TrtYr", "ClusterGroups4", "ClusterGroups4", "NMDS1", "NMDS2", "NMDS3")

scores <- scores[, col_order] # put the categorical values in order

write.csv(scores,"Data/NMDS_emerg_inverts_scores_clusters.csv") 

# for later ##

ISA.scores <- read.csv("Data/NMDS_emerg_inverts_scores_clusters.csv")
colnames(ISA.scores)

ISA.scores$ClusterGroup5 <- as.factor(ISA.scores$ClusterGroup5)
ISA.scores$ClusterGroups4 <- as.factor(ISA.scores$ClusterGroups4)
ISA.scores$ClusterGroup3 <- as.factor(ISA.scores$ClusterGroup3)

ISA.scores <- ISA.scores %>% unite("HabYr", Habitat,Year, remove = FALSE)

# selecting the ISA with p < 0.05 for the NMDS

colnames(invert.sp)
colnames(ISA.clust4)

ISA.clust4 <- read.csv("Data/ISA_fourclusters.csv")

ISA.spp <- ISA.clust4 %>% filter(pvalue < 0.001) # should be 39
target <- ISA.spp$Species # string of species names

ISA <- invert.sp %>% select(all_of(target)) # make a matrix of just those


colnames(ISA)
dim(ISA) # 39 columns and 54 sites

# calculate vectors for Indicator Species ####

(ISA.vector.12 <- envfit(nms.invert$points, ISA,
                         permutations = 999, choices = c(1,2)))                        


ISA.12 <- as.data.frame(ISA.vector.12$vectors$arrows*sqrt(ISA.vector.12$vectors$r)) #scaling vectors
ISA.12$species <- rownames(ISA.12) # add Family as a column


(ISA.vector.13 <- envfit(nms.invert$points, ISA,
                         permutations = 999, choices = c(1,3)))   


ISA.13 <- as.data.frame(ISA.vector.13$vectors$arrows*sqrt(ISA.vector.13$vectors$r)) #scaling vectors
ISA.13$species <- rownames(ISA.13)

###### NMDS with Cluster/ISA groups ###

ISA.invert.12 <- ggplot(data = ISA.scores,
                        aes(x = NMDS1, y = NMDS2)) +
  geom_point(data = ISA.scores, aes(x = NMDS1, y = NMDS2, 
                                colour = ClusterGroups4, shape = HabYr), 
             size = 4, stroke = 1.5) + # sites as points
  stat_ellipse(data = ISA.scores, aes(x = NMDS1,y = NMDS2,
                                  linetype = ClusterGroups4, colour = ClusterGroups4), 
               size = 1, type = "norm") + 
  geom_segment(data = ISA.12, aes(x = 0, xend = MDS1,
                                  y = 0, yend = MDS2), # adding in the vectors, c
               arrow = arrow(length = unit(0.5, "cm")), colour = "black") + # can add in geom_label or geom_text for labels
  theme_classic() + # no background
  theme(panel.border = element_rect(fill = NA)) + # full square around figure
  xlab("NMDS 1") +
  ylab("NMDS 2") +
  geom_label(data = ISA.12,aes(x=MDS1,y=MDS2,label=species),size=5) +
  scale_colour_viridis(discrete = TRUE) +
  scale_shape_manual(values = c(15,0, 16, 1, 17, 2, 18, 5))

ISA.invert.12

## NMDS Axis 1, 3
# same as above

ISA.invert.13 <- ggplot(data = ISA.scores,
                        aes(x = NMDS1, y = NMDS3)) +
  geom_point(data = ISA.scores , aes(x = NMDS1, y = NMDS3, 
                                colour = ClusterGroups4, shape = HabYr),
             size = 4, stroke = 1.5) + # sites as points
  stat_ellipse(data = ISA.scores, aes(x = NMDS1,y = NMDS3,
                                  colour = ClusterGroups4, linetype = ClusterGroups4), 
               size = 1, level = 0.9) + 
  geom_segment(data = ISA.13, aes(x = 0, xend = MDS1,
                                  y = 0, yend = MDS3), # adding in the vectors, c
               arrow = arrow(length = unit(0.5, "cm")), colour = "black") + # can add in geom_label or geom_text for labels
  theme_classic() + # no background
  theme(panel.border = element_rect(fill = NA)) + # full square around figure
  xlab("NMDS 1") +
  ylab("NMDS 3") +
  geom_label(data = ISA.13,aes(x = MDS1,y = MDS3,label = species), size = 5) +
  scale_colour_viridis(discrete = TRUE) + 
  scale_shape_manual(values = c(15,0, 16, 1, 17, 2, 18, 5))


ISA.invert.13

NMDS.inv.ISA <- ggarrange(ISA.invert.12, ISA.invert.13, 
                          common.legend = TRUE, 
                          legend = "bottom",
                          align = "hv")
NMDS.inv.ISA

ggsave("Figures/NMDS_invertebrate_ISAgroups.jpeg", NMDS.inv.ISA) 


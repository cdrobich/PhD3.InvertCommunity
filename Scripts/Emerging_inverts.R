library(tidyverse) # for figures, pipes etc.
library(vegan) # for NMDS
library(ggpubr) # for arranging figures
library(agricolae)
library(Hmisc)
library(car)
library(viridis)


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



abundance.lm <- lm(abundance ~ Treatment * Year, data = invert.uni)
Anova(abundance.lm, type = 3)

#Response: abundance
#Sum Sq Df F value   Pr(>F)   
#(Intercept)    901867  1  1.8186 0.183813   
#Habitat       2466358  2  2.4866 0.093849 . 
#Year            78937  1  0.1592 0.691690   
#Habitat:Year  5642917  2  5.6893 0.006063 **
#Residuals    23804326 48   


## GLMM abundance #

library(lme4)
library(MuMIn)

ggplot(invert.data, aes(x = abundance)) + 
  geom_histogram(binwidth = 100,
                 color="black", fill="white")


invert.glmm <- lmer(abundance ~ Habitat + (1|Year), data = invert.uni, REML = FALSE)
summary(invert.glmm)

#Linear mixed model fit by REML ['lmerMod']
#Formula: abundance ~ Habitat + (1 | Year)
#Data: invert.uni

#REML criterion at convergence: 833.6

#Scaled residuals: 
#  Min      1Q  Median      3Q     Max 
#-2.1911 -0.5870  0.0113  0.3136  3.8420 

#Random effects:
#Groups   Name        Variance Std.Dev.
#Year     (Intercept) 265414   515.2   
#Residual             588945   767.4   
#Number of obs: 54, groups:  Year, 2

#Fixed effects:
#Estimate Std. Error t value
#(Intercept)         382.8      406.7   0.941
#HabitatRestored    1373.6      255.8   5.369
#HabitatUninvaded    132.1      255.8   0.516

#Correlation of Fixed Effects:
#(Intr) HbttRs
#HabittRstrd -0.314       
#HabttUnnvdd -0.314  0.500

VarCorr(invert.glmm)

#Groups   Name        Std.Dev.
#Year     (Intercept) 515.18  
#Residual             767.43 

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




## Richness 
# ANOVA

rich.lm <- lm(rich ~ Habitat * Year, data = invert.uni)
Anova(rich.lm, type = 3)

#Response: rich

#Sum Sq Df  F value    Pr(>F)    
#(Intercept)  4489.0  1 140.9132 6.907e-16 ***
#Habitat       379.6  2   5.9573  0.004887 ** 
#Year          102.7  1   3.2245  0.078839 .  
#Habitat:Year   87.1  2   1.3672  0.264552    
#Residuals    1529.1 48  



# GLMM 

library(lme4)
library(MuMIn)

ggplot(invert.uni, aes(x = rich)) + 
  geom_histogram(binwidth = 1,
                 color="black", fill="white")


rich.in.glmm <- lmer(rich ~ Habitat + (1|Year), data = invert.uni, REML = FALSE)
summary(rich.in.glmm)

#Linear mixed model fit by maximum likelihood  ['lmerMod']
#Formula: rich ~ Habitat + (1 | Year)
#Data: invert.uni

#AIC      BIC   logLik deviance df.resid 
#353.0    363.0   -171.5    343.0       49 

#Scaled residuals: 
#  Min      1Q  Median      3Q     Max 
# -2.2045 -0.5940  0.1993  0.7097  2.2798 

#Random effects:
#Groups   Name        Variance Std.Dev.
#Year     (Intercept)  8.185   2.861   
#Residual             31.081   5.575   
#Number of obs: 54, groups:  Year, 2

#Fixed effects:
#Estimate Std. Error t value
#(Intercept)        24.722      2.412  10.248
#HabitatRestored    -9.556      1.858  -5.142
#HabitatUninvaded   -3.111      1.858  -1.674

#Correlation of Fixed Effects:
#  (Intr) HbttRs
#HabittRstrd -0.385       
#HabttUnnvdd -0.385  0.500


VarCorr(rich.in.glmm)

#Groups   Name        Std.Dev.
#Year     (Intercept) 2.8610  
#Residual             5.5751 

r.squaredGLMM(rich.in.glmm) # marginal and conditional r2 for model 1

#       R2m       R2c
#[1,] 0.2912272 0.438974

## coefficients
coef(rich.in.glmm)

#(Intercept) HabitatRestored HabitatUninvaded
#2017    22.04341       -9.555556        -3.111111
#2018    27.40104       -9.555556        -3.111111

# examine the residuals 
plot(rich.in.glmm)
qqnorm(resid(rich.in.glmm))
qqline(resid(rich.in.glmm))



# Univariate Figures ------------------------------------------------------

# Abundance

abundance <- ggplot(data = invert.uni, aes(x = Treatment, y = abundance, 
                              shape = Year, colour = Treatment),
       size = 4) +
  geom_violin(trim = FALSE, lwd = 1, position = position_dodge(0.5)) +
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
  geom_violin(trim = FALSE, lwd = 1, position = position_dodge(0.8)) +
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

#Data:     invert.rel 
#Distance: bray 

#Dimensions: 3 
#Stress:     0.2014091 
#Stress type 1, weak ties
#Two convergent solutions found after 57 tries
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
nms.invert$iters # 129

# Goodness of fit
(g <- goodness(nms.invert)) # smaller the number the better the fit
sum(g^2)
nms.invert$stress^2  # 0.04056562

1-nms.invert$stress^2 # 0.9594344 #analogous to square correlation coefficient


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
                          (alltaxa12$vectors)$pvals) #take list and make into dataframe

write.csv(all.taxa.df, "Data/NMDS_emerg_vectors_axis12.csv") # save vector scores as csv


alltaxa12$vectors$r[alltaxa12$vectors$r > 0.25] # selecting vectors (Family) that are reasonably correlated (r2 > 0.2)

# Anthicidae   Chironomidae Dolichopodidae  Ichneumonidae 
# 0.3335932      0.2748834      0.3038629      0.2892768 


# taking out those correlated ones to add to the figure
corr.taxa12 <- invert.rel %>% select(Anthicidae,
                                    Chironomidae,
                                    Dolichopodidae,
                                    Ichneumonidae)

# recalculated it because I am lazy but you could probably pull them from all.taxa.df
corrtaxa12 <- envfit(nms.invert$points, corr.taxa12, 
                   permutations = 999, choices = c(1,2))


corrtaxa12

# make a new data frame for the figure
species.12 <- as.data.frame(corrtaxa12$vectors$arrows*sqrt(corrtaxa12$vectors$r)) #scaling vectors so they correspond with r2
species.12$species <- rownames(species.12)


#### Vectors correlated with axis 1 & 3 ###
# same as above but now for the other axes/figure

alltaxa.13 <- envfit(nms.invert, invert.rel, 
                     permutations = 999, choices = c(1,3)) 


all.taxa13.df <- data.frame((alltaxa.13$vectors)$arrows,
                            (alltaxa.13$vectors)$r,
                            (alltaxa.13$vectors)$pvals)

write.csv(all.taxa13.df, "Data/NMDS_emerging_vectors_axis13.csv")


alltaxa.13$vectors$r[alltaxa.13$vectors$r > 0.25] 

#Latridiidae Chironomidae 
#0.3044679    0.2699855 

corr.taxa.13 <- invert.rel %>% select(Latridiidae,
                                       Chironomidae)


corrtaxa.13 <- envfit(nms.invert$points, corr.taxa.13, 
                      permutations = 999, choices = c(1,3))


corrtaxa.13

species.13 <- as.data.frame(corrtaxa.13$vectors$arrows*sqrt(corrtaxa.13$vectors$r)) #scaling vectors
species.13$species <- rownames(species.13)



## plotting with base R ##
# I don't like base R but it has some nice things built in for ordinations

# axis 1, 2 

col_vec <- c("#9970ab", "#1b7837", "#2166ac")

ordiplot(nms.invert, choices = c(1,2), 
         type = "points",
         display = "sites")


ordihull(nms.invert, groups = scores$Treatment, # ellipse hull
         col = col_vec)

plot(alltaxa12, p.max = 0.013, col = "black")


#### ggPlot Figures ####

## NMDS Axis 1, 2 

scores # coordinates we extracted
species.12 # reasonably correlated vectors with axis 1,2

invert.12 <- ggplot(data = scores,
                    aes(x = NMDS1, y = NMDS2)) +
  geom_point(data = scores, 
             aes(x = NMDS1, y = NMDS2, colour = Treatment, shape = Habitat),
             size = 4) + # sites as points
  stat_ellipse(data = scores, aes(x = NMDS1,y = NMDS2,
                                  linetype = TrtYr, colour = Treatment), size = 1) + # a 95% CI ellipses
  geom_segment(data = species.12, aes(x = 0, xend = MDS1, y = 0, yend = MDS2), # adding in the vectors, c
               arrow = arrow(length = unit(0.5, "cm")), colour = "black") + # can add in geom_label or geom_text for labels
  theme_minimal() + # no background
  theme(panel.border = element_rect(fill = NA)) + # full square around figure
  xlab("NMDS 1") +
  ylab("NMDS 2") +
  #ylim(-1, 1.5) +
  #xlim(-1.45, 1) +
  #theme(legend.position = "none") +
  geom_text_repel(data = species.12, 
                  aes(x = MDS1, y = MDS2, label = species),
                  color="black",
                  size = 6) +
  scale_color_manual(values = c("#969696","#35978f", "#2166ac")) +
  scale_shape_manual(values = c(15, 16, 17, 18)) +
  guides(linetype = "none") +
  coord_fixed()

invert.12

## NMDS Axis 1, 3
# same as above

scores
species.13 # vectors correlated with axis 1, 3

invert.13 <- ggplot(data = scores,
                    aes(x = NMDS1, y = NMDS3)) +
  geom_point(data = scores, 
             aes(x = NMDS1, y = NMDS3, 
                 colour = Treatment, shape = Habitat), size = 4) +
  stat_ellipse(data = scores, 
               aes(x = NMDS1,y = NMDS3,linetype = TrtYr, 
                   colour = Treatment), size = 1) +
  geom_segment(data = species.13, 
               aes(x = 0, xend = MDS1, y = 0, yend = MDS3),
               arrow = arrow(length = unit(0.5, "cm")), colour = "black") +
  theme_minimal() +
  theme(panel.border = element_rect(fill = NA)) +
  geom_text_repel(data = species.13, 
                  aes(x = MDS1, y = MDS3, label = species),
                  color="black",
                  size = 6) +
  scale_color_manual(values = c("#969696","#35978f", "#2166ac")) +
  scale_shape_manual(values = c(15, 16, 17, 18)) +
  guides(linetype = "none") +
  coord_fixed()

invert.13


(NMS.emerging.panel <- ggarrange(invert.12, invert.13,
                                 common.legend = TRUE,
                                 legend = "bottom"))
  

ggsave("Figures/Emerging_NMDS_panel.TIFF", NMS.emerging.panel,
       dpi = 300)





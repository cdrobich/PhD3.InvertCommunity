
# Load libraries ----------------------------------------------------------

library(tidyverse) # for figures, pipes etc.
library(vegan) # for NMDS
library(ggpubr) # for arranging figures

library(ggrepel) # labels on nmds

library(viridis) # colours

# pairwise perMANOVA comparison
library(devtools)
install_github("GuillemSalazar/EcolUtils")
library(EcolUtils)


# Load data ---------------------------------------------------------------

invert <- read.csv("Data/Emerging/emerging_invertebrates.csv")

invert$Year <- as.factor(invert$Year)
invert <- invert %>% unite("TrtYr", Treatment,Year, remove = FALSE)
invert <- invert %>% unite("HabYr", Habitat,Year, remove = FALSE)

colnames(invert)

invert.data <- invert %>% select(Araneae:Crambidae)
invert.env <- invert %>% select(Site:N)

# Multivariate Analyses ---------------------------------------------------

invert.rel <- decostand(invert.data, "max", 2, na.rm = NULL) # relativize by column max
write.csv(invert.rel, "Data/Emerging/emerging_invert_relativized.csv")

invert.rel <- read.csv("Data/Emerging/emerging_invert_relativized.csv")

# Two-way perMANOVA -------------------------------------------------------

(per.inv <- adonis2(invert.rel ~ Treatment * Year,
                    data = invert,
                    permutations = 999, method = "bray"))

#Permutation test for adonis under reduced model
#Terms added sequentially (first to last)
#Permutation: free
#Number of permutations: 999

#                Df SumOfSqs      R2      F Pr(>F)    
#Treatment       2   2.4255 0.11995 3.7364  0.001 ***
#year            1   1.2082 0.05975 3.7225  0.001 ***
#Treatment:Year  2   1.0067 0.04979 1.5509  0.010 ** 
#Residual       48  15.5795 0.77050                  
#Total          53  20.2199 1.00000


# betadisper --------------------------------------------------------------
# multivariate homogeneity of groups dispersions

bugs.b <- vegdist(invert.rel, method = "bray")

groups <- factor(invert$TrtYr)

(dispersion <- betadisper(bugs.b, groups)) 

#Homogeneity of multivariate dispersions
#
#Call: betadisper(d = bugs.b, group = groups)
#
#No. of Positive Eigenvalues: 33
#No. of Negative Eigenvalues: 20
#
#Average distance to median:
# Invaded_2017   Invaded_2018   Treated_2017   Treated_2018 Uninvaded_2017 
#0.40914        0.14508        0.14913        0.07226        0.13264 
#Uninvaded_2018 
#0.10442 
#
#Eigenvalues for PCoA axes:
#  (Showing 8 of 53 eigenvalues)
#PCoA1  PCoA2  PCoA3  PCoA4  PCoA5  PCoA6  PCoA7  PCoA8 
#3.5783 1.4234 0.4978 0.3054 0.2956 0.2588 0.1874 0.1819 

anova(dispersion)

#Analysis of Variance Table, TrtYr

#Response: Distances
#Df  Sum Sq  Mean Sq F value    Pr(>F)    
#Groups     5 0.66138 0.132275  53.522 < 2.2e-16 ***
#Residuals 48 0.11863 0.002471 

plot(dispersion)
boxplot(dispersion)

bugs.b <- vegdist(invert.rel, method = "bray")

groups.trt <- factor(invert$Treatment)

(dispersion2 <- betadisper(bugs.b, groups.trt)) 

#Homogeneity of multivariate dispersions
#
#Call: betadisper(d = bugs.b, group = groups.trt)
#
#No. of Positive Eigenvalues: 33
#No. of Negative Eigenvalues: 20
#
#Average distance to median:
#  Invaded   Treated Uninvaded 
#0.4223    0.2742    0.2398 
#
#Eigenvalues for PCoA axes:
#  (Showing 8 of 53 eigenvalues)
#PCoA1  PCoA2  PCoA3  PCoA4  PCoA5  PCoA6  PCoA7  PCoA8 
#3.5783 1.4234 0.4978 0.3054 0.2956 0.2588 0.1874 0.1819 

anova(dispersion2)

#Analysis of Variance Table

#Response: Distances
#Df  Sum Sq  Mean Sq F value    Pr(>F)    
#Groups     2 0.33851 0.169257  9.5471 0.0003007 ***
#Residuals 51 0.90416 0.017729 

permutest(dispersion2, pairwise = TRUE, permutations = 99)

#Response: Distances
#         Df  Sum Sq  Mean Sq      F N.Perm Pr(>F)   
#Groups     2 0.33851 0.169257 9.5471     99   0.01 **
#Residuals 51 0.90416 0.017729  

disp2.HSD <- TukeyHSD(dispersion2)

Fit: aov(formula = distances ~ group, data = df)

#$group
#diff        lwr         upr     p adj
#Treated-Invaded   -0.14808096 -0.2552203 -0.04094161 0.0044561
#Uninvaded-Invaded -0.18250022 -0.2896396 -0.07536088 0.0004140
#Uninvaded-Treated -0.03441926 -0.1415586  0.07272008 0.7196320

plot(dispersion2)
boxplot(dispersion2)


(adonis.pair(bugs.b, groups, nper = 1000, corr.method = "bonferroni"))

#                          combination SumsOfSqs   MeanSqs  F.Model         R2     P.value P.value.corrected
#1      Invaded_2017 <-> Invaded_2018 0.8301707 0.8301707 2.634429 0.14137430 0.000999001        0.01498501
#2      Invaded_2017 <-> Treated_2017 1.0521256 1.0521256 3.141491 0.16411947 0.000999001        0.01498501
#3      Invaded_2017 <-> Treated_2018 1.3454776 1.3454776 4.388906 0.21525951 0.000999001        0.01498501
#4    Invaded_2017 <-> Uninvaded_2017 0.5291157 0.5291157 1.510852 0.08628089 0.061938062        0.92907093
#5    Invaded_2017 <-> Uninvaded_2018 0.7430752 0.7430752 2.204071 0.12107570 0.001998002        0.02997003
#6      Invaded_2018 <-> Treated_2017 1.2617401 1.2617401 4.006533 0.20026125 0.000999001        0.01498501
#7      Invaded_2018 <-> Treated_2018 1.2498756 1.2498756 4.361483 0.21420263 0.000999001        0.01498501
#8    Invaded_2018 <-> Uninvaded_2017 0.9966132 0.9966132 3.018045 0.15869378 0.000999001        0.01498501
#9    Invaded_2018 <-> Uninvaded_2018 0.4871958 0.4871958 1.536190 0.08760112 0.030969031        0.46453546
#10     Treated_2017 <-> Treated_2018 0.6511566 0.6511566 2.125459 0.11726373 0.001998002        0.02997003
#11   Treated_2017 <-> Uninvaded_2017 0.8534432 0.8534432 2.438358 0.13224378 0.001998002        0.02997003
#12   Treated_2017 <-> Uninvaded_2018 1.0586197 1.0586197 3.141913 0.16413787 0.000999001        0.01498501
#13   Treated_2018 <-> Uninvaded_2017 1.1524646 1.1524646 3.582890 0.18296021 0.000999001        0.01498501
#14   Treated_2018 <-> Uninvaded_2018 0.9765609 0.9765609 3.164638 0.16512903 0.000999001        0.01498501
#15 Uninvaded_2017 <-> Uninvaded_2018 0.7336290 0.7336290 2.082799 0.11518124 0.000999001        0.01498501



# Emerging Invert NMDS ----------------------------------------------------


## stress plot/ see how many axes make sense

k_vec <- 1:10 #dimensions 1 - 10
stress <- numeric(length(k_vec)) # stress of each model put here
dune_dij <- metaMDSdist(invert.rel)

set.seed(25)

for(i in seq_along(k_vec)) {
  sol <- metaMDSiter(dune_dij, k = i, 
                     trace = FALSE)
  stress[i] <- sol$stress
}

plot(stress) # 3D makes sense

#### NMDS analysis 

set.seed(120) 

nms.invert <- metaMDS(invert.rel, distance = "bray", # species data, bray-curtis dissimilarity
                      autotransform = FALSE,  # NMDS will do autotransformations for you
                      k = 3, trymax = 1000)   # k = number of axes
nms.invert

#Call:
#metaMDS(comm = invert.rel, distance = "bray", k = 3, trymax = 1000,      autotransform = FALSE) 

#global Multidimensional Scaling using monoMDS

#Data:     invert.rel 
#Distance: bray 

#Dimensions: 3 
#Stress:     0.2014684 
#Stress type 1, weak ties
#Two convergent solutions found after 99 tries
#Scaling: centring, PC rotation, half-change scaling 
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
nms.invert$iters # 186

# Goodness of fit
(g <- goodness(nms.invert)) # smaller the number the better the fit
sum(g^2)
nms.invert$stress^2  # 0.0405895

1-nms.invert$stress^2 # 0.9594105 #analogous to square correlation coefficient


## extract the scores for plotting 
scr <- as.data.frame(scores(nms.invert, display = "sites")) # extract NMDS scores
colnames(scr)
 
# adding categorical info to scores
invert.env$NMDS1 <- scr$NMDS1
invert.env$NMDS2 <- scr$NMDS2
invert.env$NMDS3 <- scr$NMDS3

scores <- invert.env %>% select(Site:Habitat,NMDS1:NMDS3)

write.csv(scores,"Data/Emerging/NMDS/NMDS_emerging_inverts_NMDSscores.csv") # save this as a csv

### Vectors correlated with Axis 1 & 2 

alltaxa12 <- envfit(nms.invert, invert.rel,
                  choices = c(1,2)) #produces a list with r2, p value, and NMDS coordinates

all.taxa.df <- data.frame((alltaxa12$vectors)$arrows,
                          (alltaxa12$vectors)$r,
                          (alltaxa12$vectors)$pvals) #take list and make into data frame

write.csv(all.taxa.df, "Data/Emerging/NMDS/NMDS_emerg_vectors_axis12.csv") # save vector scores as csv


#### Vectors correlated with axis 1 & 3 

alltaxa.13 <- envfit(nms.invert, invert.rel, 
                     permutations = 999, choices = c(1,3)) 


all.taxa13.df <- data.frame((alltaxa.13$vectors)$arrows,
                            (alltaxa.13$vectors)$r,
                            (alltaxa.13$vectors)$pvals)

write.csv(all.taxa13.df, "Data/Emerging/NMDS/NMDS_emerging_vectors_axis13.csv")



## Picking the vectors we want for the figure based on fit (r > 0.2)

nmds.axis12 <- read.csv("Data/Emerging/NMDS/NMDS_emerg_vectors_axis12.csv")
nmds.axis13 <- read.csv("Data/Emerging/NMDS/NMDS_emerging_vectors_axis13.csv")

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

write.csv(corr.vectors.12, "Data/Emerging/NMDS/NMDS_emerging_correlatedvectors_axis12.csv")

# axis 1, 3
(nmds.vectors.13 <- envfit(nms.invert$points, axis13.vectors,
                           permutations = 999, choices = c(1,3)))                        


corr.vectors.13 <- as.data.frame(nmds.vectors.13$vectors$arrows*sqrt(nmds.vectors.13$vectors$r)) #scaling vectors
corr.vectors.13$species <- rownames(corr.vectors.13) # add Family as a column

write.csv(corr.vectors.13, "Data/Emerging/NMDS/NMDS_emerging_correlatedvectors_axis13.csv")

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

# NMDS Emerging Invert Figures ------------------------------------------------------------

# load the data so you don't have to redo it every time

nmds.scores <- read.csv("Data/Emerging/NMDS/NMDS_emerging_inverts_NMDSscores.csv") #points
nmds.scores$Year <- as.factor(nmds.scores$Year)


#vectors r > 0.2
vectors.12 <- read.csv("Data/Emerging/NMDS/NMDS_emerging_correlatedvectors_axis12.csv") 
vectors.13 <- read.csv("Data/Emerging/NMDS/NMDS_emerging_correlatedvectors_axis13.csv")


fill = c("Invaded_2017" = "#440C53", 
            "Invaded_2018" = "#440C53",
            "Treated_2017" = "#24908C",
            "Treated_2018" = "#24908C", 
            "Uninvaded_2017" = "#FDE825",
            "Uninvaded_2018" = "#FDE825")

colour = c("Invaded" = "#440C53",
           "Treated" = "#24908C",
           "Uninvaded" = "#FDE825")

shape = c("Invaded_2017" = 21, 
          "Invaded_2018" = 24,
          "Treated_2017" = 21,
          "Treated_2018" = 24, 
          "Uninvaded_2017" = 21,
          "Uninvaded_2018" = 24)

unique(nmds.scores$TrtYr)

## NMDS Axis 1, 2 

invert.12 <- ggplot(data = nmds.scores,
                    aes(x = NMDS1, y = NMDS2)) +
  geom_point(data = nmds.scores, 
             aes(x = NMDS1, y = NMDS2, 
                 shape = TrtYr, fill = TrtYr),
             size = 4, stroke = 1.5) +
  stat_ellipse(data = nmds.scores, 
               aes(x = NMDS1,y = NMDS2,
                   colour = Treatment), 
               size = 1, level = 0.9) +
    geom_segment(data = vectors.12, 
               aes(x = 0, xend = MDS1, y = 0, yend = MDS2),
               arrow = arrow(length = unit(0.5, "cm")),
               colour = "black") +
  geom_label_repel(data = vectors.12, 
                  aes(x = MDS1, y = MDS2, label = X),
                  color="black",
                  size = 5) +
  theme_minimal() + # no background
  theme(panel.border = element_rect(fill = NA)) + # full square around figure
  xlab("NMDS 1") +
  ylab("NMDS 2") +
  scale_fill_manual(values = fill) +
  scale_shape_manual(values = shape) +
  scale_colour_manual(values = colour) +
  theme(legend.position = "none")

invert.12

## NMDS Axis 1, 3

invert.13 <- ggplot(data = nmds.scores,
                    aes(x = NMDS1, y = NMDS3)) +
  geom_point(data = nmds.scores, 
             aes(x = NMDS1, y = NMDS3, 
                 fill = TrtYr, shape = TrtYr), 
             size = 4, stroke = 1.5) +
  stat_ellipse(data = nmds.scores, 
               aes(x = NMDS1,y = NMDS3,
                   colour = Treatment), 
               size = 1, level = 0.9) +
  geom_segment(data = vectors.13, 
               aes(x = 0, xend = MDS1, y = 0, yend = MDS3),
               arrow = arrow(length = unit(0.5, "cm")), colour = "black") +
  theme_minimal() +
  theme(panel.border = element_rect(fill = NA)) +
  geom_label_repel(data = vectors.13, 
                  aes(x = MDS1, y = MDS3, label = X),
                  color="black",
                  size = 5) +
  xlab("NMDS 1") +
  ylab("NMDS 3") +
  scale_fill_manual(values = fill) +
  scale_shape_manual(values = shape) +
  scale_colour_manual(values = colour) +
  theme(legend.position = "right")


invert.13


(NMS.emerging.panel <- ggarrange(invert.12, invert.13,
                                 common.legend = TRUE,
                                 legend = "right"))
  
nmds.emerging <- annotate_figure(NMS.emerging.panel,
                                 top = "Emerging Invertebrates")

ggsave("Figures/Emerging_NMDS_panel.TIFF", NMS.emerging.panel,
       dpi = 300,
       height = 8.19,
       width = 14.5,
       units = "in")


ggarrange(nmds.emerging,nms.aquatic,
          nrow = 2)


library(patchwork)

invert.12 <- invert.12 + ggtitle("B.    Emerging Invertebrates")

NMDS.aquatic.benthic <- benthic.12 + benthic.13 + invert.12 + invert.13 

ggsave("Figures/emerg_aquatic_NMDS_panel.TIFF", NMDS.aquatic.benthic,
       dpi = 300,
       width = 14,
       height = 10,
       units = "in")



# 2018 Emerging Invertebrates ---------------------------------------------

invert.18 <- read.csv("Data/Emerging/emerging_2018.csv") # rares and relativized in PC-Ord
str(invert.18)
colnames(invert.18)

invert.18.taxa <- invert.18 %>% select(Araneae:Crambidae)
invert.18.ev <- invert.18 %>% select(Site:ZIZPALUS)
invert.18.vectors <- invert.18 %>% select(Canopy.Height:Crambidae)


# 2018 betadisper ---------------------------------------------------------
?betadisper

invert.18 <- invert.18 %>% #rename the factors
  mutate(Treatment = fct_recode(Treatment,
                                "Treated" = "Restored")) 

invert.b <- vegdist(invert.18.taxa, method = "bray")

trt <- factor(invert.18$Treatment)

(emerg.disp <- betadisper(invert.b, trt)) 


#Average distance to median:
#  Invaded  Restored Uninvaded 
#0.5059    0.4921    0.5407 

anova(emerg.disp)

#Response: Distances
#Df   Sum Sq   Mean Sq F value Pr(>F)
#Groups     2 0.011278 0.0056391  0.8156 0.4543
#Residuals 24 0.165941 0.0069142 

permutest(emerg.disp, pairwise = TRUE, permutations = 999)


#Permutation test for homogeneity of multivariate dispersions
#Permutation: free
#Number of permutations: 999
#
#Response: Distances
#          Df   Sum Sq   Mean Sq      F N.Perm Pr(>F)
#Groups     2 0.011278 0.0056391 0.8156    999  0.444
#Residuals 24 0.165941 0.0069142                     
#
#Pairwise comparisons:
#  (Observed p-value below diagonal, permuted p-value above diagonal)
#Invaded Restored Uninvaded
#Invaded            0.73100      0.46
#Restored  0.72842               0.17
#Uninvaded 0.44344  0.16809          

boxplot(emerg.disp,
        xlab = " ")


# 2018 perMANOVA ----------------------------------------------------------


(inv.18.p <- adonis2(invert.18.taxa ~ Treatment,
                    data = invert.18,
                    permutations = 999, method = "bray"))

#Terms added sequentially (first to last)
#Permutation: free
#Number of permutations: 999
#
#adonis2(formula = invert.18.taxa ~ Treatment, data = invert.18, permutations = 999, method = "bray")
#           Df SumOfSqs      R2      F Pr(>F)    
#Treatment  2   1.7793 0.19682 2.9406  0.001 ***
#Residual  24   7.2609 0.80318                  
#Total     26   9.0402 1.00000

(adonis.pair(invert.b, trt, 
             nper = 999, corr.method = "bonferroni"))

#combination SumsOfSqs   MeanSqs  F.Model        R2 P.value P.value.corrected
#1   Invaded <-> Restored 1.2464902 1.2464902 4.354111 0.2139180   0.002             0.006
#2  Invaded <-> Uninvaded 0.4483599 0.4483599 1.418583 0.0814408   0.067             0.201
#3 Restored <-> Uninvaded 0.9740588 0.9740588 3.190785 0.1662665   0.001             0.003







# 2018 LCBD ---------------------------------------------------------------
(emerg.beta <- beta.div(invert.b, method = "percentdiff",
                        sqrt.D = FALSE, samp = FALSE,
                        nperm = 999))

emerg.beta$LCBD[emerg.beta$LCBD > mean(emerg.beta$LCBD)] # LCBD > average

#2          4          5          7          9         10         13         14         15         16 
#0.04101365 0.04390262 0.04065086 0.04627327 0.04597074 0.03988721 0.04708878 0.05299037 0.04509054 0.05575606 

#17         25 
#0.04627138 0.06554662 

invert.18.ev$LCBD <- emerg.beta$LCBD 

write.csv(invert.18.ev, "Data/Emerging/LCBD_2018.csv")

emerg.LCBD <- lm(LCBD ~ Treatment, data = invert.18.ev)
anova(emerg.LCBD)

#Response: LCBD
#Df     Sum Sq    Mean Sq F value  Pr(>F)  
#Treatment  2 0.00061547 0.00030773  2.6823 0.08884 .
#Residuals 24 0.00275342 0.00011473  

mean(invert.18.ev$LCBD) #0.037

## SCBD

invert.18.taxa.h <- decostand(invert.18.taxa, "hellinger")

(emerg.SCBD <- beta.div(invert.18.taxa.h, method = "hellinger",
                        sqrt.D = FALSE, samp = FALSE))

taxa.e.SCBD <- as.data.frame(emerg.SCBD$SCBD)
write.csv(taxa.e.SCBD, "Data/Emerging/SCBD_taxa18.csv")

## SCBD by habitat type

colnames(invert.18)
unique(invert.18$Treatment)

em.inv <- invert.18 %>% filter(Treatment == "Invaded")
em.unin <- invert.18 %>% filter(Treatment == "Uninvaded")
em.trt <- invert.18 %>% filter(Treatment == "Restored")


colnames(em.trt)

em.inv.t <- em.inv %>% select(Araneae:Crambidae)
em.inv.e <- em.inv %>% select(Site:N)

em.unin.t <- em.unin %>% select(Araneae:Crambidae)
em.unin.e <- em.unin %>% select(Site:N)

em.trt.t <- em.trt %>% select(Araneae:Crambidae)
em.trt.e <- em.trt %>% select(Site:N)


em.inv.h <- decostand(em.inv.t, "hellinger")
em.unin.h <- decostand(em.unin.t, "hellinger")
em.trt.h <- decostand(em.trt.t, "hellinger")


# Invaded

(emerg.in.SCBD <- beta.div(em.inv.h, method = "hellinger",
                        sqrt.D = FALSE, samp = FALSE))

inv.e.SCBD <- as.data.frame(emerg.in.SCBD$SCBD)
write.csv(inv.e.SCBD, "Data/Emerging/SCBD_taxa18_invaded.csv")


# Uninvaded
(emerg.unin.SCBD <- beta.div(em.unin.h, method = "hellinger",
                           sqrt.D = FALSE, samp = FALSE))

uninv.e.SCBD <- as.data.frame(emerg.unin.SCBD$SCBD)
write.csv(uninv.e.SCBD, "Data/Emerging/SCBD_taxa18_uninvaded.csv")


# Treated
(emerg.trt.SCBD <- beta.div(em.trt.h , method = "hellinger",
                             sqrt.D = FALSE, samp = FALSE))

trt.e.SCBD <- as.data.frame(emerg.trt.SCBD$SCBD)
write.csv(trt.e.SCBD , "Data/Emerging/SCBD_taxa18_treated.csv")




invert.18.ev %>% group_by(Treatment) %>% 
  summarise(meanLCBD = mean(LCBD),
            sdLCBD = sd(LCBD),
            std.error(LCBD))

# A tibble: 3 x 4
#Treatment meanLCBD  sdLCBD `std.error(LCBD)`
#1 Invaded     0.0364 0.00893           0.00298
#2 Restored    0.0432 0.00840           0.00280
#3 Uninvaded   0.0315 0.0139            0.00464


ggplot(invert.18.ev, aes(x = Treatment, y = LCBD)) +
  geom_jitter(data = invert.18.ev,
              aes(fill = Treatment, shape = Treatment),
              size = 5,
              stroke = 1.5,
              width = 0.25) +
  theme_classic(14) +
  labs(x = " ",
       y = "Emerging LCBD") +
  scale_fill_viridis(discrete = TRUE) +
  scale_shape_manual(values = c(21, 24, 22)) +
  theme(panel.border = element_rect(fill = NA)) +
  theme(legend.position = "none") +
  theme(axis.text = element_text(size = 14)) +
  ylim(0, 0.1) +
  geom_hline(yintercept = 0.037,
             linetype = "dashed")


# Beta diversity ----------------------------------------------------------
library(betapart)
citation("betapart")
citation("picante")

em.inv.t 
em.unin.t
em.trt.t 

em.inv.t[em.inv.t > 0] <- 1
em.unin.t[em.unin.t > 0] <- 1
em.trt.t[em.trt.t > 0] <- 1 

inv.bd <- as.data.frame(beta.multi(em.inv.t, index.family = "sorensen"))
inv.bd$Habitat <- c("Invaded")


unin.bd <- as.data.frame(beta.multi(em.unin.t, index.family = "sorensen"))
unin.bd$Habitat <- c("Uninvaded")

trt.bd <- as.data.frame(beta.multi(em.trt.t, index.family = "sorensen"))
trt.bd$Habitat <- c("Treated")

beta.em <- rbind(trt.bd, inv.bd, unin.bd)

write.csv(beta.em, "Data/Emerging/beta_diversity_em.csv")

### Null model
em.inv.t
em.unin.t
em.trt.t


## Invaded null model

library(picante)

invaded.null <- data.frame(matrix(as.numeric(0), ncol=(3), nrow=(999)))
colnames(invaded.null) <- c("Turnover","Nestedness","Sum")


for (i in 1:999){ #for 1000 iterations
  tempm <- as.data.frame(randomizeMatrix(em.inv.t, null.model = "frequency")) #Randomize the data
  inv.bv <- beta.multi(tempm, index.family = "sorensen")
  invaded.null[i,] <- data.frame(matrix(unlist(inv.bv), nrow = length(1), byrow = T)) 
}

inv.null <- invaded.null %>% 
  summarise(n = n(),
            N.avg = mean(Nestedness),
            N.sd = sd(Nestedness),
            N.CI = qnorm(0.95)*(N.sd/sqrt(9)), # 95% CI
            T.avg = mean(Turnover),
            T.sd = sd(Turnover),
            T.CI = qnorm(0.95)*(T.sd/sqrt(9)), # specify correct sample size
            S.avg = mean(Sum),
            S.sd = sd(Sum),
            S.CI = qnorm(0.95)*(S.sd/sqrt(9)))


inv.null$Turn <- inv.bd$beta.SIM
inv.null$Nest <- inv.bd$beta.SNE
inv.null$Over <- inv.bd$beta.SOR
inv.null$Habitat <- c("Invaded")


# Uninvaded 
uninvaded.null <- data.frame(matrix(as.numeric(0), ncol=(3), nrow=(999)))
colnames(uninvaded.null) <- c("Turnover","Nestedness","Sum")


for (i in 1:999){ #for 1000 iterations
  tempu <- as.data.frame(randomizeMatrix(em.unin.t, null.model = "frequency")) #Randomize the data
  unin.bv <- beta.multi(tempu, index.family = "sorensen")
  uninvaded.null[i,] <- data.frame(matrix(unlist(unin.bv), nrow = length(1), byrow = T)) 
}

unin.null <- uninvaded.null %>% 
  summarise(n = n(),
            N.avg = mean(Nestedness),
            N.sd = sd(Nestedness),
            N.CI = qnorm(0.95)*(N.sd/sqrt(9)), # 95% CI
            T.avg = mean(Turnover),
            T.sd = sd(Turnover),
            T.CI = qnorm(0.95)*(T.sd/sqrt(9)), # specify correct sample size
            S.avg = mean(Sum),
            S.sd = sd(Sum),
            S.CI = qnorm(0.95)*(S.sd/sqrt(9)))


unin.null$Turn <- unin.bd$beta.SIM
unin.null$Nest <- unin.bd$beta.SNE
unin.null$Over <- unin.bd$beta.SOR
unin.null$Habitat <- c("Uninvaded")

# Treated

treated.null <- data.frame(matrix(as.numeric(0), ncol=(3), nrow=(999)))
colnames(treated.null) <- c("Turnover","Nestedness","Sum")


for (i in 1:999){ #for 1000 iterations
  tempt <- as.data.frame(randomizeMatrix(em.trt.t, null.model = "frequency")) #Randomize the data
  trt.bv <- beta.multi(tempt, index.family = "sorensen")
  treated.null[i,] <- data.frame(matrix(unlist(trt.bv), nrow = length(1), byrow = T)) 
}

treat.null <- treated.null %>% 
  summarise(n = n(),
            N.avg = mean(Nestedness),
            N.sd = sd(Nestedness),
            N.CI = qnorm(0.95)*(N.sd/sqrt(9)), # 95% CI
            T.avg = mean(Turnover),
            T.sd = sd(Turnover),
            T.CI = qnorm(0.95)*(T.sd/sqrt(9)), # specify correct sample size
            S.avg = mean(Sum),
            S.sd = sd(Sum),
            S.CI = qnorm(0.95)*(S.sd/sqrt(9)))


treat.null$Turn <- trt.bd$beta.SIM
treat.null$Nest <- trt.bd$beta.SNE
treat.null$Over <- trt.bd$beta.SOR
treat.null$Habitat <- c("Treated")


bd.null.em <- rbind(treat.null, unin.null,inv.null)

write.csv(bd.null.em, "Data/Emerging/beta_diversity_null_true.csv")


sumbd <- ggplot(bd.null.em, aes(x = Habitat, y = Over, 
                                fill = Habitat, shape = Habitat)) +
  geom_errorbar(aes(ymin = S.avg - S.CI, ymax = S.avg + S.CI),
                colour = "black",
                size = 0.5,
                width = 0.3,
                position = position_dodge(0.6)) +
    geom_point(size = 7,
               alpha = 0.7,
               stroke = 1.5) +
  geom_point(aes(y = S.avg),
             fill = "black",
             size = 3) +
  labs(x = " ",
       y = expression(paste("Beta Diversity Sum"))) +
  theme_classic(14)+
  theme(legend.position = "none") +
  scale_shape_manual(values = c(21, 24, 22)) +
  scale_fill_viridis(discrete = TRUE) +
  ylim(0.5, 0.75) +
  theme(axis.text = element_text(size = 14))
  



nest <- ggplot(bd.null.em, aes(x = Habitat, y = Nest, 
                           fill = Habitat, shape = Habitat)) +
  geom_errorbar(aes(ymin = N.avg - N.CI, ymax = N.avg + N.CI),
                colour = "black",
                size = 0.5,
                width = 0.3) +
  geom_point(alpha = 0.7, size = 7,
             stroke = 1.5) +
  geom_point(aes(y = N.avg),
             fill = "black",
             size = 3) +
  labs(x = " ",
       y = expression(paste("Nestedness"))) +
  theme_classic(14)+
  theme(legend.position = "none") +
  scale_shape_manual(values = c(21, 24, 22)) +
  scale_fill_viridis(discrete = TRUE) +
  ylim(0.00, 0.2)+
  theme(axis.text = element_text(size = 14))



turn <- ggplot(bd.null.em, aes(x = Habitat, y = Turn, 
                            fill = Habitat, shape = Habitat)) +
  geom_errorbar(aes(ymin = T.avg - T.CI, ymax = T.avg + T.CI),
                colour = "black",
                size = 0.5,
                width = 0.3) +
  geom_point(alpha = 0.7, size = 7,
             stroke = 1.5) +
  geom_point(aes(y = T.avg),
             fill = "black",
             size = 3) +
  labs(x = " ",
       y = expression(paste("Turnover"))) +
  theme_classic(14) +
  theme(legend.position = "none") +
  scale_shape_manual(values = c(21, 24, 22)) +
  scale_fill_viridis(discrete = TRUE) +
  ylim(0.3, 0.75) +
  theme(axis.text = element_text(size = 14))



emerg.beta <- ggarrange(sumbd, nest, turn,
                        labels = c("D","E","F"),
                        nrow = 1)

library(patchwork)

aqem.beta <- aqua.beta + emerg.beta +
  plot_layout(nrow = 2)

ggsave("Figures/aquaitc_emerging_betadiversity.jpeg",
       aqem.beta)


# 2018 NMDS ---------------------------------------------------------------


## stress plot/ see how many axes make sense

k_vec <- 1:10 #dimensions 1 - 10
stress <- numeric(length(k_vec)) # stress of each model put here
dune_dij <- metaMDSdist(invert.18.taxa)

set.seed(25)

for(i in seq_along(k_vec)) {
  sol <- metaMDSiter(dune_dij, k = i, 
                     trace = FALSE)
  stress[i] <- sol$stress
}

plot(stress) # 3D makes sense

#### NMDS analysis 

set.seed(120) 

nms.invert18 <- metaMDS(invert.18.taxa, distance = "bray", # species data, bray-curtis dissimilarity
                      autotransform = FALSE,  # NMDS will do autotransformations for you
                      k = 3, trymax = 1000)   # k = number of axes
nms.invert18

#global Multidimensional Scaling using monoMDS
#
#Data:     invert.18.taxa 
#Distance: bray 
#
#Dimensions: 3 
#Stress:     0.159081 
#Stress type 1, weak ties
#Two convergent solutions found after 20 tries

# look at points and the stress real quick
layout(matrix(1:2, ncol = 2))
plot(nms.invert18, main = "Invertebrate NMDS plot"); stressplot(nms.invert18, main = "Shepard plot")
layout(1)

ordiplot(nms.invert18, type = "n")
orditorp(nms.invert18, display = "species")
orditorp(nms.invert18, display = "sites")


ordiplot(nms.invert18, type = "n", choices = c(1,3))
orditorp(nms.invert18, display = "species", choices = c(1,3))
orditorp(nms.invert18, display = "sites", choices = c(1,3))


# how many iterations of the NMDS
nms.invert18$iters # 81

# Goodness of fit
(g <- goodness(nms.invert18)) # smaller the number the better the fit
sum(g^2)
nms.invert18$stress^2  # 0.02530676

1-nms.invert18$stress^2 # 0.9746932 #analogous to square correlation coefficient


## extract the scores for plotting 
scr18 <- as.data.frame(scores(nms.invert18, display = "sites")) # extract NMDS scores
colnames(scr18)

# adding categorical info to scores
invert.18.ev$NMDS1 <- scr18$NMDS1
invert.18.ev$NMDS2 <- scr18$NMDS2
invert.18.ev$NMDS3 <- scr18$NMDS3

scores18 <- invert.18.ev %>% select(Site,Treatment,NMDS1:NMDS3)

write.csv(scores18,"Data/Emerging/NMDS/NMDS_emerging_2018_NMDSscores.csv") # save this as a csv

### Vectors correlated with Axis 1 & 2 

alltaxa12.2018 <- envfit(nms.invert18, invert.18.vectors,
                    choices = c(1,2),
                    na.rm = TRUE) #produces a list with r2, p value, and NMDS coordinates

all.taxa.df.2018 <- data.frame((alltaxa12.2018$vectors)$arrows,
                          (alltaxa12.2018$vectors)$r,
                          (alltaxa12.2018$vectors)$pvals) #take list and make into data frame

write.csv(all.taxa.df.2018, "Data/Emerging/NMDS/NMDS_emerg_vectors_2018_axis12.csv") # save vector scores as csv


#### Vectors correlated with axis 1 & 3 

alltaxa13.2018 <- envfit(nms.invert18, invert.18.vectors, 
                     permutations = 999, choices = c(1,3),
                     na.rm = TRUE) 


all.taxa13.df.2018 <- data.frame((alltaxa13.2018$vectors)$arrows,
                            (alltaxa13.2018$vectors)$r,
                            (alltaxa13.2018$vectors)$pvals)

write.csv(all.taxa13.df.2018, "Data/Emerging/NMDS/NMDS_emerging_vectors_2018_axis13.csv")


## Picking the vectors we want for the figure based on fit (r > 0.2)

nmds.axis1218 <- read.csv("Data/Emerging/NMDS/NMDS_emerg_vectors_2018_axis12.csv")
nmds.axis1318 <- read.csv("Data/Emerging/NMDS/NMDS_emerging_vectors_2018_axis13.csv")

colnames(nmds.axis1218)

# axis 1, 2
corr.sp.18 <- nmds.axis1218 %>% filter(X.alltaxa12.2018.vectors..r > 0.3) 
target12.18 <- corr.sp.18$X # string of the Family names

axis12.vectors.18 <- invert.18.vectors %>% select(all_of(target12.18)) # make a matrix of just those

# axis 1, 3

corr.sp13.18 <- nmds.axis1318 %>% filter(X.alltaxa13.2018.vectors..r > 0.3) 
target13.18 <- corr.sp13.18$X # string of the Family names

axis13.vectors.18 <- invert.18.vectors %>% select(all_of(target13.18)) # make a matrix of just those


# fit them to the nms
# axis 1, 2
(nmds.vectors.12.18 <- envfit(nms.invert18$points, axis12.vectors.18,
                           permutations = 999, choices = c(1,2),
                           na.rm = TRUE))                        


corr.vectors.12.18 <- as.data.frame(nmds.vectors.12.18 $vectors$arrows*sqrt(nmds.vectors.12.18$vectors$r)) #scaling vectors
corr.vectors.12.18$taxa <- rownames(corr.vectors.12.18) # add Family as a column

write.csv(corr.vectors.12.18, "Data/Emerging/NMDS/NMDS_emerging_correlatedvectors_axis12_2018.csv")

# axis 1, 3
(nmds.vectors.13.18 <- envfit(nms.invert18$points, axis13.vectors.18,
                           permutations = 999, choices = c(1,3),
                           na.rm = TRUE))                        


corr.vectors.13.18 <- as.data.frame(nmds.vectors.13.18$vectors$arrows*sqrt(nmds.vectors.13.18$vectors$r)) #scaling vectors
corr.vectors.13.18$taxa <- rownames(corr.vectors.13.18) # add Family as a column

write.csv(corr.vectors.13.18, "Data/Emerging/NMDS/NMDS_emerging_correlatedvectors_axis13_2018.csv")


# 2018 NMDS Figure --------------------------------------------------------

nmds.scores18 <- read.csv("Data/Emerging/NMDS/NMDS_emerging_2018_NMDSscores.csv") #points

unique(nmds.scores18$Treatment)

nmds.scores18 <- nmds.scores18 %>% #rename the factors
  mutate(Treatment = fct_recode(Treatment,
                              "Treated" = "Restored")) 


#vectors r > 0.2
vectors.1218 <- read.csv("Data/Emerging/NMDS/NMDS_emerging_correlatedvectors_axis12_2018.csv") 
vectors.1318 <- read.csv("Data/Emerging/NMDS/NMDS_emerging_correlatedvectors_axis13_2018.csv")


fill = c("Invaded" = "#440C53",
         "Treated" = "#24908C",
         "Uninvaded" = "#FDE825")

colour = c("Invaded" = "#440C53",
           "Treated" = "#24908C",
           "Uninvaded" = "#FDE825")

shape = c("Invaded" = 21,
          "Treated" = 24,
          "Uninvaded" = 22)


## NMDS Axis 1, 2 

invert.1218 <- ggplot(data = nmds.scores18,
                    aes(x = NMDS1, y = NMDS2)) +
  geom_point(data = nmds.scores18, 
             aes(x = NMDS1, y = NMDS2, 
                 shape = Treatment, fill = Treatment),
             size = 7, stroke = 1.5,
             alpha = 0.7) +
  stat_ellipse(data = nmds.scores18, 
               aes(x = NMDS1,y = NMDS2,
                   colour = Treatment), 
               size = 1, level = 0.9) +
  geom_segment(data = vectors.1218, 
               aes(x = 0, xend = MDS1, y = 0, yend = MDS2),
               arrow = arrow(length = unit(0.5, "cm")),
               colour = "black") +
  geom_label_repel(data = vectors.1218, 
                   aes(x = MDS1, y = MDS2, label = taxa),
                   color = "black",
                   size = 5,
                   force = 2) +
  theme_minimal() + # no background
  theme(panel.border = element_rect(fill = NA)) + # full square around figure
  xlab("NMDS 1") +
  ylab("NMDS 2") +
  scale_fill_manual(values = fill) +
  scale_shape_manual(values = shape) +
  scale_colour_manual(values = colour) +
  theme(legend.position = "none")

invert.1218

## NMDS Axis 1, 3

invert.1318 <- ggplot(data = nmds.scores18,
                    aes(x = NMDS1, y = NMDS3)) +
  geom_point(data = nmds.scores18, 
             aes(x = NMDS1, y = NMDS3, 
                 fill = Treatment, shape = Treatment), 
             size = 7, stroke = 1.5,
             alpha = 0.7) +
  stat_ellipse(data = nmds.scores18, 
               aes(x = NMDS1,y = NMDS3,
                   colour = Treatment), 
               size = 1, level = 0.9) +
  geom_segment(data = vectors.1318, 
               aes(x = 0, xend = MDS1, y = 0, yend = MDS3),
               arrow = arrow(length = unit(0.5, "cm")), colour = "black") +
  theme_minimal() +
  theme(panel.border = element_rect(fill = NA)) +
  geom_label_repel(data = vectors.1318, 
                   aes(x = MDS1, y = MDS3, label = taxa),
                   color="black",
                   size = 5,
                   force = 2) +
  xlab("NMDS 1") +
  ylab("NMDS 3") +
  scale_fill_manual(values = fill) +
  scale_shape_manual(values = shape) +
  scale_colour_manual(values = colour) +
  theme(legend.position = "bottom")


invert.1318


(NMS.emerging.panel18 <- ggarrange(invert.1218, invert.1318,
                                 common.legend = TRUE,
                                 legend = "none"))

nmds.emerging <- annotate_figure(NMS.emerging.panel18,
                                 top = "Emerging Invertebrates")

ggsave("Figures/Emerging_NMDS_panel.TIFF", NMS.emerging.panel,
       dpi = 300,
       height = 8.19,
       width = 14.5,
       units = "in")


ggarrange(nmds.emerging,nms.aquatic,
          nrow = 2)

library(patchwork)

invert.1218 <- invert.1218 + ggtitle("B.    Emerging Invertebrates")


(nmds.panel <- benthic.12 + benthic.13 + invert.1218 + invert.1318)


ggsave("Figures/emergaquat_NMDS_panel.jpeg",
       nmds.panel,
       dpi = 300,
       height = 10,
       width = 13.2,
       units = "in")





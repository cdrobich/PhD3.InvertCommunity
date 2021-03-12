
# Load packages -----------------------------------------------------------


library(tidyverse) # for figures, pipes etc.
library(vegan) # for NMDS
library(ggpubr) # for arranging figures

library(car)
library(viridis) # colours


library(ggrepel)


library(devtools)
install_github("GuillemSalazar/EcolUtils")

library(EcolUtils)

# Load Data ---------------------------------------------------------------

benthic <- read.csv("Data/Aquatic/benthic_inverts_relativized.csv") # occurrences <=2 removed, relativized

unique(benthic$Habitat)
colnames(benthic)

benthic.data <- benthic %>% select(Oligochaetae:Leptoceridae)
benthic.env <- benthic %>% select(ID:Collection.date)


# Multivariate Analyses ---------------------------------------------------

# perMANOVA ---------------------------------------------------------------

#### perMANOVA #####

(per.treat <- adonis2(benthic.data ~ Habitat, data = benthic,
                     permutations = 999, method = "bray"))


#         Df SumOfSqs      R2      F Pr(>F)    
#Habitat   2   2.2329 0.28701 4.4279  0.001 ***
#Residual 22   5.5471 0.71299                  
#Total    24   7.7800 1.00000                  


# homogeneity of groups dispersion

bugs.b <- vegdist(benthic.data, method = "bray")

groupsb <- factor(benthic$Habitat)

(dispersionb <- betadisper(bugs.b, groupsb)) # spatial median default

#Homogeneity of multivariate dispersions
##
##Call: betadisper(d = bugs.b, group = groupsb)
##
##No. of Positive Eigenvalues: 19
##No. of Negative Eigenvalues: 5
##
##Average distance to median:
##  Invaded   Treated Uninvaded 
##0.4468    0.4541    0.4860 
##
##Eigenvalues for PCoA axes:
##  (Showing 8 of 24 eigenvalues)
##PCoA1  PCoA2  PCoA3  PCoA4  PCoA5  PCoA6  PCoA7  PCoA8 
#1.5648 1.2761 1.0066 0.7229 0.5963 0.5049 0.4730 0.4529 

anova(dispersionb)

#Analysis of Variance Table
#
#Response: Distances
#Df   Sum Sq   Mean Sq F value Pr(>F)
#Groups     2 0.007034 0.0035169  0.3297 0.7227
#Residuals 22 0.234694 0.0106679 

plot(dispersionb)
boxplot(dispersionb) # actually look really good!


(adonis.pair(bugs.b, groupsb, nper = 1000, corr.method = "bonferroni"))

#              combination SumsOfSqs   MeanSqs  F.Model        R2     P.value P.value.corrected
#1   Invaded <-> Restored 0.9649239 0.9649239 3.946800 0.2083096 0.000999001       0.002997003
#2  Invaded <-> Uninvaded 1.0586760 1.0586760 4.140050 0.2282270 0.000999001       0.002997003
#3 Restored <-> Uninvaded 1.3160709 1.3160709 4.999042 0.2499641 0.000999001       0.002997003


# NMDS  -------------------------------------------------------------------

## stress plot/ see how many axes make sense

k_vec <- 1:10 #dimensions 1 - 10
stress <- numeric(length(k_vec)) # stress of each model put here
dune_dij <- metaMDSdist(benthic.data)
set.seed(25)
for(i in seq_along(k_vec)) {
  sol <- metaMDSiter(dune_dij, k = i, # you can tell I lifted this from a tutorial on the dune package
                     trace = FALSE)
  stress[i] <- sol$stress
}

plot(stress) # 3D makes sense


# NMDS analysis code ------------------------------------------------------

set.seed(120) 

nms.invert <- metaMDS(benthic.data, distance = "bray", # species data, bray-curtis dissimilarity
                      autotransform = FALSE,  # NMDS will do autotransformations for you
                      k = 3, trymax = 1000)   # k = number of axes
nms.invert

#global Multidimensional Scaling using monoMDS

#Data:     benthic.data 
#Distance: bray 

#Dimensions: 3 
#Stress:     0.1450098 
#Stress type 1, weak ties
#Two convergent solutions found after 28 tries
#Scaling: centring, PC rotation, halfchange scal

# look at points and the stress real quick
layout(matrix(1:2, ncol = 2))
plot(nms.invert, main = "Invertebrate NMDS plot"); stressplot(nms.invert, main = "Shepard plot")
layout(1)

ordiplot(nms.invert, type = "n")
orditorp(nms.invert, display = "species")
orditorp(nms.invert, display = "sites")

# how many iterations of the NMDS
nms.invert$iters # 66

# Goodness of fit
(g <- goodness(nms.invert)) # smaller the number the better the fit
sum(g^2)
nms.invert$stress^2  # 0.02102784

1-nms.invert$stress^2 # 0.9789722 #analogous to square correlation coefficient


# Scores and vectors for plotting -----------------------------------------

scr <- as.data.frame(scores(nms.invert, display = "sites")) # extract NMDS scores

# adding categorical info to scores
colnames(scr)

benthic.env$NMDS1 <- scr$NMDS1
benthic.env$NMDS2 <- scr$NMDS2
benthic.env$NMDS3 <- scr$NMDS3

scores <- benthic.env

write.csv(scores,"Data/Aquatic/NMDS/NMDS_benthic_inverts_scores.csv") # save this as a csv


### Vectors correlated with Axis 1 & 2 

alltaxa <- envfit(nms.invert, benthic.data,
                 choices = c(1,2)) #produces a list with r2, p value, and NMDS coordinates

all.taxa.df <- data.frame((alltaxa$vectors)$arrows,
                          (alltaxa$vectors)$r,
                          (alltaxa$vectors)$pvals) #take list and make into dataframe


corr.spp12 <- all.taxa.df %>% filter(X.alltaxa.vectors..r > 0.25)
corr.spp12$species <- rownames(corr.spp12)

corr.species12 <- corr.spp12$species # string of the Family names

axis12.vectors <- benthic.data %>% select(all_of(corr.species12)) # make a matrix of just those


corrtaxa12 <- envfit(nms.invert$points, axis12.vectors, 
                   permutations = 999, choices = c(1,2))
corrtaxa12


species.12 <- as.data.frame(corrtaxa12$vectors$arrows*sqrt(corrtaxa12$vectors$r)) #scaling vectors so they correspond with r2
species.12$species <- rownames(species.12)


write.csv(species.12, "Data/Aquatic/NMDS/NMDS_benthic_vectors_axis12.csv") # save vector scores as csv


#### Vectors correlated with axis 1 & 3 

alltaxa.13 <- envfit(nms.invert, benthic.data, 
                     permutations = 999, choices = c(1,3)) 

all.taxa.df.13 <- data.frame((alltaxa.13$vectors)$arrows,
                          (alltaxa.13$vectors)$r,
                          (alltaxa.13$vectors)$pvals)

corr.spp13 <- all.taxa.df.13 %>% filter(X.alltaxa.13.vectors..r > 0.25)
corr.spp13$species <- rownames(corr.spp13)

corr.species13 <- corr.spp13$species # string of the Family names

axis13.vectors <- benthic.data %>% select(all_of(corr.species13)) # make a matrix of just those


corrtaxa13 <- envfit(nms.invert$points, axis13.vectors, 
                     permutations = 999, choices = c(1,3))
corrtaxa13


species.13 <- as.data.frame(corrtaxa13$vectors$arrows*sqrt(corrtaxa13$vectors$r)) #scaling vectors so they correspond with r2
species.13$species <- rownames(species.13)

write.csv(species.13, "Data/Aquatic/NMDS/NMDS_benthic_vectors_axis13.csv")


# Base R plots ------------------------------------------------------------

# I don't like base R but it has some nice things built in for ordinations

# axis 1, 2 

col_vec <- c("#9970ab", "#1b7837", "#2166ac")

ordiplot(nms.invert, choices = c(1,2), 
         type = "points",
         display = "sites")


ordihull(nms.invert, groups = benthic.env$Habitat, # ellipse hull
         col = col_vec)

ordiellipse(nms.invert, groups = benthic.env$Habitat, # st error of centroid ellipse hull
            draw = "polygon",
            col = col_vec)

plot(alltaxa, p.max = 0.013, col = "black")



# NMDS Figures ------------------------------------------------------------

# load data from above 

benth.scores <- read.csv("Data/Aquatic/NMDS/NMDS_benthic_inverts_scores.csv")
vector.12 <- read.csv("Data/Aquatic/NMDS/NMDS_benthic_vectors_axis12.csv")
vector.13 <- read.csv("Data/Aquatic/NMDS/NMDS_benthic_vectors_axis13.csv")


str(benth.scores)

## NMDS Axis 1, 2 

benthic.12 <- ggplot(data = benth.scores,
                    aes(x = NMDS1, y = NMDS2)) +
  geom_point(data = benth.scores, 
             aes(x = NMDS1, 
                 y = NMDS2, 
                 fill = Habitat, 
                 shape = Habitat,
                 stroke = 1.5),
             size = 5) + # sites as points
  stat_ellipse(data = benth.scores, aes(x = NMDS1,y = NMDS2,
                                        colour = Habitat), 
               size = 1, level = 0.90) + # a 95% CI ellipses
  geom_segment(data = vector.12, aes(x = 0, xend = MDS1, y = 0, yend = MDS2), # adding in the vectors, c
               arrow = arrow(length = unit(0.5, "cm")), colour = "black") + # can add in geom_label or geom_text for labels
  theme_minimal() + # no background
  theme(panel.border = element_rect(fill = NA)) + # full square around figure
  xlab("NMDS 1") +
  ylab("NMDS 2") +
  theme(legend.position = "none") +
  geom_label_repel(data = vector.12, 
                  aes(x = MDS1, y = MDS2, label = species),
                  color="black",
                  size = 5) +
  scale_fill_viridis(discrete = TRUE) +
  scale_colour_viridis(discrete = TRUE) +
  scale_shape_manual(values = c(21, 24, 22)) 

benthic.12

## NMDS Axis 1, 3

benthic.13 <- ggplot(data = benth.scores,
                    aes(x = NMDS1, y = NMDS3)) +
  geom_point(data = benth.scores, 
             aes(x = NMDS1, y = NMDS3, 
                 fill = Habitat, shape = Habitat), 
             size = 5,
             stroke = 1.5) +
  stat_ellipse(data = benth.scores, 
               aes(x = NMDS1,y = NMDS3,
                   colour = Habitat), size = 1, level = 0.9) +
  geom_segment(data = vector.13, 
               aes(x = 0, xend = MDS1, y = 0, yend = MDS3),
               arrow = arrow(length = unit(0.5, "cm")), colour = "black") +
  theme_minimal() +
  theme(panel.border = element_rect(fill = NA)) +
  xlab("NMDS 1") +
  ylab("NMDS 3") +
  theme(legend.position = "none") +
  geom_label_repel(data = vector.13, 
                  aes(x = MDS1, y = MDS3, label = species),
                  color="black",
                  size = 5) +
  scale_fill_viridis(discrete = TRUE) +
  scale_colour_viridis(discrete = TRUE) +
  scale_shape_manual(values = c(21, 24, 22)) 

benthic.13

# putting both figures together

benthic.12 <- benthic.12 + ggtitle("A.    Aquatic Invertebrates")


(NMS.benthic.panel <- ggarrange(benthic.12, benthic.13,
                               align = "hv",
                               common.legend = TRUE,
                               legend = "none"))



ggsave("Figures/Benthic_NMDS_panel.jpeg", NMS.benthic.panel,
       dpi = 300,
       height = 7.33,
       width = 11.9,
       units = "in")


# Beta Diversity ----------------------------------------------------------
library(adespatial)

benthic <- read.csv("Data/Aquatic/benthic_inverts_relativized.csv") # occurrences <=2 removed, relativized

# Invaded data
aq.invaded <- benthic %>% filter(Habitat == "Invaded")
aq.inv.taxa <- aq.invaded %>% select(Oligochaetae:Leptoceridae)
aq.inv.env <- aq.invaded %>% select(ID:Collection.date)

# Uninvaded data
aq.uninvaded <- benthic %>% filter(Habitat == "Uninvaded")
aq.uninv.taxa <- aq.uninvaded %>% select(Oligochaetae:Leptoceridae)
aq.uninv.env <- aq.uninvaded %>% select(ID:Collection.date)

# Treated data
aq.treated <- benthic %>% filter(Habitat == "Treated")
aq.trt.taxa <- aq.treated %>% select(Oligochaetae:Leptoceridae)
aq.trt.env <- aq.treated %>% select(ID:Collection.date)


# Dissimilarity Matrices --------------------------------------------------
colnames(benthic)
benthic.tx <- benthic %>% select(Oligochaetae:Leptoceridae)
bent.env <- benthic %>% select(ID:Collection.date)


benthic.bc <- vegdist(benthic.tx, method = "bray")
inv.bc <- vegdist(aq.inv.taxa, method = "bray")
unin.bc <- vegdist(aq.uninv.taxa, method = "bray")
trt.bc <- vegdist(aq.trt.taxa, method = "bray")


# Local Contribution to BD

(benth.beta <- beta.div(benthic.bc, method = "percentdiff",
                      sqrt.D = FALSE, samp = FALSE,
                      nperm = 999))

benth.beta$LCBD[benth.beta$LCBD > mean(benth.beta$LCBD)] # LCBD > average

#6          7          8         10         16         17         18         19 
#0.04838159 0.04357041 0.04138004 0.05279625 0.04708411 0.04976721 0.07507751 0.04368041 

#20         21 
#0.04246180 0.05326194


bent.env$LCBD <- benth.beta$LCBD
bent.env$LCBD.p <- benth.beta$p.LCBD
write.csv(bent.env, "Data/Aquatic/LCBD_data.csv")

LCBD.benth <- lm(LCBD ~ Habitat, data = bent.env)
anova(LCBD.benth)

#Analysis of Variance Table
#
#Response: LCBD
#Df     Sum Sq    Mean Sq F value Pr(>F)
#Habitat    2 0.00019935 9.9676e-05  0.9043 0.4194
#Residuals 22 0.00242494 1.1022e-04  

mean(bent.env$LCBD) # 0.04

bent.env %>% group_by(Habitat) %>% 
  summarise(meanLCBD = mean(LCBD),
            sdLCBD = sd(LCBD),
            std.error(LCBD))

#  Habitat   meanLCBD  sdLCBD `std.error(LCBD)`
#1 Invaded     0.0365 0.00777           0.00275
#2 Treated     0.0400 0.00775           0.00258
#3 Uninvaded   0.0435 0.0147            0.00521

ggplot(bent.env, aes(x = Habitat, y = LCBD)) +
  geom_jitter(data = bent.env,
              aes(fill = Habitat, shape = Habitat),
              size = 5,
              stroke = 1.5,
              width = 0.25) +
  theme_classic(14) +
  labs(x = " ",
       y = "Aquatic LCBD") +
  scale_fill_viridis(discrete = TRUE) +
  scale_shape_manual(values = c(21, 24, 22)) +
  theme(panel.border = element_rect(fill = NA)) +
  theme(legend.position = "none") +
  theme(axis.text = element_text(size = 14)) +
  ylim(0, 0.1) +
  geom_hline(yintercept = 0.04,
             linetype = "dashed")




## Each habitat alone

(inv.beta <- beta.div(inv.bc, method = "percentdiff",
                   sqrt.D = FALSE, samp = FALSE,
                   nperm = 999))


aq.inv.env$LCBD <- inv.beta$LCBD

inv.beta$LCBD[inv.beta$LCBD > mean(inv.beta$LCBD)] # LCBD > average

#1         2         3         7 
#0.1635096 0.1573112 0.1288408 0.1386385 


(unin.beta <- beta.div(unin.bc, method = "percentdiff",
                      sqrt.D = FALSE, samp = FALSE,
                      nperm = 999))


aq.uninv.env$LCBD <- unin.beta$LCBD

unin.beta$LCBD[unin.beta$LCBD > mean(unin.beta$LCBD)] # LCBD > average

#1         3         5         6         7 
#0.1388476 0.1303577 0.1395451 0.1416225 0.1307846 


(trt.beta <- beta.div(trt.bc, method = "percentdiff",
                       sqrt.D = FALSE, samp = FALSE,
                       nperm = 999))


aq.trt.env$LCBD <- trt.beta$LCBD

trt.beta$LCBD[trt.beta$LCBD > mean(trt.beta$LCBD)] # LCBD > average

#9 
#0.3740378 

LCBD.data <- rbind(aq.trt.env, aq.inv.env, aq.uninv.env)


LCBD.anova <- lm(LCBD ~ Habitat, data = LCBD.data)
Anova(LCBD.anova)

#Anova Table (Type II tests)

#Response: LCBD
#Sum Sq Df F value Pr(>F)
#Habitat   0.001111  2  0.1422 0.8683
#Residuals 0.085978 22


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
benthic.vector <- benthic %>% select(Water_Depth:Leptoceridae)



bth.tx <- as.data.frame(sapply(benthic.data, function(col) sum(col))) # how many unique values in column

# Multivariate Analyses ---------------------------------------------------

# perMANOVA ---------------------------------------------------------------

(per.treat <- adonis2(benthic.data ~ Habitat, data = benthic,
                     permutations = 999, method = "bray"))


#         Df SumOfSqs      R2      F Pr(>F)    
#Habitat   2   2.2329 0.28701 4.4279  0.001 ***
#Residual 22   5.5471 0.71299                  
#Total    24   7.7800 1.00000                  


### Betadisper; homogeneity of groups dispersion ####

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
##Invaded   Treated Uninvaded 
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


permutest(dispersionb, pairwise = TRUE, permutations = 999)

#Permutation test for homogeneity of multivariate dispersions
#Permutation: free
#Number of permutations: 999
#
#Response: Distances
#           Df   Sum Sq   Mean Sq      F N.Perm Pr(>F)
#Groups     2 0.007034 0.0035169 0.3297    999  0.766
#Residuals 22 0.234694 0.0106679                     
#
#Pairwise comparisons:
#  (Observed p-value below diagonal, permuted p-value above diagonal)
#Invaded Treated Uninvaded
#Invaded           0.92000     0.313
#Treated   0.90261             0.574
#Uninvaded 0.30650 0.55309  

plot(dispersionb)
boxplot(dispersionb,
        xlab = " ") # actually look really good!


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

alltaxa <- envfit(nms.invert, benthic.vector,
                 choices = c(1,2)) #produces a list with r2, p value, and NMDS coordinates

all.taxa.df <- data.frame((alltaxa$vectors)$arrows,
                          (alltaxa$vectors)$r,
                          (alltaxa$vectors)$pvals) #take list and make into dataframe


corr.spp12 <- all.taxa.df %>% filter(X.alltaxa.vectors..r > 0.3)
corr.spp12$species <- rownames(corr.spp12)

corr.species12 <- corr.spp12$species # string of the Family names

axis12.vectors <- benthic.vector %>% select(all_of(corr.species12)) # make a matrix of just those


corrtaxa12 <- envfit(nms.invert$points, axis12.vectors, 
                   permutations = 999, choices = c(1,2))
corrtaxa12


species.12 <- as.data.frame(corrtaxa12$vectors$arrows*sqrt(corrtaxa12$vectors$r)) #scaling vectors so they correspond with r2
species.12$species <- rownames(species.12)


write.csv(species.12, "Data/Aquatic/NMDS/NMDS_benthic_vectors_axis12.csv") # save vector scores as csv


#### Vectors correlated with axis 1 & 3 

alltaxa.13 <- envfit(nms.invert, benthic.vector, 
                     permutations = 999, choices = c(1,3)) 

all.taxa.df.13 <- data.frame((alltaxa.13$vectors)$arrows,
                          (alltaxa.13$vectors)$r,
                          (alltaxa.13$vectors)$pvals)

corr.spp13 <- all.taxa.df.13 %>% filter(X.alltaxa.13.vectors..r > 0.3)
corr.spp13$species <- rownames(corr.spp13)

corr.species13 <- corr.spp13$species # string of the Family names

axis13.vectors <- benthic.vector %>% select(all_of(corr.species13)) # make a matrix of just those


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

fig <- ordiplot(nms.invert, choices = c(1,2), 
         type = "points",
         display = "sites")

text(fig, "sites", col="blue", cex=0.9)

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



fill = c("Invaded" = "#440C53",
         "Treated" = "#24908C",
         "Remnant" = "#FDE825")

colour = c("Invaded" = "#440C53",
           "Treated" = "#24908C",
           "Remnant" = "#FDE825")

shape = c("Invaded" = 21,
          "Treated" = 24,
          "Remnant" = 22)

## NMDS Axis 1, 2 

benthic.12 <- ggplot(data = benth.scores,
                    aes(x = NMDS1, y = NMDS2)) +
  geom_point(data = benth.scores, 
             aes(x = NMDS1, 
                 y = NMDS2, 
                 fill = Habitat, 
                 shape = Habitat,
                 stroke = 1.5),
             size = 7,
             alpha = 0.7) + # sites as points
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
                  size = 5,
                  force = 2) +
  scale_fill_manual(values = fill) +
  scale_colour_manual(values = colour) +
  scale_shape_manual(values = shape) 

benthic.12

## NMDS Axis 1, 3

benthic.13 <- ggplot(data = benth.scores,
                    aes(x = NMDS1, y = NMDS3)) +
  geom_point(data = benth.scores, 
             aes(x = NMDS1, y = NMDS3, 
                 fill = Habitat, shape = Habitat), 
             size = 7,
             stroke = 1.5,
             alpha = 0.7) +
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
                  size = 5, 
                  force = 2) +
  scale_fill_manual(values = fill) +
  scale_colour_manual(values = colour) +
  scale_shape_manual(values = shape) 

benthic.13

# putting both figures together

benthic.12 <- benthic.12 + ggtitle("A.   Aquatic Invertebrates")

benthic.12 + benthic.13



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


inv.h <- decostand(aq.inv.taxa, "hellinger")
unin.h <- decostand(aq.uninv.taxa, "hellinger")
trt.h <- decostand(aq.trt.taxa, "hellinger")

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


### Species contribution
# degree of variation of individual species

benthic.raw <- read.csv("Data/Aquatic/aquatic_inverts_rares2.csv") # occurrences <=2 removed
colnames(benthic.raw)

benth.r.taxa <- benthic.raw %>% select(Oligochaeta:Leptoceridae)
benth.r.env <- benthic.raw %>% select(Site.ID:Collection.date)

benth.r.h <- decostand(benth.r.taxa, "hellinger")

(benth.SCBD <- beta.div(benth.r.h, method = "hellinger",
                        sqrt.D = FALSE, samp = FALSE))

taxa.SCBD <- as.data.frame(benth.SCBD$SCBD)
write.csv(taxa.SCBD, "Data/Aquatic/SCBD_taxa.csv")

# Aquatic inverts by habitat ----------------------------------------------


## Each habitat alone

# Invaded 

(inv.beta <- beta.div(inv.bc, method = "percentdiff",
                   sqrt.D = FALSE, samp = FALSE,
                   nperm = 999))


aq.inv.env$LCBD <- inv.beta$LCBD

inv.beta$LCBD[inv.beta$LCBD > mean(inv.beta$LCBD)] # LCBD > average

#1         2         3         7 
#0.1635096 0.1573112 0.1288408 0.1386385 

(inv.SCBD <- beta.div(inv.h, method = "hellinger",
                      sqrt.D = FALSE, samp = FALSE))


inv.d.SCBD <- as.data.frame(inv.SCBD$SCBD)

write.csv(inv.d.SCBD, "Data/Aquatic/SCBD_taxa_invaded.csv")

## Uninvaded

(unin.beta <- beta.div(unin.bc, method = "percentdiff",
                      sqrt.D = FALSE, samp = FALSE,
                      nperm = 999))


aq.uninv.env$LCBD <- unin.beta$LCBD

unin.beta$LCBD[unin.beta$LCBD > mean(unin.beta$LCBD)] # LCBD > average

#1         3         5         6         7 
#0.1388476 0.1303577 0.1395451 0.1416225 0.1307846 

(unin.SCBD <- beta.div(unin.h, method = "hellinger",
                      sqrt.D = FALSE, samp = FALSE))


unin.d.SCBD <- as.data.frame(unin.SCBD$SCBD)

write.csv(unin.d.SCBD, "Data/Aquatic/SCBD_taxa_uninvaded.csv")



## Treated

(trt.beta <- beta.div(trt.bc, method = "percentdiff",
                       sqrt.D = FALSE, samp = FALSE,
                       nperm = 999))


aq.trt.env$LCBD <- trt.beta$LCBD

trt.beta$LCBD[trt.beta$LCBD > mean(trt.beta$LCBD)] # LCBD > average


# SCBD

(trt.SCBD <- beta.div(trt.h, method = "hellinger",
                       sqrt.D = FALSE, samp = FALSE))


trt.d.SCBD <- as.data.frame(trt.SCBD$SCBD)

write.csv(trt.d.SCBD, "Data/Aquatic/SCBD_taxa_treated.csv")


# Combine data

LCBD.data <- rbind(aq.trt.env, aq.inv.env, aq.uninv.env)


LCBD.anova <- lm(LCBD ~ Habitat, data = LCBD.data)
Anova(LCBD.anova)

#Anova Table (Type II tests)

#Response: LCBD
#Sum Sq Df F value Pr(>F)
#Habitat   0.001111  2  0.1422 0.8683
#Residuals 0.085978 22

#### Beta diversity calculations ####
library(betapart)

aq.inv.taxa
aq.uninv.taxa
aq.trt.taxa


aq.inv.taxa[aq.inv.taxa > 0] <- 1
aq.uninv.taxa[aq.uninv.taxa > 0] <- 1
aq.trt.taxa[aq.trt.taxa > 0] <- 1 

inv.bd.aq <- as.data.frame(beta.multi(aq.inv.taxa, index.family = "sorensen"))
inv.bd.aq$Habitat <- c("Invaded")

unin.bd.aq <- as.data.frame(beta.multi(aq.uninv.taxa, index.family = "sorensen"))
unin.bd.aq$Habitat <- c("Uninvaded")

trt.bd.aq <- as.data.frame(beta.multi(aq.trt.taxa, index.family = "sorensen"))
trt.bd.aq$Habitat <- c("Treated")

beta.aq <- rbind(trt.bd.aq, inv.bd.aq, unin.bd.aq)
write.csv(beta.aq, "Data/Emerging/beta_diversity_aquatic.csv")


### Null model
aq.inv.taxa
aq.uninv.taxa
aq.trt.taxa

## Invaded null model

library(picante)

invaded.null.aq <- data.frame(matrix(as.numeric(0), ncol=(3), nrow=(999)))
colnames(invaded.null.aq) <- c("Turnover","Nestedness","Sum")


for (i in 1:999){ #for 1000 iterations
  tempma <- as.data.frame(randomizeMatrix(aq.inv.taxa, null.model = "frequency")) #Randomize the data
  inv.bv.a <- beta.multi(tempma, index.family = "sorensen")
  invaded.null.aq[i,] <- data.frame(matrix(unlist(inv.bv.a), nrow = length(1), byrow = T)) 
}

inv.null.aq <- invaded.null.aq %>% 
  summarise(n = n(),
            N.avg = mean(Nestedness),
            N.sd = sd(Nestedness),
            N.CI = qnorm(0.95)*(N.sd/sqrt(8)), # 95% CI
            T.avg = mean(Turnover),
            T.sd = sd(Turnover),
            T.CI = qnorm(0.95)*(T.sd/sqrt(8)), # specify correct sample size
            S.avg = mean(Sum),
            S.sd = sd(Sum),
            S.CI = qnorm(0.95)*(S.sd/sqrt(8)))


inv.null.aq$Turn <- inv.bd.aq$beta.SIM
inv.null.aq$Nest <- inv.bd.aq$beta.SNE
inv.null.aq$Over <- inv.bd.aq$beta.SOR
inv.null.aq$Habitat <- c("Invaded")


# Uninvaded 


uninvaded.null.aq <- data.frame(matrix(as.numeric(0), ncol=(3), nrow=(999)))
colnames(uninvaded.null.aq) <- c("Turnover","Nestedness","Sum")


for (i in 1:999){ #for 1000 iterations
  tempua <- as.data.frame(randomizeMatrix(aq.uninv.taxa, null.model = "frequency")) #Randomize the data
  unin.bv.a <- beta.multi(tempua, index.family = "sorensen")
  uninvaded.null.aq[i,] <- data.frame(matrix(unlist(unin.bv.a), nrow = length(1), byrow = T)) 
}

unin.null.aq <- uninvaded.null.aq %>% 
  summarise(n = n(),
            N.avg = mean(Nestedness),
            N.sd = sd(Nestedness),
            N.CI = qnorm(0.95)*(N.sd/sqrt(8)), # 95% CI
            T.avg = mean(Turnover),
            T.sd = sd(Turnover),
            T.CI = qnorm(0.95)*(T.sd/sqrt(8)), # specify correct sample size
            S.avg = mean(Sum),
            S.sd = sd(Sum),
            S.CI = qnorm(0.95)*(S.sd/sqrt(8)))


unin.null.aq$Turn <- unin.bd.aq$beta.SIM
unin.null.aq$Nest <- unin.bd.aq$beta.SNE
unin.null.aq$Over <- unin.bd.aq$beta.SOR
unin.null.aq$Habitat <- c("Uninvaded")




# Treated

treated.null.aq <- data.frame(matrix(as.numeric(0), ncol=(3), nrow=(999)))
colnames(treated.null.aq) <- c("Turnover","Nestedness","Sum")


for (i in 1:999){ #for 1000 iterations
  tempta <- as.data.frame(randomizeMatrix(aq.trt.taxa, null.model = "frequency")) #Randomize the data
  trt.bv.aq <- beta.multi(tempta, index.family = "sorensen")
  treated.null.aq[i,] <- data.frame(matrix(unlist(trt.bv.aq), nrow = length(1), byrow = T)) 
}

treat.null.aq <- treated.null.aq %>% 
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


treat.null.aq$Turn <- trt.bd.aq$beta.SIM
treat.null.aq$Nest <- trt.bd.aq$beta.SNE
treat.null.aq$Over <- trt.bd.aq$beta.SOR
treat.null.aq$Habitat <- c("Treated")


bd.null.aq <- rbind(treat.null.aq, unin.null.aq,inv.null.aq)

write.csv(bd.null.aq, "Data/Aquatic/beta_diversity_null_true.csv")

colnames(bd.null.aq)

sumbdaq <- ggplot(bd.null.aq, aes(x = Habitat, y = Over, 
                             fill = Habitat, shape = Habitat)) +
  geom_errorbar(aes(ymin = S.avg - S.CI, ymax = S.avg + S.CI),
                colour = "black",
                size = 0.5,
                width = 0.3,
                position = position_dodge(0.6)) +
  geom_point(alpha = 0.7, size = 7,
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




nestaq <- ggplot(bd.null.aq, aes(x = Habitat, y = Nest, 
                            fill = Habitat, shape = Habitat)) +
  geom_errorbar(aes(ymin = N.avg - N.CI, ymax = N.avg + N.CI),
                colour = "black",
                size = 0.5,
                width = 0.3,
                position = position_dodge(0.6)) +
  geom_point(aes(y = N.avg),
             fill = "black",
             size = 3) +
  geom_point(alpha = 0.7, size = 7,
             stroke = 1.5) +
  labs(x = " ",
       y = expression(paste("Nestedness"))) +
  theme_classic(14)+
  theme(legend.position = "none") +
  scale_shape_manual(values = c(21, 24, 22)) +
  scale_fill_viridis(discrete = TRUE) +
  ylim(0.0, 0.2) +
  theme(axis.text = element_text(size = 14))



turnaq <- ggplot(bd.null.aq, aes(x = Habitat, y = Turn, 
                            fill = Habitat, shape = Habitat)) +
  geom_errorbar(aes(ymin = T.avg - T.CI, ymax = T.avg + T.CI),
                colour = "black",
                size = 0.5,
                width = 0.3,
                position = position_dodge(0.6)) +
  geom_point(aes(y = T.avg),
             fill = "black",
             size = 3) +
  geom_point(alpha = 0.7, size = 7,
             stroke = 1.5) +
  labs(x = " ",
       y = expression(paste("Turnover"))) +
  theme_classic(14) +
  theme(legend.position = "none") +
  scale_shape_manual(values = c(21, 24, 22)) +
  scale_fill_viridis(discrete = TRUE) +
  ylim(0.3, 0.75)+
  theme(axis.text = element_text(size = 14))



aqua.beta <- ggarrange(sumbdaq, nestaq, turnaq,
                        labels = "AUTO",
                        nrow = 1)

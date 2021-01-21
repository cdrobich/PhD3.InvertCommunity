
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

benthic <- read.csv("Data/Benthic_vegetation_QCC.csv") # occurrences = 1 removed
str(benthic)

benthic.data <- benthic %>% select(Oligochaeta:Hydroptilidae)
benthic.env <- benthic %>% select(Site.ID:Collection.date)


# Multivariate Analyses ---------------------------------------------------

# perMANOVA ---------------------------------------------------------------

benthic.rel <- decostand(benthic.data, "max", 2, na.rm = NULL) # divide by column max

#### perMANOVA #####

(per.treat <- adonis2(benthic.rel ~ Habitat, data = benthic,
                     permutations = 999, method = "bray"))


#Permutation test for adonis under reduced model
#Terms added sequentially (first to last)
#Permutation: free
#Number of permutations: 999

#         Df SumOfSqs      R2      F Pr(>F)    
#Habitat   2   2.2541 0.28114 4.3021  0.001 ***
#Residual 22   5.7636 0.71886                  
#Total    24   8.0177 1.00000                  


# homogeneity of groups dispersion

bugs.b <- vegdist(benthic.rel, method = "bray")

groups <- factor(benthic$Habitat)

(dispersion <- betadisper(bugs.b, groups)) # spatial median default

#Homogeneity of multivariate dispersions

#Call: betadisper(d = bugs.b, group = groups)

#No. of Positive Eigenvalues: 20
#No. of Negative Eigenvalues: 4

#Average distance to median:
#Invaded   Treated Uninvaded 
#0.4496    0.4673    0.4983 


plot(dispersion)
boxplot(dispersion) # actually look really good!


(adonis.pair(bugs.b, groups, nper = 1000, corr.method = "bonferroni"))

#      combination   SumsOfSqs   MeanSqs  F.Model        R2    P.value  P.value.corrected
#1 Invaded <-> Treated 0.9702695 0.9702695 3.869514 0.2050670 0.000999001       0.002997003
#2 Invaded <-> Uninvaded 1.0981493 1.0981493 4.193700 0.2305029 0.000999001       0.002997003
#3 Treated <-> Uninvaded 1.3110559 1.3110559 4.796596 0.2422940 0.000999001       0.002997003



# NMDS  -------------------------------------------------------------------

## stress plot/ see how many axes make sense

k_vec <- 1:10 #dimensions 1 - 10
stress <- numeric(length(k_vec)) # stress of each model put here
dune_dij <- metaMDSdist(benthic.rel)
set.seed(25)
for(i in seq_along(k_vec)) {
  sol <- metaMDSiter(dune_dij, k = i, # you can tell I lifted this from a tutorial on the dune package
                     trace = FALSE)
  stress[i] <- sol$stress
}

plot(stress) # 3D makes sense

#### NMDS analysis ####
set.seed(120) 

nms.invert <- metaMDS(benthic.rel, distance = "bray", # species data, bray-curtis dissimilarity
                      autotransform = FALSE,  # NMDS will do autotransformations for you
                      k = 3, trymax = 1000)   # k = number of axes
nms.invert

#Call:
#  metaMDS(comm = benthic.rel, distance = "bray", k = 3, trymax = 1000,      autotransform = FALSE) 

#global Multidimensional Scaling using monoMDS

#Data:     benthic.rel 
#Distance: bray 

#Dimensions: 3 
#Stress:     0.1437523
#Stress type 1, weak ties
#Two convergent solutions found after 20 tries
#3Scaling: centring, PC rotation, halfchange scaling 
#Species: expanded scores based on ‘benthic.rel’ 

# look at points and the stress real quick
layout(matrix(1:2, ncol = 2))
plot(nms.invert, main = "Invertebrate NMDS plot"); stressplot(nms.invert, main = "Shepard plot")
layout(1)

ordiplot(nms.invert, type = "n")
orditorp(nms.invert, display = "species")
orditorp(nms.invert, display = "sites")

# how many iterations of the NMDS
nms.invert$iters # 56

# Goodness of fit
(g <- goodness(nms.invert)) # smaller the number the better the fit
sum(g^2)
nms.invert$stress^2  # 0.02066487

1-nms.invert$stress^2 #0.9793351 #analogous to square correlation coefficient


## extract the scores for plotting 
scr <- as.data.frame(scores(nms.invert, display = "sites")) # extract NMDS scores

# adding categorical info to scores
scr$Site <- rownames(scr)
scr$Habitat <- benthic$Habitat

colnames(scr)

col_order <- c("Site", "Habitat", "NMDS1", "NMDS2", "NMDS3")

scores <- scr[, col_order] # put the categorical values in order


write.csv(scores,"Data/NMDS_benthic_inverts_scores.csv") # save this as a csv


### Vectors correlated with Axis 1 & 2 ####

alltaxa <- envfit(nms.invert, benthic.rel,
                 choices = c(1,2)) #produces a list with r2, p value, and NMDS coordinates

all.taxa.df <- data.frame((alltaxa$vectors)$arrows,
                          (alltaxa$vectors)$r,
                          (alltaxa$vectors)$pvals) #take list and make into dataframe

write.csv(all.taxa.df, "Data/NMDS_benthic_vectors_axis12.csv") # save vector scores as csv


corr.spp12 <- all.taxa.df %>% filter(X.alltaxa.vectors..r > 0.2)
corr.spp12$species <- rownames(corr.spp)

corr.species12 <- corr.spp12$species # string of the Family names

axis12.vectors <- benthic.rel %>% select(all_of(corr.species12)) # make a matrix of just those


corrtaxa12 <- envfit(nms.invert$points, axis12.vectors, 
                   permutations = 999, choices = c(1,2))
corrtaxa12

# make a new data frame for the figure
species.12 <- as.data.frame(corrtaxa12$vectors$arrows*sqrt(corrtaxa12$vectors$r)) #scaling vectors so they correspond with r2
species.12$species <- rownames(species.12)


#### Vectors correlated with axis 1 & 3 ###
# same as above but now for the other axes/figure

alltaxa.13 <- envfit(nms.invert, benthic.rel, 
                     permutations = 999, choices = c(1,3)) 


all.taxa13.df <- data.frame((alltaxa.13$vectors)$arrows,
                            (alltaxa.13$vectors)$r,
                            (alltaxa.13$vectors)$pvals)

write.csv(all.taxa13.df, "Data/NMDS_vectors_axis13.csv")



corr.spp13 <- all.taxa13.df %>% filter(X.alltaxa.13.vectors..r > 0.2)
corr.spp13$species <- rownames(corr.spp13)

corr.species13 <- corr.spp13$species # string of the Family names

axis13.vectors <- benthic.rel %>% select(all_of(corr.species13)) # make a matrix of just those


corrtaxa13 <- envfit(nms.invert$points, axis13.vectors, 
                     permutations = 999, choices = c(1,3))
corrtaxa13

# make a new data frame for the figure
species.13 <- as.data.frame(corrtaxa13$vectors$arrows*sqrt(corrtaxa13$vectors$r)) #scaling vectors so they correspond with r2
species.13$species <- rownames(species.13)


## plotting with base R ##
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


#### ggPlot Figures ####

## NMDS Axis 1, 2 

scores # coordinates we extracted
species.12 # reasonably correlated vectors with axis 1,2

benthic.12 <- ggplot(data = scores,
                    aes(x = NMDS1, y = NMDS2)) +
  geom_point(data = scores, 
             aes(x = NMDS1, 
                 y = NMDS2, 
                 colour = Habitat, 
                 shape = Habitat),
             size = 4) + # sites as points
  stat_ellipse(data = scores, aes(x = NMDS1,y = NMDS2,
                                  linetype = Habitat, colour = Habitat), 
               size = 1, level = 0.90) + # a 95% CI ellipses
  geom_segment(data = species.12, aes(x = 0, xend = MDS1, y = 0, yend = MDS2), # adding in the vectors, c
               arrow = arrow(length = unit(0.5, "cm")), colour = "black") + # can add in geom_label or geom_text for labels
  theme_minimal() + # no background
  theme(panel.border = element_rect(fill = NA)) + # full square around figure
  xlab("NMDS 1") +
  ylab("NMDS 2") +
  theme(legend.position = "none") +
  geom_text_repel(data = species.12, 
                  aes(x = MDS1, y = MDS2, label = species),
                  color="black",
                  size = 5) +
  scale_color_manual(values = c("#969696","#35978f", "#2166ac")) +
  scale_shape_manual(values = c(17, 18, 15)) +
  coord_fixed()


benthic.12

## NMDS Axis 1, 3
# same as above

scores
species.13 # vectors correlated with axis 1, 3

benthic.13 <- ggplot(data = scores,
                    aes(x = NMDS1, y = NMDS3)) +
  geom_point(data = scores, 
             aes(x = NMDS1, y = NMDS3, 
                 colour = Habitat, shape = Habitat), size = 4) +
  stat_ellipse(data = scores, 
               aes(x = NMDS1,y = NMDS3,linetype = Habitat, 
                   colour = Habitat), size = 1, level = 0.9) +
  geom_segment(data = species.13, 
               aes(x = 0, xend = MDS1, y = 0, yend = MDS3),
               arrow = arrow(length = unit(0.5, "cm")), colour = "black") +
  theme_minimal() +
  theme(panel.border = element_rect(fill = NA)) +
  xlab("NMDS 1") +
  ylab("NMDS 3") +
  theme(legend.position = c(0.85, 0.9)) +
  geom_text_repel(data = species.13, 
                  aes(x = MDS1, y = MDS3, label = species),
                  color="black",
                  size = 5) +
  scale_color_manual(values = c("#969696","#35978f", "#2166ac")) +
  scale_shape_manual(values = c(17, 18, 15)) +
  coord_fixed()

benthic.13


(NMS.benthic.panel <- ggarrange(benthic.12, benthic.13,
                               align = "hv"))

ggsave("Figures/Benthic_NMDS_panel.TIFF", NMS.benthic.panel,
       dpi = 300,
       height = 7.33,
       width = 11.9,
       units = "in")





library(tidyverse) # for figures, pipes etc.
library(vegan) # for NMDS
library(ggpubr) # for arranging figures
library(agricolae)
library(Hmisc)
library(car)
library(viridis) # colours


benthic <- read.csv("Data/Benthic_vegetation_QCC.csv") # occurrences = 1 removed
str(benthic)

benthic.data <- benthic %>% select(Oligochaeta:Hydroptilidae)
benthic.env <- benthic %>% select(Site.ID:Collection.date)


###### Univariate Analyses ###############

richness <- rowSums(benthic.data > 0) # species richness
abundance <- rowSums(benthic.data) # abundance

# make new little data frame
benthic.uni<- benthic.env
benthic.uni$rich <- richness
benthic.uni$abundance <- abundance

colnames(benthic.uni)


abundance.lm <- lm(abundance ~ Habitat, data = benthic.uni)
Anova(abundance.lm, type = 3)

abundance.hsd <- HSD.test(abundance.lm, "Habitat")

#Anova Table (Type III tests)

#Response: abundance
#               Sum Sq Df F value   Pr(>F)   
#(Intercept)  2726112  1  0.8434 0.368396   
#Habitat     47207700  2  7.3022 0.003696 **
#Residuals   71113597 22  

#          abundance       std r Min  Max  Q25    Q50     Q75
#Invaded      583.75  622.2551 8  73 1767  233  320.5  661.25
#Treated   3375.00 2920.3388 9  42 9851 2194 2472.0 3180.00
#Uninvaded    445.25  158.6368 8 234  781  381  406.5  476.00


#           abundance groups
#Treated   3375.00      a
#Invaded      583.75      b
#Uninvaded    445.25      b


rich.lm <- lm(rich ~ Habitat, data = benthic.uni)
Anova(rich.lm, type = 3)

#Anova Table (Type III tests)

#Response: rich
#             Sum Sq Df  F value   Pr(>F)    
#(Intercept) 1800.00  1 155.3703 1.91e-11 ***
#Habitat       74.57  2   3.2181  0.05944 .  
#Residuals    254.87 22 


benthic.uni %>% 
  group_by(Habitat) %>% 
  mutate(mean.s = mean(rich),
         median.s = median(rich),
         N = length(abundance),
         sd.s = sd(rich),
         sterr.s = (sd.s/(sqrt(N))),
         mean.ab = mean(abundance),
         median.ab = median(abundance),
         sd.ab = sd(abundance),
         sterr.ab = (sd.ab/sqrt(N))) -> benthic.uni


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
  scale_colour_viridis(discrete = TRUE)


richness <- ggplot(data = benthic.uni, 
                   aes(x = rich, group = Habitat, fill = Habitat)) +
  geom_density(adjust = 1.5, alpha = 0.6) +
  xlim(0, 42) +
  theme_classic() +
  xlab("Species Richness") +
  ylab("Density") +
  theme(legend.position = c(0.8, 0.8)) +
  geom_vline(aes(xintercept = mean.s, colour = Habitat),
             size = 1,
             show.legend = FALSE) +
  scale_fill_viridis(discrete = TRUE) +
  scale_colour_viridis(discrete = TRUE)



ggarrange(abundance, richness,
          labels = "AUTO",
          hjust = c(-5.5),
          vjust = 2)

# Violin Plots

violin.s <- ggplot(data = benthic.uni, 
       aes(x = Habitat, y = rich, group = Habitat, fill = Habitat)) +
  geom_violin(trim = FALSE) +
  theme_classic(base_size = 14) +
  xlab(" ") +
  ylab("Species Richness") +
  theme(legend.position = "none",
        axis.text = element_text(size = 14)) +
  ylim(0, 30) +
  scale_fill_viridis(discrete = TRUE)



violin.ab <- ggplot(data = benthic.uni, 
       aes(x = Habitat, y = abundance, group = Habitat, fill = Habitat)) +
  geom_violin(trim = FALSE, scale = "count") +
  theme_classic(base_size = 14) +
  xlab(" ") +
  ylab("Abundance") +
  theme(legend.position = "none",
        axis.text = element_text(size = 14)) +
  scale_fill_viridis(discrete = TRUE) +
  annotate("text", x = 1:3, y = c(4600, 12000, 4600),
                     label = c("a", "b", "a"),
           size = 4.5)

(violin.benthic <- ggarrange(violin.ab, violin.s,
          labels = "AUTO",
          hjust = c(-8, -5.5),
          vjust = 2))

ggsave("Figures/benthic_violin.TIFF", violin.benthic,
       dpi = 300)

# Ridge Plots

library(ggridges)

ggplot(benthic.uni, aes(x = abundance, y = Habitat, fill = Habitat)) +
  geom_density_ridges() +
  theme_ridges() + 
  theme(legend.position = "none") +
  scale_fill_viridis(discrete = TRUE)
  

ggplot(benthic.uni, aes(x = rich, y = Habitat, fill = Habitat)) +
  geom_density_ridges() +
  theme_ridges() + 
  theme(legend.position = "none") +
  scale_fill_viridis(discrete = TRUE)


############ Multivariate Analyses ################

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



library(devtools)
install_github("GuillemSalazar/EcolUtils")
library(EcolUtils)


(adonis.pair(bugs.b, groups, nper = 1000, corr.method = "bonferroni"))

#      combination   SumsOfSqs   MeanSqs  F.Model        R2    P.value  P.value.corrected
#1 Invaded <-> Treated 0.9702695 0.9702695 3.869514 0.2050670 0.000999001       0.002997003
#2 Invaded <-> Uninvaded 1.0981493 1.0981493 4.193700 0.2305029 0.000999001       0.002997003
#3 Treated <-> Uninvaded 1.3110559 1.3110559 4.796596 0.2422940 0.000999001       0.002997003

############ NMDS ############

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
#Stress:     0.143753 
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
nms.invert$stress^2  # 0.02066492

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

all.taxa.df <- data.frame((alltaxa$vectors)$arrows, (alltaxa$vectors)$r, (alltaxa$vectors)$pvals) #take list and make into dataframe
write.csv(all.taxa.df, "Data/NMDS_benthic_vectors_axis12.csv") # save vector scores as csv


alltaxa$vectors$r[alltaxa$vectors$r > 0.3] # selecting vectors (Family) that are reasonably correlated (r2 > 0.2)

#Hydrazoa    Gastropoda       Araneae    Collembola  Staphylinidae Hydrophilidae  Chironomidae Cecidomyiidae 
#0.3915531     0.3767479     0.4266004     0.6149432     0.5007713     0.4577720     0.6401936     0.3256395 

#Caenidae  Leptoceridae 
#0.3846770     0.3539771 


# taking out those correlated ones to add to the figure
corr.taxa <- benthic.rel %>% select(Hydrazoa ,
                                    Gastropoda,
                                    Araneae,
                                    Collembola,
                                    Staphylinidae,
                                    Hydrophilidae,
                                    Chironomidae,
                                    Cecidomyiidae,
                                    Caenidae,
                                    Leptoceridae)
                    
# recalculated it because I am lazy but you could probably pull them from all.taxa.df
corrtaxa <- envfit(nms.invert$points, corr.taxa, 
                   permutations = 999, choices = c(1,2))


corrtaxa

# make a new data frame for the figure
species.12 <- as.data.frame(corrtaxa$vectors$arrows*sqrt(corrtaxa$vectors$r)) #scaling vectors so they correspond with r2
species.12$species <- rownames(species.12)


#### Vectors correlated with axis 1 & 3 ###
# same as above but now for the other axes/figure

alltaxa.13 <- envfit(nms.invert, benthic.rel, 
                     permutations = 999, choices = c(1,3)) 


all.taxa13.df <- data.frame((alltaxa.13$vectors)$arrows, (alltaxa.13$vectors)$r, (alltaxa.13$vectors)$pvals)
write.csv(all.taxa13.df, "Data/NMDS_vectors_axis13.csv")


alltaxa.13$vectors$r[alltaxa.13$vectors$r > 0.3] 

#Araneae     Collembola  Staphylinidae  Hydrophilidae  Curculionidae Dolichopodidae 
#0.4264369      0.6460189      0.5405840      0.4493061      0.3753543      0.4723451 

#Cecidomyiidae   Thysanoptera Coenagrionidae 
#0.3392635      0.5573573      0.4977675

corr.taxa.13 <- benthic.rel %>% select(Araneae,
                                Collembola,
                                Staphylinidae,
                                Hydrophilidae,
                                Curculionidae,
                                Dolichopodidae,
                                Cecidomyiidae,
                                Thysanoptera,
                                Coenagrionidae)


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

invert.12 <- ggplot(data = scores,
                    aes(x = NMDS1, y = NMDS2)) +
  geom_point(data = scores, aes(x = NMDS1, y = NMDS2, colour = Habitat, shape = Habitat),
             size = 4) + # sites as points
  stat_ellipse(data = scores, aes(x = NMDS1,y = NMDS2,
                                  linetype = Habitat, colour = Habitat), size = 1) + # a 95% CI ellipses
  geom_segment(data = species.12, aes(x = 0, xend = MDS1, y = 0, yend = MDS2), # adding in the vectors, c
               arrow = arrow(length = unit(0.5, "cm")), colour = "black") + # can add in geom_label or geom_text for labels
  #xlim (-1, 1) + # setting the limits so they're symmetrical
  #ylim (-1, 1) +
  theme_minimal() + # no background
  theme(panel.border = element_rect(fill = NA)) + # full square around figure
  xlab("NMDS 1") +
  ylab("NMDS 2") +
  scale_fill_viridis(discrete = TRUE) +
  scale_colour_viridis(discrete = TRUE) +
  geom_label(data = species.12,aes(x=MDS1,y=MDS2,label=species),size=5) +
  ylim(-1, 1.5) +
  xlim(-1.45, 1) +
  theme(legend.position = "none")


## NMDS Axis 1, 3
# same as above

scores
species.13 # vectors correlated with axis 1, 3

invert.13 <- ggplot(data = scores,
                    aes(x = NMDS1, y = NMDS3)) +
  geom_point(data = scores, aes(x = NMDS1, y = NMDS3, colour = Habitat, shape = Habitat), size = 4) +
  stat_ellipse(data = scores, aes(x = NMDS1,y = NMDS3,linetype = Habitat, colour = Habitat), size = 1) +
  geom_segment(data = species.13, aes(x = 0, xend = MDS1, y = 0, yend = MDS3),
               arrow = arrow(length = unit(0.5, "cm")), colour = "black") +
  theme_minimal() +
  theme(panel.border = element_rect(fill = NA)) +
  xlab("NMDS 1") +
  ylab("NMDS 3") +
geom_label(data = species.13,aes(x=MDS1,y=MDS3,label=species),size=5) +
  scale_fill_viridis(discrete = TRUE) +
  scale_colour_viridis(discrete = TRUE) +
  theme(legend.position = c(0.9, 0.9)) +
  xlim(-1.5, 1.0)

invert.13


NMS.benthic.panel <- ggarrange(invert.12, invert.13)

ggsave("Figures/Benthic_NMDS_panel.TIFF", NMS.benthic.panel,
       dpi = 300)





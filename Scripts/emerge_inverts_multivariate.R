
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

##### perMANVOA ####

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

bugs.b <- vegdist(invert.rel, method = "bray")

groups <- factor(invert$TrtYr)

(dispersion <- betadisper(bugs.b, groups)) 

#Average distance to median:
#Invaded_2017   Invaded_2018   Treated_2017   Treated_2018 Uninvaded_2017 Uninvaded_2018 
#0.5381         0.5030         0.5366         0.4939         0.5658         0.5443 

plot(dispersion)
boxplot(dispersion)



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
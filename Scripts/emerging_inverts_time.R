
# Load packages -----------------------------------------------------------

library(tidyverse) # for figures, pipes etc.
library(vegan) # for NMDS
library(ggpubr) # for arranging figures

library(ggrepel) # labels on nmds
library(lubridate)

library(viridis) # colours

# pairwise perMANOVA comparison
library(devtools)
install_github("GuillemSalazar/EcolUtils")
library(EcolUtils)


# Load data ---------------------------------------------------------------

data <- read.csv("Data/Emerging/emerging_invert_dates.csv")
head(data)
str(data)

unique(data$Treatment)


taxa <- data %>% select(Araneae:Hesperiidae)
env <- data %>% select(ID:YrCol)


## Remove rares

taxa.occur <- as.data.frame(sapply(taxa, function(col) length(unique(col)))) # how many unique values in column


rares <- taxa[, sapply(taxa, function(col) length(unique(col))) > 2] # remove those with fewer than 2 occurrences 
rare.removes <- as.data.frame(sapply(rares, function(col) length(unique(col)))) # check to make sure that worked
min(rare.removes)


write.csv(rares, "Data/Emerging/inverts_dates_rares.csv")

# Collection periods -----------------------------------------------

env <- env %>% separate(Date, c("Day","Month", "Year"), remove = FALSE)

unique(env$Date)

#                        "19-Jun-17" "28-Jun-17" "08-Jul-17" "21-Jul-17" 
#"20-May-18" "05-Jun-18" "16-Jun-18" "25-Jun-18" "04-Jul-18" "23-Jul-18"


inverts.date <- read.csv("Data/Emerging/inverts_dates_rares_env.csv")


# NMDS Collection June 19 2017 and June 16 2018 Collection -----------------------------------------

col.1 <- c("19-Jun-17","16-Jun-18")

invert.col.1 <- inverts.date %>% filter(Date %in% col.1)

write.csv(invert.col.1, "Data/Emerging/NMDS/Col.1/inverts_collection1.csv")


invert.col1.rares <- read.csv("Data/Emerging/NMDS/Col.1/inverts_collection1_zerosremoved.csv")


# just taxa and env 
taxa.col1 <- invert.col1.rares %>% select(Araneae:Crambidae)
env.col1 <- invert.col1.rares %>% select(ID:YrCol)

# Relativize by column max

taxa.col1rel <- decostand(taxa.col1, "max", 2, na.rm = NULL)

write.csv(taxa.col1rel, "Data/Emerging/NMDS/Col.1/inverts_col1_raresrel.csv")

# NMDS Col.1  -------------------------------------------------------------

k_vec <- 1:10 #dimensions 1 - 10
stress <- numeric(length(k_vec)) # stress of each model put here
dune_dij <- metaMDSdist(taxa.col1rel)

set.seed(25)

for(i in seq_along(k_vec)) {
  sol <- metaMDSiter(dune_dij, k = i, 
                     trace = FALSE)
  stress[i] <- sol$stress
}

plot(stress) # 3D makes sense

#### NMDS analysis 

set.seed(120) 

nms.col1 <- metaMDS(taxa.col1rel, distance = "bray", # species data, bray-curtis dissimilarity
                      autotransform = FALSE,  # NMDS will do autotransformations for you
                      k = 3, trymax = 1000)   # k = number of axes
nms.col1

#Call:
#metaMDS(comm = taxa.col1rel, distance = "bray", k = 3, trymax = 1000, autotransform = FALSE) 

#global Multidimensional Scaling using monoMDS
#Data:     taxa.col1rel 
#Distance: bray 

#Dimensions: 3 
#Stress:     0.1892989 
#Stress type 1, weak ties
#Two convergent solutions found after 33 tries
#Scaling: centring, PC rotation, halfchange scaling 
#Species: expanded scores based on ‘taxa.col1rel’ 



# look at points and the stress real quick
layout(matrix(1:2, ncol = 2))
plot(nms.col1, main = "Invertebrate NMDS plot"); stressplot(nms.col1, main = "Shepard plot")
layout(1)

ordiplot(nms.col1, type = "n")
orditorp(nms.col1, display = "species")
orditorp(nms.col1, display = "sites")


ordiplot(nms.col1, type = "n", choices = c(1,3))
orditorp(nms.col1, display = "species", choices = c(1,3))
orditorp(nms.col1, display = "sites", choices = c(1,3))


# how many iterations of the NMDS
nms.col1$iters # 127

# Goodness of fit
(g <- goodness(nms.col1)) # smaller the number the better the fit
sum(g^2) 0.03583409
nms.col1$stress^2  

1-nms.col1$stress^2 # 0.9641659 #analogous to square correlation coefficient


## extract the scores for plotting 
scr <- as.data.frame(scores(nms.col1, display = "sites")) # extract NMDS scores
colnames(scr)

# adding categorical info to scores
env.col1$NMDS1 <- scr$NMDS1
env.col1$NMDS2 <- scr$NMDS2
env.col1$NMDS3 <- scr$NMDS3

write.csv(env.col1,"Data/Emerging/NMDS/Col.1/emerging_inverts_col1_NMDSscores.csv") # save this as a csv

### Taxa correlated with Axis 1 & 2 

col1.12 <- envfit(nms.col1, taxa.col1rel,
                    choices = c(1,2)) #produces a list with r2, p value, and NMDS coordinates

col1.df12 <- data.frame((col1.12$vectors)$arrows,
                          (col1.12$vectors)$r,
                          (col1.12$vectors)$pvals) #take list and make into data frame

col1.df12 <- tibble::rownames_to_column(col1.df12, "Taxa")

write.csv(col1.df12, "Data/Emerging/NMDS/Col.1/NMDS_emerg_col1_vector12.csv") # save vector scores as csv


## Taxa correlated with Axis 1 & 3

col1.13 <- envfit(nms.col1, taxa.col1rel, 
                     permutations = 999, choices = c(1,3)) 


col1.df13 <- data.frame((col1.13$vectors)$arrows,
                            (col1.13$vectors)$r,
                            (col1.13$vectors)$pvals)

col1.df13 <- tibble::rownames_to_column(col1.df13, "Taxa")

write.csv(col1.df13, "Data/Emerging/NMDS/Col.1/NMDS_emerg_col1_vector13.csv")

# Selecting taxa for vectors on 1 & 2

colnames(col1.df12)
col1.corrsp.12 <- col1.df12 %>% filter(X.col1.12.vectors..r > 0.2) #those r>0.2
target12.c1 <- col1.corrsp.12$Taxa # string of the Family names

axis12.vectors.c1 <- taxa.col1rel %>% select(all_of(target12.c1)) # make a matrix of just those


# Selecting taxa for vectors on axis 1, 3
colnames(col1.df13)
col1.corrsp.13 <- col1.df13 %>% filter(X.col1.13.vectors..r > 0.2) 
target13.c1 <- col1.corrsp.13$Taxa # string of the Family names

axis13.vectors.c1 <- taxa.col1rel %>% select(all_of(target13.c1)) # make a matrix of just those


# Fit vectors to NMDS

# axis 1, 2
(nmds.c1.vectors.12 <- envfit(nms.col1$points, axis12.vectors.c1,
                           permutations = 999, choices = c(1,2)))                        


corr.c1.vectors.12 <- as.data.frame(nmds.c1.vectors.12$vectors$arrows*sqrt(nmds.c1.vectors.12$vectors$r)) #scaling vectors
corr.c1.vectors.12$Taxa <- rownames(corr.c1.vectors.12)


write.csv(corr.c1.vectors.12, "Data/Emerging/NMDS/Col.1/NMDS_emerg_col1_vector12_final.csv")

# axis 1, 3
(nmds.c1.vectors.13 <- envfit(nms.col1$points, axis13.vectors.c1,
                           permutations = 999, choices = c(1,3)))                        


corr.c1.vectors.13 <- as.data.frame(nmds.c1.vectors.13$vectors$arrows*sqrt(nmds.c1.vectors.13$vectors$r)) #scaling vectors
corr.c1.vectors.13$Taxa <- rownames(corr.c1.vectors.13) 

write.csv(corr.c1.vectors.13, "Data/Emerging/NMDS/Col.1/NMDS_emerg_col1_vector13_final.csv")

# NMDS Emerging Invert Figure Collection 1 ------------------------------------------------------------

col1.scores <- read.csv("Data/Emerging/NMDS/Col.1/emerging_inverts_col1_NMDSscores.csv")
col1.scores$Year <- as.factor(col1.scores$Year)

col1.axis12 <- read.csv("Data/Emerging/NMDS/Col.1/NMDS_emerg_col1_vector12_final.csv")
col1.axis13 <- read.csv("Data/Emerging/NMDS/Col.1/NMDS_emerg_col1_vector13_final.csv")


## NMDS Axis 1, 2 


invert.12.c1 <- ggplot(data = col1.scores,
                    aes(x = NMDS1, y = NMDS2)) +
  geom_point(data = col1.scores, 
             aes(x = NMDS1, y = NMDS2, 
                 colour = Treatment, shape = Year),
             size = 4, stroke = 1.5) + # sites as points
  stat_ellipse(data = col1.scores, 
               aes(x = NMDS1,y = NMDS2,
                   linetype = Treatment, colour = Treatment), 
               size = 1, level = 0.9) + 
  geom_segment(data = col1.axis12, 
               aes(x = 0, xend = MDS1, y = 0, yend = MDS2), # adding in the vectors, c
               arrow = arrow(length = unit(0.5, "cm")), colour = "black") + # can add in geom_label or geom_text for labels
  theme_minimal() + # no background
  theme(panel.border = element_rect(fill = NA)) + # full square around figure
  xlab("NMDS 1") +
  ylab("NMDS 2") +
  #ylim(-1, 1.5) +
  #xlim(-1.45, 1) +
  #theme(legend.position = "none") +
  geom_text_repel(data = col1.axis12, 
                  aes(x = MDS1, y = MDS2, label = Taxa),
                  color="black",
                  size = 5) +
  scale_color_manual(values = c("#969696","#35978f", "#2166ac")) +
  scale_shape_manual(values = c(17, 16, 1, 17, 2, 18, 5)) 

invert.12.c1

## NMDS Axis 1, 3

invert.13.c1 <- ggplot(data = col1.scores,
                       aes(x = NMDS1, y = NMDS3)) +
  geom_point(data = col1.scores, 
             aes(x = NMDS1, y = NMDS3, 
                 colour = Treatment, shape = Year),
             size = 4, stroke = 1.5) + # sites as points
  stat_ellipse(data = col1.scores, 
               aes(x = NMDS1,y = NMDS3,
                   linetype = Treatment, colour = Treatment), 
               size = 1, level = 0.9) + 
  geom_segment(data = col1.axis13, 
               aes(x = 0, xend = MDS1, y = 0, yend = MDS3), # adding in the vectors, c
               arrow = arrow(length = unit(0.5, "cm")), colour = "black") + # can add in geom_label or geom_text for labels
  theme_minimal() + # no background
  theme(panel.border = element_rect(fill = NA)) + # full square around figure
  xlab("NMDS 1") +
  ylab("NMDS 3") +
  #ylim(-1, 1.5) +
  #xlim(-1.45, 1) +
  #theme(legend.position = "none") +
  geom_text_repel(data = col1.axis13, 
                  aes(x = MDS1, y = MDS3, label = Taxa),
                  color="black",
                  size = 5) +
  scale_color_manual(values = c("#969696","#35978f", "#2166ac")) +
  scale_shape_manual(values = c(17, 16, 1, 17, 2, 18, 5)) 

invert.13.c1

(NMS.emerging.panel.c1 <- ggarrange(invert.12.c1, invert.13.c1,
                                 common.legend = TRUE,
                                 legend = "none",
                                 labels = c("A", ""),
                                 align = "hv"))

NMDS.col1 <- annotate_figure(NMS.emerging.panel.c1,
                top = text_grob("Collection 19-Jun-2017 and 16-Jun-2018"))

ggsave("Figures/NMDS_emerging_19Jun17_16Jun18.jpeg", NMDS.col1)



# NMDS Collection 28-Jun-17 and 25-JUN-18 ---------------------------------

col.2 <- c("28-Jun-17","25-Jun-18")

invert.col.2 <- inverts.date %>% filter(Date %in% col.2)

write.csv(invert.col.2, "Data/Emerging/NMDS/Col.2/inverts_collection2.csv")

# remove zeros in excel because i dont have the time
invert.col2.rares <- read.csv("Data/Emerging/NMDS/Col.2/inverts_collection2_zeros.csv")
colnames(invert.col2.rares)

# just taxa and env 
taxa.col2 <- invert.col2.rares %>% select(Araneae:Noctuidae)
env.col2 <- invert.col2.rares %>% select(ID:YrCol)

# Relativize by column max

taxa.col1re2 <- decostand(taxa.col2, "max", 2, na.rm = NULL)

write.csv(taxa.col1re2, "Data/Emerging/NMDS/Col.2/inverts_col2_raresrel.csv") # use for NMDS

### NMDS "Col 2" 

data.col2 <- read.csv("Data/Emerging/NMDS/Col.2/inverts_col2_raresrel.csv")

k_vec <- 1:10 #dimensions 1 - 10
stress <- numeric(length(k_vec)) # stress of each model put here
dune_dij <- metaMDSdist(data.col2)

set.seed(25)

for(i in seq_along(k_vec)) {
  sol <- metaMDSiter(dune_dij, k = i, 
                     trace = FALSE)
  stress[i] <- sol$stress
}
plot(stress) # 3D




set.seed(120) 

nms.col2 <- metaMDS(data.col2, distance = "bray", # species data, bray-curtis dissimilarity
                    autotransform = FALSE,  # NMDS will do autotransformations for you
                    k = 3, trymax = 1000)   # k = number of axes
nms.col2

#Dimensions: 3 
#Stress:     0.02987563 
#Stress type 1, weak ties
#Two convergent solutions found after 718 tries
#Scaling: centring, PC rotation, halfchange scaling 
#Species: expanded scores based on ‘data.col2’

nms.col2$iters #200

nms.col2$stress^2   #0.0008925531
1-nms.col2$stress^2 #0.9991074

scr2 <- as.data.frame(scores(nms.col2, display = "sites")) # extract NMDS scores
colnames(scr2)

env.col2$NMDS1 <- scr2$NMDS1
env.col2$NMDS2 <- scr2$NMDS2
env.col2$NMDS3 <- scr2$NMDS3

write.csv(env.col2,"Data/Emerging/NMDS/Col.2/emerging_col2_NMDSscores.csv") # save this as a csv

## Taxa for vectors

# Axis 1 and 2
taxa.col2.axis12 <- envfit(nms.col2, taxa.col2,
                           choices = c(1,2))

taxa.col2.axis12df <- data.frame((taxa.col2.axis12$vectors)$arrows,
                                 (taxa.col2.axis12$vectors)$r,
                                 (taxa.col2.axis12$vectors)$pvals)

taxa.col2.axis12df <- tibble::rownames_to_column(taxa.col2.axis12df, "Taxa")

write.csv(taxa.col2.axis12df, "Data/Emerging/NMDS/Col.2/col2_allvectors_axis12.csv") # save vector scores as csv

# Axis 1 and 3

taxa.col2.axis13 <- envfit(nms.col2, taxa.col2,
                           choices = c(1,3))

taxa.col2.axis13df <- data.frame((taxa.col2.axis13$vectors)$arrows,
                                 (taxa.col2.axis13$vectors)$r,
                                 (taxa.col2.axis13$vectors)$pvals)

taxa.col2.axis13df <- tibble::rownames_to_column(taxa.col2.axis13df, "Taxa")

write.csv(taxa.col2.axis13df, "Data/Emerging/NMDS/Col.2/col2_allvectors_axis13.csv") # save vector scores as csv

# Correlated taxa

# Axis 1 and 2

colnames(taxa.col2.axis12df)

corrspp.col2.axis12 <- taxa.col2.axis12df %>% filter(X.taxa.col2.axis12.vectors..r > 0.2)
target12.c2 <- corrspp.col2.axis12$Taxa # string of the Family names


axis12.vectors.c2 <- taxa.col2 %>% select(all_of(target12.c2)) # make a matrix of just those

(nmds.c2.vectors.12 <- envfit(nms.col2$points, axis12.vectors.c2,
                              permutations = 999, choices = c(1,2)))                        


corr.c2.vectors.12 <- as.data.frame(nmds.c2.vectors.12$vectors$arrows*sqrt(nmds.c2.vectors.12$vectors$r)) #scaling vectors
corr.c2.vectors.12$Taxa <- rownames(corr.c2.vectors.12)


write.csv(corr.c2.vectors.12, "Data/Emerging/NMDS/Col.2/emerging_correlated_vector12.csv")

## Axis 1 and 3

colnames(taxa.col2.axis13df)

corrspp.col2.axis13 <- taxa.col2.axis13df %>% filter(X.taxa.col2.axis13.vectors..r > 0.2)
target13.c2 <- corrspp.col2.axis13$Taxa # string of the Family names


axis13.vectors.c2 <- taxa.col2 %>% select(all_of(target13.c2)) # make a matrix of just those

(nmds.c2.vectors.13 <- envfit(nms.col2$points, axis13.vectors.c2,
                              permutations = 999, choices = c(1,3)))                        


corr.c2.vectors.13 <- as.data.frame(nmds.c2.vectors.13$vectors$arrows*sqrt(nmds.c2.vectors.13$vectors$r)) #scaling vectors
corr.c2.vectors.13$Taxa <- rownames(corr.c2.vectors.13)


write.csv(corr.c2.vectors.13, "Data/Emerging/NMDS/Col.2/emerging_correlated_vector13.csv")


## Actual figure
nmds.col2.scores <- read.csv("Data/Emerging/NMDS/Col.2/emerging_col2_NMDSscores.csv")
nmds.col2.scores$Year <- as.factor(nmds.col2.scores$Year)

col2.axis12 <- read.csv("Data/Emerging/NMDS/Col.2/emerging_correlated_vector12.csv")
col2.axis13 <- read.csv("Data/Emerging/NMDS/Col.2/emerging_correlated_vector13.csv")



invert.12.c2 <- ggplot(data = nmds.col2.scores,
                       aes(x = NMDS1, y = NMDS2)) +
  geom_point(data = nmds.col2.scores, 
             aes(x = NMDS1, y = NMDS2, 
                 colour = Treatment, shape = Year),
             size = 4, stroke = 1.5) + # sites as points
  stat_ellipse(data = nmds.col2.scores, 
               aes(x = NMDS1,y = NMDS2,
                   linetype = Treatment, colour = Treatment), 
               size = 1, level = 0.9) + 
  geom_segment(data = col2.axis12, 
               aes(x = 0, xend = MDS1, y = 0, yend = MDS2), # adding in the vectors, c
               arrow = arrow(length = unit(0.5, "cm")), colour = "black") + # can add in geom_label or geom_text for labels
  theme_minimal() + # no background
  theme(panel.border = element_rect(fill = NA)) + # full square around figure
  xlab("NMDS 1") +
  ylab("NMDS 2") +
  #ylim(-1, 1.5) +
  #xlim(-1.45, 1) +
  #theme(legend.position = "none") +
  geom_text_repel(data = col2.axis12, 
                  aes(x = MDS1, y = MDS2, label = Taxa),
                  color="black",
                  size = 5) +
  scale_color_manual(values = c("#969696","#35978f", "#2166ac")) +
  scale_shape_manual(values = c(17, 16, 1, 17, 2, 18, 5)) 

invert.12.c2


invert.13.c2 <- ggplot(data = nmds.col2.scores,
                       aes(x = NMDS1, y = NMDS3)) +
  geom_point(data = nmds.col2.scores, 
             aes(x = NMDS1, y = NMDS3, 
                 colour = Treatment, shape = Year),
             size = 4, stroke = 1.5) + # sites as points
  stat_ellipse(data = nmds.col2.scores, 
               aes(x = NMDS1,y = NMDS3,
                   linetype = Treatment, colour = Treatment), 
               size = 1, level = 0.9) + 
  geom_segment(data = col2.axis13, 
               aes(x = 0, xend = MDS1, y = 0, yend = MDS3), # adding in the vectors, c
               arrow = arrow(length = unit(0.5, "cm")), colour = "black") + # can add in geom_label or geom_text for labels
  theme_minimal() + # no background
  theme(panel.border = element_rect(fill = NA)) + # full square around figure
  xlab("NMDS 1") +
  ylab("NMDS 2") +
  #ylim(-1, 1.5) +
  #xlim(-1.45, 1) +
  #theme(legend.position = "none") +
  geom_text_repel(data = col2.axis13, 
                  aes(x = MDS1, y = MDS3, label = Taxa),
                  color="black",
                  size = 5) +
  scale_color_manual(values = c("#969696","#35978f", "#2166ac")) +
  scale_shape_manual(values = c(17, 16, 1, 17, 2, 18, 5)) 

invert.13.c2


(NMS.emerging.panel.c2 <- ggarrange(invert.12.c2, invert.13.c2,
                                    common.legend = TRUE,
                                    legend = "none",
                                    labels = c("B", ""),
                                    align = "hv"))

NMDS.col2 <- annotate_figure(NMS.emerging.panel.c2,
                top = text_grob("Collection 28-Jun-17 and 25-JUN-18"))


ggsave("Figures/NMDS_emerging_28Jun17_15Jun18.jpeg", NMDS.col2)









# NMDS Collection 08-Jul-17 and 04-Jul-18 -----------------------------------

col.3 <- c("08-Jul-17","04-Jul-18")

invert.col.3 <- inverts.date %>% filter(Date %in% col.3)
write.csv(invert.col.3, "Data/Emerging/NMDS/Col.3/inverts_collection3.csv")

invert.col3.rares <- read.csv("Data/Emerging/NMDS/Col.3/inverts_collection3_zeros.csv")

# just taxa and env 
colnames(invert.col3.rares)
taxa.col3 <- invert.col3.rares %>% select(Araneae:Crambidae)
env.col3 <- invert.col3.rares %>% select(ID:YrCol)

# Relativize by column max

taxa.col1re3 <- decostand(taxa.col3, "max", 2, na.rm = NULL)

write.csv(taxa.col1re3, "Data/Emerging/NMDS/Col.3/inverts_col3_raresrel.csv")

# NMDS collection 3

data.col3 <- read.csv("Data/Emerging/NMDS/Col.3/inverts_col3_raresrel.csv")

k_vec <- 1:10 #dimensions 1 - 10
stress <- numeric(length(k_vec)) # stress of each model put here
dune_dij <- metaMDSdist(data.col3)

set.seed(25)

for(i in seq_along(k_vec)) {
  sol <- metaMDSiter(dune_dij, k = i, 
                     trace = FALSE)
  stress[i] <- sol$stress
}
plot(stress) # 3D


# Actual NMDS

set.seed(120) 

nms.col3 <- metaMDS(data.col3, distance = "bray", # species data, bray-curtis dissimilarity
                    autotransform = FALSE,  # NMDS will do autotransformations for you
                    k = 3, trymax = 1000)   # k = number of axes
nms.col3

#Dimensions: 3 
#Stress:     0.03830475 
#Stress type 1, weak ties
#Two convergent solutions found after 33 tries
#Scaling: centring, PC rotation, halfchange scaling 
#Species: expanded scores based on ‘data.col3’ 

nms.col3$iters #200

nms.col3$stress^2   #0.001467254
1-nms.col3$stress^2 #0.9985327

scr3 <- as.data.frame(scores(nms.col3, display = "sites")) # extract NMDS scores
colnames(scr3)

env.col3$NMDS1 <- scr3$NMDS1
env.col3$NMDS2 <- scr3$NMDS2
env.col3$NMDS3 <- scr3$NMDS3

write.csv(env.col3,"Data/Emerging/NMDS/Col.3/emerging_col3_NMDSscores.csv") # save this as a csv

# Vector Taxa

# Axis 1 and 2
taxa.col3.axis12 <- envfit(nms.col3, taxa.col3,
                           choices = c(1,2))

taxa.col3.axis12df <- data.frame((taxa.col3.axis12$vectors)$arrows,
                                 (taxa.col3.axis12$vectors)$r,
                                 (taxa.col3.axis12$vectors)$pvals)

taxa.col3.axis12df <- tibble::rownames_to_column(taxa.col3.axis12df, "Taxa")

write.csv(taxa.col3.axis12df, "Data/Emerging/NMDS/Col.3/col3_allvectors_axis12.csv") # save vecto


# Axis 1 and 3


taxa.col3.axis13 <- envfit(nms.col3, taxa.col3,
                           choices = c(1,3))

taxa.col3.axis13df <- data.frame((taxa.col3.axis13$vectors)$arrows,
                                 (taxa.col3.axis13$vectors)$r,
                                 (taxa.col3.axis13$vectors)$pvals)

taxa.col3.axis13df <- tibble::rownames_to_column(taxa.col3.axis13df, "Taxa")

write.csv(taxa.col3.axis13df, "Data/Emerging/NMDS/Col.3/col3_allvectors_axis13.csv") # save vector scores as csv

# Correlated taxa

# Axis 1 and 2
colnames(taxa.col3.axis12df)

corrspp.col3.axis12 <- taxa.col3.axis12df %>% filter(X.taxa.col3.axis12.vectors..r > 0.2)
target12.c3 <- corrspp.col3.axis12$Taxa # string of the Family names


axis12.vectors.c3 <- taxa.col3 %>% select(all_of(target12.c3)) # make a matrix of just those

(nmds.c3.vectors.12 <- envfit(nms.col3$points, axis12.vectors.c3,
                              permutations = 999, choices = c(1,2)))                        


corr.c3.vectors.12 <- as.data.frame(nmds.c3.vectors.12$vectors$arrows*sqrt(nmds.c3.vectors.12$vectors$r)) #scaling vectors
corr.c3.vectors.12$Taxa <- rownames(corr.c3.vectors.12)

write.csv(corr.c3.vectors.12, "Data/Emerging/NMDS/Col.3/emerging_correlated_vector12.csv")

# Axis 1 and 3

colnames(taxa.col3.axis13df)

corrspp.col3.axis13 <- taxa.col3.axis13df %>% filter(X.taxa.col3.axis13.vectors..r > 0.2)
target13.c3 <- corrspp.col3.axis13$Taxa # string of the Family names


axis13.vectors.c3 <- taxa.col3 %>% select(all_of(target13.c3)) # make a matrix of just those

(nmds.c3.vectors.13 <- envfit(nms.col3$points, axis13.vectors.c3,
                              permutations = 999, choices = c(1,3)))                        


corr.c3.vectors.13 <- as.data.frame(nmds.c3.vectors.13$vectors$arrows*sqrt(nmds.c3.vectors.13$vectors$r)) #scaling vectors
corr.c3.vectors.13$Taxa <- rownames(corr.c3.vectors.13)

write.csv(corr.c3.vectors.13, "Data/Emerging/NMDS/Col.3/emerging_correlated_vector13.csv")

## Figure 
nmds.col3.scores <- read.csv("Data/Emerging/NMDS/Col.3/emerging_col3_NMDSscores.csv")
nmds.col3.scores$Year <- as.factor(nmds.col3.scores$Year)

col3.axis12 <- read.csv("Data/Emerging/NMDS/Col.3/emerging_correlated_vector12.csv")
col3.axis13 <- read.csv("Data/Emerging/NMDS/Col.3/emerging_correlated_vector13.csv")



invert.12.c3 <- ggplot(data = nmds.col3.scores,
                       aes(x = NMDS1, y = NMDS2)) +
  geom_point(data = nmds.col3.scores, 
             aes(x = NMDS1, y = NMDS2, 
                 colour = Treatment, shape = Year),
             size = 4, stroke = 1.5) + # sites as points
  stat_ellipse(data = nmds.col3.scores, 
               aes(x = NMDS1,y = NMDS2,
                   linetype = Treatment, colour = Treatment), 
               size = 1, level = 0.9) + 
  geom_segment(data = col3.axis12, 
               aes(x = 0, xend = MDS1, y = 0, yend = MDS2), # adding in the vectors, c
               arrow = arrow(length = unit(0.5, "cm")), colour = "black") + # can add in geom_label or geom_text for labels
  theme_minimal() + # no background
  theme(panel.border = element_rect(fill = NA)) + # full square around figure
  xlab("NMDS 1") +
  ylab("NMDS 2") +
  #ylim(-1, 1.5) +
  #xlim(-1.45, 1) +
  #theme(legend.position = "none") +
  geom_text_repel(data = col3.axis12, 
                  aes(x = MDS1, y = MDS2, label = Taxa),
                  color="black",
                  size = 5) +
  scale_color_manual(values = c("#969696","#35978f", "#2166ac")) +
  scale_shape_manual(values = c(17, 16, 1, 17, 2, 18, 5)) 

invert.12.c3



invert.13.c3 <- ggplot(data = nmds.col3.scores,
                       aes(x = NMDS1, y = NMDS3)) +
  geom_point(data = nmds.col3.scores, 
             aes(x = NMDS1, y = NMDS3, 
                 colour = Treatment, shape = Year),
             size = 4, stroke = 1.5) + # sites as points
  stat_ellipse(data = nmds.col3.scores, 
               aes(x = NMDS1,y = NMDS3,
                   linetype = Treatment, colour = Treatment), 
               size = 1, level = 0.9) + 
  geom_segment(data = col3.axis13, 
               aes(x = 0, xend = MDS1, y = 0, yend = MDS3), # adding in the vectors, c
               arrow = arrow(length = unit(0.5, "cm")), colour = "black") + # can add in geom_label or geom_text for labels
  theme_minimal() + # no background
  theme(panel.border = element_rect(fill = NA)) + # full square around figure
  xlab("NMDS 1") +
  ylab("NMDS 2") +
  #ylim(-1, 1.5) +
  #xlim(-1.45, 1) +
  #theme(legend.position = "none") +
  geom_text_repel(data = col3.axis13, 
                  aes(x = MDS1, y = MDS3, label = Taxa),
                  color="black",
                  size = 5) +
  scale_color_manual(values = c("#969696","#35978f", "#2166ac")) +
  scale_shape_manual(values = c(17, 16, 1, 17, 2, 18, 5)) 

invert.13.c3


(NMS.emerging.panel.c3 <- ggarrange(invert.12.c3, invert.13.c3,
                                    common.legend = TRUE,
                                    legend = "none",
                                    labels = c("C", ""),
                                    align = "hv"))

NMDS.col3 <- annotate_figure(NMS.emerging.panel.c3,
                             top = text_grob("Collection 08-Jul-17 and 04-JUL-18"))


ggsave("Figures/NMDS_emerging_08Jul17_04Jul18.jpeg", NMDS.col3)



# NMDS Collection 21-Jul-17 and 23-Jul-18 -----------------------------------

col.4 <- c("21-Jul-17","23-Jul-18")

invert.col.4 <- inverts.date %>% filter(Date %in% col.4)

write.csv(invert.col.4, "Data/Emerging/NMDS/Col.4/inverts_collection4.csv")

invert.col4.rares <- read.csv("Data/Emerging/NMDS/Col.4/inverts_collection4_zeros.csv")


# just taxa and env 
colnames(invert.col4.rares)
taxa.col4 <- invert.col4.rares %>% select(Araneae:Crambidae)
env.col4 <- invert.col4.rares %>% select(ID:YrCol)

# Relativize by column max

taxa.col4rel <- decostand(taxa.col4, "max", 2, na.rm = NULL)

write.csv(taxa.col4rel, "Data/Emerging/NMDS/Col.4/inverts_col4_raresrel.csv")


# NMDS collection 4

data.col4 <- read.csv("Data/Emerging/NMDS/Col.4/inverts_col4_raresrel.csv")

k_vec <- 1:10 #dimensions 1 - 10
stress <- numeric(length(k_vec)) # stress of each model put here
dune_dij <- metaMDSdist(data.col4)

set.seed(25)

for(i in seq_along(k_vec)) {
  sol <- metaMDSiter(dune_dij, k = i, 
                     trace = FALSE)
  stress[i] <- sol$stress
}
plot(stress) # 3D


# actual NMDS 

set.seed(120) 

nms.col4 <- metaMDS(data.col4, distance = "bray", # species data, bray-curtis dissimilarity
                    autotransform = FALSE,  # NMDS will do autotransformations for you
                    k = 3, trymax = 1000)   # k = number of axes
nms.col4

#global Multidimensional Scaling using monoMDS

#Data:     data.col4 
#Distance: bray 

#Dimensions: 3 
#Stress:     0.03259949 
#Stress type 1, weak ties
#Two convergent solutions found after 260 tries
#Scaling: centring, PC rotation, halfchange scaling 
#Species: expanded scores based on ‘data.col4’

nms.col4$iters #200

nms.col4$stress^2   #0.001062727
1-nms.col4$stress^2 #0.9989373


scr4 <- as.data.frame(scores(nms.col4, display = "sites")) # extract NMDS scores
colnames(scr4)

env.col4$NMDS1 <- scr4$NMDS1
env.col4$NMDS2 <- scr4$NMDS2
env.col4$NMDS3 <- scr4$NMDS3

write.csv(env.col4,"Data/Emerging/NMDS/Col.4/emerging_col4_NMDSscores.csv") # save this as a csv

## Vectors 

# Axis 1 and 2
taxa.col4.axis12 <- envfit(nms.col4, taxa.col4,
                           choices = c(1,2))

taxa.col4.axis12df <- data.frame((taxa.col4.axis12$vectors)$arrows,
                                 (taxa.col4.axis12$vectors)$r,
                                 (taxa.col4.axis12$vectors)$pvals)

taxa.col4.axis12df <- tibble::rownames_to_column(taxa.col4.axis12df, "Taxa")

write.csv(taxa.col4.axis12df, "Data/Emerging/NMDS/Col.4/col4_allvectors_axis12.csv") # save vecto


# Axis 1 and 3

taxa.col4.axis13 <- envfit(nms.col4, taxa.col4,
                           choices = c(1,3))

taxa.col4.axis13df <- data.frame((taxa.col4.axis13$vectors)$arrows,
                                 (taxa.col4.axis13$vectors)$r,
                                 (taxa.col4.axis13$vectors)$pvals)

taxa.col4.axis13df <- tibble::rownames_to_column(taxa.col4.axis13df, "Taxa")

write.csv(taxa.col4.axis13df, "Data/Emerging/NMDS/Col.4/col4_allvectors_axis13.csv") # save vecto


# Axis 1 and 2
colnames(taxa.col4.axis12df)

corrspp.col4.axis12 <- taxa.col4.axis12df %>% filter(X.taxa.col4.axis12.vectors..r > 0.2)
target12.c4 <- corrspp.col4.axis12$Taxa # string of the Family names


axis12.vectors.c4 <- taxa.col4 %>% select(all_of(target12.c4)) # make a matrix of just those

(nmds.c4.vectors.12 <- envfit(nms.col4$points, axis12.vectors.c4,
                              permutations = 999, choices = c(1,2)))                        


corr.c4.vectors.12 <- as.data.frame(nmds.c4.vectors.12$vectors$arrows*sqrt(nmds.c4.vectors.12$vectors$r)) #scaling vectors
corr.c4.vectors.12$Taxa <- rownames(corr.c4.vectors.12)

write.csv(corr.c4.vectors.12, "Data/Emerging/NMDS/Col.4/emerging_correlated_vector12.csv")

# axis 1 and 3

colnames(taxa.col4.axis13df)

corrspp.col4.axis13 <- taxa.col4.axis13df %>% filter(X.taxa.col4.axis13.vectors..r > 0.2)
target13.c4 <- corrspp.col4.axis13$Taxa # string of the Family names


axis13.vectors.c4 <- taxa.col4 %>% select(all_of(target13.c4)) # make a matrix of just those

(nmds.c4.vectors.13 <- envfit(nms.col4$points, axis13.vectors.c4,
                              permutations = 999, choices = c(1,3)))                        


corr.c4.vectors.13 <- as.data.frame(nmds.c4.vectors.13$vectors$arrows*sqrt(nmds.c4.vectors.13$vectors$r)) #scaling vectors
corr.c4.vectors.13$Taxa <- rownames(corr.c4.vectors.13)

write.csv(corr.c4.vectors.13, "Data/Emerging/NMDS/Col.4/emerging_correlated_vector13.csv")


## NMDS figure

nmds.col4.scores <- read.csv("Data/Emerging/NMDS/Col.4/emerging_col4_NMDSscores.csv")
nmds.col4.scores$Year <- as.factor(nmds.col4.scores$Year)

col4.axis12 <- read.csv("Data/Emerging/NMDS/Col.4/emerging_correlated_vector12.csv")
col4.axis13 <- read.csv("Data/Emerging/NMDS/Col.4/emerging_correlated_vector13.csv")



invert.12.c4 <- ggplot(data = nmds.col4.scores,
                       aes(x = NMDS1, y = NMDS2)) +
  geom_point(data = nmds.col4.scores, 
             aes(x = NMDS1, y = NMDS2, 
                 colour = Treatment, shape = Year),
             size = 4, stroke = 1.5) + # sites as points
  stat_ellipse(data = nmds.col4.scores, 
               aes(x = NMDS1,y = NMDS2,
                   linetype = Treatment, colour = Treatment), 
               size = 1, level = 0.9) + 
  geom_segment(data = col4.axis12, 
               aes(x = 0, xend = MDS1, y = 0, yend = MDS2), # adding in the vectors, c
               arrow = arrow(length = unit(0.5, "cm")), colour = "black") + # can add in geom_label or geom_text for labels
  theme_minimal() + # no background
  theme(panel.border = element_rect(fill = NA)) + # full square around figure
  xlab("NMDS 1") +
  ylab("NMDS 2") +
  #ylim(-1, 1.5) +
  #xlim(-1.45, 1) +
  #theme(legend.position = "none") +
  geom_text_repel(data = col4.axis12, 
                  aes(x = MDS1, y = MDS2, label = Taxa),
                  color="black",
                  size = 5) +
  scale_color_manual(values = c("#969696","#35978f", "#2166ac")) +
  scale_shape_manual(values = c(17, 16, 1, 17, 2, 18, 5)) 

invert.12.c4


invert.13.c4 <- ggplot(data = nmds.col4.scores,
                       aes(x = NMDS1, y = NMDS3)) +
  geom_point(data = nmds.col4.scores, 
             aes(x = NMDS1, y = NMDS3, 
                 colour = Treatment, shape = Year),
             size = 4, stroke = 1.5) + # sites as points
  stat_ellipse(data = nmds.col4.scores, 
               aes(x = NMDS1,y = NMDS3,
                   linetype = Treatment, colour = Treatment), 
               size = 1, level = 0.9) + 
  geom_segment(data = col4.axis13, 
               aes(x = 0, xend = MDS1, y = 0, yend = MDS3), # adding in the vectors, c
               arrow = arrow(length = unit(0.5, "cm")), colour = "black") + # can add in geom_label or geom_text for labels
  theme_minimal() + # no background
  theme(panel.border = element_rect(fill = NA)) + # full square around figure
  xlab("NMDS 1") +
  ylab("NMDS 3") +
  #ylim(-1, 1.5) +
  #xlim(-1.45, 1) +
  #theme(legend.position = "none") +
  geom_text_repel(data = col4.axis13, 
                  aes(x = MDS1, y = MDS3, label = Taxa),
                  color="black",
                  size = 5) +
  scale_color_manual(values = c("#969696","#35978f", "#2166ac")) +
  scale_shape_manual(values = c(17, 16, 1, 17, 2, 18, 5)) 

invert.13.c4



(NMS.emerging.panel.c4 <- ggarrange(invert.12.c4, invert.13.c4,
                                    common.legend = TRUE,
                                    legend = "bottom",
                                    labels = c("D", ""),
                                    align = "hv"))

NMDS.col4 <- annotate_figure(NMS.emerging.panel.c4,
                             top = text_grob("Collection 21-Jul-17 and 23-JUL-18"))


ggsave("Figures/NMDS_emerging_21Jul17_23Jul18.jpeg", NMDS.col3)


ggarrange(NMDS.col1, NMDS.col2,
          NMDS.col3, NMDS.col4)

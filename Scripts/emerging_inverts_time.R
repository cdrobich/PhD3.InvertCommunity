
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

library(Hmisc)

# Load data ---------------------------------------------------------------

inverts.date <- read.csv("Data/Emerging/ermerging_time_rares.csv") # Occurences <= 2 removed
head(inverts.date)
str(inverts.date)

unique(inverts.date$Treatment)
colnames(inverts.date)

taxa <- inverts.date %>% select(Araneae:Crambidae)
env <- inverts.date %>% select(ID:YrCol)

# Collection periods -----------------------------------------------

env <- env %>% separate(Date, c("Day","Month", "Year"), remove = FALSE)

unique(env$Date)

#                        "19-Jun-17" "28-Jun-17" "08-Jul-17" "21-Jul-17" 
#"20-May-18" "05-Jun-18" "16-Jun-18" "25-Jun-18" "04-Jul-18" "23-Jul-18"


# NMDS Collection June 19 2017 and June 16 2018 Collection -----------------------------------------

col.1 <- c("19-Jun-17","16-Jun-18")

invert.col.1 <- inverts.date %>% filter(Date %in% col.1)

write.csv(invert.col.1, "Data/Emerging/NMDS/Col.1/inverts_collection1.csv")


invert.col1.rares <- read.csv("Data/Emerging/NMDS/Col.1/inverts_collection1_zerosremoved.csv")


# Empty columns removed, relativized by col max. 

col1.data <- read.csv("Data/Emerging/NMDS/Col.1/col1_zeroes_relativized.csv")
colnames(col1.data)

taxa.col1 <- col1.data %>% select(Araneae:Crambidae)
env.col1 <- col1.data %>% select(Sites:YrCol)

taxa.col1 <- taxa.col1[1:51,]
env.col1 <- env.col1[1:51,]
# NMDS Col.1 

#### NMDS analysis 

set.seed(120) 

nms.col1 <- metaMDS(taxa.col1, distance = "bray", # species data, bray-curtis dissimilarity
                      autotransform = FALSE,  # NMDS will do autotransformations for you
                      k = 3, trymax = 1000)   # k = number of axes
nms.col1

#Data:     taxa.col1 
#Distance: bray 

#Dimensions: 3 
#Stress:     0.1933805 
#Stress type 1, weak ties
#Two convergent solutions found after 26 tries
#Scaling: centring, PC rotation, halfchange scaling 
#Species: expanded scores based on ‘taxa.col1’ 

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
sum(g^2) #0.037396
nms.col1$stress^2  

1-nms.col1$stress^2 # 0.962604


## extract the scores for plotting 
scr <- as.data.frame(scores(nms.col1, display = "sites")) # extract NMDS scores
colnames(scr)

# adding categorical info to scores
env.col1$NMDS1 <- scr$NMDS1
env.col1$NMDS2 <- scr$NMDS2
env.col1$NMDS3 <- scr$NMDS3

write.csv(env.col1,"Data/Emerging/NMDS/Col.1/emerging_inverts_col1_NMDSscores.csv") # save this as a csv

### Taxa correlated with Axis 1 & 2 

col1.12 <- envfit(nms.col1, taxa.col1,
                    choices = c(1,2)) #produces a list with r2, p value, and NMDS coordinates

col1.df12 <- data.frame((col1.12$vectors)$arrows,
                          (col1.12$vectors)$r,
                          (col1.12$vectors)$pvals) #take list and make into data frame

col1.df12 <- tibble::rownames_to_column(col1.df12, "Taxa")

write.csv(col1.df12, "Data/Emerging/NMDS/Col.1/NMDS_emerg_col1_vector12.csv") # save vector scores as csv


## Taxa correlated with Axis 1 & 3

col1.13 <- envfit(nms.col1, taxa.col1, 
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

axis12.vectors.c1 <- taxa.col1 %>% select(all_of(target12.c1)) # make a matrix of just those


# Selecting taxa for vectors on axis 1, 3
colnames(col1.df13)
col1.corrsp.13 <- col1.df13 %>% filter(X.col1.13.vectors..r > 0.2) 
target13.c1 <- col1.corrsp.13$Taxa # string of the Family names

axis13.vectors.c1 <- taxa.col1 %>% select(all_of(target13.c1)) # make a matrix of just those


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

# NMDS Emerging Invert Figure Collection 1 
col1.scores <- read.csv("Data/Emerging/NMDS/Col.1/emerging_inverts_col1_NMDSscores.csv")
col1.scores$Year <- as.factor(col1.scores$Year)

col1.axis12 <- read.csv("Data/Emerging/NMDS/Col.1/NMDS_emerg_col1_vector12_final.csv")
col1.axis13 <- read.csv("Data/Emerging/NMDS/Col.1/NMDS_emerg_col1_vector13_final.csv")


## NMDS Axis 1, 2 


invert.12.c1 <- ggplot(data = col1.scores,
                    aes(x = NMDS1, y = NMDS2)) +
  geom_point(data = col1.scores, 
             aes(x = NMDS1, y = NMDS2, 
                 fill = Treatment, shape = Year),
             size = 5, stroke = 1.5,
             alpha = 0.8) + # sites as points
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
  geom_label_repel(data = col1.axis12, 
                  aes(x = MDS1, y = MDS2, label = Taxa),
                  color="black",
                  size = 6)  +
  scale_fill_viridis(discrete = TRUE) +
  scale_colour_viridis(discrete = TRUE) +
  scale_shape_manual(values = c(21, 24, 22)) +
  theme(legend.position = "none") +
  ylim(-0.8, 0.8) +
  xlim(-0.8, 0.8)

invert.12.c1

## NMDS Axis 1, 3

invert.13.c1 <- ggplot(data = col1.scores,
                       aes(x = NMDS1, y = NMDS3)) +
  geom_point(data = col1.scores, 
             aes(x = NMDS1, y = NMDS3, 
                 fill = Treatment, shape = Year),
             size = 5, stroke = 1.5,
             alpha = 0.7) + # sites as points
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
  geom_label_repel(data = col1.axis13, 
                  aes(x = MDS1, y = MDS3, label = Taxa),
                  color="black",
                  size = 6) +
  scale_fill_viridis(discrete = TRUE) +
  scale_colour_viridis(discrete = TRUE) +
  scale_shape_manual(values = c(21, 24, 22)) +
  theme(legend.position = "none") +
  ylim(-0.8, 0.8) +
  xlim(-0.8, 0.8)

invert.13.c1

(NMS.emerging.panel.c1 <- ggarrange(invert.12.c1, invert.13.c1,
                                 common.legend = TRUE,
                                 legend = "none",
                                 labels = c("C", ""),
                                 align = "hv"))

NMDS.col1 <- annotate_figure(NMS.emerging.panel.c1,
                top = text_grob("Collection 19-Jun-2017 and 16-Jun-2018"))

ggsave("Figures/NMDS_emerging_19Jun17_16Jun18.jpeg", NMDS.col1)



# NMDS Collection 28-Jun-17 and 25-JUN-18 ---------------------------------

col.2 <- c("28-Jun-17","25-Jun-18")

invert.col.2 <- inverts.date %>% filter(Date %in% col.2)

write.csv(invert.col.2, "Data/Emerging/NMDS/Col.2/inverts_collection2.csv")

# relativized by col max and removed empty columns 

col2.data <- read.csv("Data/Emerging/NMDS/Col.2/col2_zeroes_relativized.csv")

# just taxa and env 
colnames(col2.data)

taxa.col2 <- col2.data %>% select(Araneae:Noctuidae)
env.col2 <- col2.data %>% select(Site:YrCol)


### NMDS "Col 2" 

k_vec <- 1:10 #dimensions 1 - 10
stress <- numeric(length(k_vec)) # stress of each model put here
dune_dij <- metaMDSdist(taxa.col2)

set.seed(25)

for(i in seq_along(k_vec)) {
  sol <- metaMDSiter(dune_dij, k = i, 
                     trace = FALSE)
  stress[i] <- sol$stress
}
plot(stress) # 3D




set.seed(120) 

nms.col2 <- metaMDS(taxa.col2, distance = "bray", # species data, bray-curtis dissimilarity
                    autotransform = FALSE,  # NMDS will do autotransformations for you
                    k = 3, trymax = 1000)   # k = number of axes
nms.col2


#global Multidimensional Scaling using monoMDS

#Data:     taxa.col2 
#Distance: bray 

#Dimensions: 3 
#Stress:     0.2153085 
#Stress type 1, weak ties
#Two convergent solutions found after 20 tries

# look at points and the stress real quick
layout(matrix(1:2, ncol = 2))
plot(nms.col2, main = "Invertebrate NMDS plot"); stressplot(nms.col2, main = "Shepard plot")
layout(1)

ordiplot(nms.col2, type = "n")
orditorp(nms.col2, display = "species")
orditorp(nms.col2, display = "sites")



nms.col2$iters #97

nms.col2$stress^2   #0.04635775
1-nms.col2$stress^2 #0.9536422

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
                 fill = Treatment, shape = Year),
             size = 4, stroke = 1.5,
             alpha = 0.7) + # sites as points
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
  geom_label_repel(data = col2.axis12, 
                  aes(x = MDS1, y = MDS2, label = Taxa),
                  color="black",
                  size = 5) +
  scale_fill_viridis(discrete = TRUE) +
  scale_colour_viridis(discrete = TRUE) +
  scale_shape_manual(values = c(21, 24, 22)) +
  theme(legend.position = "none") +
  ylim(-0.8, 0.8) +
  xlim(-0.8, 0.8)

invert.12.c2


invert.13.c2 <- ggplot(data = nmds.col2.scores,
                       aes(x = NMDS1, y = NMDS3)) +
  geom_point(data = nmds.col2.scores, 
             aes(x = NMDS1, y = NMDS3, 
                 fill = Treatment, shape = Year),
             size = 4, stroke = 1.5,
             alpha = 0.7) + # sites as points
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
  ylab("NMDS 3") +
  #ylim(-1, 1.5) +
  #xlim(-1.45, 1) +
  #theme(legend.position = "none") +
  geom_label_repel(data = col2.axis13, 
                  aes(x = MDS1, y = MDS3, label = Taxa),
                  color="black",
                  size = 5) +
  scale_fill_viridis(discrete = TRUE) +
  scale_colour_viridis(discrete = TRUE) +
  scale_shape_manual(values = c(21, 24, 22)) +
  theme(legend.position = "none") +
  ylim(-0.8, 0.8) +
  xlim(-0.8, 0.8)

invert.13.c2


(NMS.emerging.panel.c2 <- ggarrange(invert.12.c2, invert.13.c2,
                                    common.legend = TRUE,
                                    legend = "none",
                                    labels = c("D", ""),
                                    align = "hv"))

(NMDS.col2 <- annotate_figure(NMS.emerging.panel.c2,
                top = text_grob("Collection 28-Jun-17 and 25-JUN-18")))


ggsave("Figures/NMDS_emerging_28Jun17_15Jun18.jpeg", NMDS.col2)









# NMDS Collection 08-Jul-17 and 04-Jul-18 -----------------------------------

col.3 <- c("08-Jul-17","04-Jul-18")

invert.col.3 <- inverts.date %>% filter(Date %in% col.3)
write.csv(invert.col.3, "Data/Emerging/NMDS/Col.3/inverts_collection3.csv")


col3.data <- read.csv("Data/Emerging/NMDS/Col.3/col3_zeros_rares.csv")

# just taxa and env 
colnames(col3.data)
taxa.col3 <- col3.data %>% select(Araneae:Crambidae)
env.col3 <- col3.data %>% select(ID:YrCol)

# NMDS collection 3

k_vec <- 1:10 #dimensions 1 - 10
stress <- numeric(length(k_vec)) # stress of each model put here
dune_dij <- metaMDSdist(taxa.col3)

set.seed(25)

for(i in seq_along(k_vec)) {
  sol <- metaMDSiter(dune_dij, k = i, 
                     trace = FALSE)
  stress[i] <- sol$stress
}
plot(stress) # 3D


# Actual NMDS

set.seed(120) 

nms.col3 <- metaMDS(taxa.col3, distance = "bray", # species data, bray-curtis dissimilarity
                    autotransform = FALSE,  # NMDS will do autotransformations for you
                    k = 3, trymax = 1000)   # k = number of axes
nms.col3

#global Multidimensional Scaling using monoMDS

#Data:     taxa.col3 
#Distance: bray 

#Dimensions: 3 
#Stress:     0.2088858 
#Stress type 1, weak ties
#Two convergent solutions found after 43 tries 

nms.col3$iters #126

nms.col3$stress^2   #0.04363328
1-nms.col3$stress^2 #0.9563667

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
                 fill = Treatment, shape = Year),
             size = 5, stroke = 1.5,
             alpha = 0.7) + # sites as points
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
  geom_label_repel(data = col3.axis12, 
                  aes(x = MDS1, y = MDS2, label = Taxa),
                  color="black",
                  size = 5) +
  scale_fill_viridis(discrete = TRUE) +
  scale_colour_viridis(discrete = TRUE) +
  scale_shape_manual(values = c(21, 24, 22)) +
  theme(legend.position = "none") +
  ylim(-0.8, 0.8) +
  xlim(-0.8, 0.8)

invert.12.c3



invert.13.c3 <- ggplot(data = nmds.col3.scores,
                       aes(x = NMDS1, y = NMDS3)) +
  geom_point(data = nmds.col3.scores, 
             aes(x = NMDS1, y = NMDS3, 
                 fill = Treatment, shape = Year),
             size = 5, stroke = 1.5,
             alpha = 0.7) + # sites as points
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
  ylab("NMDS 3") +
  #ylim(-1, 1.5) +
  #xlim(-1.45, 1) +
  #theme(legend.position = "none") +
  geom_label_repel(data = col3.axis13, 
                  aes(x = MDS1, y = MDS3, label = Taxa),
                  color="black",
                  size = 5) +
  scale_fill_viridis(discrete = TRUE) +
  scale_colour_viridis(discrete = TRUE) +
  scale_shape_manual(values = c(21, 24, 22)) +
  theme(legend.position = "none") +
  ylim(-0.8, 0.8) +
  xlim(-0.8, 0.8)


invert.13.c3


(NMS.emerging.panel.c3 <- ggarrange(invert.12.c3, invert.13.c3,
                                    common.legend = TRUE,
                                    legend = "none",
                                    labels = c("E", ""),
                                    align = "hv"))

NMDS.col3 <- annotate_figure(NMS.emerging.panel.c3,
                             top = text_grob("Collection 08-Jul-17 and 04-JUL-18"))


ggsave("Figures/NMDS_emerging_08Jul17_04Jul18.jpeg", NMDS.col3)



# NMDS Collection 21-Jul-17 and 23-Jul-18 -----------------------------------

col.4 <- c("21-Jul-17","23-Jul-18")

invert.col.4 <- inverts.date %>% filter(Date %in% col.4)

write.csv(invert.col.4, "Data/Emerging/NMDS/Col.4/inverts_collection4.csv")

col4.data <- read.csv("Data/Emerging/NMDS/Col.4/col4_zeros_rares.csv")


# just taxa and env 
colnames(col4.data)
taxa.col4 <- col4.data  %>% select(Araneae:Crambidae)
env.col4 <- col4.data  %>% select(ID:YrCol)


# NMDS collection 4
k_vec <- 1:10 #dimensions 1 - 10
stress <- numeric(length(k_vec)) # stress of each model put here
dune_dij <- metaMDSdist(taxa.col4)

set.seed(25)

for(i in seq_along(k_vec)) {
  sol <- metaMDSiter(dune_dij, k = i, 
                     trace = FALSE)
  stress[i] <- sol$stress
}
plot(stress) # 3D


# actual NMDS 

set.seed(120) 

nms.col4 <- metaMDS(taxa.col4, distance = "bray", # species data, bray-curtis dissimilarity
                    autotransform = FALSE,  # NMDS will do autotransformations for you
                    k = 3, trymax = 1000)   # k = number of axes
nms.col4

#global Multidimensional Scaling using monoMDS

#Data:     taxa.col4 
#Distance: bray 

#Dimensions: 3 
#Stress:     0.2204784 
#Stress type 1, weak ties
#Two convergent solutions found after 20 tries

nms.col4$iters #192

nms.col4$stress^2   #0.04861074
1-nms.col4$stress^2 #0.9513893


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
                 fill = Treatment, shape = Year),
             size = 5, stroke = 1.5,
             alpha = 0.7) + # sites as points
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
  geom_label_repel(data = col4.axis12, 
                  aes(x = MDS1, y = MDS2, label = Taxa),
                  color="black",
                  size = 5) +
  scale_fill_viridis(discrete = TRUE) +
  scale_colour_viridis(discrete = TRUE) +
  scale_shape_manual(values = c(21, 24)) +
  theme(legend.position = "none") +
  ylim(-0.6, 0.6) +
  xlim(-0.6, 0.6)

invert.12.c4


invert.13.c4 <- ggplot(data = nmds.col4.scores,
                       aes(x = NMDS1, y = NMDS3)) +
  geom_point(data = nmds.col4.scores, 
             aes(x = NMDS1, y = NMDS3, 
                 fill = Treatment, shape = Year),
             size = 4, stroke = 1.5,
             alpha = 0.7) + # sites as points
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
  geom_label_repel(data = col4.axis13, 
                  aes(x = MDS1, y = MDS3, label = Taxa),
                  color="black",
                  size = 5) +
  scale_fill_viridis(discrete = TRUE) +
  scale_colour_viridis(discrete = TRUE) +
  scale_shape_manual(values = c(21, 24)) +
  theme(legend.position = "none") +
  ylim(-0.6, 0.6) +
  xlim(-0.6, 0.6)

invert.13.c4



(NMS.emerging.panel.c4 <- ggarrange(invert.12.c4, invert.13.c4,
                                    labels = c("F", ""),
                                    align = "hv"))

NMDS.col4 <- annotate_figure(NMS.emerging.panel.c4,
                             top = text_grob("Collection 21-Jul-17 and 23-JUL-18"))


ggsave("Figures/NMDS_emerging_21Jul17_23Jul18.jpeg", NMDS.col3)



# May 20 2018 -------------------------------------------------------------

invert.May <- inverts.date %>% filter(Date == "20-May-18")

write.csv(invert.May, "Data/Emerging/NMDS/inverts_20-May-18.csv")



may.data <- read.csv("Data/Emerging/NMDS/may_zeros_rares.csv")
colnames(may.data)

taxa.may <- may.data %>% select(Araneae:Limnephilida)
env.may <- may.data %>% select(Sites:YrCol)




k_vec <- 1:10 #dimensions 1 - 10
stress <- numeric(length(k_vec)) # stress of each model put here
dune_dij <- metaMDSdist(taxa.may)

set.seed(25)

for(i in seq_along(k_vec)) {
  sol <- metaMDSiter(dune_dij, k = i, 
                     trace = FALSE)
  stress[i] <- sol$stress
}
plot(stress) # Going with two based on PCORD



set.seed(120) 

nms.may <- metaMDS(taxa.may, distance = "bray", # species data, bray-curtis dissimilarity
                    autotransform = FALSE,  # NMDS will do autotransformations for you
                    k = 2, trymax = 1000)   # k = number of axes
nms.may

#Data:     taxa.may 
#Distance: bray 

#Dimensions: 2 
#Stress:     0.2425735 
#Stress type 1, weak ties
#Two convergent solutions found after 20 tries
#Scaling: centring, PC rotation, halfchange scaling 
#Species: expanded scores based on ‘taxa.may’

layout(matrix(1:2, ncol = 2))
plot(nms.may, main = "Invertebrate NMDS plot"); stressplot(nms.may, main = "Shepard plot")
layout(1)

ordiplot(nms.may, type = "n")
orditorp(nms.may, display = "species")
orditorp(nms.may, display = "sites")


nms.may$iters #82

nms.may$stress^2   # 0.0588419
1-nms.may$stress^2 # 0.9411581

scr.may <- as.data.frame(scores(nms.may, display = "sites")) # extract NMDS scores
colnames(scr.may)

env.may$NMDS1 <- scr.may$NMDS1
env.may$NMDS2 <- scr.may$NMDS2
env.may$NMDS3 <- scr.may$NMDS3


write.csv(env.may,"Data/Emerging/NMDS/emerging_May18_NMDSscores.csv") # save this as a csv



## Taxa for vectors

# Axis 1 and 2
taxa.may.axis12 <- envfit(nms.may, taxa.may,
                           choices = c(1,2))

taxa.may.axis12df <- data.frame((taxa.may.axis12$vectors)$arrows,
                                 (taxa.may.axis12$vectors)$r,
                                 (taxa.may.axis12$vectors)$pvals)

taxa.may.axis12df <- tibble::rownames_to_column(taxa.may.axis12df, "Taxa")

write.csv(taxa.may.axis12df, "Data/Emerging/NMDS/may_allvectors_axis12.csv") # save vector scores as csv

# Axis 1 and 3

taxa.may.axis13 <- envfit(nms.may, taxa.may,
                           choices = c(1,3))

taxa.may.axis13df <- data.frame((taxa.may.axis13$vectors)$arrows,
                                 (taxa.may.axis13$vectors)$r,
                                 (taxa.may.axis13$vectors)$pvals)

taxa.may.axis13df <- tibble::rownames_to_column(taxa.may.axis13df, "Taxa")

write.csv(taxa.may.axis13df, "Data/Emerging/NMDS/may_allvectors_axis13.csv") # save vector scores as csv


# Correlated taxa

# Axis 1 and 2

colnames(taxa.may.axis12df)

corrspp.may.axis12 <- taxa.may.axis12df %>% filter(X.taxa.may.axis12.vectors..r > 0.2)
target12.may <- corrspp.may.axis12$Taxa # string of the Family names


axis12.vectors.may <- taxa.may %>% select(all_of(target12.may)) # make a matrix of just those

(nmds.may.vectors.12 <- envfit(nms.may$points, axis12.vectors.may,
                              permutations = 999, choices = c(1,2)))                        


corr.may.vectors.12 <- as.data.frame(nmds.may.vectors.12$vectors$arrows*sqrt(nmds.may.vectors.12$vectors$r)) #scaling vectors
corr.may.vectors.12$Taxa <- rownames(corr.may.vectors.12)


write.csv(corr.may.vectors.12, "Data/Emerging/NMDS/May_emerging_correlated_vector12.csv")

## Axis 1 and 3

colnames(taxa.may.axis13df)

corrspp.may.axis13 <- taxa.may.axis13df %>% filter(X.taxa.may.axis13.vectors..r > 0.2)
target13.may <- corrspp.may.axis13$Taxa # string of the Family names


axis13.vectors.may <- taxa.may %>% select(all_of(target13.may)) # make a matrix of just those

(nmds.may.vectors.13 <- envfit(nms.may$points, axis13.vectors.may,
                              permutations = 999, choices = c(1,3)))                        


corr.may.vectors.13 <- as.data.frame(nmds.may.vectors.13$vectors$arrows*sqrt(nmds.may.vectors.13$vectors$r)) #scaling vectors
corr.may.vectors.13$Taxa <- rownames(corr.may.vectors.13)


write.csv(corr.may.vectors.13, "Data/Emerging/NMDS/May_emerging_correlated_vector13.csv")


## Actual figure
nmds.may.scores <- read.csv("Data/Emerging/NMDS/emerging_May18_NMDSscores.csv")
nmds.may.scores$Year <- as.factor(nmds.may.scores$Year)

may.axis12 <- read.csv("Data/Emerging/NMDS/May_emerging_correlated_vector12.csv")
may.axis13 <- read.csv("Data/Emerging/NMDS/May_emerging_correlated_vector13.csv")



invert.12.may <- ggplot(data = nmds.may.scores,
                       aes(x = NMDS1, y = NMDS2)) +
  geom_point(data = nmds.may.scores, 
             aes(x = NMDS1, y = NMDS2, 
                 fill = Treatment, shape = Year),
             size = 5, stroke = 1.5,
             alpha = 0.7) + # sites as points
  stat_ellipse(data = nmds.may.scores, 
               aes(x = NMDS1,y = NMDS2,
                   linetype = Treatment, colour = Treatment), 
               size = 1, level = 0.9) + 
  geom_segment(data = may.axis12, 
               aes(x = 0, xend = MDS1, y = 0, yend = MDS2), # adding in the vectors, c
               arrow = arrow(length = unit(0.5, "cm")), colour = "black") + # can add in geom_label or geom_text for labels
  theme_minimal() + # no background
  theme(panel.border = element_rect(fill = NA)) + # full square around figure
  xlab("NMDS 1") +
  ylab("NMDS 2") +
  #ylim(-1, 1.5) +
  #xlim(-1.45, 1) +
  #theme(legend.position = "none") +
  geom_label_repel(data = may.axis12, 
                  aes(x = MDS1, y = MDS2, label = Taxa),
                  color="black",
                  size = 5) +
  scale_fill_viridis(discrete = TRUE) +
  scale_colour_viridis(discrete = TRUE) +
  scale_shape_manual(values = c(24)) +
  theme(legend.position = "none") +
  ylim(-0.8, 0.8) +
  xlim(-0.8, 0.8) 

invert.12.may 

nmds.may <- ggarrange(legends, invert.12.may,
          labels = "A",
          widths = c(0.25,1))

NMDS.may <- annotate_figure(nmds.may,
                top = text_grob("Collection 20-May-18"))

ggsave("Figures/May_emerging_NMDS.jpeg")

# 5 June 2018 -------------------------------------------------------------

invert.June <- inverts.date %>% filter(Date == "05-Jun-18")

write.csv(invert.June, "Data/Emerging/NMDS/inverts_05-Jun-18.csv")



june.data <- read.csv("Data/Emerging/NMDS/june18_zeros_rel.csv")
colnames(june.data)

taxa.june <- june.data %>% select(Araneae:Crambidae)
env.june <- june.data %>% select(Sites:YrCol)



set.seed(120) 

nms.june <- metaMDS(taxa.june, distance = "bray", # species data, bray-curtis dissimilarity
                    autotransform = FALSE,  # NMDS will do autotransformations for you
                    k = 3, trymax = 1000)   # k = number of axes
nms.june

#global Multidimensional Scaling using monoMDS

#Data:     taxa.june 
#Distance: bray 

#Dimensions: 3 
#Stress:     0.1581744 
#Stress type 1, weak ties
#Two convergent solutions found after 48 tries
#Scaling: centring, PC rotation, halfchange scalin

layout(matrix(1:2, ncol = 2))
plot(nms.june , main = "Invertebrate NMDS plot"); stressplot(nms.june , main = "Shepard plot")
layout(1)

nms.june$iters #118

nms.june$stress^2   #0.02500559
1-nms.june$stress^2 #0.9749944


scr.june <- as.data.frame(scores(nms.june, display = "sites")) # extract NMDS scores
colnames(scr.june)

env.june$NMDS1 <- scr.june$NMDS1
env.june$NMDS2 <- scr.june$NMDS2
env.june$NMDS3 <- scr.june$NMDS3

write.csv(env.june,"Data/Emerging/NMDS/emerging_june_NMDSscores.csv") # save this as a csv



## Taxa for vectors

# Axis 1 and 2
taxa.june.axis12 <- envfit(nms.june, taxa.june,
                           choices = c(1,2))

taxa.june.axis12df <- data.frame((taxa.june.axis12$vectors)$arrows,
                                 (taxa.june.axis12$vectors)$r,
                                 (taxa.june.axis12$vectors)$pvals)

taxa.june.axis12df <- tibble::rownames_to_column(taxa.june.axis12df, "Taxa")

write.csv(taxa.june.axis12df, "Data/Emerging/NMDS/june_allvectors_axis12.csv") # save vector scores as csv

# Axis 1 and 3

taxa.june.axis13 <- envfit(nms.june, taxa.june,
                           choices = c(1,3))

taxa.june.axis13df <- data.frame((taxa.june.axis13$vectors)$arrows,
                                 (taxa.june.axis13$vectors)$r,
                                 (taxa.june.axis13$vectors)$pvals)

taxa.june.axis13df <- tibble::rownames_to_column(taxa.june.axis13df, "Taxa")

write.csv(taxa.june.axis13df, "Data/Emerging/NMDS/june_allvectors_axis13.csv") # save vector scores as csv



# Correlated taxa

# Axis 1 and 2

colnames(taxa.june.axis12df)

corrspp.june.axis12 <- taxa.june.axis12df %>% filter(X.taxa.june.axis12.vectors..r > 0.2)
target12.june <- corrspp.june.axis12$Taxa # string of the Family names


axis12.vectors.june <- taxa.june %>% select(all_of(target12.june)) # make a matrix of just those

(nmds.june.vectors.12 <- envfit(nms.june$points, axis12.vectors.june,
                              permutations = 999, choices = c(1,2)))                        


corr.june.vectors.12 <- as.data.frame(nmds.june.vectors.12$vectors$arrows*sqrt(nmds.june.vectors.12$vectors$r)) #scaling vectors
corr.june.vectors.12$Taxa <- rownames(corr.june.vectors.12)


write.csv(corr.june.vectors.12, "Data/Emerging/NMDS/june_emerging_correlated_vector12.csv")

## Axis 1 and 3

colnames(taxa.june.axis13df)

corrspp.june.axis13 <- taxa.june.axis13df %>% filter(X.taxa.june.axis13.vectors..r > 0.2)
target13.june <- corrspp.june.axis13$Taxa # string of the Family names


axis13.vectors.june <- taxa.june %>% select(all_of(target13.june)) # make a matrix of just those

(nmds.june.vectors.13 <- envfit(nms.june$points, axis13.vectors.june,
                              permutations = 999, choices = c(1,3)))                        


corr.june.vectors.13 <- as.data.frame(nmds.june.vectors.13$vectors$arrows*sqrt(nmds.june.vectors.13$vectors$r)) #scaling vectors
corr.june.vectors.13$Taxa <- rownames(corr.june.vectors.13)


write.csv(corr.june.vectors.13, "Data/Emerging/NMDS/june_emerging_correlated_vector13.csv")



## Actual figure
nmds.june.scores <- read.csv("Data/Emerging/NMDS/emerging_june_NMDSscores.csv")
nmds.june.scores$Year <- as.factor(nmds.june.scores$Year)

june.axis12 <- read.csv("Data/Emerging/NMDS/june_emerging_correlated_vector12.csv")
june.axis13 <- read.csv("Data/Emerging/NMDS/june_emerging_correlated_vector13.csv")



invert.12.june <- ggplot(data = nmds.june.scores,
                       aes(x = NMDS1, y = NMDS2)) +
  geom_point(data = nmds.june.scores, 
             aes(x = NMDS1, y = NMDS2, 
                 fill = Treatment, shape = Year),
             size = 4, stroke = 1.5,
             alpha = 0.7) + # sites as points
  stat_ellipse(data = nmds.june.scores, 
               aes(x = NMDS1,y = NMDS2,
                   linetype = Treatment, colour = Treatment), 
               size = 1, level = 0.9) + 
  geom_segment(data = june.axis12, 
               aes(x = 0, xend = MDS1, y = 0, yend = MDS2), # adding in the vectors, c
               arrow = arrow(length = unit(0.5, "cm")), colour = "black") + # can add in geom_label or geom_text for labels
  theme_minimal() + # no background
  theme(panel.border = element_rect(fill = NA)) + # full square around figure
  xlab("NMDS 1") +
  ylab("NMDS 2") +
  #ylim(-1, 1.5) +
  #xlim(-1.45, 1) +
  #theme(legend.position = "none") +
  geom_label_repel(data = june.axis12, 
                  aes(x = MDS1, y = MDS2, label = Taxa),
                  color="black",
                  size = 5) +
  scale_fill_viridis(discrete = TRUE) +
  scale_colour_viridis(discrete = TRUE) +
  scale_shape_manual(values = c(24)) +
  theme(legend.position = "none") 

invert.12.june


invert.13.june <- ggplot(data = nmds.june.scores,
                       aes(x = NMDS1, y = NMDS3)) +
  geom_point(data = nmds.june.scores, 
             aes(x = NMDS1, y = NMDS3, 
                 fill = Treatment, shape = Year),
             size = 4, stroke = 1.5,
             alpha = 0.7) + # sites as points
  stat_ellipse(data = nmds.june.scores, 
               aes(x = NMDS1,y = NMDS3,
                   linetype = Treatment, colour = Treatment), 
               size = 1, level = 0.9) + 
  geom_segment(data = june.axis13, 
               aes(x = 0, xend = MDS1, y = 0, yend = MDS3), # adding in the vectors, c
               arrow = arrow(length = unit(0.5, "cm")), colour = "black") + # can add in geom_label or geom_text for labels
  theme_minimal() + # no background
  theme(panel.border = element_rect(fill = NA)) + # full square around figure
  xlab("NMDS 1") +
  ylab("NMDS 3") +
  #ylim(-1, 1.5) +
  #xlim(-1.45, 1) +
  #theme(legend.position = "none") +
  geom_label_repel(data = june.axis13, 
                  aes(x = MDS1, y = MDS3, label = Taxa),
                  color="black",
                  size = 5) +
  scale_fill_viridis(discrete = TRUE) +
  scale_colour_viridis(discrete = TRUE) +
  scale_shape_manual(values = c(24)) +
  theme(legend.position = "none") 

invert.13.june

nmds.june <- ggarrange(invert.12.june, invert.13.june,
                      labels = c("B", ""))


NMDS.june <- annotate_figure(nmds.june,
                top = text_grob("Collection 05-June-18"))

ggsave("Figures/NMDS_emerging_5Jun18.jpeg")


# Panel all together ------------------------------------------------------

NMDS.panel <- ggarrange(NMDS.may, NMDS.june, NMDS.col1, 
                        NMDS.col2, NMDS.col3, NMDS.col4,
                        widths = c(0.7, 1))


ggsave("Figures/NMDS.panel.jpeg", NMDS.panel,
       width = 20,
       height = 9.5,
       units = "in")



### get a legend 
test <- ggplot(data = nmds.col4.scores,
                       aes(x = NMDS1, y = NMDS3)) +
  geom_point(data = nmds.col4.scores, 
             aes(x = NMDS1, y = NMDS3, 
                 colour = Treatment, shape = Year),
             size = 4) + # sites as points
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
  geom_label_repel(data = col4.axis13, 
                   aes(x = MDS1, y = MDS3, label = Taxa),
                   color="black",
                   size = 5) +
  scale_fill_viridis(discrete = TRUE) +
  scale_colour_viridis(discrete = TRUE) +
  scale_shape_manual(values = c(19, 17)) +
  ylim(-0.6, 0.6) +
  xlim(-0.6, 0.6)



legend <- get_legend(test)
legends <- as_ggplot(legend)






# try patchwork -----------------------------------------------------------
library(patchwork)

(NMDS.may | NMDS.june | NMDS.col1) /
  (NMDS.col2 | NMDS.col4)











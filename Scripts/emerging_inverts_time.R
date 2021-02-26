
# Load packages -----------------------------------------------------------

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

# NMDS Col. 1 Data Prep ---------------------------------------------------

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


### Taxa correlated with Axis 1 & 3

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

## NMDS Axis 1, 2 
env.col1$Year <- as.factor(env.col1$Year)

invert.12.c1 <- ggplot(data = env.col1,
                    aes(x = NMDS1, y = NMDS2)) +
  geom_point(data = env.col1, 
             aes(x = NMDS1, y = NMDS2, 
                 colour = Treatment, shape = Year),
             size = 4, stroke = 1.5) + # sites as points
  stat_ellipse(data = env.col1, 
               aes(x = NMDS1,y = NMDS2,
                   linetype = Treatment, colour = Treatment), 
               size = 1, level = 0.9) + 
  geom_segment(data = corr.c1.vectors.12, 
               aes(x = 0, xend = MDS1, y = 0, yend = MDS2), # adding in the vectors, c
               arrow = arrow(length = unit(0.5, "cm")), colour = "black") + # can add in geom_label or geom_text for labels
  theme_minimal() + # no background
  theme(panel.border = element_rect(fill = NA)) + # full square around figure
  xlab("NMDS 1") +
  ylab("NMDS 2") +
  #ylim(-1, 1.5) +
  #xlim(-1.45, 1) +
  #theme(legend.position = "none") +
  geom_text_repel(data = corr.c1.vectors.12, 
                  aes(x = MDS1, y = MDS2, label = Taxa),
                  color="black",
                  size = 5) +
  scale_color_manual(values = c("#969696","#35978f", "#2166ac")) +
  scale_shape_manual(values = c(17, 16, 1, 17, 2, 18, 5)) 

invert.12.c1

## NMDS Axis 1, 3

invert.13.c1 <- ggplot(data = env.col1,
                       aes(x = NMDS1, y = NMDS3)) +
  geom_point(data = env.col1, 
             aes(x = NMDS1, y = NMDS3, 
                 colour = Treatment, shape = Year),
             size = 4, stroke = 1.5) + # sites as points
  stat_ellipse(data = env.col1, 
               aes(x = NMDS1,y = NMDS3,
                   linetype = Treatment, colour = Treatment), 
               size = 1, level = 0.9) + 
  geom_segment(data = corr.c1.vectors.13, 
               aes(x = 0, xend = MDS1, y = 0, yend = MDS3), # adding in the vectors, c
               arrow = arrow(length = unit(0.5, "cm")), colour = "black") + # can add in geom_label or geom_text for labels
  theme_minimal() + # no background
  theme(panel.border = element_rect(fill = NA)) + # full square around figure
  xlab("NMDS 1") +
  ylab("NMDS 3") +
  #ylim(-1, 1.5) +
  #xlim(-1.45, 1) +
  #theme(legend.position = "none") +
  geom_text_repel(data = corr.c1.vectors.13, 
                  aes(x = MDS1, y = MDS3, label = Taxa),
                  color="black",
                  size = 5) +
  scale_color_manual(values = c("#969696","#35978f", "#2166ac")) +
  scale_shape_manual(values = c(17, 16, 1, 17, 2, 18, 5)) 

invert.13.c1

(NMS.emerging.panel.c1 <- ggarrange(invert.12.c1, invert.13.c1,
                                 common.legend = TRUE,
                                 legend = "bottom",
                                 labels = c("A", "")))

annotate_figure(NMS.emerging.panel.c1,
                top = text_grob("Collection 19-Jun-2017 and 16-Jun-2018"))

# NMDS Collection 28-Jun-17 and 25-JUN-18 ---------------------------------

col.2 <- c("28-Jun-17","25-Jun-18")

invert.col.2 <- inverts.date %>% filter(Date %in% col.2)



# NMDS Collection 08-Jul-17 and 04-Jul-18 -----------------------------------

col.3 <- c("08-Jul-17","04-Jul-18")

invert.col.3 <- inverts.date %>% filter(Date %in% col.3)


# NMDS Collection 21-Jul-17 and 23-Jul-18 -----------------------------------

col.4 <- c("21-Jul-17","23-Jul-18")

invert.col.4 <- inverts.date %>% filter(Date %in% col.4)




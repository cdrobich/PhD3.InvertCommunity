
library(tidyverse) # for figures, pipes etc.
library(vegan) # for NMDS
library(ggpubr) # for arranging figures

sol <- read.csv("Data/Inverts_relmax_raresn2.csv")

bugs <- sol[,5:50] # just the Families
env <- sol[,1:4] # Categorical variables

env$TrtYr <- as.factor(env$TrtYr) 
env$Treatment <- as.factor(env$Treatment)

####### NMDS ########

## stress plot/ see how many axes make sense

k_vec <- 1:10 #dimensions 1 - 10
stress <- numeric(length(k_vec)) # stress of each model put here
dune_dij <- metaMDSdist(bugs)
set.seed(25)
for(i in seq_along(k_vec)) {
  sol <- metaMDSiter(dune_dij, k = i, # you can tell I lifted this from a tutorial on the dune package
                     trace = FALSE)
  stress[i] <- sol$stress
}

plot(stress) # 3D makes sense

#### NMDS analysis ####
set.seed(120) # set this so you can re-run it since its iterative

nms.invert <- metaMDS(bugs, distance = "bray", # species data, bray-curtis dissimilarity
                       autotransform = FALSE,  # NMDS will do autotransformations for you
                       k = 3, trymax = 1000)   # k = number of axes
nms.invert

# Output
#metaMDS(comm = bugs, distance = "bray", k = 3, trymax = 1000, autotransform = FALSE) 

#global Multidimensional Scaling using monoMDS

#Data:     bugs 
#Distance: bray 

#Dimensions: 3 
#Stress:     0.2042451 
#Stress type 1, weak ties
#Two convergent solutions found after 27 tries
#Scaling: centring, PC rotation, halfchange scaling 
#Species: expanded scores based on ‘bugs’


# look at points and the stress real quick
layout(matrix(1:2, ncol = 2))
plot(nms.invert, main = "Invertebrate NMDS plot"); stressplot(nms.invert, main = "Shepard plot")
layout(1)

# how many iterations of the NMDS
nms.invert$iters # 85

# Goodness of fit
(g <- goodness(nms.invert)) # smaller the number the better the fit
sum(g^2)
nms.invert$stress^2  # 0.0415


## extract the scores for plotting 
scr <- as.data.frame(scores(nms.invert, display = "sites")) # extract NMDS scores

# adding categorical info to scores
scr$site <- rownames(scr)
scr$Treatment <- env$Treatment
scr$Year <- env$Year
scr$TrtYr <- env$TrtYr

colnames(scr)
col_order <- c("site", "Treatment", "Year", "TrtYr", "NMDS1", "NMDS2", "NMDS3")

scores <- scr[, col_order] # put the categorical values in order
scores$Year <- as.factor(scores$Year)

write.csv(scores,"Data/NMDS_inverts_scores.csv") # save this as a csv

### Vectors correlated with Axis 1 & 2 ####

alltaxa <- envfit(nms.invert, bugs, 
                  permutations = 999, choices = c(1,2)) #produces a list with r2, p value, and NMDS coordinates

all.taxa.df <- data.frame((alltaxa$vectors)$arrows, (alltaxa$vectors)$r, (alltaxa$vectors)$pvals) #take list and make into dataframe
write.csv(all.taxa.df, "Data/NMDS_vectors_axis12.csv") # save vector scores as csv


alltaxa$vectors$r[alltaxa$vectors$r > 0.2] # selecting vectors (Family) that are reasonably correlated (r2 > 0.2)

#Anthicidae   Erotylidae Cecidomyiida    Sciaridae Chironomidae 
#0.2527976    0.2553976    0.2508132    0.2029796    0.3025265 

#Dolichopodid     Carnidae  Chloropidae Anthomyiidae Calliphorida 
#0.2821084    0.2510678    0.2039437    0.2640309    0.2126135

# taking out those correlated ones to add to the figure
corr.taxa <- bugs %>% select(Anthicidae,
                             Erotylidae,
                             Cecidomyiida,
                             Sciaridae,
                             Chironomidae,
                             Dolichopodid,
                             Carnidae,
                             Chloropidae,
                             Anthomyiidae,
                             Calliphorida)

# recalculated it because I am lazy but you could probably pull them from all.taxa.df
corrtaxa <- envfit(nms.invert$points, corr.taxa, 
                  permutations = 999, choices = c(1,2))


corrtaxa

# make a new data frame for the figure
species.12 <- as.data.frame(corrtaxa$vectors$arrows*sqrt(corrtaxa$vectors$r)) #scaling vectors so they correspond with r2
species.12$species <- rownames(species.12)


#### Vectors correlated with axis 1 & 3 ###
# same as above but now for the other axes/figure

alltaxa.13 <- envfit(nms.invert, bugs, 
                  permutations = 999, choices = c(1,3)) 


all.taxa13.df <- data.frame((alltaxa.13$vectors)$arrows, (alltaxa.13$vectors)$r, (alltaxa.13$vectors)$pvals)
write.csv(all.taxa13.df, "Data/NMDS_vectors_axis13.csv")


alltaxa.13$vectors$r[alltaxa.13$vectors$r > 0.2] 

#Carabidae   Erotylidae   Lampyridae Cecidomyiida Chironomidae Dolichopodid 
#0.3472899    0.2545980    0.2309457    0.2695790    0.4049056    0.2324400 

#Carnidae 
#0.2484244 

corr.taxa.13 <- bugs %>% select(Carabidae,
                             Erotylidae,
                             Cecidomyiida,
                             Lampyridae,
                             Chironomidae,
                             Dolichopodid,
                             Carnidae)



corrtaxa.13 <- envfit(nms.invert$points, corr.taxa.13, 
                   permutations = 999, choices = c(1,3))


corrtaxa.13

species.13 <- as.data.frame(corrtaxa.13$vectors$arrows*sqrt(corrtaxa.13$vectors$r)) #scaling vectors
species.13$species <- rownames(species.13)



## plotting with base R ##
# I don't like base R but it has some nice things built in for ordinations

# axis 1, 2 

col_vecs <- c("#762a83", "#9970ab", "#5aae61", "#1b7837","#4393c3","#2166ac") 
col_vec <- c("#9970ab", "#1b7837", "#2166ac")

ordiplot(nms.invert, choices = c(1,2), 
         type = "points",
         display = "sites")


ordihull(nms.invert, groups = env$Treatment, # ellipse hull
            col = col_vec)

ordiellipse(nms.invert, groups = env$Treatment, # st error of centroid ellipse hull
            draw = "polygon",
            col = col_vec)

lvl <- with(env, levels(Treatment))

legend("topright", legend = lvl,
       bty = "n", col = col_vec, # create the legend by hand and map colours to levels
       pch = 19)

plot(alltaxa, p.max = 0.05, col = "black")


#### ggPlot Figures ####

## NMDS Axis 1, 2 

scores # coordinates we extracted
species.12 # reasonably correlated vectors with axis 1,2

invert.12 <- ggplot(data = scores,
                    aes(x = NMDS1, y = NMDS2)) +
  geom_point(data = scores, aes(x = NMDS1, y = NMDS2, colour = Treatment, shape = Treatment), size = 4) + # sites as points
  stat_ellipse(data = scores, aes(x = NMDS1,y = NMDS2,linetype = Treatment, colour = Treatment), size = 1) + # a 95% CI ellipses
  geom_segment(data = species.12, aes(x = 0, xend = MDS1, y = 0, yend = MDS2), # adding in the vectors, c
               arrow = arrow(length = unit(0.5, "cm")), colour = "black") + # can add in geom_label or geom_text for labels
  xlim (-1, 1) + # setting the limits so they're symmetrical
  ylim (-1, 1) +
  theme_classic() + # no background
  scale_color_manual(values = c("#9970ab", "#1b7837", "#2166ac")) + # adding colours
  theme(panel.border = element_rect(fill = NA)) + # full square around figure
  xlab("NMDS 1") +
  ylab("NMDS 2") 

invert.12

## NMDS Axis 1, 3
# same as above
  
scores
species.13 # vectors correlated with axis 1, 3

invert.13 <- ggplot(data = scores,
                    aes(x = NMDS1, y = NMDS3)) +
  geom_point(data = scores, aes(x = NMDS1, y = NMDS3, colour = Treatment, shape = Treatment), size = 4) +
  stat_ellipse(data = scores, aes(x = NMDS1,y = NMDS3,linetype = Treatment, colour = Treatment), size = 1) +
  geom_segment(data = species.13, aes(x = 0, xend = MDS1, y = 0, yend = MDS3),
               arrow = arrow(length = unit(0.5, "cm")), colour = "black") +
  xlim (-1.1, 1.1) +
  ylim (-1.1, 1.1) +
  theme_classic() +
  scale_color_manual(values = c("#9970ab", "#1b7837", "#2166ac")) +
  theme(panel.border = element_rect(fill = NA)) +
  xlab("NMDS 1") +
  ylab("NMDS 3") 

invert.13


## this is from ggpubr
# arranging the two figures
NMDS.inv <- ggarrange(invert.13, invert.12, # put the plot items in
          nrow = 2, # I want them on top of each other
          common.legend = TRUE, # they have the same legend
          legend = "bottom") # I want it on the bottom 
NMDS.inv

ggsave("Figures/NMDS_invertebrate.jpeg", NMDS.inv) # save that figure to my folder

library(tidyverse)
library(lubridate)
library(vegan)
library(viridis)



# Emerging inverts 5-June-18 to 23-July-18
emerg.raw <- read.csv("Data/Emerging/Procrustes/emerging_invert_junjuly.csv")

colnames(emerg.raw)

# Sum over the visits 

df2 <- emerg.raw %>% 
  mutate(Date = dmy(Date)) %>% 
  group_by(ID, Treatment) %>%
  summarise_at(vars(Araneae:Hesperiidae), sum, na.rm = TRUE)


# output to remove OW and select taxa
write.csv(df2, "Data/Emerging/Procrustes/emerging_sum_junjuly.csv")


# Rarefied data -----------------------------------------------------------

emerging <- read.csv("Data/Emerging/Procrustes/ermerging_procrustes.csv")

aquatic <- read.csv("Data/Emerging/Procrustes/aquatic_procrustes.csv")

# Aquatic invertebrate data prep ------------------------------------------

# Split into three matrices

colnames(aquatic)

# filter the habitat type
aquatic.inv <- aquatic %>% filter(Habitat == "Invaded")
aquatic.unin <- aquatic %>% filter(Habitat == "Uninvaded")
aquatic.trt <- aquatic %>% filter(Habitat == "Treated")

# Remove empty columns
aquatic.inv <- aquatic.inv[, colSums(aquatic.inv != 0) > 0]
aquatic.unin <- aquatic.unin[, colSums(aquatic.unin  != 0) > 0]
aquatic.trt <- aquatic.trt[, colSums(aquatic.trt != 0) > 0]

# Just taxa
aquatic.inv.taxa <- aquatic.inv %>% select(Chironomidae:Leptoceridae)
aquatic.unin.taxa <- aquatic.unin %>% select(Lampyridae:Leptoceridae)
aquatic.trt.taxa <- aquatic.trt %>% select(Chironomidae:Leptoceridae)


#make proportions
colnames(aquatic.inv)

aquatic.inv.prop <- aquatic.inv.taxa/rowSums(aquatic.inv.taxa)
aquatic.unin.prop <- aquatic.unin.taxa/rowSums(aquatic.unin.taxa)
aquatic.trt.prop <- aquatic.trt.taxa/rowSums(aquatic.trt.taxa)

# uh puts them in the wrong spot
myList <- list(aquatic.inv.prop = aquatic.inv.prop, 
               aquatic.unin.prop = aquatic.unin.prop,
               aquatic.trt.prop = aquatic.trt.prop)
for(i in names(myList)){
  write.csv(myList[[i]], paste0(i,".csv"))
}


# Emerging invertebrate data prep -----------------------------------------

colnames(emerging)

# filter the habitat type
emerg.inv <- emerging %>% filter(Habitat == "Invaded")
emerg.unin <- emerging %>% filter(Habitat == "Uninvaded")
emerg.trt <- emerging %>% filter(Habitat == "Treated")

# Remove empty columns
emerg.inv <- emerg.inv[, colSums(emerg.inv != 0) > 0]
emerg.unin <- emerg.unin[, colSums(emerg.unin  != 0) > 0]
emerg.trt <- emerg.trt[, colSums(emerg.trt != 0) > 0]

# Just taxa
emerg.inv.taxa <- emerg.inv %>% select(Lampyridae:Crambidae)
emerg.unin.taxa <- emerg.unin %>% select(Lampyridae:Crambidae)
emerg.trt.taxa <- emerg.trt %>% select(Phoridae:Leptoceridae)


#make proportions - with Hellinger this part isnt required
colnames(emerg.inv)

emerg.inv.prop <- emerg.inv.taxa/rowSums(emerg.inv.taxa)
emerg.unin.prop <- emerg.unin.taxa/rowSums(emerg.unin.taxa)
emerg.trt.prop <- emerg.trt.taxa/rowSums(emerg.trt.taxa)

# or just decostand(emerg.inv.taxa, "total")

# uh puts them in the wrong spot
myList2 <- list(emerg.inv.prop = emerg.inv.prop, 
               emerg.unin.prop = emerg.unin.prop,
               emerg.trt.prop = emerg.trt.prop)
for(i in names(myList2)){
  write.csv(myList2[[i]], paste0(i,".csv"))
}

# Transform the matrices

emerg.inv.h <- decostand(emerg.inv.taxa, "hellinger")
emerg.unin.h <- decostand(emerg.unin.taxa, "hellinger")
emerg.trt.h <- decostand(emerg.trt.taxa, "hellinger")

aq.inv.h <- decostand(aquatic.inv.prop, "hellinger")
aq.unin.h <- decostand(aquatic.unin.prop, "hellinger")
aq.trt.h <- decostand(aquatic.trt.prop, "hellinger")

# Procrustes --------------------------------------------------------------

## Invaded habitat
inv.pro <- protest(emerg.inv.h, aq.inv.h)
summary(inv.pro)

inv.pro$t0 # 0.8937901 Procrustes correlation from non-permuted solution
inv.pro$signif # 0.006 #Significant of permutations (999)

plot(inv.pro)

# just proportions
inv.pro1 <- protest(emerg.inv.prop, aquatic.inv.prop)
summary(inv.pro1)
inv.pro1$t0 # 0.7285
inv.pro1$signif #0.033


# Uninvaded habitat

unin.pro <- protest(emerg.unin.h, aq.unin.h)

summary(unin.pro)

unin.pro$t0 # 0.8648073
unin.pro$signif # 0.049

# just proportions
unin.pro1 <- protest(emerg.unin.prop, aquatic.unin.prop)
summary(unin.pro1)
unin.pro1$t0 # 0.5978492
unin.pro1$signif #0.235


# Treated habitat

trt.pro <- protest(emerg.trt.h, aq.trt.h)
summary(trt.pro)

trt.pro$t0 #0.4854401
trt.pro$signif # 0.647

trt.pro1 <- protest(emerg.trt.prop, aquatic.trt.prop)
summary(trt.pro1)
trt.pro1$t0 # 0.3179
trt.pro1$signif #0.605


# NMDS -----------------------------------------------


# Invaded NMDS ------------------------------------------------------------
invaded.both <- full_join(emerg.inv, aquatic.inv)

invaded.both[is.na(invaded.both)] <- 0

write.csv(invaded.both, "Data/Emerging/Procrustes/invaded_matrix.csv")

invaded.taxa <- invaded.both %>% select(Lampyridae:Crambidae)
invaded.ev <- invaded.both %>% select(ID:Type)

invaded.taxa <- decostand(invaded.taxa, "hellinger")



k_vec <- 1:10 #dimensions 1 - 10
stress <- numeric(length(k_vec)) # stress of each model put here
dune_dij <- metaMDSdist(invaded.taxa)
set.seed(25)
for(i in seq_along(k_vec)) {
  sol <- metaMDSiter(dune_dij, k = i, # you can tell I lifted this from a tutorial on the dune package
                     trace = FALSE)
  stress[i] <- sol$stress
}

plot(stress) # Big drop after 1, lets do 2

inv.nms <- metaMDS(invaded.taxa, distance = "bray", 
                    autotransform = FALSE,  # NMDS will do autotransformations for you
                    k = 2, trymax = 1000)   # k = number of axes
inv.nms

#Data:     invaded.taxa 
#Distance: bray 

#Dimensions: 2 
#Stress:     0.1379622 
#Stress type 1, weak ties
#Two convergent solutions found after 20 tries

layout(matrix(1:2, ncol = 2))
plot(inv.nms, main = "Invertebrate NMDS plot"); stressplot(inv.nms, main = "Shepard plot")
layout(1)

inv.nms$iters #48

1-inv.nms$stress^2 #0.981

# Scores for points
scores.inv <- as.data.frame(scores(inv.nms, display = "sites"))

invaded.ev$NMDS1 <- scores.inv$NMDS1
invaded.ev$NMDS2 <- scores.inv$NMDS2

write.csv(invaded.ev, "Data/Emerging/Procrustes/invaded_NMDS_scores.csv")

# Scores for segments

emerg.inv.sc <- invaded.ev %>% filter(Type == "Emerging")

aqua.inv.sc <- invaded.ev %>% filter(Type == "Aquatic")


emerg.inv.sc <- rename(emerg.inv.sc, 
                        NMDS1e = NMDS1,
                        NMDS2e = NMDS2)

aqua.inv.sc  <- rename(aqua.inv.sc, 
                        NMDS1a = NMDS1,
                        NMDS2a = NMDS2)


emerging.inv.sc <- emerg.inv.sc %>% select(ID, Habitat, NMDS1e, NMDS2e)
aquatic.inv.sc <- aqua.inv.sc %>% select(ID, Habitat, NMDS1a, NMDS2a)

pro.inv.scores <- full_join(aquatic.inv.sc, emerging.inv.sc)

write.csv(pro.inv.scores, "Data/Emerging/Procrustes/NMDS_inv_segment_score.csv")

## Correalted vectors

taxa.axis12.inv <- envfit(inv.nms, invaded.taxa,
                           choices = c(1,2))

taxa.axis12df.inv <- data.frame((taxa.axis12.inv$vectors)$arrows,
                                 (taxa.axis12.inv$vectors)$r,
                                 (taxa.axis12.inv$vectors)$pvals)

taxa.axis12df.inv <- tibble::rownames_to_column(taxa.axis12df.inv, "Taxa")


colnames(taxa.axis12df.inv)

corrspp.axis12.inv <- taxa.axis12df.inv %>% filter(X.taxa.axis12.inv.vectors..r > 0.3)
target12.taxa.inv <- corrspp.axis12.inv$Taxa # string of the Family names


axis12.vectors.inv <- invaded.taxa %>% select(all_of(target12.taxa.inv)) # make a matrix of just those

(nmds.inv.vectors.12 <- envfit(inv.nms$points, axis12.vectors.inv,
                                permutations = 999, choices = c(1,2)))                        


corr.inv.vectors.12 <- as.data.frame(nmds.inv.vectors.12$vectors$arrows*sqrt(nmds.inv.vectors.12$vectors$r)) #scaling vectors
corr.inv.vectors.12$Taxa <- rownames(corr.inv.vectors.12)

write.csv(corr.inv.vectors.12, "Data/Emerging/Procrustes/NMDS_correlated_vector12_invaded.csv")

# Uninvaded NMDS ----------------------------------------------------------

uninvaded.both <- full_join(emerg.unin, aquatic.unin)

uninvaded.both[is.na(uninvaded.both)] <- 0

write.csv(uninvaded.both, "Data/Emerging/Procrustes/uninvaded_matrix.csv")

uninvaded.taxa <- uninvaded.both %>% select(Lampyridae:Corduliidae)
uninvaded.ev <- uninvaded.both %>% select(ID:Type)

unin.taxa <- decostand(uninvaded.taxa, "hellinger")



k_vec <- 1:10 #dimensions 1 - 10
stress <- numeric(length(k_vec)) # stress of each model put here
dune_dij <- metaMDSdist(unin.taxa )
set.seed(25)
for(i in seq_along(k_vec)) {
  sol <- metaMDSiter(dune_dij, k = i, # you can tell I lifted this from a tutorial on the dune package
                     trace = FALSE)
  stress[i] <- sol$stress
}

plot(stress) # Big drop after 1, lets do 2

unin.nms <- metaMDS(unin.taxa , distance = "bray", 
                   autotransform = FALSE,  # NMDS will do autotransformations for you
                   k = 2, trymax = 1000)   # k = number of axes
unin.nms

#Data:     unin.taxa 
#Distance: bray 

#Dimensions: 2 
#Stress:     0.1609472 
#Stress type 1, weak ties
#Two convergent solutions found after 20 tries

layout(matrix(1:2, ncol = 2))
plot(unin.nms, main = "Invertebrate NMDS plot"); stressplot(unin.nms, main = "Shepard plot")
layout(1)

unin.nms$iters #88

1-unin.nms$stress^2 #0.974

# Scores for points
scores.unin <- as.data.frame(scores(unin.nms, display = "sites"))

uninvaded.ev$NMDS1 <- scores.unin$NMDS1
uninvaded.ev$NMDS2 <- scores.unin$NMDS2

write.csv(uninvaded.ev, "Data/Emerging/Procrustes/uninvaded_NMDS_scores.csv")

# Scores for segments

emerg.unin.sc <- uninvaded.ev %>% filter(Type == "Emerging")

aqua.unin.sc <- uninvaded.ev %>% filter(Type == "Aquatic")


emerg.unin.sc <- rename(emerg.unin.sc, 
                       NMDS1e = NMDS1,
                       NMDS2e = NMDS2)

aqua.unin.sc  <- rename(aqua.unin.sc, 
                       NMDS1a = NMDS1,
                       NMDS2a = NMDS2)


emerging.unin.sc <- emerg.unin.sc %>% select(ID, Habitat, NMDS1e, NMDS2e)
aquatic.unin.sc <- aqua.unin.sc %>% select(ID, Habitat, NMDS1a, NMDS2a)

pro.unin.scores <- full_join(aquatic.unin.sc, emerging.unin.sc)

write.csv(pro.unin.scores, "Data/Emerging/Procrustes/NMDS_unin_segment_score.csv")

## Correalted vectors

taxa.axis12.unin <- envfit(unin.nms, unin.taxa,
                          choices = c(1,2))

taxa.axis12df.unin <- data.frame((taxa.axis12.unin$vectors)$arrows,
                                (taxa.axis12.unin$vectors)$r,
                                (taxa.axis12.unin$vectors)$pvals)

taxa.axis12df.unin <- tibble::rownames_to_column(taxa.axis12df.unin, "Taxa")


colnames(taxa.axis12df.unin)

corrspp.axis12.unin <- taxa.axis12df.unin %>% filter(X.taxa.axis12.unin.vectors..r > 0.3)
target12.taxa.unin <- corrspp.axis12.unin$Taxa # string of the Family names


axis12.vectors.unin <- unin.taxa %>% select(all_of(target12.taxa.unin)) # make a matrix of just those

(nmds.unin.vectors.12 <- envfit(unin.nms$points, axis12.vectors.unin,
                               permutations = 999, choices = c(1,2)))                        


corr.unin.vectors.12 <- as.data.frame(nmds.unin.vectors.12$vectors$arrows*sqrt(nmds.unin.vectors.12$vectors$r)) #scaling vectors
corr.unin.vectors.12$Taxa <- rownames(corr.unin.vectors.12)

write.csv(corr.unin.vectors.12, "Data/Emerging/Procrustes/NMDS_correlated_vector12_uninvaded.csv")


# Treated NMDS -------------------------------------------------------------

treated.both <- full_join(emerg.trt, aquatic.trt)

treated.both[is.na(treated.both)] <- 0

write.csv(treated.both, "Data/Emerging/Procrustes/treated_matrix.csv")

treated.taxa <- treated.both %>% select(Phoridae:Corduliidae)
treated.ev <- treated.both %>% select(ID:Type)

trt.taxa <- decostand(treated.taxa, "hellinger")



k_vec <- 1:10 #dimensions 1 - 10
stress <- numeric(length(k_vec)) # stress of each model put here
dune_dij <- metaMDSdist(trt.taxa)
set.seed(25)
for(i in seq_along(k_vec)) {
  sol <- metaMDSiter(dune_dij, k = i, # you can tell I lifted this from a tutorial on the dune package
                     trace = FALSE)
  stress[i] <- sol$stress
}

plot(stress) # Big drop after 1, lets do 2

trt.nms <- metaMDS(trt.taxa, distance = "bray", 
                     autotransform = FALSE,  # NMDS will do autotransformations for you
                     k = 2, trymax = 1000)   # k = number of axes
trt.nms

#Data:     trt.taxa 
#Distance: bray 

#Dimensions: 2 
#Stress:     0.07999074 
#Stress type 1, weak ties
#Two convergent solutions found after 20 tries

layout(matrix(1:2, ncol = 2))
plot(trt.nms, main = "Invertebrate NMDS plot"); stressplot(trt.nms, main = "Shepard plot")
layout(1)

trt.nms$iters #49

1-trt.nms$stress^2 #0.993

# Scores for points
scores.trt <- as.data.frame(scores(trt.nms, display = "sites"))

treated.ev$NMDS1 <- scores.trt$NMDS1
treated.ev$NMDS2 <- scores.trt$NMDS2

write.csv(treated.ev, "Data/Emerging/Procrustes/treated_NMDS_scores.csv")

# Scores for segments

emerg.trt.sc <- treated.ev %>% filter(Type == "Emerging")

aqua.trt.sc <- treated.ev %>% filter(Type == "Aquatic")


emerg.trt.sc <- rename(emerg.trt.sc, 
                          NMDS1e = NMDS1,
                          NMDS2e = NMDS2)

aqua.trt.sc  <- rename(aqua.trt.sc, 
                          NMDS1a = NMDS1,
                          NMDS2a = NMDS2)


emerging.trt.sc <- emerg.trt.sc %>% select(ID, Habitat, NMDS1e, NMDS2e)
aquatic.trt.sc <- aqua.trt.sc %>% select(ID, Habitat, NMDS1a, NMDS2a)

pro.trt.scores <- full_join(aquatic.trt.sc, emerging.trt.sc)

write.csv(pro.trt.scores, "Data/Emerging/Procrustes/NMDS_trt_segment_score.csv")

## Correalted vectors

taxa.axis12.trt <- envfit(trt.nms, trt.taxa,
                      choices = c(1,2))

taxa.axis12df.trt <- data.frame((taxa.axis12.trt$vectors)$arrows,
                            (taxa.axis12.trt$vectors)$r,
                            (taxa.axis12.trt$vectors)$pvals)

taxa.axis12df.trt <- tibble::rownames_to_column(taxa.axis12df.trt, "Taxa")


colnames(taxa.axis12df.trt)

corrspp.axis12.trt <- taxa.axis12df.trt %>% filter(X.taxa.axis12.trt.vectors..r > 0.3)
target12.taxa.trt <- corrspp.axis12.trt$Taxa # string of the Family names


axis12.vectors.trt <- trt.taxa%>% select(all_of(target12.taxa.trt )) # make a matrix of just those

(nmds.trt.vectors.12 <- envfit(trt.nms$points, axis12.vectors.trt,
                              permutations = 999, choices = c(1,2)))                        


corr.trt.vectors.12 <- as.data.frame(nmds.trt.vectors.12$vectors$arrows*sqrt(nmds.trt.vectors.12$vectors$r)) #scaling vectors
corr.trt.vectors.12$Taxa <- rownames(corr.trt.vectors.12)

write.csv(corr.trt.vectors.12, "Data/Emerging/Procrustes/NMDS_correlated_vector12_treated.csv")


# NMDS all data together --------------------------------------------------


# put the data together

invert.taxa <- full_join(emerging, aquatic)

invert.taxa[is.na(invert.taxa)] <- 0

write.csv(invert.taxa, "Data/Emerging/Procrustes/invertebrates_combined.csv")

colnames(invert.taxa)

taxa.in <- invert.taxa %>% select(Lampyridae:Limnephilidae)
env.in <- invert.taxa  %>% select(ID:Type)

## NMDS 

invert.h <- decostand(taxa.in, "hellinger")


k_vec <- 1:10 #dimensions 1 - 10
stress <- numeric(length(k_vec)) # stress of each model put here
dune_dij <- metaMDSdist(invert.h)
set.seed(25)
for(i in seq_along(k_vec)) {
  sol <- metaMDSiter(dune_dij, k = i, # you can tell I lifted this from a tutorial on the dune package
                     trace = FALSE)
  stress[i] <- sol$stress
}

plot(stress) # Big drop after 1, lets do 2

emerg.nms <- metaMDS(invert.h, distance = "bray", 
                    autotransform = FALSE,  # NMDS will do autotransformations for you
                    k = 2, trymax = 1000)   # k = number of axes
emerg.nms

#Data:     invert.h 
#Distance: bray 

#Dimensions: 2 
#Stress:     0.1439528 
#Stress type 1, weak ties
#Two convergent solutions found after 50 tries

layout(matrix(1:2, ncol = 2))
plot(emerg.nms, main = "Invertebrate NMDS plot"); stressplot(emerg.nms, main = "Shepard plot")
layout(1)

# Scores for points
scores <- as.data.frame(scores(emerg.nms, display = "sites"))

env.in$NMDS1 <- scores$NMDS1
env.in$NMDS2 <- scores$NMDS2

write.csv(env.in, "Data/Emerging/Procrustes/both_NMDS.csv")


# Scores for segments

emerg.invert.sc <- env.in %>% filter(Type == "Emerging")

aqua.invert.sc <- env.in %>% filter(Type == "Aquatic")

emerg.invert.sc <- rename(emerg.invert.sc, 
       NMDS1e = NMDS1,
       NMDS2e = NMDS2)

aqua.invert.sc  <- rename(aqua.invert.sc , 
                          NMDS1a = NMDS1,
                          NMDS2a = NMDS2)



emerging.sc <- emerg.invert.sc %>% select(ID, Habitat, NMDS1e, NMDS2e)
aquatic.sc <- aqua.invert.sc %>% select(ID, Habitat, NMDS1a, NMDS2a)

pro.scores <- full_join(aquatic.sc, emerging.sc)

write.csv(pro.scores, "Data/Emerging/Procrustes/NMDS_scores_segment.csv")


## Correlated vectors

taxa.axis12 <- envfit(emerg.nms, taxa.in,
                           choices = c(1,2))

taxa.axis12df <- data.frame((taxa.axis12$vectors)$arrows,
                                 (taxa.axis12$vectors)$r,
                                 (taxa.axis12$vectors)$pvals)

taxa.axis12df <- tibble::rownames_to_column(taxa.axis12df, "Taxa")


colnames(taxa.axis12df)

corrspp.axis12 <- taxa.axis12df %>% filter(X.taxa.axis12.vectors..r > 0.3)
target12.taxa <- corrspp.axis12$Taxa # string of the Family names


axis12.vectors.c3 <- taxa.in %>% select(all_of(target12.taxa)) # make a matrix of just those

(nmds.c3.vectors.12 <- envfit(emerg.nms$points, axis12.vectors.c3,
                              permutations = 999, choices = c(1,2)))                        


corr.c3.vectors.12 <- as.data.frame(nmds.c3.vectors.12$vectors$arrows*sqrt(nmds.c3.vectors.12$vectors$r)) #scaling vectors
corr.c3.vectors.12$Taxa <- rownames(corr.c3.vectors.12)

write.csv(corr.c3.vectors.12, "Data/Emerging/Procrustes/NMDS_correlated_vector12.csv")

# Figures -----------------------------------------------------------------


# Invaded NMDS figure -----------------------------------------------------

inv.scores <- read.csv("Data/Emerging/Procrustes/invaded_NMDS_scores.csv")
inv.segment <- read.csv("Data/Emerging/Procrustes/NMDS_inv_segment_score.csv")
inv.vector <- read.csv("Data/Emerging/Procrustes/NMDS_correlated_vector12_invaded.csv")


inv.nms <- ggplot(data = inv.scores,
                   aes(x = NMDS1, y = NMDS2)) +
  geom_point(data = inv.scores, 
             aes(x = NMDS1, 
                 y = NMDS2,
                 shape = Type,
                 stroke = 1.5),
             size = 7,
             fill = "#440C53",
             alpha = 0.7) +
  geom_segment(data = inv.segment, 
               aes(x = NMDS1a, xend = NMDS1e, y = NMDS2a, yend = NMDS2e), 
               arrow = arrow(length = unit(0.5, "cm")),
               size = 1, linetype = 2,
               colour = "#440C53") + 
  geom_segment(data = inv.vector, 
               aes(x = 0, xend = MDS1, y = 0, yend = MDS2), 
               arrow = arrow(length = unit(0.5, "cm")), colour = "black",
               size = 1) +
  geom_label_repel(data = inv.vector, 
                   aes(x = MDS1, y = MDS2, label = Taxa),
                   color = "black",
                   size = 5,
                   force = 2) +
  theme_minimal() + # no background
  theme(panel.border = element_rect(fill = NA)) + 
  xlab("NMDS 1") +
  ylab("NMDS 2") +
  scale_shape_manual(name = " ",
                     labels = c("Aquatic",
                                "Emerging"),
                     values = c(21, 24, 22)) +
  ggtitle("Invaded sites") +
  ylim(-1, 1) +
  xlim(-1, 1) +
  theme(legend.position = "bottom")

inv.nms 


# Uninvaded NMDS figure ---------------------------------------------------

unin.scores <- read.csv("Data/Emerging/Procrustes/uninvaded_NMDS_scores.csv")
unin.segment <- read.csv("Data/Emerging/Procrustes/NMDS_unin_segment_score.csv")
unin.vector <- read.csv("Data/Emerging/Procrustes/NMDS_correlated_vector12_uninvaded.csv")

unin.nms <- ggplot(data = unin.scores,
                  aes(x = NMDS1, y = NMDS2)) +
  geom_point(data = unin.scores, 
             aes(x = NMDS1, 
                 y = NMDS2,
                 shape = Type,
                 stroke = 1.5),
             size = 7,
             fill = "#FDE825",
             alpha = 0.7) +
  geom_segment(data = unin.segment, 
               aes(x = NMDS1a, xend = NMDS1e, y = NMDS2a, yend = NMDS2e), 
               arrow = arrow(length = unit(0.5, "cm")),
               size = 1, linetype = 2,
               colour = "#FDE825") + 
  geom_segment(data = unin.vector, 
               aes(x = 0, xend = MDS1, y = 0, yend = MDS2), 
               arrow = arrow(length = unit(0.5, "cm")), colour = "black",
               size = 1) +
  geom_label_repel(data = unin.vector, 
                   aes(x = MDS1, y = MDS2, label = Taxa),
                   color = "black",
                   size = 5, force = 2) +
  theme_minimal() + # no background
  theme(panel.border = element_rect(fill = NA)) + 
  xlab("NMDS 1") +
  ylab("NMDS 2") +
  scale_shape_manual(name = " ",
                     labels = c("Aquatic",
                                "Emerging"),
                     values = c(21, 24, 22)) +
  ggtitle("Remnant marsh sites") +
  ylim(-1, 1) +
  xlim(-1, 1) +
  theme(legend.position = "bottom")

unin.nms 



# Treated NMDS figure -----------------------------------------------------

trt.scores <- read.csv("Data/Emerging/Procrustes/treated_NMDS_scores.csv")
trt.segment <- read.csv("Data/Emerging/Procrustes/NMDS_trt_segment_score.csv")
trt.vector <- read.csv("Data/Emerging/Procrustes/NMDS_correlated_vector12_treated.csv")


trt.nms <- ggplot(data = trt.scores,
                     aes(x = NMDS1, y = NMDS2)) +
  geom_point(data = trt.scores, 
             aes(x = NMDS1, 
                 y = NMDS2,
                 shape = Type,
                 stroke = 1.5),
             size = 7,
             fill = "#24908C",
             alpha = 0.7) +
  geom_segment(data = trt.segment, 
               aes(x = NMDS1a, xend = NMDS1e, y = NMDS2a, yend = NMDS2e), 
               arrow = arrow(length = unit(0.5, "cm")),
               size = 1, linetype = 2,
               colour = "#24908C") + 
  geom_segment(data = trt.vector, 
               aes(x = 0, xend = MDS1, y = 0, yend = MDS2), 
               arrow = arrow(length = unit(0.5, "cm")), colour = "black",
               size = 1) +
  geom_label_repel(data = trt.vector, 
                   aes(x = MDS1, y = MDS2, label = Taxa),
                   color = "black",
                   size = 5,
                   force = 2) +
  theme_minimal() + # no background
  theme(panel.border = element_rect(fill = NA)) + 
  xlab("NMDS 1") +
  ylab("NMDS 2") +
  scale_shape_manual(name = " ",
                     labels = c("Aquatic",
                                "Emerging"),
                     values = c(21, 24, 22)) +
  ggtitle("Herbicide-treated sites") +
  ylim(-1, 1) +
  xlim(-1, 1) +
  theme(legend.position = "bottom")

trt.nms 


# Patchwork panel ---------------------------------------------------------

library(patchwork)

unin.nms + inv.nms + trt.nms

nmds.panle <- ggarrange(unin.nms, inv.nms, trt.nms,
          ncol = 3,
          labels = "AUTO")

ggsave("Figures/NMDS_emerg_aquatic_panel.jpeg")


# All NMDS together -------------------------------------------------------

colours = c("Invaded" = "#440C53",
            "Treated" = "#24908C",
            "Uninvaded" = "#FDE825")


invert.sc <- read.csv("Data/Emerging/Procrustes/both_NMDS.csv")
invert.seg <- read.csv("Data/Emerging/Procrustes/NMDS_scores_segment.csv")
vectors <- read.csv("Data/Emerging/Procrustes/NMDS_correlated_vector12.csv")

colnames(invert.seg)

pro.figure <- ggplot(data = invert.sc,
                     aes(x = NMDS1, y = NMDS2)) +
  geom_point(data = invert.sc, 
             aes(x = NMDS1, 
                 y = NMDS2, 
                 fill = Habitat, 
                 shape = Type,
                 stroke = 1.5),
             size = 5) +
  geom_segment(data = invert.seg, 
               aes(x = NMDS1a, xend = NMDS1e, y = NMDS2a, yend = NMDS2e,
                   colour = Habitat), 
               arrow = arrow(length = unit(0.5, "cm")),
               size = 1, linetype = 2) + 
  geom_segment(data = vectors, 
               aes(x = 0, xend = MDS1, y = 0, yend = MDS2), 
               arrow = arrow(length = unit(0.5, "cm")), colour = "black",
               size = 0.5) +
  geom_label_repel(data = vectors, 
                   aes(x = MDS1, y = MDS2, label = Taxa),
                   color = "black",
                   size = 5) +
  theme_minimal() + # no background
  theme(panel.border = element_rect(fill = NA)) + 
  xlab("NMDS 1") +
  ylab("NMDS 2") +
  scale_fill_manual(values = colours) +
  scale_colour_manual(values = colours) +
  scale_shape_manual(values = c(21, 24, 22)) +
  guides(fill = guide_legend(override.aes = list(colour = colours)))

pro.figure

ggsave("Figures/procruestes_NMDS.jpeg")

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



# Presence absence
aquatic.invert[aquatic.invert > 0] <- 1

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
inv.pro$signif # 0.007 #Significant of permutations (999)

plot(inv.pro)

# Uninvaded habitat

unin.pro <- protest(emerg.unin.h, aq.unin.h)

summary(unin.pro)

unin.pro$t0 # 0.8648073
unin.pro$signif # 0.05

# Treated habitat

trt.pro <- protest(emerg.trt.h, aq.trt.h)
summary(trt.pro)

trt.pro$t0 #0.4854401
trt.pro$signif # 0.675

# NMDS -----------------------------------------------

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

corrspp.axis12 <- taxa.axis12df %>% filter(X.taxa.axis12.vectors..r > 0.25)
target12.taxa <- corrspp.axis12$Taxa # string of the Family names


axis12.vectors.c3 <- taxa.in %>% select(all_of(target12.taxa)) # make a matrix of just those

(nmds.c3.vectors.12 <- envfit(emerg.nms$points, axis12.vectors.c3,
                              permutations = 999, choices = c(1,2)))                        


corr.c3.vectors.12 <- as.data.frame(nmds.c3.vectors.12$vectors$arrows*sqrt(nmds.c3.vectors.12$vectors$r)) #scaling vectors
corr.c3.vectors.12$Taxa <- rownames(corr.c3.vectors.12)

write.csv(corr.c3.vectors.12, "Data/Emerging/Procrustes/NMDS_correlated_vector12.csv")

# Figures -----------------------------------------------------------------

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

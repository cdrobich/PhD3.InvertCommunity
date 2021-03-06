library(tidyverse)
library(lubridate)
library(vegan)
library(janitor) # remove_empty

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
emerg.inv.taxa <- emerg.inv %>% select(Chironomidae:Leptoceridae)
emerg.unin.taxa <- emerg.unin %>% select(Lampyridae:Leptoceridae)
emerg.trt.taxa <- emerg.trt %>% select(Chironomidae:Leptoceridae)


#make proportions
colnames(emerg.inv)

emerg.inv.prop <- emerg.inv.taxa/rowSums(emerg.inv.taxa)
emerg.unin.prop <- emerg.unin.taxa/rowSums(emerg.unin.taxa)
emerg.trt.prop <- emerg.trt.taxa/rowSums(emerg.trt.taxa)

# uh puts them in the wrong spot
myList2 <- list(emerg.inv.prop = emerg.inv.prop, 
               emerg.unin.prop = emerg.unin.prop,
               emerg.trt.prop = emerg.trt.prop)
for(i in names(myList2)){
  write.csv(myList2[[i]], paste0(i,".csv"))
}












# Distance Measure and NMDS -----------------------------------------------


emerg.jaccard <- vegdist(emerg.invert, "jaccard", binary = TRUE)

aquat.jaccard <- vegdist(aquatic.invert, "jaccard", binary = TRUE)


### Emergent NMDS

k_vec <- 1:10 #dimensions 1 - 10
stress <- numeric(length(k_vec)) # stress of each model put here
dune_dij <- metaMDSdist(emerg.invert)
set.seed(25)
for(i in seq_along(k_vec)) {
  sol <- metaMDSiter(dune_dij, k = i, # you can tell I lifted this from a tutorial on the dune package
                     trace = FALSE)
  stress[i] <- sol$stress
}

plot(stress) # 3D makes sense


emerg.nms <- metaMDS(emerg.invert, distance = "jaccard", # species data, bray-curtis dissimilarity
                    autotransform = FALSE,  # NMDS will do autotransformations for you
                    k = 3, trymax = 1000)   # k = number of axes
emerg.nms

#Dimensions: 3 
#Stress:     0.1482871 
#Stress type 1, weak ties
#Two convergent solutions found after 52 tries

layout(matrix(1:2, ncol = 2))
plot(emerg.nms, main = "Invertebrate NMDS plot"); stressplot(emerg.nms, main = "Shepard plot")
layout(1)

##### Aquatic NMDS


k_vec <- 1:10 #dimensions 1 - 10
stress <- numeric(length(k_vec)) # stress of each model put here
dune_dij <- metaMDSdist(aquatic.invert)
set.seed(25)
for(i in seq_along(k_vec)) {
  sol <- metaMDSiter(dune_dij, k = i, # you can tell I lifted this from a tutorial on the dune package
                     trace = FALSE)
  stress[i] <- sol$stress
}

plot(stress) # 3D makes sense


aquatic.nms <- metaMDS(aquatic.invert, distance = "jaccard", # species data, bray-curtis dissimilarity
                     autotransform = FALSE,  # NMDS will do autotransformations for you
                     k = 3, trymax = 1000)   # k = number of axes
aquatic.nms

#Data:     aquatic.invert 
#Distance: jaccard 

#Dimensions: 3 
#Stress:     0.1188899 
#Stress type 1, weak ties
#Two convergent solutions found after 20 tries

layout(matrix(1:2, ncol = 2))
plot(aquatic.nms, main = "Invertebrate NMDS plot"); stressplot(aquatic.nms, main = "Shepard plot")
layout(1)


# Procrustes --------------------------------------------------------------



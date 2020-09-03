library(tidyverse) # for figures, pipes etc.
library(vegan) # for NMDS
library(ggpubr) # for arranging figures


### Mantel test ####

# Compare benthic inverts to emerging inverts in each community

benthic <- read.csv("Data/benthics_rares_rel.csv") # veg QCC samples, rares (<=2) removed and rel. by col max
emerge <- read.csv("Data/inverts_rare_rel_2018.csv") # 2018 emergent community, rares and rel. same as above

str(benthic)
str(emerge)

##### Mantel tests ######
?mantel

## uninvaded Mantel first first

uninvaded.b <- subset(benthic,
                      Habitat =="Uninvaded")

uninvaded.b.inverts <- uninvaded.b[,4:25]
uninvaded.b.env <- uninvaded.b[,1:3]


uninvaded.e <- subset(emerge,
                      Habitat =="Uninvaded")

uninvaded.e.inverts <- uninvaded.e[,4:93]
uninvaded.e.env <- uninvaded.e[,1:3]


# B-C dissimilarity matrix 

uninvaded.benthic.bc <- vegdist(uninvaded.b.inverts, method = "bray")
uninvaded.emerge.bc <- vegdist(uninvaded.e.inverts, method = "bray")

(uninvaded.mantel <- mantel(uninvaded.emerge.bc, uninvaded.benthic.bc,
                            method = "pearson",
                            permutations = 999))


#Mantel statistic r: 0.07563 
#Significance: 0.346 

# Upper quantiles of permutations (null model):
#  90%   95%  97.5%   99% 
#  0.241  0.341 0.406 0.456 
# Permutation: free
# Number of permutations: 999


### invaded ###

## invaded Mantel first first

invaded.b <- subset(benthic,
                      Habitat =="Invaded")

invaded.b.inverts <- invaded.b[,4:25]
invaded.b.env <- invaded.b[,1:3]


invaded.e <- subset(emerge,
                      Habitat =="Invaded")

invaded.e.inverts <- invaded.e[,4:93]
invaded.e.env <- invaded.e[,1:3]


# B-C dissimilarity matrix 

invaded.benthic.bc <- vegdist(invaded.b.inverts, method = "bray")
invaded.emerge.bc <- vegdist(invaded.e.inverts, method = "bray")

(invaded.mantel <- mantel(invaded.emerge.bc, invaded.benthic.bc,
                            method = "pearson",
                            permutations = 999))

#Mantel statistic based on Pearson's product-moment correlation 
#Mantel statistic r: -0.3187 
#      Significance: 0.925 

#Upper quantiles of permutations (null model):
#  90%   95% 97.5%   99% 
#0.287 0.350 0.393 0.454 
#Permutation: free
#Number of permutations: 999


###  Restored ###

## Restored Mantel first first

restored.b <- subset(benthic,
                    Habitat =="Restored")

Restored.b.inverts <- restored.b[,4:25]
Restored.b.env <- restored.b[,1:3]


Restored.e <- subset(emerge,
                    Habitat =="Restored")

Restored.e.inverts <- Restored.e[,4:93]
Restored.e.env <- Restored.e[,1:3]


# B-C dissimilarity matrix 

Restored.benthic.bc <- vegdist(Restored.b.inverts, method = "bray")
Restored.emerge.bc <- vegdist(Restored.e.inverts, method = "bray")

(Restored.mantel <- mantel(Restored.emerge.bc, Restored.benthic.bc,
                          method = "pearson",
                          permutations = 999))

#Mantel statistic based on Pearson's product-moment correlation 
#Mantel statistic r: 0.1133 
#      Significance: 0.319 

#Upper quantiles of permutations (null model):
#  90%   95% 97.5%   99% 
# 0.260 0.322 0.362 0.398 
#Permutation: free
#Number of permutations: 999

############ Procrustes #############

(res.pro <- procrustes(Restored.emerge.bc, Restored.benthic.bc))

summary(res.pro)
plot(res.pro)








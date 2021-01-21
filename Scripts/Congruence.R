
# Libraries ---------------------------------------------------------------


library(tidyverse) # for figures, pipes etc.
library(vegan) # for NMDS


# Load data ---------------------------------------------------------------

# deleted INOW and UNOW so same number of rows (25) in both 

benthic <- read.csv("Data/benthic_congruence.csv") # should be only flying taxa
emerge <- read.csv("Data/emerging_congruence.csv") # 2018 emergent community, only flying

str(benthic)
str(emerge)

# just the environmental variables
benthic.env <- benthic %>% select(Site.ID:Sample.Type)
emerge.env <- emerge %>% select(Site:Habitat)

benthic.env$ID <- rownames(benthic.env)
emerge.env$ID <- rownames(emerge.env)

# Just the species
benthic.spp <- benthic %>% select(Lampyridae:Hydroptilidae)
emerge.spp <- emerge %>% select(Phalacridae:Crambidae)


# make them presence/absence

benthic.spp[benthic.spp > 0] <- 1
emerge.spp[emerge.spp >0] <- 1

benthic.spp$ID <- rownames(benthic.spp)
emerge.spp$ID <- rownames(emerge.spp)

# put env variables back in

benthic.pa <- full_join(benthic.env, benthic.spp)
emerge.pa <- full_join(emerge.env, emerge.spp)

##### Mantel tests ######
?mantel

## uninvaded Mantel first first

uninvaded.b <- subset(benthic.pa,
                      Habitat =="Uninvaded")

uninvaded.b.inverts <- uninvaded.b %>% select(Lampyridae:Hydroptilidae)



uninvaded.e <- subset(emerge.pa,
                      Treatment =="Uninvaded")

uninvaded.e.inverts <- uninvaded.e %>% select(Phalacridae:Crambidae)




# B-C dissimilarity matrix 

uninvaded.benthic.bc <- vegdist(uninvaded.b.inverts, method = "bray")
uninvaded.emerge.bc <- vegdist(uninvaded.e.inverts, method = "bray")

(uninvaded.mantel <- mantel(uninvaded.emerge.bc, uninvaded.benthic.bc,
                            method = "pearson",
                            permutations = 999))

# Jaccard

#Mantel statistic r: 0.09216 
#Significance: 0.298 

#Upper quantiles of permutations (null model):
#  90%   95% 97.5%   99% 
#  0.255 0.321 0.383 0.472 
#Permutation: free
#Number of permutations: 999


# Bray Curtis

#Mantel statistic r: 0.1223 
#Significance: 0.248 

#Upper quantiles of permutations (null model):
#  90%   95% 97.5%   99% 
#  0.244 0.306 0.380 0.425 
#Permutation: free
#Number of permutations: 999







### invaded ###

## invaded Mantel first first

invaded.b <- subset(benthic,
                      Habitat =="Invaded")

invaded.b.inverts <- invaded.b %>% select(Lampyridae:Hydroptilidae)


invaded.e <- subset(emerge,
                      Treatment =="Invaded")

invaded.e.inverts <- invaded.e %>% select(Phalacridae:Crambidae)



# B-C dissimilarity matrix 

invaded.benthic.bc <- vegdist(invaded.b.inverts, method = "bray")
invaded.emerge.bc <- vegdist(invaded.e.inverts, method = "bray")

(invaded.mantel <- mantel(invaded.emerge.bc, invaded.benthic.bc,
                            method = "pearson",
                            permutations = 999))

# Jaccard

#Mantel statistic r: -0.0658 
#Significance: 0.61 

#Upper quantiles of permutations (null model):
#  90%   95% 97.5%   99% 
#  0.203 0.259 0.289 0.359 
#Permutation: free
#Number of permutations: 999


# Bray-Curtis

#Mantel statistic r: -0.04851 
#Significance: 0.561 

#Upper quantiles of permutations (null model):
#  90%   95% 97.5%   99% 
#  0.259 0.303 0.340 0.372 
#Permutation: free
#Number of permutations: 999






### Treated ###

## Restored Mantel first first

restored.b <- subset(benthic,
                    Habitat =="Treated")

Restored.b.inverts <- restored.b %>% select(Lampyridae:Hydroptilidae)


Restored.e <- subset(emerge,
                    Treatment =="Treated")

Restored.e.inverts <- Restored.e %>% select(Phalacridae:Crambidae)



# dissimilarity matrix 

Restored.benthic.bc <- vegdist(Restored.b.inverts, method = "bray")
Restored.emerge.bc <- vegdist(Restored.e.inverts, method = "bray")

(Restored.mantel <- mantel(Restored.emerge.bc, Restored.benthic.bc,
                          method = "pearson",
                          permutations = 999))
# Jaccard Distance

#Mantel statistic r: -0.3519 
#Significance: 0.925 

#Upper quantiles of permutations (null model):
#  90%   95% 97.5%   99% 
#  0.439 0.568 0.644 0.687 
#Permutation: free
#Number of permutations: 999


# Bray Curtis

#Mantel statistic r: -0.3498 
#Significance: 0.936 

#Upper quantiles of permutations (null model):
#  90%   95% 97.5%   99% 
#  0.567 0.666 0.743 0.788 
#Permutation: free
#Number of permutations: 999






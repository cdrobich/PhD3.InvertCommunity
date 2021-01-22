
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

#Mantel statistic r: -0.217 
#Significance: 0.878 

#Upper quantiles of permutations (null model):
#  90%   95% 97.5%   99% 
#  0.246 0.326 0.389 0.426 
#Permutation: free
#Number of permutations: 999


# Bray Curtis


#Mantel statistic r: -0.2184 
#Significance: 0.884 

#Upper quantiles of permutations (null model):
#  90%   95% 97.5%   99% 
#  0.249 0.308 0.366 0.421 
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

#Mantel statistic r: 0.1327 
#Significance: 0.197 

#Upper quantiles of permutations (null model):
#  90%   95% 97.5%   99% 
#  0.197 0.261 0.296 0.345 
#Permutation: free
#Number of permutations: 999


# Bray-Curtis

#Mantel statistic r: 0.197 
#Significance: 0.171 

#Upper quantiles of permutations (null model):
#  90%   95% 97.5%   99% 
#  0.268 0.314 0.341 0.393 
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

#Mantel statistic r: 0.5635 
#Significance: 0.06 

#Upper quantiles of permutations (null model):
#  90%   95% 97.5%   99% 
#  0.458 0.586 0.667 0.706 
#Permutation: free
#Number of permutations: 999


# Bray Curtis

#Mantel statistic r: 0.7026 
#Significance: 0.027 

#Upper quantiles of permutations (null model):
#  90%   95% 97.5%   99% 
#  0.532 0.632 0.709 0.782 
#Permutation: free
#Number of permutations: 999


library(ecodist)

mantel(Restored.emerge.bc ~ Restored.benthic.bc)

# mantelr      pval1      pval2      pval3  llim.2.5% ulim.97.5% 
# 0.7025840  0.0310000  0.9700000  0.0310000  0.1321817  0.8329853



mantel(invaded.emerge.bc ~ invaded.benthic.bc)

#mantelr      pval1      pval2      pval3  llim.2.5% ulim.97.5% 
#0.1970081  0.1430000  0.8580000  0.2550000 -0.1159272  0.6508261




mantel(uninvaded.emerge.bc ~ uninvaded.benthic.bc)

#mantelr       pval1       pval2       pval3   llim.2.5%  ulim.97.5% 
#-0.21843955  0.87700000  0.12400000  0.26100000 -0.39572492 -0.03993152 




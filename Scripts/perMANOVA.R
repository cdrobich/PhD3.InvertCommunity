
library(tidyverse) # for figures, pipes etc.
library(vegan) # for NMDS
library(ggpubr) # for arranging figures
library(cluster) # for agnes (agglomerative hierarchical clustering)

sol <- read.csv("Data/Inverts_relmax_raresn2.csv")

bugs <- sol[,5:94] # just the Families, data is rel. by col max and rares (<= 2) removed
env <- sol[,1:4] # Categorical variables

env$TrtYr <- as.factor(env$TrtYr) 
env$Treatment <- as.factor(env$Treatment)

## perMANOVA

(per.treat <- adonis(bugs ~ Treatment, data = env,
                    permutations = 999, method = "bray"))

#adonis(formula = bugs ~ Treatment, data = env, 
#   permutations = 999,      method = "bray") 

#Permutation: free
#Number of permutations: 999

#Terms added sequentially (first to last)

#             Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
#  Treatment  2    2.4255 1.21273  3.4758 0.11995  0.001 ***
#  Residuals 51   17.7944 0.34891         0.88005           
# Total     53   20.2199                 1.00000 



# homogeneity of groups dispersion

bugs.b <- vegdist(bugs, method = "bray")

groups <- factor(env$Treatment)

(dispersion <- betadisper(bugs.b, groups)) # spatial median default

plot(dispersion)
boxplot(dispersion) # actually look really good!

install.packages("devtools")
library(devtools)

?install_github
install_github("GuillemSalazar/EcolUtils")

library(EcolUtils)
citation("EcolUtils")

(adonis.pair(bugs.b, groups, nper = 1000, corr.method = "bonferroni"))

#            combination   SumsOfSqs   MeanSqs  F.Model       R2     P.value  P.value.corrected
#1   Invaded <-> Restored 1.7139459 1.7139459 5.100551 0.13044704 0.000999001       0.002997003
#2  Invaded <-> Uninvaded 0.5961001 0.5961001 1.655632 0.04643395 0.015984016       0.047952048
#3 Restored <-> Uninvaded 1.3281515 1.3281515 3.787631 0.10023467 0.000999001       0.002997003

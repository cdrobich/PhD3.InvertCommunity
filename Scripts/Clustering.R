
library(tidyverse) # for figures, pipes etc.
library(vegan) # for NMDS
library(ggpubr) # for arranging figures
library(cluster) # for agnes (agglomerative hierarchical clustering)

sol <- read.csv("Data/Inverts_relmax_raresn2.csv")

bugs <- sol[,5:94] # just the Families, data is rel. by col max and rares (<= 2) removed
env <- sol[,1:4] # Categorical variables

env$TrtYr <- as.factor(env$TrtYr) 
env$Treatment <- as.factor(env$Treatment)

# Agglomeration Hierarchical Clustering Analysis
# use sorensen distance, flexible beta linkage (-0.25)

library(labdsv) # Dufrene-Legendre Indicator Species Analysis

cluster <- env$Treatment

bugs.b <- vegdist(bugs, method = "bray") # Bray Curtis dissimilarity matrix

?agnes

# flexible beta

(bugs.cluster <- agnes(bugs.b, diss=TRUE, method = "flexible",
                       par.method = c(0.625, 0.625, -0.25))) 

##https://stat.ethz.ch/R-manual/R-devel/library/cluster/html/agnes.html
#alpha = 0.625 ==> beta = -1/4  is "recommended" by some

bugs.cluster$ac # 0.598717 agglomerative coefficient (1 = strongest cluster structure)
plot(bugs.cluster)

# all birds
ncluster.bd.isa <-vector("list",1)
for (i in 2:20){
  ncluster.bd.isa[[i]] <- cutree(as.hclust(allbirds.agnes.flex),k=i)
}


### Did some stuff in PCORD
# found that a group of 3 is best

groups <- c("10","9","8","7","6","5","4","3","2")
pvalue <- c("0.2878","0.319","0.292","0.2613","0.216","0.1757","0.1636","0.1847","0.239")
sig.ind <- c("20","15","17","19","34","35","40","42","33")
scree.plot <- data.frame(groups, pvalue, sig.ind)
scree.plot
plot(scree.plot)

str(scree.plot)
scree.plot$groups <- as.factor(scree.plot$groups)

scree.plot <- scree.plot %>% 
  mutate(groups = fct_relevel(groups,
                              "2","3","4","5","6","7","8","9","10"))


p.clust <- ggplot(scree.plot, aes(x = groups, y = pvalue)) +
  geom_point(size = 4) +
  theme_classic(base_size = 16) +
  theme(panel.border = element_rect(fill = NA)) +
  xlab("Number of Groups") +
  ylab("Average P-value") 

p.clust

ind.clust <- ggplot(scree.plot, aes(x = groups, y = sig.ind)) +
  geom_point(size = 4) +
  theme_classic(base_size = 16) +
  theme(panel.border = element_rect(fill = NA)) +
  xlab("Number of Groups") +
  ylab("Significant Indicators") 

ind.clust

assessment <- ggarrange(p.clust, ind.clust,
                        labels = "AUTO",
                        hjust = c(-9, -6),
                        vjust = 2.5)
assessment

ggsave("Figures/Screeplot_clusters.jpeg", assessment)

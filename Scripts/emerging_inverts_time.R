
# Load packages -----------------------------------------------------------

library(tidyverse) # for figures, pipes etc.
library(vegan) # for NMDS
library(ggpubr) # for arranging figures

library(ggrepel) # labels on nmds
library(lubridate)

library(viridis) # colours

# pairwise perMANOVA comparison
library(devtools)
install_github("GuillemSalazar/EcolUtils")
library(EcolUtils)

library(Hmisc)

# Load data ---------------------------------------------------------------

inverts.date <- read.csv("Data/Emerging/ermerging_time_rares.csv") # Occurences <= 2 removed
head(inverts.date)
str(inverts.date)

unique(inverts.date$Treatment)
colnames(inverts.date)

taxa <- inverts.date %>% select(Araneae:Crambidae)
env <- inverts.date %>% select(ID:YrCol)

# Collection periods -----------------------------------------------

env <- env %>% separate(Date, c("Day","Month", "Year"), remove = FALSE)

unique(env$Date)

#                        "19-Jun-17" "28-Jun-17" "08-Jul-17" "21-Jul-17" 
#"20-May-18" "05-Jun-18" "16-Jun-18" "25-Jun-18" "04-Jul-18" "23-Jul-18"


# Collection 1:Str June 19 2017 and June 16 2018  -----------------------------------------

col.1 <- c("19-Jun-17","16-Jun-18")

invert.col.1 <- inverts.date %>% filter(Date %in% col.1)

write.csv(invert.col.1, "Data/Emerging/NMDS/Col.1/inverts_collection1.csv")

# Empty columns removed, relativized by col max. in PCORD

col1.data <- read.csv("Data/Emerging/NMDS/Col.1/col1_zeroes_relativized.csv")
colnames(col1.data)

taxa.col1 <- col1.data %>% select(Araneae:Crambidae)
env.col1 <- col1.data %>% select(Sites:YrCol)


## perMANOVA
taxa.col1

(per.col1 <- adonis2(taxa.col1 ~ Treatment * Year, data = col1.data,
                      permutations = 999, method = "bray"))


#               Df SumOfSqs      R2      F Pr(>F)    
#Treatment       2   1.9848 0.09578 2.6608  0.001 ***
#Year            1   0.8817 0.04255 2.3641  0.001 ***
#Treatment:Year  2   1.0734 0.05180 1.4390  0.025 *  
#Residual       45  16.7839 0.80988                  
#Total          50  20.7239 1.00000



#### NMDS analysis 

set.seed(120) 

nms.col1 <- metaMDS(taxa.col1, distance = "bray", # species data, bray-curtis dissimilarity
                      autotransform = FALSE,  # NMDS will do autotransformations for you
                      k = 3, trymax = 1000)   # k = number of axes
nms.col1

#Data:     taxa.col1 
#Distance: bray 

#Dimensions: 3 
#Stress:     0.1933805 
#Stress type 1, weak ties
#Two convergent solutions found after 26 tries
#Scaling: centring, PC rotation, halfchange scaling 
#Species: expanded scores based on ‘taxa.col1’ 

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
sum(g^2) #0.037396
nms.col1$stress^2  

1-nms.col1$stress^2 # 0.962604


## extract the scores for plotting 
scr <- as.data.frame(scores(nms.col1, display = "sites")) # extract NMDS scores
colnames(scr)

# adding categorical info to scores
env.col1$NMDS1 <- scr$NMDS1
env.col1$NMDS2 <- scr$NMDS2
env.col1$NMDS3 <- scr$NMDS3

write.csv(env.col1,"Data/Emerging/NMDS/Col.1/emerging_inverts_col1_NMDSscores.csv") # save this as a csv

### Taxa correlated with Axis 1 & 2 

col1.12 <- envfit(nms.col1, taxa.col1,
                    choices = c(1,2)) #produces a list with r2, p value, and NMDS coordinates

col1.df12 <- data.frame((col1.12$vectors)$arrows,
                          (col1.12$vectors)$r,
                          (col1.12$vectors)$pvals) #take list and make into data frame

col1.df12 <- tibble::rownames_to_column(col1.df12, "Taxa")

write.csv(col1.df12, "Data/Emerging/NMDS/Col.1/NMDS_emerg_col1_vector12.csv") # save vector scores as csv


## Taxa correlated with Axis 1 & 3

col1.13 <- envfit(nms.col1, taxa.col1, 
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

axis12.vectors.c1 <- taxa.col1 %>% select(all_of(target12.c1)) # make a matrix of just those


# Selecting taxa for vectors on axis 1, 3
colnames(col1.df13)
col1.corrsp.13 <- col1.df13 %>% filter(X.col1.13.vectors..r > 0.2) 
target13.c1 <- col1.corrsp.13$Taxa # string of the Family names

axis13.vectors.c1 <- taxa.col1 %>% select(all_of(target13.c1)) # make a matrix of just those


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



# Collection 2: 28-Jun-17 and 25-JUN-18 ---------------------------------

col.2 <- c("28-Jun-17","25-Jun-18")

invert.col.2 <- inverts.date %>% filter(Date %in% col.2)

write.csv(invert.col.2, "Data/Emerging/NMDS/Col.2/inverts_collection2.csv")

# relativized by col max and removed empty columns 

col2.data <- read.csv("Data/Emerging/NMDS/Col.2/col2_zeroes_relativized.csv")

# just taxa and env 
colnames(col2.data)

taxa.col2 <- col2.data %>% select(Araneae:Noctuidae)
env.col2 <- col2.data %>% select(Site:YrCol)

## perMANOVA

(per.col2 <- adonis2(taxa.col2 ~ Treatment * Year, data = col2.data,
                     permutations = 999, method = "bray"))

#adonis2(formula = taxa.col2 ~ Treatment * Year, data = col2.data, permutations = 999, method = "bray")
#               Df SumOfSqs      R2      F Pr(>F)    
#Treatment       2   1.7543 0.08037 2.2080  0.001 ***
#Year            1   0.7473 0.03423 1.8811  0.003 ** 
#Treatment:Year  2   1.0525 0.04822 1.3247  0.058 .  
#Residual       46  18.2736 0.83718                  
#Total          51  21.8276 1.00000



col2.b <- vegdist(taxa.col2, method = "bray")

trt.col2 <- factor(env.col2$Treatment)
yr.col2 <- factor(env.col2$Year)

(dispersion <- betadisper(col2.b, trt.col2))
plot(dispersion)

(adonis.pair(col2.b, trt.col2, 
             nper = 1000, corr.method = "bonferroni"))


#            combination SumsOfSqs   MeanSqs  F.Model         R2     P.value P.value.corrected
#1   Invaded <-> Treated 1.0171130 1.0171130 2.567345 0.07427083 0.000999001       0.002997003
#2 Invaded <-> Uninvaded 0.6503846 0.6503846 1.520858 0.04281592 0.031968032       0.095904096
#3 Treated <-> Uninvaded 0.9772179 0.9772179 2.418616 0.07027058 0.001998002       0.005994006




### NMDS 

k_vec <- 1:10 #dimensions 1 - 10
stress <- numeric(length(k_vec)) # stress of each model put here
dune_dij <- metaMDSdist(taxa.col2)

set.seed(25)

for(i in seq_along(k_vec)) {
  sol <- metaMDSiter(dune_dij, k = i, 
                     trace = FALSE)
  stress[i] <- sol$stress
}
plot(stress) # 3D



set.seed(120) 

nms.col2 <- metaMDS(taxa.col2, distance = "bray", # species data, bray-curtis dissimilarity
                    autotransform = FALSE,  # NMDS will do autotransformations for you
                    k = 3, trymax = 1000)   # k = number of axes
nms.col2


#global Multidimensional Scaling using monoMDS

#Data:     taxa.col2 
#Distance: bray 

#Dimensions: 3 
#Stress:     0.2153085 
#Stress type 1, weak ties
#Two convergent solutions found after 20 tries

# look at points and the stress real quick
layout(matrix(1:2, ncol = 2))
plot(nms.col2, main = "Invertebrate NMDS plot"); stressplot(nms.col2, main = "Shepard plot")
layout(1)

ordiplot(nms.col2, type = "n")
orditorp(nms.col2, display = "species")
orditorp(nms.col2, display = "sites")



nms.col2$iters #97

nms.col2$stress^2   #0.04635775
1-nms.col2$stress^2 #0.9536422

scr2 <- as.data.frame(scores(nms.col2, display = "sites")) # extract NMDS scores
colnames(scr2)

env.col2$NMDS1 <- scr2$NMDS1
env.col2$NMDS2 <- scr2$NMDS2
env.col2$NMDS3 <- scr2$NMDS3

write.csv(env.col2,"Data/Emerging/NMDS/Col.2/emerging_col2_NMDSscores.csv") # save this as a csv

## Taxa for vectors

# Axis 1 and 2
taxa.col2.axis12 <- envfit(nms.col2, taxa.col2,
                           choices = c(1,2))

taxa.col2.axis12df <- data.frame((taxa.col2.axis12$vectors)$arrows,
                                 (taxa.col2.axis12$vectors)$r,
                                 (taxa.col2.axis12$vectors)$pvals)

taxa.col2.axis12df <- tibble::rownames_to_column(taxa.col2.axis12df, "Taxa")

write.csv(taxa.col2.axis12df, "Data/Emerging/NMDS/Col.2/col2_allvectors_axis12.csv") # save vector scores as csv

# Axis 1 and 3

taxa.col2.axis13 <- envfit(nms.col2, taxa.col2,
                           choices = c(1,3))

taxa.col2.axis13df <- data.frame((taxa.col2.axis13$vectors)$arrows,
                                 (taxa.col2.axis13$vectors)$r,
                                 (taxa.col2.axis13$vectors)$pvals)

taxa.col2.axis13df <- tibble::rownames_to_column(taxa.col2.axis13df, "Taxa")

write.csv(taxa.col2.axis13df, "Data/Emerging/NMDS/Col.2/col2_allvectors_axis13.csv") # save vector scores as csv

# Correlated taxa

# Axis 1 and 2

colnames(taxa.col2.axis12df)

corrspp.col2.axis12 <- taxa.col2.axis12df %>% filter(X.taxa.col2.axis12.vectors..r > 0.2)
target12.c2 <- corrspp.col2.axis12$Taxa # string of the Family names


axis12.vectors.c2 <- taxa.col2 %>% select(all_of(target12.c2)) # make a matrix of just those

(nmds.c2.vectors.12 <- envfit(nms.col2$points, axis12.vectors.c2,
                              permutations = 999, choices = c(1,2)))                        


corr.c2.vectors.12 <- as.data.frame(nmds.c2.vectors.12$vectors$arrows*sqrt(nmds.c2.vectors.12$vectors$r)) #scaling vectors
corr.c2.vectors.12$Taxa <- rownames(corr.c2.vectors.12)


write.csv(corr.c2.vectors.12, "Data/Emerging/NMDS/Col.2/emerging_correlated_vector12.csv")

## Axis 1 and 3

colnames(taxa.col2.axis13df)

corrspp.col2.axis13 <- taxa.col2.axis13df %>% filter(X.taxa.col2.axis13.vectors..r > 0.2)
target13.c2 <- corrspp.col2.axis13$Taxa # string of the Family names


axis13.vectors.c2 <- taxa.col2 %>% select(all_of(target13.c2)) # make a matrix of just those

(nmds.c2.vectors.13 <- envfit(nms.col2$points, axis13.vectors.c2,
                              permutations = 999, choices = c(1,3)))                        


corr.c2.vectors.13 <- as.data.frame(nmds.c2.vectors.13$vectors$arrows*sqrt(nmds.c2.vectors.13$vectors$r)) #scaling vectors
corr.c2.vectors.13$Taxa <- rownames(corr.c2.vectors.13)


write.csv(corr.c2.vectors.13, "Data/Emerging/NMDS/Col.2/emerging_correlated_vector13.csv")










# Collection 3: 08-Jul-17 and 04-Jul-18 -----------------------------------

col.3 <- c("08-Jul-17","04-Jul-18")

invert.col.3 <- inverts.date %>% filter(Date %in% col.3)
write.csv(invert.col.3, "Data/Emerging/NMDS/Col.3/inverts_collection3.csv")


col3.data <- read.csv("Data/Emerging/NMDS/Col.3/col3_zeros_rares.csv")

# just taxa and env 
colnames(col3.data)
taxa.col3 <- col3.data %>% select(Araneae:Crambidae)
env.col3 <- col3.data %>% select(ID:YrCol)

### perMANOVA

(per.col3 <- adonis2(taxa.col3 ~ Treatment * Year, data = col3.data,
                     permutations = 999, method = "bray"))

#               Df SumOfSqs      R2      F Pr(>F)    
#Treatment       2   2.0057 0.09061 2.6004  0.001 ***
#Year            1   0.6452 0.02915 1.6729  0.018 *  
#Treatment:Year  2   0.9722 0.04392 1.2605  0.071 .  
#Residual       48  18.5116 0.83632                  
#Total          53  22.1348 1.00000    

col3.b <- vegdist(taxa.col3, method = "bray")

trt.col3 <- factor(env.col3$Treatment)


(dispersion <- betadisper(col3.b, trt.col3))
plot(dispersion)

(adonis.pair(col3.b, trt.col3, 
             nper = 1000, corr.method = "bonferroni"))

#             combination SumsOfSqs   MeanSqs  F.Model         R2     P.value P.value.corrected
#1   Invaded <-> Treated 1.2138287 1.2138287 3.240205 0.08700825 0.000999001       0.002997003
#2 Invaded <-> Uninvaded 0.8156821 0.8156821 1.996619 0.05546685 0.000999001       0.002997003
#3 Treated <-> Uninvaded 0.9790501 0.9790501 2.442039 0.06701160 0.000999001       0.002997003

# NMDS collection 3

k_vec <- 1:10 #dimensions 1 - 10
stress <- numeric(length(k_vec)) # stress of each model put here
dune_dij <- metaMDSdist(taxa.col3)

set.seed(25)

for(i in seq_along(k_vec)) {
  sol <- metaMDSiter(dune_dij, k = i, 
                     trace = FALSE)
  stress[i] <- sol$stress
}
plot(stress) # 3D


# Actual NMDS

set.seed(120) 

nms.col3 <- metaMDS(taxa.col3, distance = "bray", # species data, bray-curtis dissimilarity
                    autotransform = FALSE,  # NMDS will do autotransformations for you
                    k = 3, trymax = 1000)   # k = number of axes
nms.col3

#global Multidimensional Scaling using monoMDS

#Data:     taxa.col3 
#Distance: bray 

#Dimensions: 3 
#Stress:     0.2088858 
#Stress type 1, weak ties
#Two convergent solutions found after 43 tries 

nms.col3$iters #126

nms.col3$stress^2   #0.04363328
1-nms.col3$stress^2 #0.9563667

scr3 <- as.data.frame(scores(nms.col3, display = "sites")) # extract NMDS scores
colnames(scr3)

env.col3$NMDS1 <- scr3$NMDS1
env.col3$NMDS2 <- scr3$NMDS2
env.col3$NMDS3 <- scr3$NMDS3

write.csv(env.col3,"Data/Emerging/NMDS/Col.3/emerging_col3_NMDSscores.csv") # save this as a csv

# Vector Taxa

# Axis 1 and 2
taxa.col3.axis12 <- envfit(nms.col3, taxa.col3,
                           choices = c(1,2))

taxa.col3.axis12df <- data.frame((taxa.col3.axis12$vectors)$arrows,
                                 (taxa.col3.axis12$vectors)$r,
                                 (taxa.col3.axis12$vectors)$pvals)

taxa.col3.axis12df <- tibble::rownames_to_column(taxa.col3.axis12df, "Taxa")

write.csv(taxa.col3.axis12df, "Data/Emerging/NMDS/Col.3/col3_allvectors_axis12.csv") # save vecto


# Axis 1 and 3

taxa.col3.axis13 <- envfit(nms.col3, taxa.col3,
                           choices = c(1,3))

taxa.col3.axis13df <- data.frame((taxa.col3.axis13$vectors)$arrows,
                                 (taxa.col3.axis13$vectors)$r,
                                 (taxa.col3.axis13$vectors)$pvals)

taxa.col3.axis13df <- tibble::rownames_to_column(taxa.col3.axis13df, "Taxa")

write.csv(taxa.col3.axis13df, "Data/Emerging/NMDS/Col.3/col3_allvectors_axis13.csv") # save vector scores as csv

# Correlated taxa

# Axis 1 and 2
colnames(taxa.col3.axis12df)

corrspp.col3.axis12 <- taxa.col3.axis12df %>% filter(X.taxa.col3.axis12.vectors..r > 0.2)
target12.c3 <- corrspp.col3.axis12$Taxa # string of the Family names


axis12.vectors.c3 <- taxa.col3 %>% select(all_of(target12.c3)) # make a matrix of just those

(nmds.c3.vectors.12 <- envfit(nms.col3$points, axis12.vectors.c3,
                              permutations = 999, choices = c(1,2)))                        


corr.c3.vectors.12 <- as.data.frame(nmds.c3.vectors.12$vectors$arrows*sqrt(nmds.c3.vectors.12$vectors$r)) #scaling vectors
corr.c3.vectors.12$Taxa <- rownames(corr.c3.vectors.12)

write.csv(corr.c3.vectors.12, "Data/Emerging/NMDS/Col.3/emerging_correlated_vector12.csv")

# Axis 1 and 3

colnames(taxa.col3.axis13df)

corrspp.col3.axis13 <- taxa.col3.axis13df %>% filter(X.taxa.col3.axis13.vectors..r > 0.2)
target13.c3 <- corrspp.col3.axis13$Taxa # string of the Family names


axis13.vectors.c3 <- taxa.col3 %>% select(all_of(target13.c3)) # make a matrix of just those

(nmds.c3.vectors.13 <- envfit(nms.col3$points, axis13.vectors.c3,
                              permutations = 999, choices = c(1,3)))                        


corr.c3.vectors.13 <- as.data.frame(nmds.c3.vectors.13$vectors$arrows*sqrt(nmds.c3.vectors.13$vectors$r)) #scaling vectors
corr.c3.vectors.13$Taxa <- rownames(corr.c3.vectors.13)

write.csv(corr.c3.vectors.13, "Data/Emerging/NMDS/Col.3/emerging_correlated_vector13.csv")


# Collection 4: 21-Jul-17 and 23-Jul-18 -----------------------------------

col.4 <- c("21-Jul-17","23-Jul-18")

invert.col.4 <- inverts.date %>% filter(Date %in% col.4)

write.csv(invert.col.4, "Data/Emerging/NMDS/Col.4/inverts_collection4.csv")

col4.data <- read.csv("Data/Emerging/NMDS/Col.4/col4_zeros_rares.csv")

# just taxa and env 
colnames(col4.data)
taxa.col4 <- col4.data  %>% select(Araneae:Crambidae)
env.col4 <- col4.data  %>% select(ID:YrCol)


## perMANOVA
(per.col4 <- adonis2(taxa.col4 ~ Treatment * Year, data = col4.data,
                     permutations = 999, method = "bray"))


#               Df SumOfSqs      R2      F Pr(>F)    
#Treatment       2   1.7464 0.08366 2.3216  0.001 ***
#Year            1   0.8085 0.03873 2.1496  0.001 ***
#Treatment:Year  2   1.3951 0.06683 1.8545  0.001 ***
#Residual       45  16.9259 0.81078                  
#Total          50  20.8759 1.00000 



# NMDS collection 4
k_vec <- 1:10 #dimensions 1 - 10
stress <- numeric(length(k_vec)) # stress of each model put here
dune_dij <- metaMDSdist(taxa.col4)

set.seed(25)

for(i in seq_along(k_vec)) {
  sol <- metaMDSiter(dune_dij, k = i, 
                     trace = FALSE)
  stress[i] <- sol$stress
}
plot(stress) # 3D


# actual NMDS 

set.seed(120) 

nms.col4 <- metaMDS(taxa.col4, distance = "bray", # species data, bray-curtis dissimilarity
                    autotransform = FALSE,  # NMDS will do autotransformations for you
                    k = 3, trymax = 1000)   # k = number of axes
nms.col4

#global Multidimensional Scaling using monoMDS

#Data:     taxa.col4 
#Distance: bray 

#Dimensions: 3 
#Stress:     0.2204784 
#Stress type 1, weak ties
#Two convergent solutions found after 20 tries

nms.col4$iters #192

nms.col4$stress^2   #0.04861074
1-nms.col4$stress^2 #0.9513893


scr4 <- as.data.frame(scores(nms.col4, display = "sites")) # extract NMDS scores
colnames(scr4)

env.col4$NMDS1 <- scr4$NMDS1
env.col4$NMDS2 <- scr4$NMDS2
env.col4$NMDS3 <- scr4$NMDS3

write.csv(env.col4,"Data/Emerging/NMDS/Col.4/emerging_col4_NMDSscores.csv") # save this as a csv

## Vectors 

# Axis 1 and 2
taxa.col4.axis12 <- envfit(nms.col4, taxa.col4,
                           choices = c(1,2))

taxa.col4.axis12df <- data.frame((taxa.col4.axis12$vectors)$arrows,
                                 (taxa.col4.axis12$vectors)$r,
                                 (taxa.col4.axis12$vectors)$pvals)

taxa.col4.axis12df <- tibble::rownames_to_column(taxa.col4.axis12df, "Taxa")

write.csv(taxa.col4.axis12df, "Data/Emerging/NMDS/Col.4/col4_allvectors_axis12.csv") # save vecto


# Axis 1 and 3

taxa.col4.axis13 <- envfit(nms.col4, taxa.col4,
                           choices = c(1,3))

taxa.col4.axis13df <- data.frame((taxa.col4.axis13$vectors)$arrows,
                                 (taxa.col4.axis13$vectors)$r,
                                 (taxa.col4.axis13$vectors)$pvals)

taxa.col4.axis13df <- tibble::rownames_to_column(taxa.col4.axis13df, "Taxa")

write.csv(taxa.col4.axis13df, "Data/Emerging/NMDS/Col.4/col4_allvectors_axis13.csv") # save vecto


# Axis 1 and 2
colnames(taxa.col4.axis12df)

corrspp.col4.axis12 <- taxa.col4.axis12df %>% filter(X.taxa.col4.axis12.vectors..r > 0.2)
target12.c4 <- corrspp.col4.axis12$Taxa # string of the Family names


axis12.vectors.c4 <- taxa.col4 %>% select(all_of(target12.c4)) # make a matrix of just those

(nmds.c4.vectors.12 <- envfit(nms.col4$points, axis12.vectors.c4,
                              permutations = 999, choices = c(1,2)))                        


corr.c4.vectors.12 <- as.data.frame(nmds.c4.vectors.12$vectors$arrows*sqrt(nmds.c4.vectors.12$vectors$r)) #scaling vectors
corr.c4.vectors.12$Taxa <- rownames(corr.c4.vectors.12)

write.csv(corr.c4.vectors.12, "Data/Emerging/NMDS/Col.4/emerging_correlated_vector12.csv")

# axis 1 and 3

colnames(taxa.col4.axis13df)

corrspp.col4.axis13 <- taxa.col4.axis13df %>% filter(X.taxa.col4.axis13.vectors..r > 0.2)
target13.c4 <- corrspp.col4.axis13$Taxa # string of the Family names


axis13.vectors.c4 <- taxa.col4 %>% select(all_of(target13.c4)) # make a matrix of just those

(nmds.c4.vectors.13 <- envfit(nms.col4$points, axis13.vectors.c4,
                              permutations = 999, choices = c(1,3)))                        


corr.c4.vectors.13 <- as.data.frame(nmds.c4.vectors.13$vectors$arrows*sqrt(nmds.c4.vectors.13$vectors$r)) #scaling vectors
corr.c4.vectors.13$Taxa <- rownames(corr.c4.vectors.13)

write.csv(corr.c4.vectors.13, "Data/Emerging/NMDS/Col.4/emerging_correlated_vector13.csv")




# May 20 2018 -------------------------------------------------------------

invert.May <- inverts.date %>% filter(Date == "20-May-18")

write.csv(invert.May, "Data/Emerging/NMDS/inverts_20-May-18.csv")



may.data <- read.csv("Data/Emerging/NMDS/may_zeros_rares.csv")
colnames(may.data)

taxa.may <- may.data %>% select(Araneae:Limnephilida)
env.may <- may.data %>% select(Sites:YrCol)


## perMANOVA
(per.may <- adonis2(taxa.may ~ Treatment * Year, data = may.data,
                     permutations = 999, method = "bray"))


#           Df SumOfSqs      R2      F Pr(>F)    
#Treatment  2   2.1306 0.23086 3.1517  0.001 ***
#Residual  21   7.0984 0.76914                  
#Total     23   9.2290 1.00000 



may.b <- vegdist(taxa.may, method = "bray")

trt.may <- factor(env.may$Treatment)


(dispersion <- betadisper(may.b,trt.may))
plot(dispersion)

(adonis.pair(may.b, trt.may, 
             nper = 1000, corr.method = "bonferroni"))


#           combination SumsOfSqs   MeanSqs  F.Model        R2     P.value P.value.corrected
#1   Invaded <-> Treated 1.3545809 1.3545809 4.062205 0.2249008 0.000999001       0.002997003
#2 Invaded <-> Uninvaded 0.6701022 0.6701022 1.946093 0.1220420 0.032967033       0.098901099
#3 Treated <-> Uninvaded 1.1712700 1.1712700 3.483203 0.1992314 0.000999001       0.002997003

## NMDS

k_vec <- 1:10 #dimensions 1 - 10
stress <- numeric(length(k_vec)) # stress of each model put here
dune_dij <- metaMDSdist(taxa.may)

set.seed(25)

for(i in seq_along(k_vec)) {
  sol <- metaMDSiter(dune_dij, k = i, 
                     trace = FALSE)
  stress[i] <- sol$stress
}
plot(stress) # Going with two based on PCORD



set.seed(120) 

nms.may <- metaMDS(taxa.may, distance = "bray", # species data, bray-curtis dissimilarity
                    autotransform = FALSE,  # NMDS will do autotransformations for you
                    k = 2, trymax = 1000)   # k = number of axes
nms.may

#Data:     taxa.may 
#Distance: bray 

#Dimensions: 2 
#Stress:     0.2425735 
#Stress type 1, weak ties
#Two convergent solutions found after 20 tries
#Scaling: centring, PC rotation, halfchange scaling 
#Species: expanded scores based on ‘taxa.may’

layout(matrix(1:2, ncol = 2))
plot(nms.may, main = "Invertebrate NMDS plot"); stressplot(nms.may, main = "Shepard plot")
layout(1)

ordiplot(nms.may, type = "n")
orditorp(nms.may, display = "species")
orditorp(nms.may, display = "sites")


nms.may$iters #82

nms.may$stress^2   # 0.0588419
1-nms.may$stress^2 # 0.9411581

scr.may <- as.data.frame(scores(nms.may, display = "sites")) # extract NMDS scores
colnames(scr.may)

env.may$NMDS1 <- scr.may$NMDS1
env.may$NMDS2 <- scr.may$NMDS2
env.may$NMDS3 <- scr.may$NMDS3


write.csv(env.may,"Data/Emerging/NMDS/emerging_May18_NMDSscores.csv") # save this as a csv



## Taxa for vectors

# Axis 1 and 2
taxa.may.axis12 <- envfit(nms.may, taxa.may,
                           choices = c(1,2))

taxa.may.axis12df <- data.frame((taxa.may.axis12$vectors)$arrows,
                                 (taxa.may.axis12$vectors)$r,
                                 (taxa.may.axis12$vectors)$pvals)

taxa.may.axis12df <- tibble::rownames_to_column(taxa.may.axis12df, "Taxa")

write.csv(taxa.may.axis12df, "Data/Emerging/NMDS/may_allvectors_axis12.csv") # save vector scores as csv

# Axis 1 and 3

taxa.may.axis13 <- envfit(nms.may, taxa.may,
                           choices = c(1,3))

taxa.may.axis13df <- data.frame((taxa.may.axis13$vectors)$arrows,
                                 (taxa.may.axis13$vectors)$r,
                                 (taxa.may.axis13$vectors)$pvals)

taxa.may.axis13df <- tibble::rownames_to_column(taxa.may.axis13df, "Taxa")

write.csv(taxa.may.axis13df, "Data/Emerging/NMDS/may_allvectors_axis13.csv") # save vector scores as csv


# Correlated taxa

# Axis 1 and 2

colnames(taxa.may.axis12df)

corrspp.may.axis12 <- taxa.may.axis12df %>% filter(X.taxa.may.axis12.vectors..r > 0.2)
target12.may <- corrspp.may.axis12$Taxa # string of the Family names


axis12.vectors.may <- taxa.may %>% select(all_of(target12.may)) # make a matrix of just those

(nmds.may.vectors.12 <- envfit(nms.may$points, axis12.vectors.may,
                              permutations = 999, choices = c(1,2)))                        


corr.may.vectors.12 <- as.data.frame(nmds.may.vectors.12$vectors$arrows*sqrt(nmds.may.vectors.12$vectors$r)) #scaling vectors
corr.may.vectors.12$Taxa <- rownames(corr.may.vectors.12)


write.csv(corr.may.vectors.12, "Data/Emerging/NMDS/May_emerging_correlated_vector12.csv")

## Axis 1 and 3

colnames(taxa.may.axis13df)

corrspp.may.axis13 <- taxa.may.axis13df %>% filter(X.taxa.may.axis13.vectors..r > 0.2)
target13.may <- corrspp.may.axis13$Taxa # string of the Family names


axis13.vectors.may <- taxa.may %>% select(all_of(target13.may)) # make a matrix of just those

(nmds.may.vectors.13 <- envfit(nms.may$points, axis13.vectors.may,
                              permutations = 999, choices = c(1,3)))                        


corr.may.vectors.13 <- as.data.frame(nmds.may.vectors.13$vectors$arrows*sqrt(nmds.may.vectors.13$vectors$r)) #scaling vectors
corr.may.vectors.13$Taxa <- rownames(corr.may.vectors.13)


write.csv(corr.may.vectors.13, "Data/Emerging/NMDS/May_emerging_correlated_vector13.csv")


# 5 June 2018 -------------------------------------------------------------

invert.June <- inverts.date %>% filter(Date == "05-Jun-18")

write.csv(invert.June, "Data/Emerging/NMDS/inverts_05-Jun-18.csv")



june.data <- read.csv("Data/Emerging/NMDS/june18_zeros_rel.csv")
colnames(june.data)

taxa.june <- june.data %>% select(Araneae:Crambidae)
env.june <- june.data %>% select(Sites:YrCol)

## perMANOVA
(per.june <- adonis2(taxa.june ~ Treatment * Year, data = june.data,
                    permutations = 999, method = "bray"))

#         Df SumOfSqs      R2      F Pr(>F)    
#Treatment  2   1.9636 0.20714 2.8739  0.001 ***
#Residual  22   7.5157 0.79286                  
#Total     24   9.4793 1.00000 

june.b <- vegdist(taxa.june, method = "bray")

trt.june <- factor(env.june$Treatment)


(dispersion <- betadisper(june.b,trt.june))
plot(dispersion)

(adonis.pair(june.b,trt.june, 
             nper = 1000, corr.method = "bonferroni"))

#            combination SumsOfSqs   MeanSqs  F.Model         R2     P.value P.value.corrected
#1  Invaded <-> Treated 1.4088944 1.4088944 4.417348 0.22749492 0.001998002       0.005994006
#2 Invaded <-> Uninvaded 0.5323191 0.5323191 1.474066 0.08947798 0.076923077       0.230769231
#3 Treated <-> Uninvaded 1.0055694 1.0055694 2.914424 0.17230408 0.001998002       0.005994006



## NMDS
set.seed(120) 

nms.june <- metaMDS(taxa.june, distance = "bray", # species data, bray-curtis dissimilarity
                    autotransform = FALSE,  # NMDS will do autotransformations for you
                    k = 3, trymax = 1000)   # k = number of axes
nms.june

#global Multidimensional Scaling using monoMDS

#Data:     taxa.june 
#Distance: bray 

#Dimensions: 3 
#Stress:     0.1581744 
#Stress type 1, weak ties
#Two convergent solutions found after 48 tries
#Scaling: centring, PC rotation, halfchange scalin

layout(matrix(1:2, ncol = 2))
plot(nms.june , main = "Invertebrate NMDS plot"); stressplot(nms.june , main = "Shepard plot")
layout(1)

nms.june$iters #118

nms.june$stress^2   #0.02500559
1-nms.june$stress^2 #0.9749944


scr.june <- as.data.frame(scores(nms.june, display = "sites")) # extract NMDS scores
colnames(scr.june)

env.june$NMDS1 <- scr.june$NMDS1
env.june$NMDS2 <- scr.june$NMDS2
env.june$NMDS3 <- scr.june$NMDS3

write.csv(env.june,"Data/Emerging/NMDS/emerging_june_NMDSscores.csv") # save this as a csv



## Taxa for vectors

# Axis 1 and 2
taxa.june.axis12 <- envfit(nms.june, taxa.june,
                           choices = c(1,2))

taxa.june.axis12df <- data.frame((taxa.june.axis12$vectors)$arrows,
                                 (taxa.june.axis12$vectors)$r,
                                 (taxa.june.axis12$vectors)$pvals)

taxa.june.axis12df <- tibble::rownames_to_column(taxa.june.axis12df, "Taxa")

write.csv(taxa.june.axis12df, "Data/Emerging/NMDS/june_allvectors_axis12.csv") # save vector scores as csv

# Axis 1 and 3

taxa.june.axis13 <- envfit(nms.june, taxa.june,
                           choices = c(1,3))

taxa.june.axis13df <- data.frame((taxa.june.axis13$vectors)$arrows,
                                 (taxa.june.axis13$vectors)$r,
                                 (taxa.june.axis13$vectors)$pvals)

taxa.june.axis13df <- tibble::rownames_to_column(taxa.june.axis13df, "Taxa")

write.csv(taxa.june.axis13df, "Data/Emerging/NMDS/june_allvectors_axis13.csv") # save vector scores as csv



# Correlated taxa

# Axis 1 and 2

colnames(taxa.june.axis12df)

corrspp.june.axis12 <- taxa.june.axis12df %>% filter(X.taxa.june.axis12.vectors..r > 0.2)
target12.june <- corrspp.june.axis12$Taxa # string of the Family names


axis12.vectors.june <- taxa.june %>% select(all_of(target12.june)) # make a matrix of just those

(nmds.june.vectors.12 <- envfit(nms.june$points, axis12.vectors.june,
                              permutations = 999, choices = c(1,2)))                        


corr.june.vectors.12 <- as.data.frame(nmds.june.vectors.12$vectors$arrows*sqrt(nmds.june.vectors.12$vectors$r)) #scaling vectors
corr.june.vectors.12$Taxa <- rownames(corr.june.vectors.12)


write.csv(corr.june.vectors.12, "Data/Emerging/NMDS/june_emerging_correlated_vector12.csv")

## Axis 1 and 3

colnames(taxa.june.axis13df)

corrspp.june.axis13 <- taxa.june.axis13df %>% filter(X.taxa.june.axis13.vectors..r > 0.2)
target13.june <- corrspp.june.axis13$Taxa # string of the Family names


axis13.vectors.june <- taxa.june %>% select(all_of(target13.june)) # make a matrix of just those

(nmds.june.vectors.13 <- envfit(nms.june$points, axis13.vectors.june,
                              permutations = 999, choices = c(1,3)))                        


corr.june.vectors.13 <- as.data.frame(nmds.june.vectors.13$vectors$arrows*sqrt(nmds.june.vectors.13$vectors$r)) #scaling vectors
corr.june.vectors.13$Taxa <- rownames(corr.june.vectors.13)


write.csv(corr.june.vectors.13, "Data/Emerging/NMDS/june_emerging_correlated_vector13.csv")










library(tidyverse)
library(lubridate)

# Emerging inverts 5-June-18 to 23-July-18
emerg.raw <- read.csv("Data/Emerging/Procrustes/emerging_invert_junjuly.csv")

colnames(emerg.raw)

# Sum over the visits 

df2 <- emerg.raw %>% 
  mutate(Date = dmy(Date)) %>% 
  group_by(ID, Treatment) %>%
  summarise_at(vars(Araneae:Hesperiidae), sum, na.rm = TRUE)

count <- emerg.raw %>% 
  mutate(Date = dmy(Date)) %>% 
  group_by(ID, Treatment) %>%
  mutate(N = length(ID))

colnames(count)

count <- count %>% filter(YrCol == "2018Col._2")

n <- count$N


df2$N <- n

# output to remove OW sites and add counts
write.csv(df2, "Data/Emerging/Procrustes/emerging_sum_junjuly.csv")

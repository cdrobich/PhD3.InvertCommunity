library(tidyverse)


env.site <- read.csv("Data/Emerging/site_characteristics.csv")
str(env.site)
colnames(env.site)

df2.en <- env.site %>% 
  group_by(Site) %>%
  summarise_at(vars(Avg_Can:St.Dead), mean, na.rm = TRUE)

write.csv(df2.en, "Data/Emerging/site_characteristics_avg.csv")



veg.site <- read.csv("Data/Emerging/site_vegetation.csv")
str(veg.site)
colnames(veg.site)

df2.veg <- veg.site %>% 
  group_by(Site) %>%
  summarise_at(vars(Litter:ZIZPALUS), mean, na.rm = TRUE)

write.csv(df2.veg, "Data/Emerging/site_vegetation_avg.csv")

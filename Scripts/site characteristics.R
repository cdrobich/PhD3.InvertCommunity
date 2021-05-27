library(tidyverse)


# 2018 sites --------------------------------------------------------------

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



# 2017 sites --------------------------------------------------------------

veg.site17 <- read.csv("Data/Emerging/2017_Site_Characteristics.csv")


str(veg.site17)
colnames(veg.site17)


df7.veg <- veg.site17 %>% 
  group_by(Site) %>%
  summarise_at(vars(Avg_Can:UTRINTER), mean, na.rm = TRUE)

write.csv(df7.veg, "Data/Emerging/site_vegetation_avg_2017.csv")

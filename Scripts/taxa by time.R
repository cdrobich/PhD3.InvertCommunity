# Load libraries ----------------------------------------------------------

library(tidyverse) # for figures, pipes etc.
library(vegan) # for NMDS
library(ggpubr) # for arranging figures
library(patchwork)



invert.data <- read.csv("Data/Emerging/emerging_invertebrates.csv")


invert.data <- invert.data %>% filter(Year == "2018")
colnames(invert.data)

taxa.18 <- invert.data %>% select(Araneae:Crambidae)

sum.tx <- as.data.frame(sapply(taxa.18, function(col) sum(col))) # how many unique values in column

Chironomidae
Ceratopogoni
Hydroptilida
Caenidae
Coenagrionid
Araneae


inverts.time <- read.csv("Data/Emerging/ermerging_time_rares.csv")

inverts.time <- inverts.time %>% filter(Year == "2018")



st_err <- function(x, na.rm=FALSE) {
  if (na.rm) x <- na.omit(x)
  sqrt(var(x)/length(x))
}


chiro <- inverts.time %>% select(ID:Date, Chironomidae)

chiro.sum <- chiro %>% 
  mutate(Date = dmy(Date)) %>% 
  group_by(Treatment, Date) %>% 
  summarise(avg = mean(Chironomidae),
            sd = sd(Chironomidae),
            N = length(Chironomidae),
            str.2 = st_err(Chironomidae)) 


fill = c("Invaded" = "#440C53",
         "Treated" = "#24908C",
         "Remnant" = "#FDE825")

colour = c("Invaded" = "#440C53",
           "Treated" = "#24908C",
           "Remnant" = "#FDE825")

shape = c("Invaded" = 21,
          "Treated" = 24,
          "Remnant" = 22)


ggplot(chiro.sum, aes(x =Date, y = avg,
                      shape = Treatment,
                      group = Treatment,
                      fill = Treatment)) +
  geom_errorbar(aes(ymin = avg - str.2, ymax = avg + str.2), 
                width = 0.5) +
  geom_line(color = "black") +
  geom_point(size = 7,
             stroke = 1) + 
  scale_x_date(date_labels = "%d-%b-%y",
             date_breaks = "1 week") +
  theme_classic() +
  labs(y = "Chironomidae abundance",
       x = "") + 
  theme(panel.border = element_rect(fill = NA)) +
  theme(text = element_text(size = 14),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 14)) +
  scale_fill_manual(values = fill) +
  scale_shape_manual(values = shape) +
  scale_colour_manual(values = colour) 


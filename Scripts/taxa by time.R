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



inverts.time <- read.csv("Data/Emerging/ermerging_time_rares.csv")

inverts.time <- inverts.time %>% filter(Year == "2018")



st_err <- function(x, na.rm=FALSE) {
  if (na.rm) x <- na.omit(x)
  sqrt(var(x)/length(x))
}



fill = c("Invaded" = "#440C53",
         "Treated" = "#24908C",
         "Remnant" = "#FDE825")

colour = c("Invaded" = "#440C53",
           "Treated" = "#24908C",
           "Remnant" = "#FDE825")

shape = c("Invaded" = 21,
          "Treated" = 24,
          "Remnant" = 22)


### Chironomids 

chiro <- inverts.time %>% select(ID:Date, Chironomidae)

chiro.sum <- chiro %>% 
  mutate(Date = dmy(Date)) %>% 
  group_by(Treatment, Date) %>% 
  summarise(avg = mean(Chironomidae),
            sd = sd(Chironomidae),
            N = length(Chironomidae),
            str.2 = st_err(Chironomidae)) 


chiro.plot <- ggplot(chiro.sum, aes(x = Date, y = avg,
                      shape = Treatment,
                      group = Treatment,
                      fill = Treatment)) +
  geom_errorbar(aes(ymin = avg - str.2, ymax = avg + str.2), 
                width = 0.5) +
  geom_line(color = "black") +
  geom_point(size = 7,
             stroke = 1,
             alpha = 0.85) + 
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
  scale_colour_manual(values = colour)+
  theme(axis.text.x = element_text(angle = 45, 
                                   vjust = 0.9,
                                   hjust = 1)) 

chiro.plot

## Coenagrionid
coen <- inverts.time %>% select(ID:Date, Coenagrionidae)

coen.sum <- coen %>% 
  mutate(Date = dmy(Date)) %>% 
  group_by(Treatment, Date) %>% 
  summarise(avg = mean(Coenagrionidae),
            sd = sd(Coenagrionidae),
            N = length(Coenagrionidae),
            str.2 = st_err(Coenagrionidae)) 


coen.plot <- ggplot(coen.sum, aes(x = Date, y = avg,
                      shape = Treatment,
                      group = Treatment,
                      fill = Treatment)) +
  geom_errorbar(aes(ymin = avg - str.2, ymax = avg + str.2), 
                width = 0.5) +
  geom_line(color = "black") +
  geom_point(size = 7,
             stroke = 1,
             alpha = 0.85) + 
  scale_x_date(date_labels = "%d-%b-%y",
               date_breaks = "1 week") +
  theme_classic() +
  labs(y = "Coenagrionidae abundance",
       x = "") + 
  theme(panel.border = element_rect(fill = NA)) +
  theme(text = element_text(size = 14),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 14)) +
  scale_fill_manual(values = fill) +
  scale_shape_manual(values = shape) +
  scale_colour_manual(values = colour) +
  theme(axis.text.x = element_text(angle = 45, 
                                   vjust = 0.9,
                                   hjust = 1))

coen.plot


## Araneae

spider <- inverts.time %>% select(ID:Date, Araneae)

spid.sum <- spider %>% 
  mutate(Date = dmy(Date)) %>% 
  group_by(Treatment, Date) %>% 
  summarise(avg = mean(Araneae),
            sd = sd(Araneae),
            N = length(Araneae),
            str.2 = st_err(Araneae)) 


spider.plot <- ggplot(spid.sum, aes(x = Date, y = avg,
                      shape = Treatment,
                      group = Treatment,
                      fill = Treatment)) +
  geom_errorbar(aes(ymin = avg - str.2, ymax = avg + str.2), 
                width = 0.5) +
  geom_line(color = "black") +
  geom_point(size = 7,
             stroke = 1,
             alpha = 0.85) + 
  scale_x_date(date_labels = "%d-%b-%y",
               date_breaks = "1 week") +
  theme_classic() +
  labs(y = "Araneae abundance",
       x = "") + 
  theme(panel.border = element_rect(fill = NA)) +
  theme(text = element_text(size = 14),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 14)) +
  scale_fill_manual(values = fill) +
  scale_shape_manual(values = shape) +
  scale_colour_manual(values = colour) +
  theme(axis.text.x = element_text(angle = 45, 
                                   vjust = 0.9,
                                   hjust = 1))

spider.plot

## Caddisfly 

lepto <- inverts.time %>% select(ID:Date, Leptoceridae)

lepto.sum <- lepto %>% 
  mutate(Date = dmy(Date)) %>% 
  group_by(Treatment, Date) %>% 
  summarise(avg = mean(Leptoceridae),
            sd = sd(Leptoceridae),
            N = length(Leptoceridae),
            str.2 = st_err(Leptoceridae)) 


lepto.plot <- ggplot(lepto.sum, aes(x = Date, y = avg,
                                    shape = Treatment,
                                    group = Treatment,
                                    fill = Treatment)) +
  geom_errorbar(aes(ymin = avg - str.2, ymax = avg + str.2), 
                width = 0.5) +
  geom_line(color = "black") +
  geom_point(size = 7,
             stroke = 1,
             alpha = 0.85) + 
  scale_x_date(date_labels = "%d-%b-%y",
               date_breaks = "1 week") +
  theme_classic() +
  labs(y = "Leptoceridae abundance",
       x = "") + 
  theme(panel.border = element_rect(fill = NA)) +
  theme(text = element_text(size = 14),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 14)) +
  scale_fill_manual(values = fill) +
  scale_shape_manual(values = shape) +
  scale_colour_manual(values = colour) +
  theme(axis.text.x = element_text(angle = 45, 
                                   vjust = 0.9,
                                   hjust = 1)) 

lepto.plot

## Culicidae

cucli <- inverts.time %>% select(ID:Date, Culicidae)

cucli.sum <- cucli %>% 
  mutate(Date = dmy(Date)) %>% 
  group_by(Treatment, Date) %>% 
  summarise(avg = mean(Culicidae),
            sd = sd(Culicidae),
            N = length(Culicidae),
            str.2 = st_err(Culicidae)) 


cucli.plot <- ggplot(cucli.sum, aes(x = Date, y = avg,
                                    shape = Treatment,
                                    group = Treatment,
                                    fill = Treatment)) +
  geom_errorbar(aes(ymin = avg - str.2, ymax = avg + str.2), 
                width = 0.5) +
  geom_line(color = "black") +
  geom_point(size = 7,
             stroke = 1,
             alpha = 0.85) + 
  scale_x_date(date_labels = "%d-%b-%y",
               date_breaks = "1 week") +
  theme_classic() +
  labs(y = "Culicidae abundance",
       x = "") + 
  theme(panel.border = element_rect(fill = NA)) +
  theme(text = element_text(size = 14),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 14)) +
  scale_fill_manual(values = fill) +
  scale_shape_manual(values = shape) +
  scale_colour_manual(values = colour) +
  theme(axis.text.x = element_text(angle = 45, 
                                   vjust = 0.9,
                                   hjust = 1))

cucli.plot


## Hydroptilida

hydro <- inverts.time %>% select(ID:Date, Hydroptilidae)

hydro.sum <- hydro %>% 
  mutate(Date = dmy(Date)) %>% 
  group_by(Treatment, Date) %>% 
  summarise(avg = mean(Hydroptilidae),
            sd = sd(Hydroptilidae),
            N = length(Hydroptilidae),
            str.2 = st_err(Hydroptilidae)) 


hydro.plot <- ggplot(hydro.sum, aes(x = Date, y = avg,
                                    shape = Treatment,
                                    group = Treatment,
                                    fill = Treatment)) +
  geom_errorbar(aes(ymin = avg - str.2, ymax = avg + str.2), 
                width = 0.5) +
  geom_line(color = "black") +
  geom_point(size = 7,
             stroke = 1,
             alpha = 0.85) + 
  scale_x_date(date_labels = "%d-%b-%y",
               date_breaks = "1 week") +
  theme_classic() +
  labs(y = "Hydroptilidae abundance",
       x = "") + 
  theme(panel.border = element_rect(fill = NA)) +
  theme(text = element_text(size = 14),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 14)) +
  scale_fill_manual(values = fill) +
  scale_shape_manual(values = shape) +
  scale_colour_manual(values = colour) +
  theme(axis.text.x = element_text(angle = 45, 
                                   vjust = 0.9,
                                   hjust = 1))

hydro.plot


taxa.panel <- ggarrange(chiro.plot, coen.plot, spider.plot,
          lepto.plot, cucli.plot, hydro.plot,
          common.legend = TRUE,
          legend = "right")

ggsave("Figures/taxa_panel_time.jpeg",
       taxa.panel)


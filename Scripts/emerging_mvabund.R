
# Library -----------------------------------------------------------------

library(tidyverse)
library(vegan)
library(mvabund)
library(patchwork)
library(lubridate)
library(plotrix)
library(patchwork)

# Data --------------------------------------------------------------------
inverts <- read.csv("Data/Emerging/emerging_time_long.csv")

inverts.time <- read.csv("Data/Emerging/ermerging_time_rares.csv")


## For line plots 

inverts17 <- inverts %>% 
  filter(Year == "2017") %>% 
  mutate(Date = dmy(Date)) %>% 
  filter(count > 0)


inverts17s <- inverts17 %>% 
  group_by(Date, Time, Treatment) %>% 
  summarise(mean = mean(count),
            std = sd(count),
            str = std.error(count))

inverts18 <- inverts %>% 
  filter(Year == "2018") %>% 
  mutate(Date = dmy(Date)) %>% 
  filter(count > 0)

inverts18s <- inverts18 %>% 
  group_by(Date, Time, Treatment) %>% 
  summarise(mean = mean(count),
            std = sd(count),
            str = std.error(count))


## Chironomid line plots

chiro17s <- inverts17 %>% 
  filter(Family == "Chironomidae") %>% 
  group_by(Date, Time, Treatment) %>% 
  summarise(mean = mean(count),
            std = sd(count),
            str = std.error(count))


chiro18s <- inverts18 %>% 
  filter(Family == "Chironomidae") %>% 
  group_by(Date, Time, Treatment) %>% 
  summarise(mean = mean(count),
            std = sd(count),
            str = std.error(count))
  

# multivariate GLM --------------------------------------------------------

inv.2017 <- inverts.time %>% 
  filter(Year == "2017") %>% 
  mutate(Date = dmy(Date))

inv.2018 <- inverts.time %>% 
  filter(Year == "2018") %>% 
  mutate(Date = dmy(Date))


# 2017 mvabund ------------------------------------------------------------

fam.17 <- inv.2017 %>% select(Araneae:Crambidae)
env.17 <- inv.2017 %>% select(ID:Time)

inv.17.mv <- mvabund(fam.17)

inv7.mod <- manyglm(inv.17.mv ~ Treatment * Time,
                    data = env.17, family = "negativebinomial")


plot(inv7.mod) # looks good

output <- anova(inv7.mod, p.uni = "adjusted")

#Time elapsed: 0 hr 6 min 52 sec
#Analysis of Deviance Table

#Model: inv.17.mv ~ Treatment * Time

#Multivariate test:
#                 Res.Df Df.diff Dev    Pr(>Dev)    
#(Intercept)       102                           
#Treatment         100       2 487.7    0.001 ***
#Time               99       1 225.1    0.001 ***
#Treatment:Time     97       2 158.9    0.001 ***

mv.output <- as.data.frame(output$uni.p)

# 2018 mvabund ------------------------------------------------------------

fam.18 <- inv.2018 %>% select(Araneae:Crambidae)
env.18 <- inv.2018 %>% select(ID:Time)

inv.18.mv <- mvabund(fam.18)

inv8.mod <- manyglm(inv.18.mv ~ Treatment * Time,
                    data = env.18, family = "negativebinomial")

plot(inv8.mod) #looks good

output8 <- anova(inv8.mod, p.uni = "adjusted")

mv.output8 <- as.data.frame(output8$uni.p)


#Multivariate test:
#                  Res.Df Df.diff   Dev Pr(>Dev)    
#(Intercept)       155                           
#Treatment         153       2 679.1    0.001 ***
#Time              152       1 551.6    0.001 ***
#Treatment:Time    150       2 214.9    0.001 ***


# Figures -----------------------------------------------------------------
fill = c("Invaded" = "#440C53",
         "Treated" = "#24908C",
         "Remnant" = "#3A518B")

colour = c("Invaded" = "#440C53",
           "Treated" = "#24908C",
           "Remnant" = "#3A518B")

shape = c("Invaded" = 21,
          "Treated" = 24,
          "Remnant" = 22)



# geom_smooth figures -----------------------------------------------------

str(inverts)
unique(inverts$Treatment)

chiro17 <- ggplot(inv.2017, aes(x = Date, y = Chironomidae + 1, 
                    fill = Treatment,
                    shape = Treatment)) +
  geom_jitter(size = 3, alpha = 0.8) +
  scale_y_log10() +
  scale_x_date(date_labels = "%d-%b-%y",
               date_breaks = "1 week") +
  labs(y = (expression(paste("Chironomidae density per 1 ","", m^2))),
       x = ' ') +
  theme_bw() +
  geom_smooth(aes(colour = Treatment)) +
  scale_fill_manual(values = fill) +
  scale_colour_manual(values = colour) +
  scale_shape_manual(values = shape) +
  theme(axis.text = element_text(size = 12),
        strip.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 11),
        legend.position = "none") +
  ggtitle("2017")


chiro18 <- ggplot(inv.2018, aes(x = Date, y = Chironomidae + 1, 
                     fill = Treatment,
                     shape = Treatment)) +
  geom_jitter(size = 3, alpha = 0.8) +
  scale_y_log10() +
  scale_x_date(date_labels = "%d-%b-%y",
               date_breaks = "2 week") +
  labs(y = (expression(paste("Chironomidae density per 1 ","", m^2))),
       x = ' ') +
  theme_bw() +
  geom_smooth(aes(colour = Treatment)) +
  scale_fill_manual(values = fill) +
  scale_colour_manual(values = colour) +
  scale_shape_manual(values = shape) +
  theme(axis.text = element_text(size = 12),
        strip.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 11),
        legend.position = "none") +
  ggtitle("2018") 
theme(axis.title.y = element_text(colour = "white"),
      legend.title = element_blank(),
      axis.text.y = element_text(colour = "white"),
      axis.text.x = element_text(colour = "white"),
      axis.line = element_line(colour = "white")) +
  theme(panel.background = element_rect(
    fill = 'black'),
    plot.background = element_rect(fill = 'black'))

chiros <- chiro17 + chiro18

ggsave("Figures/Chironomidae_bothyears.TIFF",
       chiros)

### All counts

total7 <- ggplot(inverts17, aes(x = Date, y = count, 
                     fill = Treatment,
                     shape = Treatment)) +
  geom_jitter(size = 3, alpha = 0.8) +
  geom_smooth(aes(colour = Treatment)) +
  scale_y_log10() +
  scale_x_date(date_labels = "%d-%b-%y",
               date_breaks = "2 week") +
  labs(y = (expression(paste("Density per 1 ","", m^2))),
       x = ' ') +
  theme_bw() +
  scale_fill_manual(values = fill) +
  scale_colour_manual(values = colour) +
  scale_shape_manual(values = shape) +
  theme(axis.text = element_text(size = 12),
        strip.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 11),
        legend.position = "none") +
  ggtitle("2017")


total8 <- ggplot(inverts18, aes(x = Date, y = count, 
                      fill = Treatment,
                      shape = Treatment)) +
  geom_jitter(size = 3, alpha = 0.8) +
  geom_smooth(aes(colour = Treatment)) +
  scale_y_log10() +
  scale_x_date(date_labels = "%d-%b-%y",
               date_breaks = "2 week") +
  labs(y = (expression(paste("Density per 1 ","", m^2))),
       x = ' ') +
  theme_bw() +
  scale_fill_manual(values = fill) +
  scale_colour_manual(values = colour) +
  scale_shape_manual(values = shape) +
  theme(axis.text = element_text(size = 12),
        strip.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 11)) +
  ggtitle("2018")

total <- total7 + total8

ggsave("Figures/total_smooth_bothyears.TIFF",
       total)



panel <- total7 + total8 + chiro17 + chiro18 +
  plot_annotation(tag_levels = "A")

ggsave("Figures/total_chiro_bothyears.TIFF",
       panel)



# Line plots --------------------------------------------------------------


in17 <- ggplot(inverts17s, aes(x = Date, y = mean,
                               shape = Treatment,
                               group = Treatment,
                               fill = Treatment)) +
  geom_errorbar(aes(ymin = mean - str, ymax = mean + str,
                    colour = Treatment), 
                width = 0.5) +
  geom_line(aes(colour = Treatment)) +
  geom_point(size = 4, stroke = 1.5) +
  scale_x_date(date_labels = "%d-%b-%y",
               date_breaks = "1 week") +
  labs(y = (expression(paste("Average density per 1 ","", m^2))),
       x = ' ') +
  theme_bw() +
  scale_fill_manual(values = fill) +
  scale_colour_manual(values = colour) +
  scale_shape_manual(values = shape) +
  theme(axis.text = element_text(size = 12),
        strip.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        legend.text = element_text(size = 11),
        legend.position = "none") 



in18 <- ggplot(inverts18s, aes(x = Date, y = mean,
                               shape = Treatment,
                               group = Treatment,
                               fill = Treatment)) +
  geom_errorbar(aes(ymin = mean - str, ymax = mean + str,
                    colour = Treatment), 
                width = 0.5) +
  geom_line(aes(colour = Treatment)) +
  ylim(0, 300) +
  geom_point(size = 4, stroke = 1.5) +
  scale_x_date(date_labels = "%d-%b-%y",
               date_breaks = "2 week") +
  labs(y = (expression(paste("Average density per 1 ","", m^2))),
       x = ' ') +
  theme_bw() +
  scale_fill_manual(values = fill) +
  scale_colour_manual(values = colour) +
  scale_shape_manual(values = shape) +
  theme(axis.text = element_text(size = 12),
        strip.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        legend.text = element_text(size = 11)) +
  labs(colour = "Habitat",
       fill = "Habitat",
       shape = "Habitat")

inv178 <- in17 + in18 +
  plot_annotation(tag_levels = "A")

ggsave("Figures/invert_density_bothyears.TIFF",
       inv178)



chirol17 <- ggplot(chiro17s, aes(x = Date, y = mean,
                       shape = Treatment,
                       group = Treatment,
                       fill = Treatment)) +
  geom_errorbar(aes(ymin = mean - str, ymax = mean + str,
                    colour = Treatment), 
                width = 0.5) +
  geom_line(aes(colour = Treatment)) +
  geom_point(size = 4, stroke = 1.5) +
  scale_x_date(date_labels = "%d-%b-%y",
               date_breaks = "1 week") +
  labs(y = (expression(paste("Chironomidae density per 1 ","", m^2))),
       x = ' ') +
  theme_bw() +
  scale_fill_manual(values = fill) +
  scale_colour_manual(values = colour) +
  scale_shape_manual(values = shape) +
  theme(axis.text = element_text(size = 12),
        strip.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        legend.text = element_text(size = 11),
        legend.position = "none") 


chirol18 <- ggplot(chiro18s, aes(x = Date, y = mean,
                     shape = Treatment,
                     group = Treatment,
                     fill = Treatment)) +
  geom_errorbar(aes(ymin = mean - str, ymax = mean + str,
                    colour = Treatment), 
                width = 0.5) +
  geom_line(aes(colour = Treatment)) +
  geom_point(size = 4, stroke = 1.5) +
  scale_x_date(date_labels = "%d-%b-%y",
               date_breaks = "2 week") +
  labs(y = (expression(paste("Chironomidae density per 1 ","", m^2))),
       x = ' ') +
  theme_bw() +
  scale_fill_manual(values = fill) +
  scale_colour_manual(values = colour) +
  scale_shape_manual(values = shape) +
  theme(axis.text = element_text(size = 12),
        strip.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        legend.text = element_text(size = 11),
        legend.position = "none") 




line.panel <- in17 + in18 +
  chirol17 + chirol18 +
  plot_annotation(tag_levels = "A")

ggsave("Figures/lineplots_panel_20172018.tiff",
       line.panel)


# National Phrag figures (black background) -------------------------------



ggplot(inv.2018, aes(x = Date, y = Chironomidae + 1, 
                     fill = Treatment,
                     shape = Treatment)) +
  scale_y_log10() +
  scale_x_date(date_labels = "%d-%b-%y",
               date_breaks = "2 week") +
  labs(y = (expression(paste("Chironomidae density per 1 ","", m^2))),
       x = ' ') +
  theme_classic() +
  geom_smooth(aes(colour = Treatment,
                  size = 0.8)) +  
  scale_fill_manual(values = fill) +
  scale_colour_manual(values = colour) +
  scale_shape_manual(values = shape) +
  theme(axis.text = element_text(size = 12),
        strip.text = element_text(size = 14),
        axis.title = element_text(size = 14),
        legend.text = element_text(size = 11),
        legend.position = "none") +
theme(axis.title.y = element_text(colour = "white"),
      legend.title = element_blank(),
      axis.text.y = element_text(colour = "white"),
      axis.text.x = element_text(colour = "white"),
      axis.line = element_line(colour = "white")) +
  theme(panel.background = element_rect(
    fill = 'black'),
    plot.background = element_rect(fill = 'black'))

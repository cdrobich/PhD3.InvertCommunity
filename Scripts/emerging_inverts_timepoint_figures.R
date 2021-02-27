
# Load packages -----------------------------------------------------------

library(tidyverse) # for figures, pipes etc.
library(vegan) # for NMDS
library(ggpubr) # for arranging figures

library(ggrepel) # labels on nmds
library(lubridate)

library(viridis) # colours

# NMDS Collection 19-Jun-17 and 16-Jun-18 ---------------------------------------------------------------

# NMDS Emerging Invert Figure Collection 1 
col1.scores <- read.csv("Data/Emerging/NMDS/Col.1/emerging_inverts_col1_NMDSscores.csv")
col1.scores$Year <- as.factor(col1.scores$Year)

col1.axis12 <- read.csv("Data/Emerging/NMDS/Col.1/NMDS_emerg_col1_vector12_final.csv")
col1.axis13 <- read.csv("Data/Emerging/NMDS/Col.1/NMDS_emerg_col1_vector13_final.csv")


## NMDS Axis 1, 2 


invert.12.c1 <- ggplot(data = col1.scores,
                       aes(x = NMDS1, y = NMDS2)) +
  geom_point(data = col1.scores, 
             aes(x = NMDS1, y = NMDS2, 
                 fill = Treatment, shape = Year),
             size = 5, stroke = 1.5,
             alpha = 0.8) + # sites as points
  stat_ellipse(data = col1.scores, 
               aes(x = NMDS1,y = NMDS2,
                   linetype = Treatment, colour = Treatment), 
               size = 1, level = 0.9) + 
  geom_segment(data = col1.axis12, 
               aes(x = 0, xend = MDS1, y = 0, yend = MDS2), # adding in the vectors, c
               arrow = arrow(length = unit(0.5, "cm")), colour = "black") + # can add in geom_label or geom_text for labels
  theme_minimal() + # no background
  theme(panel.border = element_rect(fill = NA)) + # full square around figure
  xlab("NMDS 1") +
  ylab("NMDS 2") +
  #ylim(-1, 1.5) +
  #xlim(-1.45, 1) +
  #theme(legend.position = "none") +
  geom_label_repel(data = col1.axis12, 
                   aes(x = MDS1, y = MDS2, label = Taxa),
                   color="black",
                   size = 6)  +
  scale_fill_viridis(discrete = TRUE) +
  scale_colour_viridis(discrete = TRUE) +
  scale_shape_manual(values = c(21, 24, 22)) +
  theme(legend.position = "none") +
  ylim(-0.8, 0.8) +
  xlim(-0.8, 0.8)

invert.12.c1

## NMDS Axis 1, 3

invert.13.c1 <- ggplot(data = col1.scores,
                       aes(x = NMDS1, y = NMDS3)) +
  geom_point(data = col1.scores, 
             aes(x = NMDS1, y = NMDS3, 
                 fill = Treatment, shape = Year),
             size = 5, stroke = 1.5,
             alpha = 0.7) + # sites as points
  stat_ellipse(data = col1.scores, 
               aes(x = NMDS1,y = NMDS3,
                   linetype = Treatment, colour = Treatment), 
               size = 1, level = 0.9) + 
  geom_segment(data = col1.axis13, 
               aes(x = 0, xend = MDS1, y = 0, yend = MDS3), # adding in the vectors, c
               arrow = arrow(length = unit(0.5, "cm")), colour = "black") + # can add in geom_label or geom_text for labels
  theme_minimal() + # no background
  theme(panel.border = element_rect(fill = NA)) + # full square around figure
  xlab("NMDS 1") +
  ylab("NMDS 3") +
  #ylim(-1, 1.5) +
  #xlim(-1.45, 1) +
  #theme(legend.position = "none") +
  geom_label_repel(data = col1.axis13, 
                   aes(x = MDS1, y = MDS3, label = Taxa),
                   color="black",
                   size = 6) +
  scale_fill_viridis(discrete = TRUE) +
  scale_colour_viridis(discrete = TRUE) +
  scale_shape_manual(values = c(21, 24, 22)) +
  theme(legend.position = "none") +
  ylim(-0.8, 0.8) +
  xlim(-0.8, 0.8)

invert.13.c1

(NMS.emerging.panel.c1 <- ggarrange(invert.12.c1, invert.13.c1,
                                    common.legend = TRUE,
                                    legend = "none",
                                    labels = c("C", ""),
                                    align = "hv"))

NMDS.col1 <- annotate_figure(NMS.emerging.panel.c1,
                             top = text_grob("Collection 19-Jun-2017 and 16-Jun-2018"))

ggsave("Figures/NMDS_emerging_19Jun17_16Jun18.jpeg", NMDS.col1)



# NMDS Collection 28-Jun-17 and 25-JUN-18 ---------------------------------


## Actual figure
nmds.col2.scores <- read.csv("Data/Emerging/NMDS/Col.2/emerging_col2_NMDSscores.csv")
nmds.col2.scores$Year <- as.factor(nmds.col2.scores$Year)

col2.axis12 <- read.csv("Data/Emerging/NMDS/Col.2/emerging_correlated_vector12.csv")
col2.axis13 <- read.csv("Data/Emerging/NMDS/Col.2/emerging_correlated_vector13.csv")



invert.12.c2 <- ggplot(data = nmds.col2.scores,
                       aes(x = NMDS1, y = NMDS2)) +
  geom_point(data = nmds.col2.scores, 
             aes(x = NMDS1, y = NMDS2, 
                 fill = Treatment, shape = Year),
             size = 4, stroke = 1.5,
             alpha = 0.7) + # sites as points
  stat_ellipse(data = nmds.col2.scores, 
               aes(x = NMDS1,y = NMDS2,
                   linetype = Treatment, colour = Treatment), 
               size = 1, level = 0.9) + 
  geom_segment(data = col2.axis12, 
               aes(x = 0, xend = MDS1, y = 0, yend = MDS2), # adding in the vectors, c
               arrow = arrow(length = unit(0.5, "cm")), colour = "black") + # can add in geom_label or geom_text for labels
  theme_minimal() + # no background
  theme(panel.border = element_rect(fill = NA)) + # full square around figure
  xlab("NMDS 1") +
  ylab("NMDS 2") +
  #ylim(-1, 1.5) +
  #xlim(-1.45, 1) +
  #theme(legend.position = "none") +
  geom_label_repel(data = col2.axis12, 
                   aes(x = MDS1, y = MDS2, label = Taxa),
                   color="black",
                   size = 5) +
  scale_fill_viridis(discrete = TRUE) +
  scale_colour_viridis(discrete = TRUE) +
  scale_shape_manual(values = c(21, 24, 22)) +
  theme(legend.position = "none") +
  ylim(-0.8, 0.8) +
  xlim(-0.8, 0.8)

invert.12.c2


invert.13.c2 <- ggplot(data = nmds.col2.scores,
                       aes(x = NMDS1, y = NMDS3)) +
  geom_point(data = nmds.col2.scores, 
             aes(x = NMDS1, y = NMDS3, 
                 fill = Treatment, shape = Year),
             size = 4, stroke = 1.5,
             alpha = 0.7) + # sites as points
  stat_ellipse(data = nmds.col2.scores, 
               aes(x = NMDS1,y = NMDS3,
                   linetype = Treatment, colour = Treatment), 
               size = 1, level = 0.9) + 
  geom_segment(data = col2.axis13, 
               aes(x = 0, xend = MDS1, y = 0, yend = MDS3), # adding in the vectors, c
               arrow = arrow(length = unit(0.5, "cm")), colour = "black") + # can add in geom_label or geom_text for labels
  theme_minimal() + # no background
  theme(panel.border = element_rect(fill = NA)) + # full square around figure
  xlab("NMDS 1") +
  ylab("NMDS 3") +
  #ylim(-1, 1.5) +
  #xlim(-1.45, 1) +
  #theme(legend.position = "none") +
  geom_label_repel(data = col2.axis13, 
                   aes(x = MDS1, y = MDS3, label = Taxa),
                   color="black",
                   size = 5) +
  scale_fill_viridis(discrete = TRUE) +
  scale_colour_viridis(discrete = TRUE) +
  scale_shape_manual(values = c(21, 24, 22)) +
  theme(legend.position = "none") +
  ylim(-0.8, 0.8) +
  xlim(-0.8, 0.8)

invert.13.c2


(NMS.emerging.panel.c2 <- ggarrange(invert.12.c2, invert.13.c2,
                                    common.legend = TRUE,
                                    legend = "none",
                                    labels = c("D", ""),
                                    align = "hv"))

(NMDS.col2 <- annotate_figure(NMS.emerging.panel.c2,
                              top = text_grob("Collection 28-Jun-17 and 25-JUN-18")))


ggsave("Figures/NMDS_emerging_28Jun17_15Jun18.jpeg", NMDS.col2)




# NMDS Collection 08-Jul-17 and 04-Jul-18  --------------------------------


## Figure 
nmds.col3.scores <- read.csv("Data/Emerging/NMDS/Col.3/emerging_col3_NMDSscores.csv")
nmds.col3.scores$Year <- as.factor(nmds.col3.scores$Year)

col3.axis12 <- read.csv("Data/Emerging/NMDS/Col.3/emerging_correlated_vector12.csv")
col3.axis13 <- read.csv("Data/Emerging/NMDS/Col.3/emerging_correlated_vector13.csv")



invert.12.c3 <- ggplot(data = nmds.col3.scores,
                       aes(x = NMDS1, y = NMDS2)) +
  geom_point(data = nmds.col3.scores, 
             aes(x = NMDS1, y = NMDS2, 
                 fill = Treatment, shape = Year),
             size = 5, stroke = 1.5,
             alpha = 0.7) + # sites as points
  stat_ellipse(data = nmds.col3.scores, 
               aes(x = NMDS1,y = NMDS2,
                   linetype = Treatment, colour = Treatment), 
               size = 1, level = 0.9) + 
  geom_segment(data = col3.axis12, 
               aes(x = 0, xend = MDS1, y = 0, yend = MDS2), # adding in the vectors, c
               arrow = arrow(length = unit(0.5, "cm")), colour = "black") + # can add in geom_label or geom_text for labels
  theme_minimal() + # no background
  theme(panel.border = element_rect(fill = NA)) + # full square around figure
  xlab("NMDS 1") +
  ylab("NMDS 2") +
  #ylim(-1, 1.5) +
  #xlim(-1.45, 1) +
  #theme(legend.position = "none") +
  geom_label_repel(data = col3.axis12, 
                   aes(x = MDS1, y = MDS2, label = Taxa),
                   color="black",
                   size = 5) +
  scale_fill_viridis(discrete = TRUE) +
  scale_colour_viridis(discrete = TRUE) +
  scale_shape_manual(values = c(21, 24, 22)) +
  theme(legend.position = "none") +
  ylim(-0.8, 0.8) +
  xlim(-0.8, 0.8)

invert.12.c3



invert.13.c3 <- ggplot(data = nmds.col3.scores,
                       aes(x = NMDS1, y = NMDS3)) +
  geom_point(data = nmds.col3.scores, 
             aes(x = NMDS1, y = NMDS3, 
                 fill = Treatment, shape = Year),
             size = 5, stroke = 1.5,
             alpha = 0.7) + # sites as points
  stat_ellipse(data = nmds.col3.scores, 
               aes(x = NMDS1,y = NMDS3,
                   linetype = Treatment, colour = Treatment), 
               size = 1, level = 0.9) + 
  geom_segment(data = col3.axis13, 
               aes(x = 0, xend = MDS1, y = 0, yend = MDS3), # adding in the vectors, c
               arrow = arrow(length = unit(0.5, "cm")), colour = "black") + # can add in geom_label or geom_text for labels
  theme_minimal() + # no background
  theme(panel.border = element_rect(fill = NA)) + # full square around figure
  xlab("NMDS 1") +
  ylab("NMDS 3") +
  #ylim(-1, 1.5) +
  #xlim(-1.45, 1) +
  #theme(legend.position = "none") +
  geom_label_repel(data = col3.axis13, 
                   aes(x = MDS1, y = MDS3, label = Taxa),
                   color="black",
                   size = 5) +
  scale_fill_viridis(discrete = TRUE) +
  scale_colour_viridis(discrete = TRUE) +
  scale_shape_manual(values = c(21, 24, 22)) +
  theme(legend.position = "none") +
  ylim(-0.8, 0.8) +
  xlim(-0.8, 0.8)


invert.13.c3


(NMS.emerging.panel.c3 <- ggarrange(invert.12.c3, invert.13.c3,
                                    common.legend = TRUE,
                                    legend = "none",
                                    labels = c("E", ""),
                                    align = "hv"))

NMDS.col3 <- annotate_figure(NMS.emerging.panel.c3,
                             top = text_grob("Collection 08-Jul-17 and 04-JUL-18"))


ggsave("Figures/NMDS_emerging_08Jul17_04Jul18.jpeg", NMDS.col3)



# NMDS Collection 21-Jul-17 and 23-Jul-18 ---------------------------------


## NMDS figure

nmds.col4.scores <- read.csv("Data/Emerging/NMDS/Col.4/emerging_col4_NMDSscores.csv")
nmds.col4.scores$Year <- as.factor(nmds.col4.scores$Year)

col4.axis12 <- read.csv("Data/Emerging/NMDS/Col.4/emerging_correlated_vector12.csv")
col4.axis13 <- read.csv("Data/Emerging/NMDS/Col.4/emerging_correlated_vector13.csv")



invert.12.c4 <- ggplot(data = nmds.col4.scores,
                       aes(x = NMDS1, y = NMDS2)) +
  geom_point(data = nmds.col4.scores, 
             aes(x = NMDS1, y = NMDS2, 
                 fill = Treatment, shape = Year),
             size = 5, stroke = 1.5,
             alpha = 0.7) + # sites as points
  stat_ellipse(data = nmds.col4.scores, 
               aes(x = NMDS1,y = NMDS2,
                   linetype = Treatment, colour = Treatment), 
               size = 1, level = 0.9) + 
  geom_segment(data = col4.axis12, 
               aes(x = 0, xend = MDS1, y = 0, yend = MDS2), # adding in the vectors, c
               arrow = arrow(length = unit(0.5, "cm")), colour = "black") + # can add in geom_label or geom_text for labels
  theme_minimal() + # no background
  theme(panel.border = element_rect(fill = NA)) + # full square around figure
  xlab("NMDS 1") +
  ylab("NMDS 2") +
  #ylim(-1, 1.5) +
  #xlim(-1.45, 1) +
  #theme(legend.position = "none") +
  geom_label_repel(data = col4.axis12, 
                   aes(x = MDS1, y = MDS2, label = Taxa),
                   color="black",
                   size = 5) +
  scale_fill_viridis(discrete = TRUE) +
  scale_colour_viridis(discrete = TRUE) +
  scale_shape_manual(values = c(21, 24)) +
  theme(legend.position = "none") +
  ylim(-0.6, 0.6) +
  xlim(-0.6, 0.6)

invert.12.c4


invert.13.c4 <- ggplot(data = nmds.col4.scores,
                       aes(x = NMDS1, y = NMDS3)) +
  geom_point(data = nmds.col4.scores, 
             aes(x = NMDS1, y = NMDS3, 
                 fill = Treatment, shape = Year),
             size = 4, stroke = 1.5,
             alpha = 0.7) + # sites as points
  stat_ellipse(data = nmds.col4.scores, 
               aes(x = NMDS1,y = NMDS3,
                   linetype = Treatment, colour = Treatment), 
               size = 1, level = 0.9) + 
  geom_segment(data = col4.axis13, 
               aes(x = 0, xend = MDS1, y = 0, yend = MDS3), # adding in the vectors, c
               arrow = arrow(length = unit(0.5, "cm")), colour = "black") + # can add in geom_label or geom_text for labels
  theme_minimal() + # no background
  theme(panel.border = element_rect(fill = NA)) + # full square around figure
  xlab("NMDS 1") +
  ylab("NMDS 3") +
  #ylim(-1, 1.5) +
  #xlim(-1.45, 1) +
  #theme(legend.position = "none") +
  geom_label_repel(data = col4.axis13, 
                   aes(x = MDS1, y = MDS3, label = Taxa),
                   color="black",
                   size = 5) +
  scale_fill_viridis(discrete = TRUE) +
  scale_colour_viridis(discrete = TRUE) +
  scale_shape_manual(values = c(21, 24)) +
  theme(legend.position = "none") +
  ylim(-0.6, 0.6) +
  xlim(-0.6, 0.6)

invert.13.c4



(NMS.emerging.panel.c4 <- ggarrange(invert.12.c4, invert.13.c4,
                                    labels = c("F", ""),
                                    align = "hv"))

NMDS.col4 <- annotate_figure(NMS.emerging.panel.c4,
                             top = text_grob("Collection 21-Jul-17 and 23-JUL-18"))


ggsave("Figures/NMDS_emerging_21Jul17_23Jul18.jpeg", NMDS.col3)


# May 20 2018 -------------------------------------------------------------


## Actual figure
nmds.may.scores <- read.csv("Data/Emerging/NMDS/emerging_May18_NMDSscores.csv")
nmds.may.scores$Year <- as.factor(nmds.may.scores$Year)

may.axis12 <- read.csv("Data/Emerging/NMDS/May_emerging_correlated_vector12.csv")
may.axis13 <- read.csv("Data/Emerging/NMDS/May_emerging_correlated_vector13.csv")



invert.12.may <- ggplot(data = nmds.may.scores,
                        aes(x = NMDS1, y = NMDS2)) +
  geom_point(data = nmds.may.scores, 
             aes(x = NMDS1, y = NMDS2, 
                 fill = Treatment, shape = Year),
             size = 5, stroke = 1.5,
             alpha = 0.7) + # sites as points
  stat_ellipse(data = nmds.may.scores, 
               aes(x = NMDS1,y = NMDS2,
                   linetype = Treatment, colour = Treatment), 
               size = 1, level = 0.9) + 
  geom_segment(data = may.axis12, 
               aes(x = 0, xend = MDS1, y = 0, yend = MDS2), # adding in the vectors, c
               arrow = arrow(length = unit(0.5, "cm")), colour = "black") + # can add in geom_label or geom_text for labels
  theme_minimal() + # no background
  theme(panel.border = element_rect(fill = NA)) + # full square around figure
  xlab("NMDS 1") +
  ylab("NMDS 2") +
  #ylim(-1, 1.5) +
  #xlim(-1.45, 1) +
  #theme(legend.position = "none") +
  geom_label_repel(data = may.axis12, 
                   aes(x = MDS1, y = MDS2, label = Taxa),
                   color="black",
                   size = 5) +
  scale_fill_viridis(discrete = TRUE) +
  scale_colour_viridis(discrete = TRUE) +
  scale_shape_manual(values = c(24)) +
  theme(legend.position = "none") +
  ylim(-0.8, 0.8) +
  xlim(-0.8, 0.8) 

invert.12.may 

nmds.may <- ggarrange(legends, invert.12.may,
                      labels = "A",
                      widths = c(0.25,1))

NMDS.may <- annotate_figure(nmds.may,
                            top = text_grob("Collection 20-May-18"))

ggsave("Figures/May_emerging_NMDS.jpeg")


# NMDS collection 5-Jun-18  -----------------------------------------------



## Actual figure
nmds.june.scores <- read.csv("Data/Emerging/NMDS/emerging_june_NMDSscores.csv")
nmds.june.scores$Year <- as.factor(nmds.june.scores$Year)

june.axis12 <- read.csv("Data/Emerging/NMDS/june_emerging_correlated_vector12.csv")
june.axis13 <- read.csv("Data/Emerging/NMDS/june_emerging_correlated_vector13.csv")



invert.12.june <- ggplot(data = nmds.june.scores,
                         aes(x = NMDS1, y = NMDS2)) +
  geom_point(data = nmds.june.scores, 
             aes(x = NMDS1, y = NMDS2, 
                 fill = Treatment, shape = Year),
             size = 4, stroke = 1.5,
             alpha = 0.7) + # sites as points
  stat_ellipse(data = nmds.june.scores, 
               aes(x = NMDS1,y = NMDS2,
                   linetype = Treatment, colour = Treatment), 
               size = 1, level = 0.9) + 
  geom_segment(data = june.axis12, 
               aes(x = 0, xend = MDS1, y = 0, yend = MDS2), # adding in the vectors, c
               arrow = arrow(length = unit(0.5, "cm")), colour = "black") + # can add in geom_label or geom_text for labels
  theme_minimal() + # no background
  theme(panel.border = element_rect(fill = NA)) + # full square around figure
  xlab("NMDS 1") +
  ylab("NMDS 2") +
  #ylim(-1, 1.5) +
  #xlim(-1.45, 1) +
  #theme(legend.position = "none") +
  geom_label_repel(data = june.axis12, 
                   aes(x = MDS1, y = MDS2, label = Taxa),
                   color="black",
                   size = 5) +
  scale_fill_viridis(discrete = TRUE) +
  scale_colour_viridis(discrete = TRUE) +
  scale_shape_manual(values = c(24)) +
  theme(legend.position = "none") 

invert.12.june


invert.13.june <- ggplot(data = nmds.june.scores,
                         aes(x = NMDS1, y = NMDS3)) +
  geom_point(data = nmds.june.scores, 
             aes(x = NMDS1, y = NMDS3, 
                 fill = Treatment, shape = Year),
             size = 4, stroke = 1.5,
             alpha = 0.7) + # sites as points
  stat_ellipse(data = nmds.june.scores, 
               aes(x = NMDS1,y = NMDS3,
                   linetype = Treatment, colour = Treatment), 
               size = 1, level = 0.9) + 
  geom_segment(data = june.axis13, 
               aes(x = 0, xend = MDS1, y = 0, yend = MDS3), # adding in the vectors, c
               arrow = arrow(length = unit(0.5, "cm")), colour = "black") + # can add in geom_label or geom_text for labels
  theme_minimal() + # no background
  theme(panel.border = element_rect(fill = NA)) + # full square around figure
  xlab("NMDS 1") +
  ylab("NMDS 3") +
  #ylim(-1, 1.5) +
  #xlim(-1.45, 1) +
  #theme(legend.position = "none") +
  geom_label_repel(data = june.axis13, 
                   aes(x = MDS1, y = MDS3, label = Taxa),
                   color="black",
                   size = 5) +
  scale_fill_viridis(discrete = TRUE) +
  scale_colour_viridis(discrete = TRUE) +
  scale_shape_manual(values = c(24)) +
  theme(legend.position = "none") 

invert.13.june

nmds.june <- ggarrange(invert.12.june, invert.13.june,
                       labels = c("B", ""))


NMDS.june <- annotate_figure(nmds.june,
                             top = text_grob("Collection 05-June-18"))

ggsave("Figures/NMDS_emerging_5Jun18.jpeg")




# Legend ------------------------------------------------------

test <- ggplot(data = nmds.col4.scores,
               aes(x = NMDS1, y = NMDS3)) +
  geom_point(data = nmds.col4.scores, 
             aes(x = NMDS1, y = NMDS3, 
                 colour = Treatment, shape = Year),
             size = 4) + # sites as points
  stat_ellipse(data = nmds.col4.scores, 
               aes(x = NMDS1,y = NMDS3,
                   linetype = Treatment, colour = Treatment), 
               size = 1, level = 0.9) + 
  geom_segment(data = col4.axis13, 
               aes(x = 0, xend = MDS1, y = 0, yend = MDS3), # adding in the vectors, c
               arrow = arrow(length = unit(0.5, "cm")), colour = "black") + # can add in geom_label or geom_text for labels
  theme_minimal() + # no background
  theme(panel.border = element_rect(fill = NA)) + # full square around figure
  xlab("NMDS 1") +
  ylab("NMDS 3") +
  #ylim(-1, 1.5) +
  #xlim(-1.45, 1) +
  #theme(legend.position = "none") +
  geom_label_repel(data = col4.axis13, 
                   aes(x = MDS1, y = MDS3, label = Taxa),
                   color="black",
                   size = 5) +
  scale_fill_viridis(discrete = TRUE) +
  scale_colour_viridis(discrete = TRUE) +
  scale_shape_manual(values = c(19, 17)) +
  ylim(-0.6, 0.6) +
  xlim(-0.6, 0.6)



legend <- get_legend(test)
legends <- as_ggplot(legend)


# Combine panels ----------------------------------------------------------

NMDS.panel <- ggarrange(NMDS.may, NMDS.june, NMDS.col1, 
                        NMDS.col2, NMDS.col3, NMDS.col4,
                        widths = c(0.7, 1))


ggsave("Figures/NMDS.panel.jpeg", NMDS.panel,
       width = 20,
       height = 9.5,
       units = "in")


# Patchwork ---------------------------------------------------------------

library(patchwork)

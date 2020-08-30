
library(tidyverse)

invert <- read.csv("Data/Invertebrate community.csv")

invert$Year <- as.factor(invert$Year)
str(invert)



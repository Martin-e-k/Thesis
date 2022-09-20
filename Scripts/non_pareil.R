library(Nonpareil)
library(RColorBrewer)
library(tidyverse)
library(broom)

setwd("C:/Users/marti/Desktop/Speciale/nonpareil_files")

nonpareil <- read_csv("nonpareil.txt", col_names = FALSE)
attach(nonpareil)
nps <- Nonpareil.set(X1, plot.opts=list(plot.observed=FALSE))

plot.Nonpareil.Set(x = nps, plot.observed=FALSE, legend.opts = F)

detach(nonpareil)
summary(nps)
Nonpareil.curve(nonpareil$X1[1])
Nonpareil.add.curve(nonpareil$X1[2])

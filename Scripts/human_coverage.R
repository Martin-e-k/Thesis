library(readr)
library(tidyverse)
library(stats)

human_coverage <- read_delim("C:/Users/marti/Desktop/Speciale/data/human_coverage.tsv", 
                             delim = "\t", escape_double = FALSE, 
                             trim_ws = TRUE)

human_coverage$Ratio <- human_coverage$Mapped / human_coverage$Total

human_coverage$Location <- str_split(human_coverage$File, "_")[[1]]

str_split(human_coverage$File, "_")[[1]][1]

human_coverage$Location = sapply(human_coverage$File, function(x) strsplit(x, "_")[[1]][1])

human_coverage[5][human_coverage[5] == "Mid"] <- "Qajaa"


ggplot(human_coverage, aes(x = File, y = Ratio, color = Location)) +
  geom_point() +
  theme_minimal() +
  scale_color_manual(values = c("mediumaquamarine",  "coral2", "darkorchid3"), 
                     name = "Geographic location", 
                     labels = c( "Nuuk", "Qajaa", "Sermermiut")) +
  labs(y = "Ratio human DNA") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

ggplot(human_coverage, aes(x = Location, y = Ratio*100, fill = Location)) +
  geom_boxplot(alpha = 0.7, outlier.alpha = 0) +
  geom_jitter() +
  theme_minimal() +
  scale_fill_manual(values = c("mediumaquamarine",  "coral2", "darkorchid3"), 
                    name = "Geographic location", 
                    labels = c( "Nuuk", "Qajaa", "Sermermiut")) +
  labs(y = "Percentage human DNA") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

wilcox <- pairwise.wilcox.test(human_coverage$Ratio, human_coverage$Location, p.adj = "none")
wilcox_adj <- pairwise.wilcox.test(human_coverage$Ratio, human_coverage$Location, p.adj = "holm")

wilcox$p.value
wilcox_adj$p.value

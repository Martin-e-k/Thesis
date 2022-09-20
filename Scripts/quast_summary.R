library(tidyverse)
library(patchwork)

setwd("C:/Users/marti/Desktop/Speciale/scripts")

assembly_data <- read_tsv("C:/Users/marti/Desktop/Speciale/assembly_summary.tsv")

largest_contig_boxplot <- assembly_data %>% ggplot(aes(x=Largest_contig)) + 
  geom_boxplot() + 
  labs(title = "Length of the largest contig in each assembly") +
  theme_bw() +
  theme(axis.ticks.y = element_blank(),
        axis.text.y = element_blank()) +
  xlab("Largest contig (bp)")

largest_contig_boxplot

N50_boxplot <- assembly_data %>% ggplot(aes(x=N50)) + 
  geom_boxplot() +
  labs(title = "N50 value for each assembly") +
  #xlab("Sample Type") +
  xlab("N50 value") +
  theme_bw() +
  theme(axis.ticks.y = element_blank(),
        axis.text.y = element_blank())

N50_boxplot

total_length_boxplot <- assembly_data %>% ggplot(aes(x=Total_length)) +
  geom_boxplot() +
  labs(title = "Total length of all contigs in each assembly") +
  #xlab("Sample Type") +
  xlab("Total contig length (bp)") +
  theme_bw() +
  theme(axis.ticks.y = element_blank(),
        axis.text.y = element_blank())

total_length_boxplot            

nr_contigs_boxplot <- assembly_data %>% ggplot(aes(x=contigs)) +
  geom_boxplot() +
  labs(title = "Number of contigs in each assembly") +
  #xlab("Sample Type") +
  xlab("Number of contigs") +
  theme_bw() +
  theme(axis.ticks.y = element_blank(),
        axis.text.y = element_blank())

nr_contigs_boxplot                                                     

(N50_boxplot + largest_contig_boxplot) / (total_length_boxplot + nr_contigs_boxplot) + plot_layout(guides = "collect") + plot_annotation(tag_levels = "a")
#ggsave(filename = "largest_contig_boxplot.jpeg", plot = largest_contig_boxplot, device=, width = 16, height = 8)

assembly_data$contig_lengths = assembly_data$Total_length / assembly_data$contigs

contigs_length_boxplot <- assembly_data %>% ggplot(aes(x=contig_lengths)) +
  geom_boxplot() +
  labs(title = "Average contig lengths") +
  #xlab("Sample Type") +
  xlab("Contig lengths") +
  theme_bw()

contigs_length_boxplot

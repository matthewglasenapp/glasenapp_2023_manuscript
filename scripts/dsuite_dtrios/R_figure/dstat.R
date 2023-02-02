library(ggplot2)
library(scales)
library(svglite)
library(ggExtra)

setwd("/Users/matt/Documents/dissertation/dissertation_data/abba_baba/R_figure/")

data = read.csv("/Users/matt/Documents/dissertation/dissertation_data/abba_baba/R_figure/d.csv")

data$taxa <- factor(data$taxa, levels = rev(unique(data$taxa)), ordered=TRUE)

ggplot(data, aes(taxa, d, color=Significance)) +
  geom_point(size=2) +
  geom_errorbar(aes(ymin=d-error, ymax=d+error)) + 
  coord_flip() +
  scale_y_continuous(labels = label_number(accuracy = 0.001)) +
  labs(title="D Statistic Results",
       x="Taxa Triplet (P1,P2),P3",
       y="D Statistic") + 
  scale_color_manual(values=c("light blue", "dark blue")) + 
  #scale_color_brewer(palette="Paired") + 
  #scale_color_manual(values=c("#f1a340", "#998ec3")) + 
  
  theme_linedraw(base_line_size = 2) + 
  removeGrid(x = TRUE, y = FALSE) + 
  theme(
    plot.title = element_text(size=18, face="bold", hjust = 0.5),
    axis.title.x = element_text(size=12, face="bold"),
    axis.title.y = element_text(size=12, face="bold"),
    axis.text.x = element_text(size=10),
    axis.text.y = element_text(size=10, hjust = 1.0),
    legend.position = "bottom",
    legend.title = element_text(size=12),
    legend.text = element_text(size=12))
    #legend.position = c(0.91,0.15))
    #legend.direction = "vertical",
    #legend.key.height = unit(1.5, 'cm'),
    #legend.key.width = unit(1.5, 'cm'))

ggsave("D.svg", width=169, height = 150, units = "mm")


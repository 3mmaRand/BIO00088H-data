library(tidyverse)
library(ggforce)
library(viridis)

# 
col <- viridis(3, begin = 0.1, end = 0.9)

circle_data <- data.frame(x0 = c(1, 1.5, 1.25),
                          y0 = c(1, 1, 1.5),
                          r = c(0.6, 0.6, 0.6),
                          fill = col)
text_data <- data.frame(x = c(1.25, 1.25, 1.25, 0.87, 1.62, 0.75, 1.75),
                        y = c(1.25, 0.75, 1.75, 1.37, 1.37, 0.75, 0.75),
                        label = c("1", "2", "3", "4", "5", "6", "7"))

ggplot() +
  geom_circle(data = circle_data, 
              aes(x0 = x0, y0 = y0, r = r, fill = col), 
              alpha = 0.3) +
  scale_fill_manual(values = col) +
  geom_text(data = text_data, 
            aes(x = x, y = y, label = label), 
            size = 5) +
  geom_text(aes(x = 1, y = 1.8, label = "gene1\ngene2\ngene3"), 
            size = 4) +
  geom_text(aes(x = 2, y = 2, 
                label = "geneA\ngeneB\ngeneC\ngeneD\ngeneE\ngeneF\ngeneG"), 
            size = 4) +
  # pointerto geneA
  geom_curve(aes(x = 1.9, y = 2, 
                 xend = 1.25, yend = 1.3), 
             arrow = arrow(length = unit(0.1, "inches")), 
             curvature = 0.1) +
  coord_fixed() +
  theme_void() +
  theme(legend.position = "none")


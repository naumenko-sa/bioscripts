library(tidyverse)

inh <- read_tsv("omim_inheritance_dictionary.tsv")
inh$inheritance_mode <- factor(inh$inheritance_mode,
                               levels = inh$inheritance_mode[order(inh$count, 
                                                                   decreasing = FALSE)])

offset <- 200

png("omim_inheritance2025.png", 
    width=20, height=10, 
    units="cm", res=300)

inh |> ggplot(aes(x = inheritance_mode, 
                          y = count)) +
  geom_segment(aes(x = inheritance_mode, xend = inheritance_mode, y = 0, yend = count),
               color = "gray", size = 1) +
  geom_point(color = "blue", size = 5) +
  geom_text(aes(label = count), nudge_y = offset, color = "black", size = 4) +
  coord_flip() +  # Flip coordinates to make the chart horizontal
  theme_minimal() +
  labs(title = "Inheritance Mode Counts in OMIM-2025",
       x = "Inheritance Mode",
       y = "Count") + ylim(0, 4000)

dev.off()

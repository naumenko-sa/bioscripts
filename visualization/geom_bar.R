#remotes::install_github("clauswilke/colorblindr")
#install.packages("nyflights13")

library(nyflights13)
library(tidyverse)
#library(colorblindr)

covid <- NULL
covid_stat_13 <- flights %>% 
                select(year, month, day) %>% 
                group_by(year, month) %>% 
                tally(name = "death")

for(i in 1:10){
    covid_stat_buf <- covid_stat_13 
    covid_stat_buf$death <- covid_stat_buf$death + (-1) ^ sample(0:1, 1) * sample(500:(1000*i), 12)
    covid_stat_buf$year <- 2010 + i
    covid <- bind_rows(covid, covid_stat_buf)
}

covid$year <- factor(covid$year)
covid$month <- factor(covid$month)

min_death <- min(covid$death)
max_death <- max(covid$death)

png("death_statistics.png", 
    res = 300, 
    units = "mm", 
    width = 297, height = 210)
covid %>%  ggplot(aes(x = month,
                      y = death, 
                      group = year,
                      fill = year)) + 
    geom_bar(stat = "identity", position = position_dodge2(width = 0.9, preserve = "single")) + 
    scale_fill_brewer("year", 
                      palette = "Paired",
                      guide = guide_legend(nrow = 1))+    
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background =  element_blank(),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          axis.title = element_text(size = 12, face = "bold"),
          plot.title = element_text(size = 12, face = "bold"),
          legend.position = "top",
          legend.title = element_text(size = 12, face = "bold")) +
    scale_x_discrete(labels = month.name) +
    xlab("Month") +
    ylab("Number of deaths") +
    ggtitle("Deaths statistics in X region, 2011-2020")
dev.off()

# label bars with years
# geom_text(aes(label = year), 
#          position = position_dodge2(width = 0.9, preserve = "single"),
#          angle = 90,
#          hjust = -0.1,
#          size =  8 * 5 / 14) +

# use if <9 years
# scale_fill_OkabeIto() +


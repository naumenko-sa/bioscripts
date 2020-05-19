# VennDiagram package of no use?
#https://stackoverflow.com/questions/46801087/scale-circle-size-venn-diagram-by-relative-proportion

install.packages("eulerr")
library(eulerr)

venn_diagram <- euler(c("edits1" = 109124, 
                        "edits2" = 112409, 
                        "edits3" = 120275, 
                        "edits1&edits2" = 103233, 
                        "edits2&edits3" = 98430, 
                        "edits1&edits3" = 93869, 
                        "edits1&edits2&edits3" = 91880))

venn_diagram <- euler(c("edits1" = 109, 
                        "edits2" = 112, 
                        "edits3" = 120, 
                        "edits1&edits2" = 103, 
                        "edits2&edits3" = 98, 
                        "edits1&edits3" = 94, 
                        "edits1&edits2&edits3" = 92))


plot(venn_diagram, 
     counts = TRUE,
     font = 3,
     cex = 1,
     alpha = 0.5,
     fill = c ("red", "green", "blue"))
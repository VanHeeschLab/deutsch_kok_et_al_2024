library(tidyverse)

theme_set(theme_classic() + 
            theme(text = element_text(colour = "black"), 
                  axis.title = element_text(size=7, colour = "black"),
                  axis.text = element_text(size = 7, colour = "black"),
                  axis.ticks = element_line(color = "black"),
                  legend.text = element_text(size = 7, colour = "black"),
                  legend.title = element_text(size = 7, colour = "black")) )

blue_col <- c("#6C88A5", "#7B95B3", "#8BA3C1", "#9AB1CF", "#AAC0DD", "#B9CEE9", "#C9DBF7", "#D8E9FF", "#E7F6FF", "#F5FFFF")
grey_col <- c("#7F7F7F", "#8E8E8E", "#9D9D9D", "#ABABAB", "#BABABA", "#C9C9C9", "#D8D8D8", "#E6E6E6", "#F5F5F5", "#FFFFFF")
pink_col <- c("#EFA891", "#F2B09F", "#F5C8AD", "#F7E1BC", "#FAF9CA", "#FCFAD8", "#FFFBF6", "#FFFCFF", "#FFFDFF", "#FFFFFE")
green_col <- c("#7D905D", "#8D9E6C", "#9DA97B", "#ADB58A", "#BCC19A", "#CBDCA9", "#DAE8B8", "#EAF3C8", "#F9FFD7", "#FFFFE6")
purple_col <- c("#9F85A1", "#AD92AF", "#BCAFC1", "#CBBDD0", "#D9CAE0", "#E8D7EF", "#F6E5FF", "#FFF2FF", "#FFFFFE", "#FFFFFF")

red_col <- c("#B10027")

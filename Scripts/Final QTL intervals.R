library(ggpubr)
library(scales)
library(shades)
library(ggpattern)

##### Prep for Final Graph (setting level heights)

l1t <- 0.25
l1b <- 0.15

l2t <- 0.42
l2b <- 0.32

l3t <- 0.59
l3b <- 0.49

l4t <- 0.76
l4b <- 0.66

l5t <- 0.93
l5b <- 0.83

l6t <- 1.1
l6b <- 1.0

l7t <- 1.27
l7b <- 1.17

l8t <- 1.44
l8b <- 1.34

l9t <- 1.61
l9b <- 1.51

chr_thick <- 0.05
#####

##### Final Graph
## Quick note--2/28/2023--initial ggplot is to graph the markers to the chromosomes as a scaffold for the interval
## figure. Then, rectangles representing QTL intervals are drawn on top of that. Finally, to get the legend with
## all 17 colors, I plot a blank geom_point using the "intervals" object that then generates a legend that has all
## 17 colors matched to their corresponding trait names.
## Also a note: Sterility was dropped from gmap for this, but is still present as row 14 in the intervals object
## That is why the numbering of the rectangles jumps from 13 to 15.
## Sterility wasn't included in the "c16" variable--thats why it wasn't in the legend
## Now with "c17", how do I get Free Filament length to appear? --ended up being capitalized in intervals, but not
##     in c17, which meant it didn't map to the legend.

## 6/10/23 I had to shift the graph to incorporate one more interval for Free Filament length on Chr 6.
## The loops below change chromosome values greater than 6 to add 0.27 to allow room for another box

out.ehk <- scanone(gmap, pheno.col = 1, method = "ehk")
#levels(out.ehk$chr) <- c("1","2","3","4","5","6","7.5","8.5","9.5")
out.ehk$chr <- as.numeric(out.ehk$chr)
for (i in 1:length(out.ehk$chr)){
  if (out.ehk$chr[i] > 6){
    out.ehk$chr[i] <- out.ehk$chr[i] + 0.27
  }
}

for (i in 1:length(intervals$chr)){
  if (intervals$chr[i] > 6){
    intervals$chr[i] <- intervals$chr[i] + 0.27
  }
}



plot <- ggplot(out.ehk, aes(x=out.ehk$pos, y=out.ehk$chr)) +
  geom_point(pch = "|", size = 3, key_glyph = "rect", color = "darkgrey") +
  theme_minimal() +
  theme(legend.position = "none")+
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)))+
  theme(axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)))+
  scale_y_continuous(breaks=c(1,2,3,4,5,6,7.27,8.27,9.27), labels = c("1","2","3","4","5","6","7","8","9"))+
  
  annotate("segment", x = 0, xend = max[1], y = 1, yend = 1, size = chr_thick, color = "darkgrey")+
  annotate("segment", x = 0, xend = max[2], y = 2, yend = 2, size = chr_thick, color = "darkgrey")+
  annotate("segment", x = 0, xend = max[3], y = 3, yend = 3, size = chr_thick, color = "darkgrey")+
  annotate("segment", x = 0, xend = max[4], y = 4, yend = 4, size = chr_thick, color = "darkgrey")+
  annotate("segment", x = 0, xend = max[5], y = 5, yend = 5, size = chr_thick, color = "darkgrey")+
  annotate("segment", x = 0, xend = max[6], y = 6, yend = 6, size = chr_thick, color = "darkgrey")+
  annotate("segment", x = 0, xend = max[7], y = 7.27, yend = 7.27, size = chr_thick, color = "darkgrey")+
  annotate("segment", x = 0, xend = max[8], y = 8.27, yend = 8.27, size = chr_thick, color = "darkgrey")+
  annotate("segment", x = 0, xend = max[9], y = 9.27, yend = 9.27, size = chr_thick, color = "darkgrey")+
  labs(x = "Genetic Position (cM)",
       y = "Chromosome")+
  coord_cartesian(ylim = c(1,10.47))+
  
  ##Boxes showing right and left endpoints of the intervals
  geom_rect(aes(ymax = intervals[1,1]+l1t, ymin = intervals[1,1]+l1b, xmin = intervals[1,2], xmax = intervals[1,4]), color = 'black', fill = c17[1],size = 0.3) +
  geom_rect(aes(ymax = intervals[2,1]+l1t, ymin = intervals[2,1]+l1b, xmin = intervals[2,2], xmax = intervals[2,4]), color = 'black', fill = c17[2],size = 0.3) +
  geom_rect(aes(ymax = intervals[3,1]+l3t, ymin = intervals[3,1]+l3b, xmin = intervals[3,2], xmax = intervals[3,4]), color = 'black', fill = c17[3],size = 0.3) +
  geom_rect(aes(ymax = intervals[4,1]+l2t, ymin = intervals[4,1]+l2b, xmin = intervals[4,2], xmax = intervals[4,4]), color = 'black', fill = c17[4],size = 0.3) +
  geom_rect(aes(ymax = intervals[5,1]+l7t, ymin = intervals[5,1]+l7b, xmin = intervals[5,2], xmax = intervals[5,4]), color = 'black', fill = c17[5],size = 0.3) +
  geom_rect(aes(ymax = intervals[6,1]+l4t, ymin = intervals[6,1]+l4b, xmin = intervals[6,2], xmax = intervals[6,4]), color = 'black', fill = c17[6],size = 0.3) +
  geom_rect(aes(ymax = intervals[7,1]+l6t, ymin = intervals[7,1]+l6b, xmin = intervals[7,2], xmax = intervals[7,4]), color = 'black', fill = c17[7],size = 0.3) +
  geom_rect(aes(ymax = intervals[8,1]+l1t, ymin = intervals[8,1]+l1b, xmin = intervals[8,2], xmax = intervals[8,4]), color = 'black', fill = c17[8],size = 0.3) +
  geom_rect(aes(ymax = intervals[9,1]+l1t, ymin = intervals[9,1]+l1b, xmin = intervals[9,2], xmax = intervals[9,4]), color = 'black', fill = c17[9],size = 0.3) +
  geom_rect(aes(ymax = intervals[10,1]+l5t, ymin = intervals[10,1]+l5b, xmin = intervals[10,2], xmax = intervals[10,4]), color = 'black', fill = c17[10],size = 0.3) +
  geom_rect(aes(ymax = intervals[11,1]+l2t, ymin = intervals[11,1]+l2b, xmin = intervals[11,2], xmax = intervals[11,4]), color = 'black', fill = c17[11],size = 0.3) +
  geom_rect(aes(ymax = intervals[12,1]+l8t, ymin = intervals[12,1]+l8b, xmin = intervals[12,2], xmax = intervals[12,4]), color = 'black', fill = c17[12],size = 0.3) +
  geom_rect(aes(ymax = intervals[13,1]+l9t, ymin = intervals[13,1]+l9b, xmin = intervals[13,2], xmax = intervals[13,4]), color = 'black', fill = c17[13],size = 0.3) +
  geom_rect(aes(ymax = intervals[15,1]+l1t, ymin = intervals[15,1]+l1b, xmin = intervals[15,2], xmax = intervals[15,4]), color = 'black', fill = c17[14],size = 0.3) +
  geom_rect(aes(ymax = intervals[16,1]+l3t, ymin = intervals[16,1]+l3b, xmin = intervals[16,2], xmax = intervals[16,4]), color = 'black', fill = c17[15],size = 0.3) +
  geom_rect(aes(ymax = intervals[17,1]+l4t, ymin = intervals[17,1]+l4b, xmin = intervals[17,2], xmax = intervals[17,4]), color = 'black', fill = c17[16],size = 0.3) +
  geom_rect(aes(ymax = intervals[24,1]+l5t, ymin = intervals[24,1]+l5b, xmin = intervals[24,2], xmax = intervals[24,4]), color = 'black', fill = c17[17],size = 0.3) +
  
  geom_rect(aes(ymax = intervals[18,1]+l2t, ymin = intervals[18,1]+l2b, xmin = intervals[18,2], xmax = intervals[18,4]), color = 'black', fill = c17[7],size = 0.3) +
  geom_rect(aes(ymax = intervals[19,1]+l4t, ymin = intervals[19,1]+l4b, xmin = intervals[19,2], xmax = intervals[19,4]), color = 'black', fill = c17[8],size = 0.3) +
  geom_rect(aes(ymax = intervals[20,1]+l1t, ymin = intervals[20,1]+l1b, xmin = intervals[20,2], xmax = intervals[20,4]), color = 'black', fill = c17[14],size = 0.3) +
  geom_rect(aes(ymax = intervals[21,1]+l1t, ymin = intervals[21,1]+l1b, xmin = intervals[21,2], xmax = intervals[21,4]), color = 'black', fill = c17[14],size = 0.3) +
  geom_rect(aes(ymax = intervals[22,1]+l1t, ymin = intervals[22,1]+l1b, xmin = intervals[22,2], xmax = intervals[22,4]), color = 'black', fill = c17[14],size = 0.3) +
  geom_rect(aes(ymax = intervals[23,1]+l1t, ymin = intervals[23,1]+l1b, xmin = intervals[23,2], xmax = intervals[23,4]), color = 'black', fill = c17[15],size = 0.3) +

  
  
  
  ##These are for marking the top position within each interval (just changed xmin and xmax to both be the top position)
#  geom_rect(aes(ymax = intervals[1,1]+l1t, ymin = intervals[1,1]+l1b, xmin = intervals[1,3], xmax = intervals[1,3]), color = 'black', fill = c16[1],size = 0.3) +
  geom_rect(aes(ymax = intervals[2,1]+l1t, ymin = intervals[2,1]+l1b, xmin = intervals[2,3], xmax = intervals[2,3]), color = 'black', fill = c17[2],size = 0.3) +
  geom_rect(aes(ymax = intervals[3,1]+l3t, ymin = intervals[3,1]+l3b, xmin = intervals[3,3], xmax = intervals[3,3]), color = 'black', fill = c17[3],size = 0.3) +
  geom_rect(aes(ymax = intervals[4,1]+l2t, ymin = intervals[4,1]+l2b, xmin = intervals[4,3], xmax = intervals[4,3]), color = 'black', fill = c17[4],size = 0.3) +
  geom_rect(aes(ymax = intervals[5,1]+l7t, ymin = intervals[5,1]+l7b, xmin = intervals[5,3], xmax = intervals[5,3]), color = 'black', fill = c17[5],size = 0.3) +
  geom_rect(aes(ymax = intervals[6,1]+l4t, ymin = intervals[6,1]+l4b, xmin = intervals[6,3], xmax = intervals[6,3]), color = 'black', fill = c17[6],size = 0.3) +
  geom_rect(aes(ymax = intervals[7,1]+l6t, ymin = intervals[7,1]+l6b, xmin = intervals[7,3], xmax = intervals[7,3]), color = 'black', fill = c17[7],size = 0.3) +
  geom_rect(aes(ymax = intervals[8,1]+l1t, ymin = intervals[8,1]+l1b, xmin = intervals[8,3], xmax = intervals[8,3]), color = 'black', fill = c17[8],size = 0.3) +
  geom_rect(aes(ymax = intervals[9,1]+l1t, ymin = intervals[9,1]+l1b, xmin = intervals[9,3], xmax = intervals[9,3]), color = 'black', fill = c17[9],size = 0.3) +
  geom_rect(aes(ymax = intervals[10,1]+l5t, ymin = intervals[10,1]+l5b, xmin = intervals[10,3], xmax = intervals[10,3]), color = 'black', fill = c17[10],size = 0.3) +
  geom_rect(aes(ymax = intervals[11,1]+l2t, ymin = intervals[11,1]+l2b, xmin = intervals[11,3], xmax = intervals[11,3]), color = 'black', fill = c17[11],size = 0.3) +
  geom_rect(aes(ymax = intervals[12,1]+l8t, ymin = intervals[12,1]+l8b, xmin = intervals[12,3], xmax = intervals[12,3]), color = 'black', fill = c17[12],size = 0.3) +
  geom_rect(aes(ymax = intervals[13,1]+l9t, ymin = intervals[13,1]+l9b, xmin = intervals[13,3], xmax = intervals[13,3]), color = 'black', fill = c17[13],size = 0.3) +
  geom_rect(aes(ymax = intervals[15,1]+l1t, ymin = intervals[15,1]+l1b, xmin = intervals[15,3], xmax = intervals[15,3]), color = 'black', fill = c17[15],size = 0.3) +
  geom_rect(aes(ymax = intervals[16,1]+l3t, ymin = intervals[16,1]+l3b, xmin = intervals[16,3], xmax = intervals[16,3]), color = 'black', fill = c17[16],size = 0.3) +
  geom_rect(aes(ymax = intervals[17,1]+l4t, ymin = intervals[17,1]+l4b, xmin = intervals[17,3], xmax = intervals[17,3]), color = 'black', fill = c17[16],size = 0.3) +
  geom_rect(aes(ymax = intervals[24,1]+l5t, ymin = intervals[24,1]+l5b, xmin = intervals[24,3], xmax = intervals[24,3]), color = 'black', fill = c17[17],size = 0.3) +
  
  geom_rect(aes(ymax = intervals[18,1]+l2t, ymin = intervals[18,1]+l2b, xmin = intervals[18,3], xmax = intervals[18,3]), color = 'black', fill = c17[7],size = 0.3) +
  geom_rect(aes(ymax = intervals[19,1]+l4t, ymin = intervals[19,1]+l4b, xmin = intervals[19,3], xmax = intervals[19,3]), color = 'black', fill = c17[8],size = 0.3) +
  geom_rect(aes(ymax = intervals[20,1]+l1t, ymin = intervals[20,1]+l1b, xmin = intervals[20,3], xmax = intervals[20,3]), color = 'black', fill = c17[14],size = 0.3) +
  geom_rect(aes(ymax = intervals[21,1]+l1t, ymin = intervals[21,1]+l1b, xmin = intervals[21,3], xmax = intervals[21,3]), color = 'black', fill = c17[14],size = 0.3) +
  geom_rect(aes(ymax = intervals[22,1]+l1t, ymin = intervals[22,1]+l1b, xmin = intervals[22,3], xmax = intervals[22,3]), color = 'black', fill = c17[14],size = 0.3) +
  geom_rect(aes(ymax = intervals[23,1]+l1t, ymin = intervals[23,1]+l1b, xmin = intervals[23,3], xmax = intervals[23,3]), color = 'black', fill = c17[15],size = 0.3) +
  
  ##This is for setting up the legend, it's a blank geom_point
  geom_point(data = intervals, aes(x = 0, y = 1, fill=trait), colour = "grey50", key_glyph = draw_key_polygon) +
  scale_fill_manual(values = c17,
                    breaks = names(c17))+
  theme(legend.position = "right",
        legend.key.size = unit(0.75, 'cm'), #change legend key size
        legend.key.height = unit(0.75, 'cm'), #change legend key height
        legend.key.width = unit(0.75, 'cm'), #change legend key width
        legend.title = element_text(size=14), #change legend title font size
        legend.text = element_text(size=10)) + #change legend text font size)
  labs(fill = "Legend")
  
plot

pdf(file = "test8.pdf")
plot
dev.off()

## 2/28/2023 Saved interval figure as pdf with dimensions w 10" x h 8.5" manually from plot screen

##############################################################


## Code to generate lines for all rectangles at once
## Not correct colors for rectangles--intervals contains sterility, while c17 doesn't, creating a problem after index 14

for(i in 1:nrow(intervals)){
  newLine <- paste0("geom_rect(aes(ymax = intervals[",i,",1]+l1t, ymin = intervals[",i,",1]+l1b, xmin = intervals[",i,",2], xmax = intervals[",i,",4]), color = 'black', fill = c17[",i,"],size = 0.3) +")
  print(noquote(newLine))
}



###############################################################


## Code to generate blank graph, used to create legend
## Update: 2/28/2023 No longer needed, got legend working on other graph properly


#data <- data.frame(colnames(gmap$pheno),rep(0,14),rep(0,14))
#rownames(data) <- colnames(gmap$pheno)
#colnames(data) <- c("col","x","y")

#blank <- ggplot(data.frame(data), aes(x=x, y=y, fill = col))+
#  scale_fill_manual(values=colors)+
#  geom_point(key_glyph = "rect")+
  
#  labs(fill = "Legend")


#pdf(file = "test1.pdf")
#ggarrange(plot, widths = c(4,1))
#dev.off()

#################################################

#g_legend <- function(a.gplot){ 
#  tmp <- ggplot_gtable(ggplot_build(a.gplot)) 
#  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
#  legend <- tmp$grobs[[leg]] 
#  legend
#} 

#legend <- g_legend(blank) 

#grid.newpage()
#grid.draw(legend) 

###############################################################

"Anther.length"="#aaffc3"
"Major.length.diameter"="#9A6324"

## Color stuff

show_col(c("#aaffc3","#808000","#e6194B","#f58231","#ffe119","#bfef45","#3cb44b","#9A6324","#469990","#fabed4","#4363d8","#000075","#911eb4","#f032e6"))
colors <- c("#aaffc3","#808000","#e6194B","#f58231","#ffe119","#bfef45","#3cb44b","#9A6324","#469990","#fabed4","#4363d8","#000075","#911eb4","#f032e6")

colorset <- c("Anther length"="#aaffc3", "Anther width"="#808000", "Calyx length"="#e6194B", "Calyx tooth"="#f58231", "Corolla length"="#ffe119", "Lobe length"="#bfef45",
              "Lobe width"="#3cb44b", "Major length diameter"="#9A6324", "Midrib"="#469990", "Pedicel length"="#fabed4", "Stamen length"="#4363d8",
              "Stigma length"="#000075", "Style length"="#911eb4", "Throat"="#f032e6")

# old colors
#colors <- c("#e6194B","#ffe119","#3cb44b","#42d4f4","#4363d8","#000075","#911eb4","#f58231","#bfef45")
#colors <- c("#9A6324","#808000","#e6194B","#f58231","#ffe119","#bfef45","#3cb44b","#aaffc3","#469990","#fabed4","#4363d8","#000075","#911eb4","#f032e6")


#colorset <- c("Anther.length"="#7D0112", "Anther.width"="#9A4509", "Calyx.length"="#B17217", "Calyx.tooth"="#C29A38", "Corolla.length"="#CBBD5E", "Lobe.length"="#CDDB87",
#              "Lobe.width"="#CDEFB0", "Major.length.diameter"="#C3E9C7", "Midrib"="#90D6B8", "Pedicel.length"="#4FBEAF", "Stamen.length"="#00A3AA",
#              "Stigma.length"="#0084A5", "Style.length"="#005FA0", "Throat"="#1F28A2")





##### Test Graph--2/25/2023

plot <- ggplot(out.ehk, aes(x=out.ehk$pos, y=out.ehk$chr, color = 'black')) +
  geom_point(pch = "|", size = 3, key_glyph = "rect") +
  theme_minimal() +
  theme(legend.position = "right")+
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)))+
  theme(axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)))+
  
  annotate("segment", x = 0, xend = max[1], y = 1, yend = 1, size = chr_thick)+
  annotate("segment", x = 0, xend = max[2], y = 2, yend = 2, size = chr_thick)+
  annotate("segment", x = 0, xend = max[3], y = 3, yend = 3, size = chr_thick)+
  annotate("segment", x = 0, xend = max[4], y = 4, yend = 4, size = chr_thick)+
  annotate("segment", x = 0, xend = max[5], y = 5, yend = 5, size = chr_thick)+
  annotate("segment", x = 0, xend = max[6], y = 6, yend = 6, size = chr_thick)+
  annotate("segment", x = 0, xend = max[7], y = 7, yend = 7, size = chr_thick)+
  annotate("segment", x = 0, xend = max[8], y = 8, yend = 8, size = chr_thick)+
  annotate("segment", x = 0, xend = max[9], y = 9, yend = 9, size = chr_thick)+
  labs(x = "Genetic Position (cM)",
       y = "Chromosome",
       color = "Legend")+
  scale_color_manual(values = colorset)+
  coord_cartesian(ylim = c(1,10))+
  
  ##Solitary Flower Intervals 1
  #geom_rect(aes(ymax = intervals[1,1]+l1t, ymin = intervals[1,1]+l1b, xmin = intervals[1,2], xmax = intervals[1,3]), color = colors[1], fill = colors[1],size = 1) +
  #geom_rect(aes(ymax = intervals[1,1]+l1t, ymin = intervals[1,1]+l1b, xmin = intervals[1,2], xmax = intervals[1,2]), color = colors[1], fill = colors[1],size = 1) +
  
  geom_rect(aes(ymax = intervals[1,1]+l1t, ymin = intervals[1,1]+l1b, xmin = intervals[1,2], xmax = intervals[1,4]), color = 'black', fill = colors[1],size = 0.3) +
  geom_rect(aes(ymax = intervals[2,1]+l5t, ymin = intervals[2,1]+l5b, xmin = intervals[2,2], xmax = intervals[2,4]), color = 'black', fill = colors[2],size = 0.3) +
  geom_rect(aes(ymax = intervals[3,1]+l7t, ymin = intervals[3,1]+l7b, xmin = intervals[3,2], xmax = intervals[3,4]), color = 'black', fill = colors[3],size = 0.3) +
  geom_rect(aes(ymax = intervals[4,1]+l1t, ymin = intervals[4,1]+l1b, xmin = intervals[4,2], xmax = intervals[4,4]), color = 'black', fill = colors[4],size = 0.3) +
  geom_rect(aes(ymax = intervals[5,1]+l2t, ymin = intervals[5,1]+l2b, xmin = intervals[5,2], xmax = intervals[5,4]), color = 'black', fill = colors[5],size = 0.3) +
  geom_rect(aes(ymax = intervals[6,1]+l3t, ymin = intervals[6,1]+l3b, xmin = intervals[6,2], xmax = intervals[6,4]), color = 'black', fill = colors[6],size = 0.3) +
  geom_rect(aes(ymax = intervals[7,1]+l1t, ymin = intervals[7,1]+l1b, xmin = intervals[7,2], xmax = intervals[7,4]), color = 'black', fill = colors[7],size = 0.3) +
  geom_rect(aes(ymax = intervals[8,1]+l1t, ymin = intervals[8,1]+l1b, xmin = intervals[8,2], xmax = intervals[8,4]), color = 'black', fill = colors[8],size = 0.3) +
  geom_rect(aes(ymax = intervals[9,1]+l4t, ymin = intervals[9,1]+l4b, xmin = intervals[9,2], xmax = intervals[9,4]), color = 'black', fill = colors[9],size = 0.3) +
  geom_rect(aes(ymax = intervals[10,1]+l1t, ymin = intervals[10,1]+l1b, xmin = intervals[10,2], xmax = intervals[10,4]), color = 'black', fill = colors[10],size = 0.3) +
  geom_rect(aes(ymax = intervals[11,1]+l6t, ymin = intervals[11,1]+l6b, xmin = intervals[11,2], xmax = intervals[11,4]), color = 'black', fill = colors[11],size = 0.3) +
  geom_rect(aes(ymax = intervals[12,1]+l1t, ymin = intervals[12,1]+l1b, xmin = intervals[12,2], xmax = intervals[12,4]), color = 'black', fill = colors[12],size = 0.3) +
  geom_rect(aes(ymax = intervals[13,1]+l1t, ymin = intervals[13,1]+l1b, xmin = intervals[13,2], xmax = intervals[13,4]), color = 'black', fill = colors[13],size = 0.3) +
  geom_rect(aes(ymax = intervals[14,1]+l7t, ymin = intervals[14,1]+l7b, xmin = intervals[14,2], xmax = intervals[14,4]), color = 'black', fill = colors[14],size = 0.3) +
  
  geom_rect(aes(ymax = 6+l2t, ymin = 6+l2b, xmin = 37.95, xmax = 365.00), color = 'black', fill = colors[7],size = 0.3) +
  geom_rect(aes(ymax = 2+l1t, ymin = 2+l1b, xmin = 20.69, xmax = 137.00), color = 'black', fill = colors[12],size = 0.3)
#2

pdf(file = "test1.pdf")
plot
dev.off()









































##### testing

plot <- ggplot(out.ehk, aes(x=out.ehk$pos, y=out.ehk$chr, color = 'black')) +
  geom_point(pch = "|", size = 3, key_glyph = "rect") +
  theme_minimal() +
  theme(legend.position = "right")+
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)))+
  theme(axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)))+
  
  annotate("segment", x = 0, xend = max[1], y = 1, yend = 1, size = chr_thick)+
  annotate("segment", x = 0, xend = max[2], y = 2, yend = 2, size = chr_thick)+
  annotate("segment", x = 0, xend = max[3], y = 3, yend = 3, size = chr_thick)+
  annotate("segment", x = 0, xend = max[4], y = 4, yend = 4, size = chr_thick)+
  annotate("segment", x = 0, xend = max[5], y = 5, yend = 5, size = chr_thick)+
  annotate("segment", x = 0, xend = max[6], y = 6, yend = 6, size = chr_thick)+
  annotate("segment", x = 0, xend = max[7], y = 7, yend = 7, size = chr_thick)+
  annotate("segment", x = 0, xend = max[8], y = 8, yend = 8, size = chr_thick)+
  annotate("segment", x = 0, xend = max[9], y = 9, yend = 9, size = chr_thick)+
  labs(x = "Genetic Position (cM)",
       y = "Chromosome",
       color = "Legend")+
  scale_fill_manual(values = c17,
                    breaks = names(c17))+
  coord_cartesian(ylim = c(1,10))+
  
  ##Solitary Flower Intervals 1
  #geom_rect(aes(ymax = intervals[1,1]+l1t, ymin = intervals[1,1]+l1b, xmin = intervals[1,2], xmax = intervals[1,3]), color = colors[1], fill = colors[1],size = 1) +
  #geom_rect(aes(ymax = intervals[1,1]+l1t, ymin = intervals[1,1]+l1b, xmin = intervals[1,2], xmax = intervals[1,2]), color = colors[1], fill = colors[1],size = 1) +
  
  
  geom_rect(data = intervals, aes(ymax = chr+l1t, ymin = chr-l1b, xmin = lpos, xmax = rpos, fill = trait))
  
  
  
  



  geom_rect(aes(ymax = intervals[1,1]+l1t, ymin = intervals[1,1]+l1b, xmin = intervals[1,2], xmax = intervals[1,4]), color = 'black', fill = c17[1],size = 0.3) +
  geom_rect(aes(ymax = intervals[2,1]+l1t, ymin = intervals[2,1]+l1b, xmin = intervals[2,2], xmax = intervals[2,4]), color = 'black', fill = c17[2],size = 0.3) +
  geom_rect(aes(ymax = intervals[3,1]+l1t, ymin = intervals[3,1]+l1b, xmin = intervals[3,2], xmax = intervals[3,4]), color = 'black', fill = c17[3],size = 0.3) +
  geom_rect(aes(ymax = intervals[4,1]+l1t, ymin = intervals[4,1]+l1b, xmin = intervals[4,2], xmax = intervals[4,4]), color = 'black', fill = c17[4],size = 0.3) +
  geom_rect(aes(ymax = intervals[5,1]+l1t, ymin = intervals[5,1]+l1b, xmin = intervals[5,2], xmax = intervals[5,4]), color = 'black', fill = c17[5],size = 0.3) +
  geom_rect(aes(ymax = intervals[6,1]+l1t, ymin = intervals[6,1]+l1b, xmin = intervals[6,2], xmax = intervals[6,4]), color = 'black', fill = c17[6],size = 0.3) +
  geom_rect(aes(ymax = intervals[7,1]+l1t, ymin = intervals[7,1]+l1b, xmin = intervals[7,2], xmax = intervals[7,4]), color = 'black', fill = c17[7],size = 0.3) +
  geom_rect(aes(ymax = intervals[8,1]+l1t, ymin = intervals[8,1]+l1b, xmin = intervals[8,2], xmax = intervals[8,4]), color = 'black', fill = c17[8],size = 0.3) +
  geom_rect(aes(ymax = intervals[9,1]+l1t, ymin = intervals[9,1]+l1b, xmin = intervals[9,2], xmax = intervals[9,4]), color = 'black', fill = c17[9],size = 0.3) +
  geom_rect(aes(ymax = intervals[10,1]+l1t, ymin = intervals[10,1]+l1b, xmin = intervals[10,2], xmax = intervals[10,4]), color = 'black', fill = c17[10],size = 0.3) +
  geom_rect(aes(ymax = intervals[11,1]+l1t, ymin = intervals[11,1]+l1b, xmin = intervals[11,2], xmax = intervals[11,4]), color = 'black', fill = c17[11],size = 0.3) +
  geom_rect(aes(ymax = intervals[12,1]+l1t, ymin = intervals[12,1]+l1b, xmin = intervals[12,2], xmax = intervals[12,4]), color = 'black', fill = c17[12],size = 0.3) +
  geom_rect(aes(ymax = intervals[13,1]+l1t, ymin = intervals[13,1]+l1b, xmin = intervals[13,2], xmax = intervals[13,4]), color = 'black', fill = c17[13],size = 0.3) +
  geom_rect(aes(ymax = intervals[14,1]+l1t, ymin = intervals[14,1]+l1b, xmin = intervals[14,2], xmax = intervals[14,4]), color = 'black', fill = c17[14],size = 0.3) +
  geom_rect(aes(ymax = intervals[15,1]+l1t, ymin = intervals[15,1]+l1b, xmin = intervals[15,2], xmax = intervals[15,4]), color = 'black', fill = c17[15],size = 0.3) +
  geom_rect(aes(ymax = intervals[16,1]+l1t, ymin = intervals[16,1]+l1b, xmin = intervals[16,2], xmax = intervals[16,4]), color = 'black', fill = c17[16],size = 0.3) +
  geom_rect(aes(ymax = intervals[17,1]+l1t, ymin = intervals[17,1]+l1b, xmin = intervals[17,2], xmax = intervals[17,4]), color = 'black', fill = c17[17],size = 0.3) +
  
  geom_rect(aes(ymax = intervals[18,1]+l1t, ymin = intervals[18,1]+l1b, xmin = intervals[18,2], xmax = intervals[18,4]), color = 'black', fill = c17[7],size = 0.3) +
  geom_rect(aes(ymax = intervals[19,1]+l1t, ymin = intervals[19,1]+l1b, xmin = intervals[19,2], xmax = intervals[19,4]), color = 'black', fill = c17[8],size = 0.3) +
  geom_rect(aes(ymax = intervals[20,1]+l1t, ymin = intervals[20,1]+l1b, xmin = intervals[20,2], xmax = intervals[20,4]), color = 'black', fill = c17[15],size = 0.3) +
  geom_rect(aes(ymax = intervals[21,1]+l1t, ymin = intervals[21,1]+l1b, xmin = intervals[21,2], xmax = intervals[21,4]), color = 'black', fill = c17[15],size = 0.3) +
  geom_rect(aes(ymax = intervals[22,1]+l1t, ymin = intervals[22,1]+l1b, xmin = intervals[22,2], xmax = intervals[22,4]), color = 'black', fill = c17[15],size = 0.3) +
  geom_rect(aes(ymax = intervals[23,1]+l1t, ymin = intervals[23,1]+l1b, xmin = intervals[23,2], xmax = intervals[23,4]), color = 'black', fill = c17[16],size = 0.3)



################################################
  ################################################
  ################################################
  
plot <- ggplot(out.ehk, aes(x=out.ehk$pos, y=out.ehk$chr)) +
  geom_point(pch = "|", size = 3, key_glyph = "rect", color = "darkgrey") +
  theme_minimal() +
  theme(legend.position = "none")+
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)))+
  theme(axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)))+
  
  annotate("segment", x = 0, xend = max[1], y = 1, yend = 1, size = chr_thick, color = "darkgrey")+
  annotate("segment", x = 0, xend = max[2], y = 2, yend = 2, size = chr_thick, color = "darkgrey")+
  annotate("segment", x = 0, xend = max[3], y = 3, yend = 3, size = chr_thick, color = "darkgrey")+
  annotate("segment", x = 0, xend = max[4], y = 4, yend = 4, size = chr_thick, color = "darkgrey")+
  annotate("segment", x = 0, xend = max[5], y = 5, yend = 5, size = chr_thick, color = "darkgrey")+
  annotate("segment", x = 0, xend = max[6], y = 6, yend = 6, size = chr_thick, color = "darkgrey")+
  annotate("segment", x = 0, xend = max[7], y = 7, yend = 7, size = chr_thick, color = "darkgrey")+
  annotate("segment", x = 0, xend = max[8], y = 8, yend = 8, size = chr_thick, color = "darkgrey")+
  annotate("segment", x = 0, xend = max[9], y = 9, yend = 9, size = chr_thick, color = "darkgrey")+
  labs(x = "Genetic Position (cM)",
       y = "Chromosome")+
  coord_cartesian(ylim = c(1,10.2))+
  geom_rect(aes(ymax = intervals[21,1]+l1t, ymin = intervals[21,1]+l1b, xmin = intervals[21,3], xmax = intervals[21,3]), color = 'black', fill = c16[15],size = 0.3) +
  geom_rect(aes(ymax = intervals[22,1]+l1t, ymin = intervals[22,1]+l1b, xmin = intervals[22,3], xmax = intervals[22,3]), color = 'black', fill = c16[15],size = 0.3) +
  
  geom_rect(aes(ymax = intervals[23,1]+l1t, ymin = intervals[23,1]+l1b, xmin = intervals[23,3], xmax = intervals[23,3]), color = 'black', fill = c16[16],size = 0.3) +
#  geom_rect_pattern(aes(ymax = intervals[23,1]+l1t, ymin = intervals[23,1]+l1b, xmin = intervals[23,3], xmax = intervals[23,3]), color = 'black', fill = c16[16],size = 0.3) +
#  pattern_fill(pattern = "crosshatch", pattern_color = "blue", pattern_size = 0.2) +
  
  ##This is for setting up the legend, it's a blank geom_point
  geom_point(data = intervals, aes(x = lpos, y = rpos, fill=trait), colour = "grey50", key_glyph = draw_key_polygon) +
  scale_fill_manual(values = c16,
                    breaks = names(c16))+
  theme(legend.position = "right") +
  labs(fill = "Legend")

plot

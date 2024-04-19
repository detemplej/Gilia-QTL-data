## Inherits objects from "QTL_dataprep.R" script (must be run first)

## Levels to assist with evenness in plotting QTL interval boxes
l1t <- 0.25
l1b <- 0.15

l2t <- 0.42
l2b <- 0.32

l3t <- 0.59
l3b <- 0.49

l4t <- 0.76
l4b <- 0.66

##################### Graph notes ##############################################
## initial ggplot "geom_point" is to show the marker distribution along the
## chromosomes. Segments representing the entire length of each chromosome are
## then drawn to connect all markers on each chromosome. QTL intervals are 
## drawn, then a label with the trait name and a connecting segment is drawn
## for each individual QTL. Some adjustment was done in the code that generates
## the geom_rect, geom_text, geom_segment, at the bottom of the script, and some
## extra adjustment was done manually within the ggplot itself to avoid overlap
## between QTLs, labels, and connecting line segments

## Initial scanone output, used to plot marker positions across chromosomes
out.ehk <- scanone(gmap_all, pheno.col = 1, method = "ehk")
out.ehk$chr <- as.numeric(out.ehk$chr)

## Color object 
white <- rep("#FFFFFF",nrow(cim_intervals_significant))
names(white) <- cim_intervals_significant[,1]

################## Interval Plot ###############################################
plot <- ggplot(out.ehk, aes(x=out.ehk$pos, y=out.ehk$chr)) +
  
  ## Plotting marker positions across all chromosomes
  geom_point(pch = "|", size = 3, key_glyph = "rect", color = "darkgrey") +
  
  ## Plot themes
  theme_minimal() +
  theme(legend.position = "none")+
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0)))+
  theme(axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)))+
  scale_y_continuous(breaks=c(1,2,3,4,5,6,7,8,9), labels = c("1","2","3","4","5","6","7","8","9"))+
  
  ## Drawing in chromosomes
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
  
  ## Extra padding on top of graph to accomodate multiple QTLs on Chr. 9
  coord_cartesian(ylim = c(1,9.5))+
  
  ## Boxes representing QTL intervals
  geom_rect(aes(ymax = cim_intervals_significant[1,2]+l1t, ymin = cim_intervals_significant[1,2]+l1b, xmin = cim_intervals_significant[1,3], xmax = cim_intervals_significant[1,5]), color = 'black', fill = white[1],size = 0.3) + 
  geom_rect(aes(ymax = cim_intervals_significant[2,2]+l1t, ymin = cim_intervals_significant[2,2]+l1b, xmin = cim_intervals_significant[2,3], xmax = cim_intervals_significant[2,5]), color = 'black', fill = white[2],size = 0.3) + 
  geom_rect(aes(ymax = cim_intervals_significant[3,2]+l1t, ymin = cim_intervals_significant[3,2]+l1b, xmin = cim_intervals_significant[3,3], xmax = cim_intervals_significant[3,5]), color = 'black', fill = white[3],size = 0.3) + 
  geom_rect(aes(ymax = cim_intervals_significant[4,2]+l1t, ymin = cim_intervals_significant[4,2]+l1b, xmin = cim_intervals_significant[4,3], xmax = cim_intervals_significant[4,5]), color = 'black', fill = white[4],size = 0.3) + 
  geom_rect(aes(ymax = cim_intervals_significant[5,2]+l1t, ymin = cim_intervals_significant[5,2]+l1b, xmin = cim_intervals_significant[5,3], xmax = cim_intervals_significant[5,5]), color = 'black', fill = white[5],size = 0.3) + 
  geom_rect(aes(ymax = cim_intervals_significant[6,2]+l1t, ymin = cim_intervals_significant[6,2]+l1b, xmin = cim_intervals_significant[6,3], xmax = cim_intervals_significant[6,5]), color = 'black', fill = white[6],size = 0.3) + 
  geom_rect(aes(ymax = cim_intervals_significant[7,2]+l1t, ymin = cim_intervals_significant[7,2]+l1b, xmin = cim_intervals_significant[7,3], xmax = cim_intervals_significant[7,5]), color = 'black', fill = white[7],size = 0.3) + 
  geom_rect(aes(ymax = cim_intervals_significant[8,2]+l1t, ymin = cim_intervals_significant[8,2]+l1b, xmin = cim_intervals_significant[8,3], xmax = cim_intervals_significant[8,5]), color = 'black', fill = white[8],size = 0.3) + 
  geom_rect(aes(ymax = cim_intervals_significant[9,2]+l1t, ymin = cim_intervals_significant[9,2]+l1b, xmin = cim_intervals_significant[9,3], xmax = cim_intervals_significant[9,5]), color = 'black', fill = white[9],size = 0.3) + 
  geom_rect(aes(ymax = cim_intervals_significant[10,2]+l1t, ymin = cim_intervals_significant[10,2]+l1b, xmin = cim_intervals_significant[10,3], xmax = cim_intervals_significant[10,5]), color = 'black', fill = white[10],size = 0.3) + 
  geom_rect(aes(ymax = cim_intervals_significant[11,2]+l1t, ymin = cim_intervals_significant[11,2]+l1b, xmin = cim_intervals_significant[11,3], xmax = cim_intervals_significant[11,5]), color = 'black', fill = white[11],size = 0.3) + 
  geom_rect(aes(ymax = cim_intervals_significant[12,2]+l1t, ymin = cim_intervals_significant[12,2]+l1b, xmin = cim_intervals_significant[12,3], xmax = cim_intervals_significant[12,5]), color = 'black', fill = white[12],size = 0.3) + 
  geom_rect(aes(ymax = cim_intervals_significant[13,2]+l1t, ymin = cim_intervals_significant[13,2]+l1b, xmin = cim_intervals_significant[13,3], xmax = cim_intervals_significant[13,5]), color = 'black', fill = white[13],size = 0.3) + 
  geom_rect(aes(ymax = cim_intervals_significant[14,2]+l1t, ymin = cim_intervals_significant[14,2]+l1b, xmin = cim_intervals_significant[14,3], xmax = cim_intervals_significant[14,5]), color = 'black', fill = white[14],size = 0.3) + 
  geom_rect(aes(ymax = cim_intervals_significant[15,2]+l1t, ymin = cim_intervals_significant[15,2]+l1b, xmin = cim_intervals_significant[15,3], xmax = cim_intervals_significant[15,5]), color = 'black', fill = white[15],size = 0.3) + 
  geom_rect(aes(ymax = cim_intervals_significant[16,2]+l1t, ymin = cim_intervals_significant[16,2]+l1b, xmin = cim_intervals_significant[16,3], xmax = cim_intervals_significant[16,5]), color = 'black', fill = white[16],size = 0.3) + 
  geom_rect(aes(ymax = cim_intervals_significant[17,2]+l1t, ymin = cim_intervals_significant[17,2]+l1b, xmin = cim_intervals_significant[17,3], xmax = cim_intervals_significant[17,5]), color = 'black', fill = white[17],size = 0.3) + 
  geom_rect(aes(ymax = cim_intervals_significant[18,2]+l2t, ymin = cim_intervals_significant[18,2]+l2b, xmin = cim_intervals_significant[18,3], xmax = cim_intervals_significant[18,5]), color = 'black', fill = white[18],size = 0.3) + 
  geom_rect(aes(ymax = cim_intervals_significant[19,2]+l3t, ymin = cim_intervals_significant[19,2]+l3b, xmin = cim_intervals_significant[19,3], xmax = cim_intervals_significant[19,5]), color = 'black', fill = white[19],size = 0.3) + 
  geom_rect(aes(ymax = cim_intervals_significant[20,2]+l4t, ymin = cim_intervals_significant[20,2]+l4b, xmin = cim_intervals_significant[20,3], xmax = cim_intervals_significant[20,5]), color = 'black', fill = white[20],size = 0.3) + 
  
  ## Labels with trait names
  geom_text(aes(label = cim_intervals_significant[1,1], x = cim_intervals_significant[1,4]+0, y = cim_intervals_significant[1,2]+0.5))+ 
  geom_text(aes(label = cim_intervals_significant[2,1], x = cim_intervals_significant[2,4]+-85, y = cim_intervals_significant[2,2]+0.5))+ 
  geom_text(aes(label = cim_intervals_significant[3,1], x = cim_intervals_significant[3,4]+0, y = cim_intervals_significant[3,2]+0.75))+ 
  geom_text(aes(label = cim_intervals_significant[4,1], x = cim_intervals_significant[4,4]+0, y = cim_intervals_significant[4,2]+0.5))+ 
  geom_text(aes(label = cim_intervals_significant[5,1], x = cim_intervals_significant[5,4]+-125, y = cim_intervals_significant[5,2]+0.5))+ 
  geom_text(aes(label = cim_intervals_significant[6,1], x = cim_intervals_significant[6,4]+100, y = cim_intervals_significant[6,2]+0.5))+ 
  geom_text(aes(label = cim_intervals_significant[7,1], x = cim_intervals_significant[7,4]+0, y = cim_intervals_significant[7,2]+0.5))+ 
  geom_text(aes(label = cim_intervals_significant[8,1], x = cim_intervals_significant[8,4]+0, y = cim_intervals_significant[8,2]+0.5))+ 
  geom_text(aes(label = cim_intervals_significant[9,1], x = cim_intervals_significant[9,4]+0, y = cim_intervals_significant[9,2]+0.5))+ 
  geom_text(aes(label = cim_intervals_significant[10,1], x = cim_intervals_significant[10,4]+0, y = cim_intervals_significant[10,2]+0.75))+ 
  geom_text(aes(label = cim_intervals_significant[11,1], x = cim_intervals_significant[11,4]+0, y = cim_intervals_significant[11,2]+0.5))+ 
  geom_text(aes(label = cim_intervals_significant[12,1], x = cim_intervals_significant[12,4]+0, y = cim_intervals_significant[12,2]+0.5))+ 
  geom_text(aes(label = cim_intervals_significant[13,1], x = cim_intervals_significant[13,4]+50, y = cim_intervals_significant[13,2]+0.5))+ 
  geom_text(aes(label = cim_intervals_significant[14,1], x = cim_intervals_significant[14,4]+0, y = cim_intervals_significant[14,2]+0.75))+ 
  geom_text(aes(label = cim_intervals_significant[15,1], x = cim_intervals_significant[15,4]+35, y = cim_intervals_significant[15,2]+0.5))+ 
  geom_text(aes(label = cim_intervals_significant[16,1], x = cim_intervals_significant[16,4]+-100, y = cim_intervals_significant[16,2]+0.5))+ 
  geom_text(aes(label = cim_intervals_significant[17,1], x = cim_intervals_significant[17,4]+-125, y = cim_intervals_significant[17,2]+0.75))+ 
  geom_text(aes(label = cim_intervals_significant[18,1], x = cim_intervals_significant[18,4]+170, y = cim_intervals_significant[18,2]+0.2))+ 
  geom_text(aes(label = cim_intervals_significant[19,1], x = cim_intervals_significant[19,4]+180, y = cim_intervals_significant[19,2]+0.4))+ 
  geom_text(aes(label = cim_intervals_significant[20,1], x = cim_intervals_significant[20,4]+180, y = cim_intervals_significant[20,2]+0.6))+ 
  
  ## Line segments connecting trait names with QTLs
  geom_segment(aes(x = cim_intervals_significant[1,4]+0, xend = cim_intervals_significant[1,4], y = cim_intervals_significant[1,2]+0.5-0.1, yend = cim_intervals_significant[1,2]+l1t))+ 
  geom_segment(aes(x = cim_intervals_significant[2,4]+-85, xend = cim_intervals_significant[2,4], y = cim_intervals_significant[2,2]+0.5-0.1, yend = cim_intervals_significant[2,2]+l1t))+ 
  geom_segment(aes(x = cim_intervals_significant[3,4]+0, xend = cim_intervals_significant[3,4], y = cim_intervals_significant[3,2]+0.75-0.1, yend = cim_intervals_significant[3,2]+l1t))+ 
  geom_segment(aes(x = cim_intervals_significant[4,4]+0, xend = cim_intervals_significant[4,4], y = cim_intervals_significant[4,2]+0.5-0.1, yend = cim_intervals_significant[4,2]+l1t))+ 
  geom_segment(aes(x = cim_intervals_significant[5,4]+-125, xend = cim_intervals_significant[5,4], y = cim_intervals_significant[5,2]+0.5-0.1, yend = cim_intervals_significant[5,2]+l1t))+ 
  geom_segment(aes(x = cim_intervals_significant[6,4]+100, xend = cim_intervals_significant[6,4], y = cim_intervals_significant[6,2]+0.5-0.1, yend = cim_intervals_significant[6,2]+l1t))+ 
  geom_segment(aes(x = cim_intervals_significant[7,4]+0, xend = cim_intervals_significant[7,4], y = cim_intervals_significant[7,2]+0.5-0.1, yend = cim_intervals_significant[7,2]+l1t))+ 
  geom_segment(aes(x = cim_intervals_significant[8,4]+0, xend = cim_intervals_significant[8,4], y = cim_intervals_significant[8,2]+0.5-0.1, yend = cim_intervals_significant[8,2]+l1t))+ 
  geom_segment(aes(x = cim_intervals_significant[9,4]+0, xend = cim_intervals_significant[9,4], y = cim_intervals_significant[9,2]+0.5-0.1, yend = cim_intervals_significant[9,2]+l1t))+ 
  geom_segment(aes(x = cim_intervals_significant[10,4]+0, xend = cim_intervals_significant[10,4], y = cim_intervals_significant[10,2]+0.75-0.1, yend = cim_intervals_significant[10,2]+l1t))+ 
  geom_segment(aes(x = cim_intervals_significant[11,4]+0, xend = cim_intervals_significant[11,4], y = cim_intervals_significant[11,2]+0.5-0.1, yend = cim_intervals_significant[11,2]+l1t))+ 
  geom_segment(aes(x = cim_intervals_significant[12,4]+0, xend = cim_intervals_significant[12,4], y = cim_intervals_significant[12,2]+0.5-0.1, yend = cim_intervals_significant[12,2]+l1t))+ 
  geom_segment(aes(x = cim_intervals_significant[13,4]+50, xend = cim_intervals_significant[13,4], y = cim_intervals_significant[13,2]+0.5-0.1, yend = cim_intervals_significant[13,2]+l1t))+ 
  geom_segment(aes(x = cim_intervals_significant[14,4]+0, xend = cim_intervals_significant[14,4], y = cim_intervals_significant[14,2]+0.75-0.1, yend = cim_intervals_significant[14,2]+l1t))+ 
  geom_segment(aes(x = cim_intervals_significant[15,4]+35, xend = cim_intervals_significant[15,4], y = cim_intervals_significant[15,2]+0.5-0.1, yend = cim_intervals_significant[15,2]+l1t))+ 
  geom_segment(aes(x = cim_intervals_significant[16,4]+-100, xend = cim_intervals_significant[16,4], y = cim_intervals_significant[16,2]+0.5-0.1, yend = cim_intervals_significant[16,2]+l1t))+ 
  geom_segment(aes(x = cim_intervals_significant[17,4]+-125, xend = cim_intervals_significant[17,4], y = cim_intervals_significant[17,2]+0.75-0.1, yend = cim_intervals_significant[17,2]+l1t))+ 
  geom_segment(aes(x = cim_intervals_significant[18,4]+60, xend = cim_intervals_significant[18,4]+12, y = cim_intervals_significant[18,2]+(l1t+l1b)/2, yend = cim_intervals_significant[18,2]+l2t))+ 
  geom_segment(aes(x = cim_intervals_significant[19,4]+40, xend = cim_intervals_significant[19,4]+12, y = cim_intervals_significant[19,2]+(l2t+l2b)/2, yend = cim_intervals_significant[19,2]+l3t))+ 
  geom_segment(aes(x = cim_intervals_significant[20,4]+60, xend = cim_intervals_significant[20,4]+12, y = cim_intervals_significant[20,2]+l3t, yend = cim_intervals_significant[20,2]+l4t))


plot


################### Extra functions ############################################

## Manually defined positions for labels, adjusted to avoid overlaps
label_x <- c(0,-85,0,0,-125,+100,0,0,0,0,0,0,50,0,35,-100,-125,50+120,60+120,60+120)
label_y <- c(0.5,0.5,0.75,0.5,0.5,0.5,0.5,0.5,0.5,0.75,0.5,0.5,0.5,0.75,0.5,0.5,0.75,0.2,0.4,0.6)

line_x <- c(0,-85,0,0,-125,+100,0,0,0,0,0,0,50,0,35,-100,-125,50,60,60)

## Generates QTL interval boxes (geom_rect)
for(i in 1:nrow(cim_intervals_significant)){
  newLine <- paste0("geom_rect(aes(ymax = cim_intervals_significant[",i,",2]+l1t, ymin = cim_intervals_significant[",i,",2]+l1b, xmin = cim_intervals_significant[",i,",3], xmax = cim_intervals_significant[",i,",5]), color = 'black', fill = white[",i,"],size = 0.3) +")
  cat(newLine, "\n")
}

## Generates labels for QTLs (geom_text)
for(i in 1:nrow(cim_intervals_significant)){
  newLine <- paste0("geom_text(aes(label = cim_intervals_significant[",i,",1], x = cim_intervals_significant[",i,",4]+",label_x[i],", y = cim_intervals_significant[",i,",2]+",label_y[i],"))+")
  cat(newLine, "\n")
}

## Generates line segments connecting labels to QTLs (geom_segment)
for(i in 1:nrow(cim_intervals_significant)){
  newLine <- paste0("geom_segment(aes(x = cim_intervals_significant[",i,",4]+",label_x[i],", xend = cim_intervals_significant[",i,",4], y = cim_intervals_significant[",i,",2]+",label_y[i],"-0.1, yend = cim_intervals_significant[",i,",2]+l1t))+")
  cat(newLine, "\n")
}


#################### Saving as a high-quality image file #######################

## Optimized for 5000x3500, reduced to 4800x3360 to fit within 50Mb size req. ##
tiff(filename = "QTL_intervals.tif", res=500, width = 4800, height = 3360)
plot
dev.off()

png(filename = "Intervals.png", res=500, width = 4800, height = 3360)
plot
dev.off()


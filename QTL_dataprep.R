library(ggplot2)
library(qtl)
library(dplyr)

## final gmap -- removed Petal color, anther color, filament color, and solitary flower number ##################
# This has 25 alphabetized columns--the 24 trait columns plus one ID column
#gmap <- read.cross(format = "csv", file = "D:/backups/Pictures/Research (Whipple Lab)/My_paper/gmap5335_final.csv", map.function = "c-f", genotypes = c("AA", "AB", "BB"), alleles = c("A", "B"))
gmap <- read.cross(format = "csv", file = "gmap5335_final.csv", map.function = "c-f", genotypes = c("AA", "AB", "BB"), alleles = c("A", "B"))

#################################################################################################################

## Transforming Free filament length, Pedicel length, Internode
gmap$pheno$Free.filament <- log(gmap$pheno$Free.filament)
gmap$pheno$Pedicel.length <- log(gmap$pheno$Pedicel.length)
gmap$pheno$Internode <- log(gmap$pheno$Internode)


str(gmap$pheno)

# $ Anther.length        : num  0.7 0.616 0.792 0.624 0.668 ...
# $ Anther.width         : num  0.356 0.408 0.688 0.358 0.52 0.35 0.308 0.404 NA 0.602 ...
# $ Calyx.length         : num  2.9 3.63 3.67 2.69 3.56 ...
# $ Calyx.sinus          : num  1.57 2.09 1.93 1.07 2.21 ...
# $ Calyx.tooth          : num  1.23 1.49 1.66 1.32 1.36 ...
# $ Corolla.length       : num  6.91 7.28 9.11 6.59 8.76 7.33 4.91 6.62 NA 6.15 ...
# $ Death.early          : num  0 0 0 0 0 0 0 0 0 0 ...
# $ Floral.buds          : num  57 60 70 70 66 53 63 57 57 70 ...
# $ Free.filament        : num  0.238 0.353 1.894 0.338 1.3 ...
# $ Individual           : num  1 2 3 4 5 6 7 8 9 10 ...
# $ Internode            : num  0.5 2 5.5 NA 0.5 1 0.2 2 0.5 1 ...
# $ Lobe.length          : num  4.48 4.78 5.86 3.67 5.02 ...
# $ Lobe.width           : num  1.75 1.9 2.29 1.22 1.45 ...
# $ Main.flower          : num  68 75 99 94 84 75 75 75 68 88 ...
# $ Major.length.diameter: num  24 23.7 26 16 22 19.5 21.5 26 23 20 ...
# $ Midrib               : num  0.336 0.47 0.45 0.352 0.364 0.353 0.426 0.398 NA 0.392 ...
# $ Ovary                : num  0.883 NA 0.811 0.881 0.792 ...
# $ Pedicel.length       : num  1 1.5 3.1 NA 1 1 0.5 1.5 1 1.5 ...
# $ Stamen.length        : num  3.83 4.65 6.72 4.11 6.2 ...
# $ Sterility            : num  1 0 0 0 1 0 1 1 0 0 ...
# $ Stigma.length        : num  0.525 0.76 0.825 0.705 0.795 0.53 0.695 0.77 NA 0.815 ...
# $ Style.length         : num  4.64 5.04 7.26 4.23 5.97 4.55 4.59 5.08 NA 4.27 ...
# $ Throat               : num  1.308 2.013 2.063 0.963 1.585 ...
# $ Tube.length          : num  2.45 2.53 2.97 2.94 3.49 2.99 1.92 2.49 NA 1.87 ...
# $ Tube.width           : num  2.98 2.5 3.07 2.57 2.56 2.75 2.54 2.9 NA 2.7 ...

### Custom color palette (colors I can distinguish as a colorblind person) ##############
#colors <- c("#aaffc3","#808000","#e6194B","#f58231","#ffe119","#bfef45","#3cb44b","#9A6324","#469990","#fabed4","#4363d8","#000075","#911eb4","#f032e6")


#c17 <- c(
#  "dodgerblue2", "#E31A1C", # red
#  "green4",
#  "#6A3D9A", # purple
#  "#FF7F00", # orange
#  "black", "gold1",
#  "orchid1", "#FDBF6F", "blue1", "steelblue4",
#  "darkturquoise", "green1", "yellow4", "yellow3",
#  "darkorange4", "brown"
#)

blues <- c("#003396","#3373C4","#86CEFA")
reds <- c("#c50000","#ea5252","#f9aeae")
greens <- c("#063b00","#0a5d00","#089000","#00E227")
tans <- c("#664229","#987554","#d3b95f")
greys <- c("#212121","#575757")
yellow <- "#fbfb00"
orange <- "#B2560D"

#c17 <- c(blues,reds,greens,tans,greys,yellow,orange) # 1
#c17 <- c(orange,yellow,greys,tans,greens,reds,blues) # 2
#c17 <- c(orange,blues,yellow,reds,greens,greys,tans) # 3
#c17 <- c(tans,blues,yellow,greens,greys,reds,orange) # 4

#c17 <- c(blues,greens,tans,reds,yellow,orange,greys) # 5
#c17 <- c(greys,greens,tans,orange,yellow,reds,blues) # 6
#c17 <- c(reds,blues,greens,greys,tans,orange,yellow) # 7
c17 <- c(greens,blues,tans,reds,greys,yellow,orange) # 8

c18 <- c(greens,blues,tans,reds,greys,yellow,orange,"#FCFBF4")



#pie(rep(1, 17), col = c17)


##This loop gets the chr, lpos, pos, rpos, lod of the highest hit for each trait and puts them into individual vectors
##    These are then put together into a dataframe

chr <- integer()
pos <- integer()
lod <- integer()
lpos <- integer()
rpos <- integer()

gmap <- calc.genoprob(gmap)

for (i in 1:length(gmap$pheno)){
  out.ehk <- scanone(gmap, pheno.col = i, method = "ehk")
  #plot(out.ehk, main = colnames(gmap$pheno[i]))
  #abline(h = 4.16, lty = 2)
  
  chr[i] <- as.numeric(out.ehk[which.max(out.ehk$lod),1])
  pos[i] <- as.numeric(out.ehk[which.max(out.ehk$lod),2])
  lod[i] <- as.numeric(out.ehk[which.max(out.ehk$lod),3])
  
  int <- lodint(out.ehk, chr = chr[i], qtl.index = pos[i], drop = 1.5)
  lpos[i] <- int$pos[1]
  rpos[i] <- int$pos[3]
  
}

intervals <- data.frame('chr'=chr,'lpos'=lpos,'pos'=pos,'rpos'=rpos,'lod'=lod, row.names = colnames(gmap$pheno))

intervals <- intervals[intervals$lod>4.16,]

## Removing free filament--I don't think there is enough support for this interval (too narrow, doesn't make sense in context of the QTL graph)
## I added this back in after transformation--the hit on Chr 1 went away, and a large group of markers on Chr 6 became significant
#intervals <- intervals[!rownames(intervals)=="Free.filament",]

## Removing Tube.length--very small interval, little support
intervals <- intervals[!rownames(intervals)=="Tube.length",]

####################################
## Getting chromosome lengths

max <- numeric()

for (i in 1:9){
  map <- pull.map(gmap, chr = i)
  max[i] <- max(map[[1]])
}


#### Need to get second QTL intervals from data down below, integrate them into rest of the data




################################Lobe Length
k <- 12
str(gmap$pheno[k])
out.ehk <- scanone(gmap,pheno.col = k, method = "ehk")
sum <- summary(out.ehk)
sum <- sum[sum$lod>4.2,]
sum <- sum[which.min(sum$lod),]
out.ehk <- scanone(gmap, pheno.col = k, chr = sum$chr, method = "ehk")

chr <- numeric()
pos <- numeric()
lod <- numeric()
lpos <- numeric()
rpos <- numeric()

chr[1] <- as.numeric(sum$chr)
pos[1] <- as.numeric(out.ehk[which.max(out.ehk$lod),2])
lod[1] <- as.numeric(out.ehk[which.max(out.ehk$lod),3])

int <- lodint(out.ehk, chr = as.numeric(sum$chr), qtl.index = sum$pos, drop = 1.5)
lpos[1] <- int$pos[1]
rpos[1] <- int$pos[3]

second_qtl <- numeric()
second_qtl <- data.frame('chr'=chr,'lpos'=lpos,'pos'=pos,'rpos'=rpos,'lod'=lod, row.names = paste0(colnames(gmap$pheno[k]),"2"))
intervals <- rbind(intervals,second_qtl)


################################Lobe Width
k <- 13
str(gmap$pheno[k])
out.ehk <- scanone(gmap,pheno.col = k, method = "ehk")
sum <- summary(out.ehk)
sum <- sum[sum$lod>4.16,]
sum <- sum[which.min(sum$lod),]
out.ehk <- scanone(gmap, pheno.col = k, chr = sum$chr, method = "ehk")

chr <- numeric()
pos <- numeric()
lod <- numeric()
lpos <- numeric()
rpos <- numeric()

chr[1] <- as.numeric(sum$chr)
pos[1] <- as.numeric(out.ehk[which.max(out.ehk$lod),2])
lod[1] <- as.numeric(out.ehk[which.max(out.ehk$lod),3])

int <- lodint(out.ehk, chr = as.numeric(sum$chr), qtl.index = sum$pos, drop = 1.5)
lpos[1] <- int$pos[1]
rpos[1] <- int$pos[3]

second_qtl <- numeric()
second_qtl <- data.frame('chr'=chr,'lpos'=lpos,'pos'=pos,'rpos'=rpos,'lod'=lod, row.names = paste0(colnames(gmap$pheno[k]),"2"))
intervals <- rbind(intervals,second_qtl)




################################ Stigma Length
k <- 21
str(gmap$pheno[k])
out.ehk <- scanone(gmap,pheno.col = k, method = "ehk")
sum <- summary(out.ehk)
sum <- sum[sum$lod>4.16,]
sum <- sum[-which.max(sum$lod),]

chr <- numeric()
pos <- numeric()
lod <- numeric()
lpos <- numeric()
rpos <- numeric()
j <- 1

for (i in sum$chr){
  out.ehk <- scanone(gmap, pheno.col = k, chr = i, method = "ehk")


  
  chr[j] <- as.numeric(i)
  pos[j] <- as.numeric(out.ehk[which.max(out.ehk$lod),2])
  lod[j] <- as.numeric(out.ehk[which.max(out.ehk$lod),3])
  
  int <- lodint(out.ehk, chr = as.numeric(i), qtl.index = sum$pos[j], drop = 1.5)
  lpos[j] <- int$pos[1]
  rpos[j] <- int$pos[3]
  
  second_qtl[j,] <- data.frame('chr'=chr[j],'lpos'=lpos[j],'pos'=pos[j],'rpos'=rpos[j],'lod'=lod[j], row.names = paste0(colnames(gmap$pheno[k]),j+1))
  j <- j+1
}

second_qtl <- numeric()
second_qtl <- data.frame('chr'=chr,'lpos'=lpos,'pos'=pos,'rpos'=rpos,'lod'=lod)
rownames(second_qtl) <- c("Stigma.length2", "Stigma.length3", "Stigma.length4")
intervals <- rbind(intervals,second_qtl)

################################ Style Length
k <- 22
str(gmap$pheno[k])
out.ehk <- scanone(gmap,pheno.col = k, method = "ehk")
sum <- summary(out.ehk)
sum <- sum[sum$lod>4.16,]
sum <- sum[which.min(sum$lod),]
out.ehk <- scanone(gmap, pheno.col = k, chr = sum$chr, method = "ehk")

chr <- numeric()
pos <- numeric()
lod <- numeric()
lpos <- numeric()
rpos <- numeric()

chr[1] <- as.numeric(sum$chr)
pos[1] <- as.numeric(out.ehk[which.max(out.ehk$lod),2])
lod[1] <- as.numeric(out.ehk[which.max(out.ehk$lod),3])

int <- lodint(out.ehk, chr = as.numeric(sum$chr), qtl.index = sum$pos, drop = 1.5)
lpos[1] <- int$pos[1]
rpos[1] <- int$pos[3]

second_qtl <- numeric()
second_qtl <- data.frame('chr'=chr,'lpos'=lpos,'pos'=pos,'rpos'=rpos,'lod'=lod, row.names = paste0(colnames(gmap$pheno[k]),"2"))
intervals <- rbind(intervals,second_qtl)




################################ Moving Free Filament row to bottom

# Identify the row index of "Free.filament" using row names
row_index <- which(rownames(intervals) == "Free.filament")

# Store the row in a separate variable
row_to_move <- intervals[row_index, ]

# Remove the selected row from the original dataframe
intervals <- intervals[-row_index, ]

# Append the row to the end of the dataframe
intervals <- rbind(intervals, row_to_move)

# Print the modified dataframe
print(intervals)




############ Naming colors in c18 ###############################

names(c18) <- c(rownames(intervals[1:17,]),"Free.filament")


intervals$trait <- as.factor(rownames(intervals))
intervals$trait[18] <- "Lobe.length"
intervals$trait[19] <- "Lobe.width"
intervals$trait[20] <- "Stigma.length"
intervals$trait[21] <- "Stigma.length"
intervals$trait[22] <- "Stigma.length"
intervals$trait[23] <- "Style.length"



## Dropping Sterility
gmap$pheno <- subset(gmap$pheno, select = -Sterility)
c17 <- c18[-14]

## Updating names
#names(c17) <- c("Anther length","Anther width","Sepal length","Sepal tooth length","Petal length","Internode length","Petal lobe length","Petal lobe width","Vegetative rosette diameter","Sepal midrib width","Ovary shape","Pedicel length","Filament length","Stigma length","Style length","Petal throat length", "Free filament length")
names(c17) <- c("Petal length","Petal lobe length","Petal lobe width","Petal throat length","Filament length","Free filament length","Anther length","Anther width","Style length","Stigma length","Ovary shape","Sepal length","Sepal tooth length","Sepal midrib width","Pedicel length","Internode length","Vegetative rosette diameter")



intervals$trait <- c("Anther length","Anther width","Sepal length","Sepal tooth length","Petal length","Internode length","Petal lobe length","Petal lobe width","Vegetative rosette diameter","Sepal midrib width","Ovary shape","Pedicel length","Filament length","Sterility","Stigma length","Style length","Petal throat length","Petal lobe length","Petal lobe width","Stigma length","Stigma length","Stigma length","Style length", "Free filament length")





############ Naming colors in c17 ###############################
#
#names(c17) <- rownames(intervals[1:17,])
#
#
#intervals$trait <- as.factor(rownames(intervals))
#intervals$trait[18] <- "Lobe.length"
#intervals$trait[19] <- "Lobe.width"
#intervals$trait[20] <- "Stigma.length"
#intervals$trait[21] <- "Stigma.length"
#intervals$trait[22] <- "Stigma.length"
#intervals$trait[23] <- "Style.length"
#
#
#
## Dropping Sterility
#gmap$pheno <- subset(gmap$pheno, select = -Sterility)
#c16 <- c17[-14]
#
## Updating names
#names(c16) <- c("Anther length","Anther width","Sepal length","Sepal tooth length","Petal length","Internode length","Petal lobe length","Petal lobe width","Vegetative rosette diameter","Sepal midrib width","Ovary shape","Pedicel length","Filament length","Stigma length","Style length","Throat length")
#
#intervals$trait <- c("Anther length","Anther width","Sepal length","Sepal tooth length","Petal length","Internode length","Petal lobe length","Petal lobe width","Vegetative rosette diameter","Sepal midrib width","Ovary shape","Pedicel length","Filament length","Sterility","Stigma length","Style length","Throat length","Petal lobe length","Petal lobe width","Stigma length","Stigma length","Stigma length","Style length")

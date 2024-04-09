library(qtl)
library(dplyr)

gmap <- read.cross(format = "csv", file = "gmap5335_final.csv", map.function = "c-f", genotypes = c("AA", "AB", "BB"), alleles = c("A", "B"))
#################################################################################################################

## Transforming Free filament length, Pedicel length, Internode
gmap$pheno$Free.filament <- log(gmap$pheno$Free.filament)
gmap$pheno$Pedicel.length <- log(gmap$pheno$Pedicel.length)
gmap$pheno$Internode <- log(gmap$pheno$Internode)

gmap <- calc.genoprob(gmap)

gmap_qtl <- gmap_all
names(gmap_qtl$pheno) <- trait_names
gmap_qtl$pheno <- gmap_qtl$pheno %>%
  select(-c(Sterility,Individual))

# Updating trait names
trait_names <- c("Petal length", "Tube length", "Tube width", "Throat length", "Petal lobe length",
                 "Petal lobe width", "Filament length", "Free filament length", "Anther length", 
                 "Anther width", "Sepal length", "Sepal sinus length", "Sepal tooth length",
                 "Midrib width", "Ovary", "Style length", "Stigma length", "Death early",
                 "Days to flower", "Main flower date", "Solitary flower number", "Pedicel length",
                 "Internode", "Petal color", "Anther color", "Filament color",
                 "Vegetative rosette diameter", "Sterility", "Individual")

colnames <- c("Anther length","Anther width","Sepal length","Sepal sinus length","Sepal tooth length",
              "Petal length","Days to flower","Free filament length","Internode length","Petal lobe length",
              "Petal lobe width","Main.flower","Vegetative rosette diameter","Sepal midrib width","Ovary",
              "Pedicel length","Filament length","Stigma length","Style length","Throat length","Petal tube length",
              "Petal tube width","group")

new_order <- c("Petal length","Petal lobe length","Petal lobe width",
               "Tube length","Tube width","Throat length",
               "Filament length","Free filament length","Anther length",
               "Anther width","Style length","Stigma length","Ovary",
               "Sepal length","Sepal sinus length","Sepal tooth length",
               "Midrib width","Pedicel length","Internode",
               "Days to flower","Vegetative rosette diameter")

gmap_qtl$pheno <- gmap_qtl$pheno %>%
  select(new_order)

cim_thresholds_all_named_updated <- cim_thresholds_all_named
names(cim_thresholds_all_named_updated) <- trait_names

cim_thresholds_graph <- cim_thresholds_all_named_updated[new_order]

for (i in 1:length(gmap_qtl$pheno)){
  cim <- cim(gmap_qtl, pheno.col=i, method = "ehk", window = 10, n.marcovar=3)
  scan1 <- scanone(gmap_qtl, pheno.col = i, method = "ehk")
  
  png(filename = paste0("../CIM_graphs/","QTL_graphs_",i,".png"), res=500, width = 5000, height = 4000)
  plot(cim, col = "black")
  plot(scan1, add=T, col="lightgrey")
  plot(cim, add=T, col = "black")
  title(colnames(gmap_test_log$pheno)[i])
  abline(h = summary(cim_thresholds_graph[[i]])[1], lty = 2, col = "lightblue", lwd=3)
  abline(h = summary(cim_thresholds_graph[[i]])[2], lty = 2, col = "blue", lwd=3)
  abline(h = 4.16, lty = 2, col = "lightgrey", lwd=3)
  dev.off
}

## Save as 1000w x 800h?

par(mfrow = c(2, 2))

for (i in 1:5){
  index <- i*4
  
  specific_index <- index - 3
  cim1 <- cim(gmap_qtl, pheno.col=specific_index, method = "ehk", window = 10, n.marcovar=3)
  scan1 <- scanone(gmap_qtl, pheno.col = specific_index, method = "ehk")
  
  specific_index <- index - 2
  cim2 <- cim(gmap_qtl, pheno.col=specific_index, method = "ehk", window = 10, n.marcovar=3)
  scan2 <- scanone(gmap_qtl, pheno.col = specific_index, method = "ehk")
  
  specific_index <- index - 1
  cim3 <- cim(gmap_qtl, pheno.col=specific_index, method = "ehk", window = 10, n.marcovar=3)
  scan3 <- scanone(gmap_qtl, pheno.col = specific_index, method = "ehk")
  
  specific_index <- index - 0
  cim4 <- cim(gmap_qtl, pheno.col=specific_index, method = "ehk", window = 10, n.marcovar=3)
  scan4 <- scanone(gmap_qtl, pheno.col = specific_index, method = "ehk")
  
  png(filename = paste0("../CIM_graphs/","QTL_graphs_",i,".png"), res=500, width = 5000, height = 4000)
  par(mfrow = c(2, 2))
  
  plot(cim1, col = "black")
  plot(scan1, add=T, col="lightgrey")
  plot(cim1, add=T, col = "black")
  title(colnames(gmap_qtl$pheno)[index - 3])
  abline(h = summary(cim_thresholds_graph[[index - 3]])[1], lty = 2, col = "lightblue", lwd=3)
  abline(h = summary(cim_thresholds_graph[[index - 3]])[2], lty = 2, col = "blue", lwd=3)
  abline(h = 4.16, lty = 2, col = "lightgrey", lwd=3)
  
  plot(cim2, col = "black")
  plot(scan2, add=T, col="lightgrey")
  plot(cim2, add=T, col = "black")
  title(colnames(gmap_qtl$pheno)[index - 2])
  abline(h = summary(cim_thresholds_graph[[index - 2]])[1], lty = 2, col = "lightblue", lwd=3)
  abline(h = summary(cim_thresholds_graph[[index - 2]])[2], lty = 2, col = "blue", lwd=3)
  abline(h = 4.16, lty = 2, col = "lightgrey", lwd=3)
  
  plot(cim3, col = "black")
  plot(scan3, add=T, col="lightgrey")
  plot(cim3, add=T, col = "black")
  title(colnames(gmap_qtl$pheno)[index - 1])
  abline(h = summary(cim_thresholds_graph[[index - 1]])[1], lty = 2, col = "lightblue", lwd=3)
  abline(h = summary(cim_thresholds_graph[[index - 1]])[2], lty = 2, col = "blue", lwd=3)
  abline(h = 4.16, lty = 2, col = "lightgrey", lwd=3)
  
  plot(cim4, col = "black")
  plot(scan4, add=T, col="lightgrey")
  plot(cim4, add=T, col = "black")
  title(colnames(gmap_qtl$pheno)[index])
  abline(h = summary(cim_thresholds_graph[[index]])[1], lty = 2, col = "lightblue", lwd=3)
  abline(h = summary(cim_thresholds_graph[[index]])[2], lty = 2, col = "blue", lwd=3)
  abline(h = 4.16, lty = 2, col = "lightgrey", lwd=3)
  
  dev.off()
}

png(filename = paste0("../CIM_graphs/","QTL_graphs_",6,".png"), res=500, width = 5000, height = 4000)
par(mfrow = c(2, 2))

cim1 <- cim(gmap_qtl, pheno.col=21, method = "ehk", window = 10, n.marcovar=3)
scan1 <- scanone(gmap_qtl, pheno.col = 21, method = "ehk")

plot(cim1, col = "black")
plot(scan1, add=T, col="lightgrey")
plot(cim1, add=T, col = "black")
title(colnames(gmap_qtl$pheno)[21])
abline(h = summary(cim_thresholds_graph[[21]])[1], lty = 2, col = "lightblue", lwd=3)
abline(h = summary(cim_thresholds_graph[[21]])[2], lty = 2, col = "blue", lwd=3)
abline(h = 4.16, lty = 2, col = "lightgrey", lwd=3)

plot.new()
plot.new()
plot.new()
dev.off()




out.ehk <- scanone(gmap_qtl, pheno.col = 4, method = "ehk")
plot(out.ehk, main = colnames(gmap_qtl$pheno)[4], cex.main = 2, cex.lab = 1.5)
abline(h = 4.16, col = "blue", lty = 2)













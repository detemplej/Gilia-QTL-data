library(qtl)
library(dplyr)

gmap <- read.cross(format = "csv", file = "gmap5335_final.csv", map.function = "c-f", genotypes = c("AA", "AB", "BB"), alleles = c("A", "B"))
#################################################################################################################

## Transforming Free filament length, Pedicel length, Internode
gmap$pheno$Free.filament <- log(gmap$pheno$Free.filament)
gmap$pheno$Pedicel.length <- log(gmap$pheno$Pedicel.length)
gmap$pheno$Internode <- log(gmap$pheno$Internode)

gmap <- calc.genoprob(gmap)

gmap_qtl <- gmap
gmap_qtl$pheno <- gmap$pheno %>%
  select(-c(Main.flower,Sterility,Death.early,Individual))

# Updating trait names
colnames <- c("Anther length","Anther width","Sepal length","Sepal sinus length","Sepal tooth length",
              "Petal length","Days to flower","Free filament length","Internode length","Petal lobe length",
              "Petal lobe width","Vegetative rosette diameter","Sepal midrib width","Ovary shape",
              "Pedicel length","Filament length","Stigma length","Style length","Throat length","Petal tube length",
              "Petal tube width")
colnames(gmap_qtl$pheno) <- colnames

# Setting up the new order, matching to other graphs
new_order <- c("Petal length","Petal lobe length","Petal lobe width",
               "Petal tube length","Petal tube width","Throat length",
               "Filament length","Free filament length","Anther length",
               "Anther width","Style length","Stigma length","Ovary shape",
               "Sepal length","Sepal sinus length","Sepal tooth length",
               "Sepal midrib width","Pedicel length","Internode length",
               "Days to flower","Vegetative rosette diameter","group")
new_order <- new_order[-22]

gmap_qtl$pheno <- gmap_qtl$pheno %>%
  select(new_order)

for (i in 1:length(gmap_qtl$pheno)){
  out.ehk <- scanone(gmap_qtl, pheno.col = i, method = "ehk")
  jpeg(filename = paste0("QTL_graphs/",colnames(gmap_qtl$pheno)[i],".jpg"),width = 1000, height = 800)
  plot(out.ehk, main = colnames(gmap_qtl$pheno)[i], cex.main = 2, cex.lab = 1.5)
  abline(h = 4.16, col = "blue", lty = 2)
  dev.off()
}

## Save as 1000w x 800h?

par(mfrow = c(2, 2))

for (i in 1:5){
  index <- i*4
  
  specific_index <- index - 3
  out.ehk1 <- scanone(gmap_qtl, pheno.col = specific_index, method = "ehk")
  
  specific_index <- index - 2
  out.ehk2 <- scanone(gmap_qtl, pheno.col = specific_index, method = "ehk")
  
  specific_index <- index - 1
  out.ehk3 <- scanone(gmap_qtl, pheno.col = specific_index, method = "ehk")
  
  specific_index <- index - 0
  out.ehk4 <- scanone(gmap_qtl, pheno.col = specific_index, method = "ehk")
  
  jpeg(filename = paste0("QTL_graphs/","QTL_graphs_",i,".jpg"),width = 1000, height = 800)
  par(mfrow = c(2, 2))
  plot(out.ehk1, main = colnames(gmap_qtl$pheno)[specific_index-3], cex.main = 2, cex.lab = 1.5)
  abline(h = 4.16, col = "blue", lty = 2)
  plot(out.ehk2, main = colnames(gmap_qtl$pheno)[specific_index-2], cex.main = 2, cex.lab = 1.5)
  abline(h = 4.16, col = "blue", lty = 2)
  plot(out.ehk3, main = colnames(gmap_qtl$pheno)[specific_index-1], cex.main = 2, cex.lab = 1.5)
  abline(h = 4.16, col = "blue", lty = 2)
  plot(out.ehk4, main = colnames(gmap_qtl$pheno)[specific_index-0], cex.main = 2, cex.lab = 1.5)
  abline(h = 4.16, col = "blue", lty = 2)
  dev.off()
}

jpeg(filename = paste0("QTL_graphs/","QTL_graphs_",6,".jpg"),width = 1000, height = 800)
par(mfrow = c(2, 2))
plot(out.ehk1, main = colnames(gmap_qtl$pheno)[21], cex.main = 2, cex.lab = 1.5)
abline(h = 4.16, col = "blue", lty = 2)
plot.new()
plot.new()
plot.new()
dev.off()




out.ehk <- scanone(gmap_qtl, pheno.col = 4, method = "ehk")
plot(out.ehk, main = colnames(gmap_qtl$pheno)[4], cex.main = 2, cex.lab = 1.5)
abline(h = 4.16, col = "blue", lty = 2)













## Inherits objects from "QTL_dataprep.R" script (must be run first)

gmap_qtl <- gmap_all

new_order <- c("Petal Length","Petal Lobe Length","Petal Lobe Width",
               "Tube Length","Tube Width","Throat Length",
               "Filament Length","Free Filament Length","Anther Length",
               "Anther Width","Style Length","Stigma Length","Ovary",
               "Sepal Length","Sepal Sinus Length","Sepal Tooth Length",
               "Midrib Width","Pedicel Length","Internode",
               "Days To Flower","Vegetative Rosette Diameter")

## Selecting only the relevant phenotype columns for QTL plots
gmap_qtl$pheno <- gmap_qtl$pheno %>%
  select(new_order)

## Selecting thresholds (20,000 permutations/trait) for relevant phenotype columns
load("CIM_thresholds_all_20000.RData")
names(cim_thresholds_all) <- trait_names
cim_thresholds_graph <- cim_thresholds_all[new_order]


## Graphing QTL traces (within RStudio)
par(mfrow = c(2, 2))

for (i in 1:length(gmap_qtl$pheno)){
  cim <- cim(gmap_qtl, pheno.col=i, method = "ehk", window = 10, n.marcovar=3)
  scan1 <- scanone(gmap_qtl, pheno.col = i, method = "ehk")
  
#  png(filename = paste0("../CIM_graphs/","QTL_graphs_",i,".png"), res=500, width = 5000, height = 4000)
  plot(cim, col = "black")
  plot(scan1, add=T, col="lightgrey")
  plot(cim, add=T, col = "black")
  title(colnames(gmap_qtl$pheno)[i])
  abline(h = summary(cim_thresholds_graph[[i]])[1], lty = 2, col = "lightblue", lwd=3)
  abline(h = summary(cim_thresholds_graph[[i]])[2], lty = 2, col = "blue", lwd=3)
  abline(h = 4.16, lty = 2, col = "lightgrey", lwd=3)
#  dev.off
}


## Graphing QTL traces (for output files)
## More complicated loop so that it can automatically create and name files
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
  
  png(filename = paste0("QTL_graphs_",i,".png"), res=500, width = 5000, height = 4000)
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

png(filename = paste0("QTL_graphs_",6,".png"), res=500, width = 5000, height = 4000)
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


library(corrplot)

## Redefining new order, some manual manipulation involved after looking at graph
#new_order <- c(
#  "Tube.width",
#  "Anther.length",
#  "Anther.width",
#  "Lobe.width",
#  "Style.length",
#  "Stamen.length",
#  "Throat",
#  "Corolla.length",
#  "Lobe.length",
#  "Tube.length",
#  "Free.filament",
  
#  "Midrib",
#  "Calyx.length",
#  "Calyx.tooth",
#  "Calyx.sinus",
#  "Internode",
#  "Pedicel.length",
  
#  "Ovary",
#  "Stigma.length"
#)

new_order <- c("Petal length","Petal lobe length","Petal lobe width",
               "Petal tube length","Petal tube width","Throat length",
               "Filament length","Free filament length","Anther length",
               "Anther width","Style length","Stigma length","Ovary shape",
               "Sepal length","Sepal sinus length","Sepal tooth length",
               "Sepal midrib width","Pedicel length","Internode length",
               "Days to flower","Vegetative rosette diameter","group")
new_order <- new_order[-22]

## Get data from QTL_interval_dataprep_simplified.R and Figure_hists.R
f2_corr <- f2 %>% 
  select(any_of(new_order))

f2_corr_na <- f2_corr[rowSums(is.na(f2_corr))==0,]



## Creating corr object and then running corrplot
# Floral Traits
f2_corr_na <- cor(f2_corr_na[new_order])
corrplot(f2_corr_na, type = "upper", method = "circle", order = "original", diag = F, addgrid.col = 'white', tl.col = "black", tl.cex = 0.8, tl.srt = 60)

# 6/19/2023 output to a png file with dimensions 925w x 765h



library(qtl)
library(tidyverse)

## F2 cross object
## Contains 29 phenotypes, including a column that has 'individual' numbers and
## the sterility phenotype.
gmap_all <- read.cross(format = "csv", file = "gmap_all.csv", map.function = "c-f", genotypes = c("AA", "AB", "BB"), alleles = c("A", "B"))
  

#################################################################################################################

trait_names <- c("Petal Length", "Tube Length", "Tube Width", "Throat Length", "Petal Lobe Length",
                 "Petal Lobe Width", "Filament Length", "Free Filament Length*", "Anther Length", 
                 "Anther Width", "Sepal Length", "Sepal Sinus Length", "Sepal Tooth Length",
                 "Midrib Width", "Ovary", "Style Length", "Stigma Length", "Death early",
                 "Floral buds", "Main flower date", "Solitary flower number", "Pedicel Length*",
                 "Internode*", "Petal color", "Anther color", "Filament color",
                 "Vegetative Rosette Diameter", "Sterility", "Individual")

names(gmap_all$pheno) <- trait_names

str(gmap_all$pheno)

## Transforming Pedicel length
## After mapping with transformed and untransformed trait values for Pedicel length,
## Free filament length, and Internode length (the most severely skewed distributions),
## we found significant QTLs for untransformed Free filament length and Internode
## length traits, whereas we found a significant QTL for transformed Pedicel length.
## Here, we just transform Pedicel length to reflect this result.

gmap_all$pheno$`Pedicel Length*` <- log(gmap_all$pheno$`Pedicel Length`)

################# CIM Interval Table ###########################################

## Setting up a blank data.frame to store interval values

cim_intervals <- data.frame(trait = character(0),
                            chromosome = numeric(0),
                            lbound = numeric(0),
                            peak = numeric(0),
                            rbound = numeric(0),
                            lod = numeric(0),
                            Sig0.1 = logical(0),
                            Sig0.05 = logical(0))

## A loop to find the interval with the highest LOD score for each trait. For traits
## with multiple significant intervals, these are added to the data.frame in the code
## below the loop.

for (i in 1:29){
  cim <- cim(gmap_all, pheno.col=i, method = "ehk", window = 10, n.marcovar=3)
  
  highest_peak_marker <- which.max(cim$lod)
  highest_peak_chromosome <- cim$chr[highest_peak_marker]
  
  int <- lodint(cim, chr = highest_peak_chromosome)
  
  cim_intervals[i,1] <- names(gmap_all$pheno[i])
  cim_intervals[i,2] <- int[1,1]
  cim_intervals[i,3] <- int[1,2]
  cim_intervals[i,4] <- int[2,2]
  cim_intervals[i,5] <- int[3,2]
  cim_intervals[i,6] <- int[2,3]
  cim_intervals[i,7] <- int[2,3] > summary(cim_thresholds_all[[i]])[2]
  cim_intervals[i,8] <- int[2,3] > summary(cim_thresholds_all[[i]])[1]
}

## 2nd interval for Petal lobe length
cim <- cim(gmap_all, pheno.col=5, method = "ehk", window = 10, n.marcovar=3)
int <- lodint(cim, chr = 1)
cim_intervals[nrow(cim_intervals)+1,1] <- names(gmap_all$pheno[5])
cim_intervals[nrow(cim_intervals),2] <- int[1,1]
cim_intervals[nrow(cim_intervals),3] <- int[1,2]
cim_intervals[nrow(cim_intervals),4] <- int[2,2]
cim_intervals[nrow(cim_intervals),5] <- int[3,2]
cim_intervals[nrow(cim_intervals),6] <- int[2,3]
cim_intervals[nrow(cim_intervals),7] <- int[2,3] > summary(cim_thresholds_all[[5]])[2]
cim_intervals[nrow(cim_intervals),8] <- int[2,3] > summary(cim_thresholds_all[[5]])[1]

## 2nd interval for Petal lobe width
cim <- cim(gmap_all, pheno.col=6, method = "ehk", window = 10, n.marcovar=3)
int <- lodint(cim, chr = 6)
cim_intervals[nrow(cim_intervals)+1,1] <- names(gmap_all$pheno[6])
cim_intervals[nrow(cim_intervals),2] <- int[1,1]
cim_intervals[nrow(cim_intervals),3] <- int[1,2]
cim_intervals[nrow(cim_intervals),4] <- int[2,2]
cim_intervals[nrow(cim_intervals),5] <- int[3,2]
cim_intervals[nrow(cim_intervals),6] <- int[2,3]
cim_intervals[nrow(cim_intervals),7] <- int[2,3] > summary(cim_thresholds_all[[6]])[2]
cim_intervals[nrow(cim_intervals),8] <- int[2,3] > summary(cim_thresholds_all[[6]])[1]

## 2nd interval for Anther width
cim <- cim(gmap_all, pheno.col=10, method = "ehk", window = 10, n.marcovar=3)
int <- lodint(cim, chr = 2)
cim_intervals[nrow(cim_intervals)+1,1] <- names(gmap_all$pheno[10])
cim_intervals[nrow(cim_intervals),2] <- int[1,1]
cim_intervals[nrow(cim_intervals),3] <- int[1,2]
cim_intervals[nrow(cim_intervals),4] <- int[2,2]
cim_intervals[nrow(cim_intervals),5] <- int[3,2]
cim_intervals[nrow(cim_intervals),6] <- int[2,3]
cim_intervals[nrow(cim_intervals),7] <- int[2,3] > summary(cim_thresholds_all[[10]])[2]
cim_intervals[nrow(cim_intervals),8] <- int[2,3] > summary(cim_thresholds_all[[10]])[1]

## The QTL for these traits were previously reported in Jarvis et al. 2021 GBE, and
## are not included for this analysis.

cim_intervals_significant <- cim_intervals %>%
  filter(Sig0.1==TRUE) %>%
  filter(!trait %in% c("Petal color","Anther color","Solitary flower number","Filament color")) %>%
  arrange(chromosome,peak)

## Sepal sinus length had too small of an interval to show up properly on the 
## graph, so here I manually make it a bit wider so it will show up clearly on 
## the graph

cim_intervals_significant[11,3] <- cim_intervals_significant[11,3] - 1.5
cim_intervals_significant[11,5] <- cim_intervals_significant[11,5] + 1.5

## Removing the Sterility trait

cim_intervals_significant <- cim_intervals_significant[-17,]

## Removing rownames (makes it easier to use R object directly in the interval plot)

rownames(cim_intervals_significant) <- NULL


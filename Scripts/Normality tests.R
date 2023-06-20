# Load data from QTL interval dataprep, Figure_hists

############# Shapiro test

shapiro_f2 <- data.frame(Trait = character(), pvalue = numeric())

for (i in 1:(length(f2)-1)){
  shapiro_f2[i,1] <- colnames(f2[i])
  shapiro_f2[i,2] <- shapiro.test(f2[[i]])[[2]]
}

# See how many traits have normal distributions
sum(shapiro_f2[1:22,2] >= 0.01)
shapiro_f2[which(shapiro_f2[1:22,2] >= 0.01),1]

############# Kolmogorov-Smirnov test

ks_f2 <- data.frame(Trait = character(), pvalue = numeric())

for (i in 1:(length(f2)-1)){
  ks_f2[i,1] <- colnames(f2[i])
  ks_f2[i,2] <- ks.test(f2[[i]], "pnorm")[[2]]
}

# See how many traits have normal distributions
sum(ks_f2[1:22,2] >= 0.05)
ks_f2[which(ks_f2[1:22,2] >= 0.05),1]

############# Anderson-Darling test

#ad_f2 <- data.frame(Trait = character(), pvalue = numeric())

#for (i in 1:(length(f2)-1)){
#  ad_f2[i,1] <- colnames(f2[i])
#  ad_f2[i,2] <- ad.test(f2[[i]])[[2]]
#}

# See how many traits have normal distributions
#sum(ad_f2[1:22,2] >= 0.01)
#ad_f2[which(ad_f2[1:22,2] >= 0.01),1]



############# Dip (bimodality) test

#dip_f2 <- data.frame(Trait = character(), pvalue = numeric())

#for (i in 1:(length(f2)-1)){
#  dip_f2[i,1] <- colnames(f2[i])
#  dip_f2[i,2] <- dip(na.omit(f2[[i]]))
#}

# See how many traits have normal distributions
#sum(dip_f2[1:22,2] >= 0.05)
#dip_f2[which(dip_f2[1:22,2] >= 0.05),1]










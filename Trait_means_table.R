## Load data from "Histogram and Correlation Figures.R" for parent, f1, and f2
## phenotype data

## Removing 'group' variable from datasets
f2 <- f2 %>%
  select(-group)
f1 <- f1 %>%
  select(-group)
yor <- yor %>%
  select(-group)
cap <- cap %>%
  select(-group)

########### Yorkii
yor_mean <- numeric()
yor_se <- numeric()
yor_sd <- numeric()

for (i in 1:21){
#  print(colnames(yor[i]))
#  print(sd(na.omit(yor[,i]))/length(na.omit(yor[,i])))
  yor_mean[i] <- mean(na.omit(yor[,i]))
  yor_se[i] <- sd(na.omit(yor[,i]))/sqrt(length(na.omit(yor[,i])))
  yor_sd[i] <- sd(na.omit(yor[,i]))
}


########### Capitata
cap_mean <- numeric()
cap_se <- numeric()
cap_sd <- numeric()

for (i in 1:21){
  #  print(colnames(cap[i]))
  #  print(sd(na.omit(cap[,i]))/length(na.omit(cap[,i])))
  cap_mean[i] <- mean(na.omit(cap[,i]))
  cap_se[i] <- sd(na.omit(cap[,i]))/sqrt(length(na.omit(cap[,i])))
  cap_sd[i] <- sd(na.omit(cap[,i]))
}


############ F1
f1_mean <- numeric()
f1_se <- numeric()
f1_sd <- numeric()

for (i in 1:21){
  #  print(colnames(f1[i]))
  #  print(sd(na.omit(f1[,i]))/length(na.omit(f1[,i])))
  f1_mean[i] <- mean(na.omit(f1[,i]))
  f1_se[i] <- sd(na.omit(f1[,i]))/sqrt(length(na.omit(f1[,i])))
  f1_sd[i] <- sd(na.omit(f1[,i]))
}

############ F2
f2_mean <- numeric()
f2_se <- numeric()
f2_sd <- numeric()

for (i in 1:21){
  #  print(colnames(f2[i]))
  #  print(sd(na.omit(f2[,i]))/length(na.omit(f2[,i])))
  f2_mean[i] <- mean(na.omit(f2[,i]))
  f2_se[i] <- sd(na.omit(f2[,i]))/sqrt(length(na.omit(f2[,i])))
  f2_sd[i] <- sd(na.omit(f2[,i]))
}


p_values <- numeric()

for (i in 1:21){
  #  print(colnames(f2[i]))
  #  print(sd(na.omit(f2[,i]))/length(na.omit(f2[,i])))
  p_values[i] <- t.test(na.omit(yor[,i]),na.omit(cap[,i]))[[3]]
}



#trait_table <- data.frame(trait = colnames(yor)[1:22],
#                          yor_mean,yor_se,
#                          cap_mean,cap_se,
#                          f1_mean,f1_se,
#                          f2_mean,f2_se)

trait_table <- data.frame(trait = colnames(yor)[1:21],
                          yor_mean,yor_sd,
                          cap_mean,cap_sd,
                          f1_mean,f1_sd,
                          f2_mean,f2_sd,
                          p_values)

#write.csv(trait_table, file = "trait_table_6-19-23.csv")



####################### growth chamber parent populations ######################

library(readxl)

yor <- read_xlsx("C:/Users/detemplj/OneDrive - Iowa State University/Documents/Research (Whipple Lab)/Gilia_Pictures/Gilia parent flowers/GiliaY.xlsx")
yor <- as.data.frame(yor)
yor <- yor[,2:19]
cap <- read_xlsx("C:/Users/detemplj/OneDrive - Iowa State University/Documents/Research (Whipple Lab)/Gilia_Pictures/Gilia parent flowers/GiliaC.xlsx")
cap <- as.data.frame(cap)
cap <- cap[,2:19]

## Getting 'Ovary shape' trait
yor$`Ovary shape` <- yor$`Ovary width`/yor$`Ovary length`
cap$`Ovary shape` <- cap$`Ovary width`/cap$`Ovary length`
yor <- subset(yor, select = -c(`Ovary width`, `Ovary length`))
cap <- subset(cap, select = -c(`Ovary width`, `Ovary length`))

## Updating trait names
colnames <- c("Petal length","Petal tube length","Petal tube width","Throat length","Petal lobe length",
              "Petal lobe width","Filament length","Free filament length","Anther length","Anther width",
              "Sepal length","Sepal sinus length","Sepal tooth length","Sepal midrib width","Style length",
              "Stigma length","Ovary shape")
colnames(yor) <- colnames
colnames(cap) <- colnames

## Updating order of trait names to match other figures
new_order <- c("Petal length","Petal lobe length","Petal lobe width",
               "Petal tube length","Petal tube width","Throat length",
               "Filament length","Free filament length","Anther length",
               "Anther width","Style length","Stigma length","Ovary shape",
               "Sepal length","Sepal sinus length","Sepal tooth length",
               "Sepal midrib width")
yor <- yor %>%
  select(new_order)
cap <- cap %>%
  select(new_order)

########### Yorkii
yor_mean <- numeric()
yor_se <- numeric()
yor_sd <- numeric()

for (i in 1:17){
  yor_mean[i] <- mean(na.omit(yor[,i]))
  yor_se[i] <- sd(na.omit(yor[,i]))/sqrt(length(na.omit(yor[,i])))
  yor_sd[i] <- sd(na.omit(yor[,i]))
}


########### Capitata
cap_mean <- numeric()
cap_se <- numeric()
cap_sd <- numeric()

for (i in 1:17){
  cap_mean[i] <- mean(na.omit(cap[,i]))
  cap_se[i] <- sd(na.omit(cap[,i]))/sqrt(length(na.omit(cap[,i])))
  cap_sd[i] <- sd(na.omit(cap[,i]))
}


p_values <- numeric()

for (i in 1:17){
  p_values[i] <- t.test(na.omit(yor[,i]),na.omit(cap[,i]))[[3]]
}



#trait_table_gr <- data.frame(trait = colnames(yor)[1:17],
#                          yor_mean,yor_se,
#                          cap_mean,cap_se,
#                          p_values)

trait_table_gr <- data.frame(trait = colnames(yor)[1:17],
                          yor_mean,yor_sd,
                          cap_mean,cap_sd,
                          p_values)
write.csv(trait_table_gr, file = "trait_table_6-19-23_growth_room_parents.csv")


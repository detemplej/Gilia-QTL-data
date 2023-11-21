# Libraries:
library(qtl)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(reshape2)

#### Dataprep
## Reading in data
#    yor = read.csv("D:/backups/Pictures/Research (Whipple Lab)/3_23 histogram data/yor_clean.csv", header = TRUE)
#    cap = read.csv("D:/backups/Pictures/Research (Whipple Lab)/3_23 histogram data/cap_clean.csv", header = TRUE)
#    f1 = read.csv("D:/backups/Pictures/Research (Whipple Lab)/3_23 histogram data/f1_clean.csv", header = TRUE)
#    gmap <- read.cross(format = "csv", file = "D:/backups/Pictures/Research (Whipple Lab)/My_paper/gmap5335_final.csv", map.function = "c-f", genotypes = c("AA", "AB", "BB"), alleles = c("A", "B"))
#    f2 <- gmap$pheno
    
    yor = read.csv("yor_clean.csv", header = TRUE)
    cap = read.csv("cap_clean.csv", header = TRUE)
    f1 = read.csv("f1_clean.csv", header = TRUE)
    gmap <- read.cross(format = "csv", file = "gmap5335_final.csv", map.function = "c-f", genotypes = c("AA", "AB", "BB"), alleles = c("A", "B"))
    f2 <- gmap$pheno
    
## Matching columns between datasets    
yor <- yor %>%
  select(colnames(yor)[colnames(yor) %in% colnames(gmap$pheno)])
cap <- cap %>%
  select(colnames(cap)[colnames(cap) %in% colnames(gmap$pheno)])
f1 <- f1 %>%
  select(colnames(f1)[colnames(f1) %in% colnames(gmap$pheno)])
## All sets should now have 22 equivalent traits

## f2 has 24 (Sterility and Death.early only measured here)
f2_other <- data.frame(f2$Death.early,f2$Sterility)
f2 <- f2 %>%
  select(colnames(f2)[colnames(f2) %in% colnames(yor)])

## Adding 'group' column in preparation for combining into one data.frame
f1 <- f1 %>%
  mutate(group = rep("f1", length(f1[,1])))
yor <- yor %>%
  mutate(group = rep("yor", length(yor[,1])))
cap <- cap %>%
  mutate(group = rep("cap", length(cap[,1])))
f2 <- f2 %>%
  mutate(group = rep("f2", length(f2[,1])))

# Changing trait names
# Updated 6-17-23 from tables from paper draft

colnames <- c("Anther length","Anther width","Sepal length","Sepal sinus length","Sepal tooth length",
              "Petal length","Days to flower","Free filament length","Internode length","Petal lobe length",
              "Petal lobe width","Main.flower","Vegetative rosette diameter","Sepal midrib width","Ovary shape",
              "Pedicel length","Filament length","Stigma length","Style length","Throat length","Petal tube length",
              "Petal tube width","group")
colnames(f2) <- colnames
colnames(f1) <- colnames
colnames(yor) <- colnames
colnames(cap) <- colnames

# Removing the "Main.flower" column--Floral buds will be kept and named as "Days to Flower" since these
# are highly correlated measurements of flowering time
f2 <- f2 %>%
  select(-Main.flower)
f1 <- f1 %>%
  select(-Main.flower)
yor <- yor %>%
  select(-Main.flower)
cap <- cap %>%
  select(-Main.flower)

new_order <- c("Petal length","Petal lobe length","Petal lobe width",
  "Petal tube length","Petal tube width","Throat length",
  "Filament length","Free filament length","Anther length",
  "Anther width","Style length","Stigma length","Ovary shape",
  "Sepal length","Sepal sinus length","Sepal tooth length",
  "Sepal midrib width","Pedicel length","Internode length",
  "Days to flower","Vegetative rosette diameter","group")

f2 <- f2 %>%
  select(new_order)
f1 <- f1 %>%
  select(new_order)
yor <- yor %>%
  select(new_order)
cap <- cap %>%
  select(new_order)

## Combining data, and changing factor levels
all_data <- rbind(yor,cap,f1,f2)
blank <- rbind(yor,cap,f1)
blank$group <- factor(blank$group, levels = c("cap","yor","f1"))    

## Getting summarized columns for parents and f1
yor_sum <- yor %>%
  summarize(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)))
yor_sum <- t(yor_sum)
colnames(yor_sum) <- "value"
yor_sum <- data.frame(variable = as.factor(rownames(yor_sum)),  # Create data for lines
                      vline = yor_sum[,1])
rownames(yor_sum) <- NULL

cap_sum <- cap %>%
  summarize(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)))
cap_sum <- t(cap_sum)
cap_sum <- data.frame(variable = as.factor(rownames(cap_sum)),  # Create data for lines
                      vline = cap_sum[,1])
rownames(cap_sum) <- NULL

f1_sum <- f1 %>%
  summarize(across(where(is.numeric), ~ mean(.x, na.rm = TRUE)))
f1_sum <- t(f1_sum)
f1_sum <- data.frame(variable = as.factor(rownames(f1_sum)),  # Create data for lines
                      vline = f1_sum[,1])
rownames(f1_sum) <- NULL

## Getting long format for facet_wrap
f2_long <- melt(f2, id.vars = "group")

## All trait histograms together

plot <- ggplot(f2_long, aes(x = value))+
  geom_histogram(position = "identity", alpha = 0.7, bins = 15)+
  facet_wrap(variable ~ .,scales = "free") 

plot +                                                 # Add different lines to facet plot
  geom_vline(data = yor_sum,
             aes(xintercept = vline, color = "blue"), lwd = 2,
             alpha = 0.5, show.legend = F) +
  geom_vline(data = f1_sum,
             aes(xintercept = vline, color = "red"), lwd = 2,
             alpha = 0.5, show.legend = F) +
  geom_vline(data = cap_sum,
             aes(xintercept = vline, color = "yellow"), lwd = 2,
             alpha = 0.5, key_glyph = draw_key_rect, show.legend = T) +
  scale_color_manual(values = c("red","yellow","blue"),
                     labels = c("G. yorkii", "F1", "G. capitata"),
                     name = "Group")+
  guides(color = guide_legend(override.aes = list(alpha = 0.5)))+
  labs(x = "Measurement",
       y = "Count") +
  theme(plot.title = element_text(hjust = 0.5),
        text = element_text(size = 17))

## 6/19/2023 saved as 1500w x 1200h figure (All Traits)



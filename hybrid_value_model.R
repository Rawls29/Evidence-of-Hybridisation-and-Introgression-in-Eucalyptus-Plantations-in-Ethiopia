####Housekeeping####
rm(list=ls())
setwd('C:/Users/samra/Documents/My Documents/Uni/Imperial/Summer Project/R Scripts') #Usual setting of WD

####Libraries####
library(tidyverse)
library(pwr)
library(ggplot2)

####Loading in Necessary Data####
STR_out <- read.table("combined_STR_admix.txt")
metadata <- read.csv("eucalyptus_sample_data_2022.csv")
load("esth_results_b1000_n3000_0.1.R")
esth_results <- esth_results$hi

####Formatting Dataset####
for(n in 1:nrow(metadata)){
  if(!is.na(metadata$Sample[n])){metadata$Sample[n] <- paste0("Euc_", metadata$Sample[n])}
  else{metadata$Sample[n] <- metadata$Sample[n]}
} #Creating a column in the metadata that can be used to merge with the STRUCTURE output

dataset <- merge(STR_out, metadata, by='Sample') #Merging the metadata and structure output

colnames(dataset)[1] <- "ID" #Creating a column in the merged dataset that can be used to merge with the gghybird hybrid index results

dataset <- merge(dataset, esth_results, by="ID") #Merging the datasets

dataset$hybrid_value2 <- NA
for(i in 1:nrow(dataset)){
  dataset$hybrid_value2[i] <- 0.5-(sqrt((0.5-dataset$h_posterior_mode[i])^2))
} #Adding a column of hybird values (bound between 0 an 0.5, as opposed to the hybrid index which is bound between 0 and 1)

####Circumference Model####
dataset$circ_per_year <- dataset$CBH_cm/dataset$age_y
circ_model <- lm(log(circ_per_year)~hybrid_value2, 
             data=dataset)
plot(circ_model)
summary(circ_model)

####Height Model####
dataset$height_per_year <- dataset$height_cm/dataset$age_y
height_model <- lm(log(height_per_year)~hybrid_value2, data=dataset)
plot(height_model)
summary(height_model)

###Growth per Year Model
dataset$growth_per_year <- (dataset$height_cm*dataset$CBH_cm)/dataset$age_y
annual_growth_model <- lm(log(growth_per_year)~hybrid_value2, data=dataset)
plot(annual_growth_model)
summary(annual_growth_model)

#Volume per Year Model
dataset$volume <- pi*dataset$height_cm*(dataset$CBH_cm/2*pi)^2
dataset$vol_per_yr <- dataset$volume/dataset$age_y
annual_vol_model <- lm(log(vol_per_yr)~hybrid_value2, data=dataset)
plot(annual_vol_model)
summary(annual_vol_model)

####Power Analysis####
pwr.f2.test(u=1, v=626, f2=0.0001484, sig.level = 0.05)
pwr.f2.test(u=1, power=0.8, f2=0.0001484, sig.level = 0.05)

####Plotting Hybrid Value against Altitude###
for(i in 1:nrow(dataset)){
  dataset$Altitude[i] <- 1180 + (dataset$Alt._range[i]*200)
}

means <- data.frame()
for(i in unique(dataset$Altitude)){
  altitude <- i
  mean_hybrid_value <- mean(dataset$hybrid_value2[which(dataset$Altitude==i)])
  temp <- data.frame(altitude, mean_hybrid_value)
  means <- rbind(means, temp)
}

ggplot(data=means)+
  geom_col(aes(x=altitude, y=mean_hybrid_value))+
  labs(x="Altitude (m)", y="Mean Admixture Value")+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background= element_blank(),
        axis.line= element_line(colour="black"))

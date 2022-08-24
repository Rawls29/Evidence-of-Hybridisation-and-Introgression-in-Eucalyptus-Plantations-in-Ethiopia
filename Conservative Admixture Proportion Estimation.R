####Housekeeping####
rm(list=ls())
setwd('C:/Users/samra/Documents/My Documents/Uni/Imperial/Summer Project/R Scripts')

####Libraries####
library(dplyr)
library(ggplot2)

####Loading in requied documents####
k_2_1 <- read.table("K_2_Proportions/K2_Q1.txt")
k_2_2 <- read.table("K_2_Proportions/K2_Q2.txt")
k_2_3 <- read.table("K_2_Proportions/K2_Q3.txt")
k_2_4 <- read.table("K_2_Proportions/K2_Q4.txt")
k_2_5 <- read.table("K_2_Proportions/K2_Q5.txt")
k_2_6 <- read.table("K_2_Proportions/K2_Q6.txt")
k_2_7 <- read.table("K_2_Proportions/K2_Q7.txt")
k_2_8 <- read.table("K_2_Proportions/K2_Q8.txt")
k_2_9 <- read.table("K_2_Proportions/K2_Q9.txt")
k_2_10 <- read.table("K_2_Proportions/K2_Q10.txt")

####Data wrangling####
k_2_1 <- subset(k_2_1, select = -c(V1, V5))
k_2_2 <- subset(k_2_2, select = -c(V1, V5))
k_2_3 <- subset(k_2_3, select = -c(V1, V5))
k_2_4 <- subset(k_2_4, select = -c(V1, V5))
k_2_5 <- subset(k_2_5, select = -c(V1, V5))
k_2_6 <- subset(k_2_6, select = -c(V1, V5))
k_2_7 <- subset(k_2_7, select = -c(V1, V5))
k_2_8 <- subset(k_2_8, select = -c(V1, V5))
k_2_9 <- subset(k_2_9, select = -c(V1, V5))
k_2_10 <- subset(k_2_10, select = -c(V1, V5))

colnames(k_2_1) <- c("Sample", "%_Missing", "Transect", "Pop_2_1", "Pop_1_1")
colnames(k_2_2) <- c("Sample", "%_Missing", "Transect", "Pop_2_2", "Pop_1_2")
colnames(k_2_3) <- c("Sample", "%_Missing", "Transect", "Pop_1_3", "Pop_2_3")
colnames(k_2_4) <- c("Sample", "%_Missing", "Transect", "Pop_1_4", "Pop_2_4")
colnames(k_2_5) <- c("Sample", "%_Missing", "Transect", "Pop_2_5", "Pop_1_5")
colnames(k_2_6) <- c("Sample", "%_Missing", "Transect", "Pop_2_6", "Pop_1_6")
colnames(k_2_7) <- c("Sample", "%_Missing", "Transect", "Pop_1_7", "Pop_2_7")
colnames(k_2_8) <- c("Sample", "%_Missing", "Transect", "Pop_1_8", "Pop_2_8")
colnames(k_2_9) <- c("Sample", "%_Missing", "Transect", "Pop_1_9", "Pop_2_9")
colnames(k_2_10) <- c("Sample", "%_Missing", "Transect", "Pop_1_10", "Pop_2_10")

Admix_est <- merge(k_2_1, k_2_2, by=c("Sample", "%_Missing", "Transect"))
Admix_est <- merge(Admix_est, k_2_3, by=c("Sample", "%_Missing", "Transect"))
Admix_est <- merge(Admix_est, k_2_4, by=c("Sample", "%_Missing", "Transect"))
Admix_est <- merge(Admix_est, k_2_5, by=c("Sample", "%_Missing", "Transect"))
Admix_est <- merge(Admix_est, k_2_6, by=c("Sample", "%_Missing", "Transect"))
Admix_est <- merge(Admix_est, k_2_7, by=c("Sample", "%_Missing", "Transect"))
Admix_est <- merge(Admix_est, k_2_8, by=c("Sample", "%_Missing", "Transect"))
Admix_est <- merge(Admix_est, k_2_9, by=c("Sample", "%_Missing", "Transect"))
Admix_est <- merge(Admix_est, k_2_10, by=c("Sample", "%_Missing", "Transect"))

combined <- data.frame(matrix(NA,
                              nrow=nrow(Admix_est),
                              ncol=4))
colnames(combined) <- c("Sample", "Transect", "Pop_1", "Pop_2")
for(n in 1:nrow(Admix_est)){
  combined$Sample[n] <- Admix_est$Sample[n]
  combined$Transect[n] <- Admix_est$Transect[n]
  Admix_estimates_Pop_1 <- list(Admix_est$Pop_1_1[n],
                          Admix_est$Pop_1_2[n],
                          Admix_est$Pop_1_3[n],
                          Admix_est$Pop_1_4[n],
                          Admix_est$Pop_1_5[n],
                          Admix_est$Pop_1_6[n],
                          Admix_est$Pop_1_7[n],
                          Admix_est$Pop_1_8[n],
                          Admix_est$Pop_1_9[n],
                          Admix_est$Pop_1_10[n])
  Admix_estimates_Pop_2 <- list(Admix_est$Pop_2_1[n],
                                Admix_est$Pop_2_2[n],
                                Admix_est$Pop_2_3[n],
                                Admix_est$Pop_2_4[n],
                                Admix_est$Pop_2_5[n],
                                Admix_est$Pop_2_6[n],
                                Admix_est$Pop_2_7[n],
                                Admix_est$Pop_2_8[n],
                                Admix_est$Pop_2_9[n],
                                Admix_est$Pop_2_10[n])
  if(Admix_est$Pop_1_1[n]>Admix_est$Pop_2_1[n]){
    combined$Pop_1[n] <- max(unlist(Admix_estimates_Pop_1))
    combined$Pop_2[n] <- min(unlist(Admix_estimates_Pop_2))
  }
  if(Admix_est$Pop_1_1[n]<Admix_est$Pop_2_1[n]){
    combined$Pop_2[n] <- max(unlist(Admix_estimates_Pop_2))
    combined$Pop_1[n] <- min(unlist(Admix_estimates_Pop_1))
  }
}

for(n in 1:nrow(combined)){
  if(combined$Pop_1[n]+combined$Pop_2[n]!=1){
    print(paste(combined$Sample, "Doesn't equal 1"))}
}

write.table(combined, "combined_STR_admix.txt", sep="\t", quote=FALSE)

####Q-Values Plot####
combined_ordered <- arrange(combined, Pop_2)
combined_ordered$n <- NA
for(i in 1:nrow(combined_ordered)){
  combined_ordered$n[i] <- i
}
STRUCTURE_Plot <- ggplot(data = combined_ordered)+
  geom_col(aes(x=n, y=1), fill="#56B4E9")+
  geom_col(aes(x=n, y=Pop_1), fill= "#FF3366")+
  ylab("Ancestry Proportions")+
  xlab("")+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background= element_blank(),
        axis.line= element_line(colour="black"))
print(STRUCTURE_Plot)

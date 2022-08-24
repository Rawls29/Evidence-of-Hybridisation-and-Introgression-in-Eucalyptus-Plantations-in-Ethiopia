####Housekeeping####
rm(list=ls())
setwd('C:/Users/samra/Documents/My Documents/Uni/Imperial/Summer Project/R Scripts')

####Libraries####
library(tidyverse)
library(lme4)
library(lmerTest)
library(car)

####Loading Required Files####
a <- read.table("combined_STR_admix.txt")
metadata <- read.csv("eucalyptus_sample_data_2022.csv")
load("esth_results_b1000_n3000_0.1.R")

####Combining Files into One Dataframe####
for(n in 1:nrow(metadata)){
  if(!is.na(metadata$Sample[n])){metadata$Sample[n] <- paste0("Euc_", metadata$Sample[n])}
  else{metadata$Sample[n] <- metadata$Sample[n]}
}

b <- merge(a, metadata, by='Sample')

colnames(b)[1] <- "ID"
b <- merge(b, esth_results$hi, by="ID")

#Assigning population based on STRUCTURE output values
b$Population <- NA
for(i in 1:nrow(b)){
  if(!is.na(b$Pop_1[i]))
  {if(b$Pop_1[i]>0.8){b$Population[i]<-1}
    else{
      if(b$Pop_2[i]>0.8){b$Population[i]<-2}
      else{b$Population[i]<-"Admixed"}
    }}
}

####Subsetting by Transect and Site####
T1 <- subset(b, b$Transect.x==1)
T1_1 <- subset(T1, T1$Alt._range==1)
T1_2 <- subset(T1, T1$Alt._range==2)
T1_3 <- subset(T1, T1$Alt._range==3)
T1_4 <- subset(T1, T1$Alt._range==4)
T1_5 <- subset(T1, T1$Alt._range==5)
T1_6 <- subset(T1, T1$Alt._range==6)
T1_7 <- subset(T1, T1$Alt._range==7)
T1_8 <- subset(T1, T1$Alt._range==8)
T1_9 <- subset(T1, T1$Alt._range==9)
T1_10 <- subset(T1, T1$Alt._range==10)

T2 <- subset(b, b$Transect.x==2)
T2_1 <- subset(T2, T2$Alt._range==1)
T2_2 <- subset(T2, T2$Alt._range==2)
T2_3 <- subset(T2, T2$Alt._range==3)
T2_4 <- subset(T2, T2$Alt._range==4)
T2_5 <- subset(T2, T2$Alt._range==5)
T2_6 <- subset(T2, T2$Alt._range==6)
T2_7 <- subset(T2, T2$Alt._range==7)
T2_8 <- subset(T2, T2$Alt._range==8)
T2_9 <- subset(T2, T2$Alt._range==9)
T2_10 <- subset(T2, T2$Alt._range==10)

T3 <- subset(b, b$Transect.x==3)
T3_1 <- subset(T3, T3$Alt._range==1)
T3_2 <- subset(T3, T3$Alt._range==2)
T3_3 <- subset(T3, T3$Alt._range==3)
T3_4 <- subset(T3, T3$Alt._range==4)
T3_5 <- subset(T3, T3$Alt._range==5)
T3_6 <- subset(T3, T3$Alt._range==6)
T3_7 <- subset(T3, T3$Alt._range==7)
T3_8 <- subset(T3, T3$Alt._range==8)
T3_9 <- subset(T3, T3$Alt._range==9)
T3_10 <- subset(T3, T3$Alt._range==10)

T4 <- subset(b, b$Transect.x==4)
T4_1 <- subset(T4, T4$Alt._range==1)
T4_2 <- subset(T4, T4$Alt._range==2)
T4_3 <- subset(T4, T4$Alt._range==3)
T4_4 <- subset(T4, T4$Alt._range==4)
T4_5 <- subset(T4, T4$Alt._range==5)
T4_6 <- subset(T4, T4$Alt._range==6)
T4_7 <- subset(T4, T4$Alt._range==7)
T4_8 <- subset(T4, T4$Alt._range==8)
T4_9 <- subset(T4, T4$Alt._range==9)
T4_10 <- subset(T4, T4$Alt._range==10)

T5 <- subset(b, b$Transect.x==5)
T5_1 <- subset(T5, T5$Alt._range==1)
T5_2 <- subset(T5, T5$Alt._range==2)
T5_3 <- subset(T5, T5$Alt._range==3)
T5_4 <- subset(T5, T5$Alt._range==4)
T5_5 <- subset(T5, T5$Alt._range==5)
T5_6 <- subset(T5, T5$Alt._range==6)
T5_7 <- subset(T5, T5$Alt._range==7)
T5_8 <- subset(T5, T5$Alt._range==8)
T5_9 <- subset(T5, T5$Alt._range==9)
T5_10 <- subset(T5, T5$Alt._range==10)

T6 <- subset(b, b$Transect.x==6)
T6_1 <- subset(T6, T6$Alt._range==1)
T6_2 <- subset(T6, T6$Alt._range==2)
T6_3 <- subset(T6, T6$Alt._range==3)
T6_4 <- subset(T6, T6$Alt._range==4)
T6_5 <- subset(T6, T6$Alt._range==5)
T6_6 <- subset(T6, T6$Alt._range==6)
T6_7 <- subset(T6, T6$Alt._range==7)
T6_8 <- subset(T6, T6$Alt._range==8)
T6_9 <- subset(T6, T6$Alt._range==9)
T6_10 <- subset(T6, T6$Alt._range==10)

T7 <- subset(b, b$Transect.x==7)
T7_1 <- subset(T7, T7$Alt._range==1)
T7_2 <- subset(T7, T7$Alt._range==2)
T7_3 <- subset(T7, T7$Alt._range==3)
T7_4 <- subset(T7, T7$Alt._range==4)
T7_5 <- subset(T7, T7$Alt._range==5)
T7_6 <- subset(T7, T7$Alt._range==6)
T7_7 <- subset(T7, T7$Alt._range==7)
T7_8 <- subset(T7, T7$Alt._range==8)
T7_9 <- subset(T7, T7$Alt._range==9)
T7_10 <- subset(T7, T7$Alt._range==10)

T8 <- subset(b, b$Transect.x==8)
T8_1 <- subset(T8, T8$Alt._range==1)
T8_2 <- subset(T8, T8$Alt._range==2)
T8_3 <- subset(T8, T8$Alt._range==3)
T8_4 <- subset(T8, T8$Alt._range==4)
T8_5 <- subset(T8, T8$Alt._range==5)
T8_6 <- subset(T8, T8$Alt._range==6)
T8_7 <- subset(T8, T8$Alt._range==7)
T8_8 <- subset(T8, T8$Alt._range==8)
T8_9 <- subset(T8, T8$Alt._range==9)
T8_10 <- subset(T8, T8$Alt._range==10)

T9 <- subset(b, b$Transect.x==8)
T9_1 <- subset(T9, T9$Alt._range==1)
T9_2 <- subset(T9, T9$Alt._range==2)
T9_3 <- subset(T9, T9$Alt._range==3)
T9_4 <- subset(T9, T9$Alt._range==4)
T9_5 <- subset(T9, T9$Alt._range==5)
T9_6 <- subset(T9, T9$Alt._range==6)
T9_7 <- subset(T9, T9$Alt._range==7)
T9_8 <- subset(T9, T9$Alt._range==8)
T9_9 <- subset(T9, T9$Alt._range==9)
T9_10 <- subset(T9, T9$Alt._range==10)

T10 <- subset(b, b$Transect.x==1)
T10_1 <- subset(T10, T10$Alt._range==1)
T10_2 <- subset(T10, T10$Alt._range==2)
T10_3 <- subset(T10, T10$Alt._range==3)
T10_4 <- subset(T10, T10$Alt._range==4)
T10_5 <- subset(T10, T10$Alt._range==5)
T10_6 <- subset(T10, T10$Alt._range==6)
T10_7 <- subset(T10, T10$Alt._range==7)
T10_8 <- subset(T10, T10$Alt._range==8)
T10_9 <- subset(T10, T10$Alt._range==9)
T10_10 <- subset(T10, T10$Alt._range==10)

####Getting average hybird indices (from gghybrid) per site####
T1_1_av <- mean(T1_1$h_posterior_mode)
T1_2_av <- mean(T1_2$h_posterior_mode)
T1_3_av <- mean(T1_3$h_posterior_mode)
T1_4_av <- mean(T1_4$h_posterior_mode)
T1_5_av <- mean(T1_5$h_posterior_mode)
T1_6_av <- mean(T1_6$h_posterior_mode)
T1_7_av <- mean(T1_7$h_posterior_mode)
T1_8_av <- mean(T1_8$h_posterior_mode)
T1_9_av <- mean(T1_9$h_posterior_mode)
T1_10_av <- mean(T1_10$h_posterior_mode)
Transect_1 <- c(T1_1_av, T1_2_av, T1_3_av, T1_4_av, T1_5_av, T1_6_av, T1_7_av,
                T1_8_av, T1_9_av, T1_10_av)

T2_1_av <- mean(T2_1$h_posterior_mode)
T2_2_av <- mean(T2_2$h_posterior_mode)
T2_3_av <- mean(T2_3$h_posterior_mode)
T2_4_av <- mean(T2_4$h_posterior_mode)
T2_5_av <- mean(T2_5$h_posterior_mode)
T2_6_av <- mean(T2_6$h_posterior_mode)
T2_7_av <- mean(T2_7$h_posterior_mode)
T2_8_av <- mean(T2_8$h_posterior_mode)
T2_9_av <- mean(T2_9$h_posterior_mode)
T2_10_av <- mean(T2_10$h_posterior_mode)
Transect_2 <- c(T2_1_av, T2_2_av, T2_3_av, T2_4_av, T2_5_av, T2_6_av, T2_7_av,
                T2_8_av, T2_9_av, T2_10_av)

T3_1_av <- mean(T3_1$h_posterior_mode)
T3_2_av <- mean(T3_2$h_posterior_mode)
T3_3_av <- mean(T3_3$h_posterior_mode)
T3_4_av <- mean(T3_4$h_posterior_mode)
T3_5_av <- mean(T3_5$h_posterior_mode)
T3_6_av <- mean(T3_6$h_posterior_mode)
T3_7_av <- mean(T3_7$h_posterior_mode)
T3_8_av <- mean(T3_8$h_posterior_mode)
T3_9_av <- mean(T3_9$h_posterior_mode)
T3_10_av <- mean(T3_10$h_posterior_mode)
Transect_3 <- c(T3_1_av, T3_2_av, T3_3_av, T3_4_av, T3_5_av, T3_6_av, T3_7_av,
                T3_8_av, T3_9_av, T3_10_av)

T4_1_av <- mean(T4_1$h_posterior_mode)
T4_2_av <- mean(T4_2$h_posterior_mode)
T4_3_av <- mean(T4_3$h_posterior_mode)
T4_4_av <- mean(T4_4$h_posterior_mode)
T4_5_av <- mean(T4_5$h_posterior_mode)
T4_6_av <- mean(T4_6$h_posterior_mode)
T4_7_av <- mean(T4_7$h_posterior_mode)
T4_8_av <- mean(T4_8$h_posterior_mode)
T4_9_av <- mean(T4_9$h_posterior_mode)
T4_10_av <- mean(T4_10$h_posterior_mode)
Transect_4 <- c(T4_1_av, T4_2_av, T4_3_av, T4_4_av, T4_5_av, T4_6_av, T4_7_av,
                T4_8_av, T4_9_av, T4_10_av)

T5_1_av <- mean(T5_1$h_posterior_mode)
T5_2_av <- mean(T5_2$h_posterior_mode)
T5_3_av <- mean(T5_3$h_posterior_mode)
T5_4_av <- mean(T5_4$h_posterior_mode)
T5_5_av <- mean(T5_5$h_posterior_mode)
T5_6_av <- mean(T5_6$h_posterior_mode)
T5_7_av <- mean(T5_7$h_posterior_mode)
T5_8_av <- mean(T5_8$h_posterior_mode)
T5_9_av <- mean(T5_9$h_posterior_mode)
T5_10_av <- mean(T5_10$h_posterior_mode)
Transect_5 <- c(T5_1_av, T5_2_av, T5_3_av, T5_4_av, T5_5_av, T5_6_av, T5_7_av,
                T5_8_av, T5_9_av, T5_10_av)

T6_1_av <- mean(T6_1$h_posterior_mode)
T6_2_av <- mean(T6_2$h_posterior_mode)
T6_3_av <- mean(T6_3$h_posterior_mode)
T6_4_av <- mean(T6_4$h_posterior_mode)
T6_5_av <- mean(T6_5$h_posterior_mode)
T6_6_av <- mean(T6_6$h_posterior_mode)
T6_7_av <- mean(T6_7$h_posterior_mode)
T6_8_av <- mean(T6_8$h_posterior_mode)
T6_9_av <- mean(T6_9$h_posterior_mode)
T6_10_av <- mean(T6_10$h_posterior_mode)
Transect_6 <- c(T6_1_av, T6_2_av, T6_3_av, T6_4_av, T6_5_av, T6_6_av, T6_7_av,
                T6_8_av, T6_9_av, T6_10_av)

T7_1_av <- mean(T7_1$h_posterior_mode)
T7_2_av <- mean(T7_2$h_posterior_mode)
T7_3_av <- mean(T7_3$h_posterior_mode)
T7_4_av <- mean(T7_4$h_posterior_mode)
T7_5_av <- mean(T7_5$h_posterior_mode)
T7_6_av <- mean(T7_6$h_posterior_mode)
T7_7_av <- mean(T7_7$h_posterior_mode)
T7_8_av <- mean(T7_8$h_posterior_mode)
T7_9_av <- mean(T7_9$h_posterior_mode)
T7_10_av <- mean(T7_10$h_posterior_mode)
Transect_7 <- c(T7_1_av, T7_2_av, T7_3_av, T7_4_av, T7_5_av, T7_6_av, T7_7_av,
                T7_8_av, T7_9_av, T7_10_av)

T8_1_av <- mean(T8_1$h_posterior_mode)
T8_2_av <- mean(T8_2$h_posterior_mode)
T8_3_av <- mean(T8_3$h_posterior_mode)
T8_4_av <- mean(T8_4$h_posterior_mode)
T8_5_av <- mean(T8_5$h_posterior_mode)
T8_6_av <- mean(T8_6$h_posterior_mode)
T8_7_av <- mean(T8_7$h_posterior_mode)
T8_8_av <- mean(T8_8$h_posterior_mode)
T8_9_av <- mean(T8_9$h_posterior_mode)
T8_10_av <- mean(T8_10$h_posterior_mode)
Transect_8 <- c(T8_1_av, T8_2_av, T8_3_av, T8_4_av, T8_5_av, T8_6_av, T8_7_av,
                T8_8_av, T8_9_av, T8_10_av)

T9_1_av <- mean(T9_1$h_posterior_mode)
T9_2_av <- mean(T9_2$h_posterior_mode)
T9_3_av <- mean(T9_3$h_posterior_mode)
T9_4_av <- mean(T9_4$h_posterior_mode)
T9_5_av <- mean(T9_5$h_posterior_mode)
T9_6_av <- mean(T9_6$h_posterior_mode)
T9_7_av <- mean(T9_7$h_posterior_mode)
T9_8_av <- mean(T9_8$h_posterior_mode)
T9_9_av <- mean(T9_9$h_posterior_mode)
T9_10_av <- mean(T9_10$h_posterior_mode)
Transect_9 <- c(T9_1_av, T9_2_av, T9_3_av, T9_4_av, T9_5_av, T9_6_av, T9_7_av,
                T9_8_av, T9_9_av, T9_10_av)

T10_1_av <- mean(T10_1$h_posterior_mode)
T10_2_av <- mean(T10_2$h_posterior_mode)
T10_3_av <- mean(T10_3$h_posterior_mode)
T10_4_av <- mean(T10_4$h_posterior_mode)
T10_5_av <- mean(T10_5$h_posterior_mode)
T10_6_av <- mean(T10_6$h_posterior_mode)
T10_7_av <- mean(T10_7$h_posterior_mode)
T10_8_av <- mean(T10_8$h_posterior_mode)
T10_9_av <- mean(T10_9$h_posterior_mode)
T10_10_av <- mean(T10_10$h_posterior_mode)
Transect_10 <- c(T10_1_av, T10_2_av, T10_3_av, T10_4_av, T10_5_av, T10_6_av, T10_7_av,
                T10_8_av, T10_9_av, T10_10_av)

####Creating a dataframe with average site hybrid indices####
df <- data.frame(Transect_1, Transect_2, Transect_3, Transect_4, Transect_5,
                 Transect_6, Transect_7, Transect_8, Transect_9, Transect_10)

Transect <- c(1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
              2, 2, 2, 2, 2, 2, 2, 2, 2, 2,
              3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
              4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
              5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
              6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
              7, 7, 7, 7, 7, 7, 7, 7, 7, 7,
              8, 8, 8, 8, 8, 8, 8, 8, 8, 8,
              9, 9, 9, 9, 9, 9, 9, 9, 9, 9,
              10, 10, 10, 10, 10, 10, 10, 10, 10, 10)
Altitude_Range <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10,
                    1, 2, 3, 4, 5, 6, 7, 8, 9, 10,
                    1, 2, 3, 4, 5, 6, 7, 8, 9, 10,
                    1, 2, 3, 4, 5, 6, 7, 8, 9, 10,
                    1, 2, 3, 4, 5, 6, 7, 8, 9, 10,
                    1, 2, 3, 4, 5, 6, 7, 8, 9, 10,
                    1, 2, 3, 4, 5, 6, 7, 8, 9, 10,
                    1, 2, 3, 4, 5, 6, 7, 8, 9, 10,
                    1, 2, 3, 4, 5, 6, 7, 8, 9, 10,
                    1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
Av_Hybrid_Index <- c(T1_1_av, T1_2_av, T1_3_av, T1_4_av, T1_5_av,
                    T1_6_av, T1_7_av, T1_8_av, T1_9_av, T1_10_av,
                    T2_1_av, T2_2_av, T2_3_av, T2_4_av, T2_5_av,
                    T2_6_av, T2_7_av, T2_8_av, T2_9_av, T2_10_av,
                    T3_1_av, T3_2_av, T3_3_av, T3_4_av, T3_5_av,
                    T3_6_av, T3_7_av, T3_8_av, T3_9_av, T3_10_av,
                    T4_1_av, T4_2_av, T4_3_av, T4_4_av, T4_5_av,
                    T4_6_av, T4_7_av, T4_8_av, T4_9_av, T4_10_av,
                    T5_1_av, T5_2_av, T5_3_av, T5_4_av, T5_5_av,
                    T5_6_av, T5_7_av, T5_8_av, T5_9_av, T5_10_av,
                    T6_1_av, T6_2_av, T6_3_av, T6_4_av, T6_5_av,
                    T6_6_av, T6_7_av, T6_8_av, T6_9_av, T6_10_av,
                    T7_1_av, T7_2_av, T7_3_av, T7_4_av, T7_5_av,
                    T7_6_av, T7_7_av, T7_8_av, T7_9_av, T7_10_av,
                    T8_1_av, T8_2_av, T8_3_av, T8_4_av, T8_5_av,
                    T8_6_av, T8_7_av, T8_8_av, T8_9_av, T8_10_av,
                    T9_1_av, T9_2_av, T9_3_av, T9_4_av, T9_5_av,
                    T9_6_av, T9_7_av, T9_8_av, T9_9_av, T9_10_av,
                    T10_1_av, T10_2_av, T10_3_av, T10_4_av, T10_5_av,
                    T10_6_av, T10_7_av, T10_8_av, T10_9_av, T10_10_av)
d2 <- data.frame(Transect, Altitude_Range, Av_Hybrid_Index)

####Model of Hybrid Index against Altitude####
b$Pop_1 <- as.numeric(b$Pop_1)
b$Pop_2 <- as.numeric(b$Pop_2)
b$Alt._range <- as.numeric(b$Alt._range)
for(i in 1:nrow(b)){
  b$Altitude[i] <- 1180 + (b$Alt._range[i]*200)
}
d2$Altitude <- NA
for(i in 1:nrow(d2)){
  d2$Altitude[i] <- 1180 + (d2$Altitude_Range[i]*200)
}


b$Transect.x <- as.character(b$Transect.x)
d2$Transect <- as.character(d2$Transect)
no_pure <- b[which(b$h_posterior_mode!=0&b$h_posterior_mode!=1),]

hist(b$Pop_1)
hist(logit(b$Pop_1))
hist(b$Alt._range)

model <- lmer(h_posterior_mode~Altitude+(1|Transect.x), data=no_pure)
qqnorm(resid(model))
qqline(resid(model))
plot(model)
summary(model)

fixed <- as.data.frame(coef(summary(model)))

ggplot()+
  geom_point(data=d2, aes(x=Altitude, y=Av_Hybrid_Index), colour="black", size=2)+
  geom_abline(aes(intercept=fixed[1,1] , slope=fixed[2,1]))+
  labs(x="Altitude (m)", y="Hybrid Index")+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background= element_blank(),
        axis.line= element_line(colour="black"))

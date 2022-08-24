####Housekeeping####
rm(list=ls())
setwd('C:/Users/samra/Documents/My Documents/Uni/Imperial/Summer Project') #Usual setting of WD

####Loading in Necessary Data####
STR_out <- read.table("combined_STR_admix.txt")
metadata <- read.csv("eucalyptus_sample_data_2022.csv")
load("esth_results_b1000_n3000_0.1.R")
esth_results <- esth_results$hi

####Formatting Dataset####
for(n in 1:nrow(metadata)){
  if(!is.na(metadata$Sample[n])){metadata$Sample[n] <- paste0("Euc_", metadata$Sample[n])}
  else{metadata$Sample[n] <- metadata$Sample[n]}
}

dataset <- merge(STR_out, metadata, by='Sample')

dataset$Population <- NA
for(i in 1:nrow(dataset)){
  if(!is.na(dataset$Pop_1[i]))
  {if(dataset$Pop_1[i]>0.8){dataset$Population[i]<-"Camaldulensis"}
    else{
      if(dataset$Pop_2[i]>0.8){dataset$Population[i]<-"Globulus"}
      else{dataset$Population[i]<-"Admixed"}
    }}
}

colnames(dataset)[1] <- "ID"

dataset <- merge(dataset, esth_results, by="ID")

dataset$hybrid_value2 <- NA
for(i in 1:nrow(dataset)){
  dataset$hybrid_value2[i] <- 0.5-(sqrt((0.5-dataset$h_posterior_mode[i])^2))
}

dataset$Population_Hybrid <- NA
for(i in 1:nrow(dataset)){
  if(dataset$h_posterior_mode[i]<0.2){dataset$Population_Hybrid[i] <- "Cam"}else
    {if(dataset$h_posterior_mode[i]>0.8){dataset$Population_Hybrid[i]<-"Glob"}else
      {dataset$Population_Hybrid[i] <- "Admix"}}
}

####This section is getting sub-datasets for each site, and making tables to see which ones have both admixed and parental individuals####
T1 <- dataset[which(dataset$Transect.x==1),]
T2 <- dataset[which(dataset$Transect.x==2),]
T3 <- dataset[which(dataset$Transect.x==3),]
T4 <- dataset[which(dataset$Transect.x==4),]
T5 <- dataset[which(dataset$Transect.x==5),]
T6 <- dataset[which(dataset$Transect.x==6),]
T7 <- dataset[which(dataset$Transect.x==7),]
T8 <- dataset[which(dataset$Transect.x==8),]

T1_1 <- T1[which(T1$Site==1),]
table(T1_1$Population_Hybrid)
T1_2 <- T1[which(T1$Site==2),]
table(T1_2$Population_Hybrid)
T1_3 <- T1[which(T1$Site==3),]
table(T1_3$Population_Hybrid)
T1_4 <- T1[which(T1$Site==4),]
table(T1_4$Population_Hybrid)
T1_5 <- T1[which(T1$Site==5),]
table(T1_5$Population_Hybrid)
T1_6 <- T1[which(T1$Site==6),]
table(T1_6$Population_Hybrid)
T1_7 <- T1[which(T1$Site==7),]
table(T1_7$Population_Hybrid)
T1_8 <- T1[which(T1$Site==8),]
table(T1_8$Population_Hybrid)
T1_9 <- T1[which(T1$Site==9),]
table(T1_9$Population_Hybrid)
T1_10 <- T1[which(T1$Site==10),]
table(T1_10$Population_Hybrid)

T4_1 <- T4[which(T4$Site==1),]
table(T4_1$Population_Hybrid)
T4_2 <- T4[which(T4$Site==2),]
table(T4_2$Population_Hybrid)
T4_3 <- T4[which(T4$Site==3),]
table(T4_3$Population_Hybrid)
T4_4 <- T4[which(T4$Site==4),]
table(T4_4$Population_Hybrid)
T4_5 <- T4[which(T4$Site==5),]
table(T4_5$Population_Hybrid)
T4_6 <- T4[which(T4$Site==6),]
table(T4_6$Population_Hybrid)
T4_7 <- T4[which(T4$Site==7),]
table(T4_7$Population_Hybrid)
T4_8 <- T4[which(T4$Site==8),]
table(T4_8$Population_Hybrid)
T4_9 <- T4[which(T4$Site==9),]
table(T4_9$Population_Hybrid)
T4_10 <- T4[which(T4$Site==10),]
table(T4_10$Population_Hybrid)

T6_1 <- T6[which(T6$Site==1),]
table(T6_1$Population_Hybrid)
T6_2 <- T6[which(T6$Site==2),]
table(T6_2$Population_Hybrid)
T6_3 <- T6[which(T6$Site==3),]
table(T6_3$Population_Hybrid)
T6_4 <- T6[which(T6$Site==4),]
table(T6_4$Population_Hybrid)
T6_5 <- T6[which(T6$Site==5),]
table(T6_5$Population_Hybrid)
T6_6 <- T6[which(T6$Site==6),]
table(T6_6$Population_Hybrid)
T6_7 <- T6[which(T6$Site==7),]
table(T6_7$Population_Hybrid)
T6_8 <- T6[which(T6$Site==8),]
table(T6_8$Population_Hybrid)
T6_9 <- T6[which(T6$Site==9),]
table(T6_9$Population_Hybrid)
T6_10 <- T6[which(T6$Site==10),]
table(T6_10$Population_Hybrid)

####This bit is me mostly manually pulling out the entries I want for the interesting sites to get pairs####
T1_1_cam <- T1_1[which(T1_1$h_posterior_mode==min(T1_1$h_posterior_mode)),]
T1_1_most_cam <- T1_1_cam[which(T1_1_cam$h_cred_int_upper==min(T1_1_cam$h_cred_int_upper)),]
T1_1_ad <- T1_1[which(T1_1$hybrid_value2==max(T1_1$hybrid_value2)),]
T1_1_ad_hybrid <- T1_1_ad$hybrid_value2
T1_1_cam_value <- (pi*T1_1_most_cam$height_cm*(T1_1_most_cam$CBH_cm/2*pi)^2)/T1_1_most_cam$age_y #New grwoth rate formula, looks more horrible but actually gives cm^3 per year, which is neat
T1_1_ad_value <- (pi*T1_1_ad$height_cm*(T1_1_ad$CBH_cm/2*pi)^2)/T1_1_ad$age_y

pairs <- data.frame(T1_1_cam_value, NA, T1_1_ad_value, T1_1_ad_hybrid)
colnames(pairs) <- c("Cam_Value", "Glob_Value", "Admix_Value", "Hybrid_Value")

T1_3_cam <- T1_3[which(T1_3$h_posterior_mode==min(T1_3$h_posterior_mode)),]
T1_3_most_cam <- T1_3_cam[which(T1_3_cam$h_cred_int_upper==min(T1_3_cam$h_cred_int_upper)),]
T1_3_glob <- T1_3[which(T1_3$h_posterior_mode==max(T1_3$h_posterior_mode)),]
T1_3_most_glob <- T1_3_glob[which(T1_3_glob$h_cred_int_lower==max(T1_3_glob$h_cred_int_lower)),]
T1_3_ad <- T1_3[which(T1_3$hybrid_value2==max(T1_3$hybrid_value2)),]
T1_3_ad_hybrid <- T1_3_ad$hybrid_value2
T1_3_cam_value <- (pi*T1_3_most_cam$height_cm*(T1_3_most_cam$CBH_cm/2*pi)^2)/T1_3_most_cam$age_y
T1_3_glob_value <- (pi*T1_3_most_glob$height_cm*(T1_3_most_glob$CBH_cm/2*pi)^2)/T1_3_most_glob$age_y
T1_3_ad_value <- (pi*T1_3_ad$height_cm*(T1_3_ad$CBH_cm/2*pi)^2)/T1_3_ad$age_y

pairs_temp <- data.frame(T1_3_cam_value, T1_3_glob_value, T1_3_ad_value, T1_3_ad_hybrid)
colnames(pairs_temp) <- c("Cam_Value", "Glob_Value", "Admix_Value", "Hybrid_Value")
pairs <- rbind(pairs, pairs_temp)

T1_9_glob <- T1_9[which(T1_9$h_posterior_mode==max(T1_9$h_posterior_mode)),]
T1_9_most_glob <- T1_9_glob[which(T1_9_glob$h_cred_int_lower==max(T1_9_glob$h_cred_int_lower)),]
T1_9_ad <- T1_9[which(T1_9$hybrid_value2==max(T1_9$hybrid_value2)),]
T1_9_ad_hybrid <- T1_9_ad$hybrid_value2
T1_9_glob_value <- (pi*T1_9_most_glob$height_cm*(T1_9_most_glob$CBH_cm/2*pi)^2)
T1_9_ad_value <- (pi*T1_9_ad$height_cm*(T1_9_ad$CBH_cm/2*pi)^2)

pairs_temp <- data.frame(NA, T1_9_glob_value, T1_9_ad_value, T1_9_ad_hybrid)
colnames(pairs_temp) <- c("Cam_Value", "Glob_Value", "Admix_Value", "Hybrid_Value")
pairs <- rbind(pairs, pairs_temp)

T1_10_glob <- T1_10[which(T1_10$h_posterior_mode==max(T1_10$h_posterior_mode)),]
T1_10_most_glob <- T1_10_glob[which(T1_10_glob$h_cred_int_lower==max(T1_10_glob$h_cred_int_lower)),]
T1_10_ad <- T1_10[which(T1_10$hybrid_value2==max(T1_10$hybrid_value2)),]
T1_10_ad_hybrid <- T1_10_ad$hybrid_value2
T1_10_glob_value <- (pi*T1_10_most_glob$height_cm*(T1_10_most_glob$CBH_cm/2*pi)^2)/T1_10_most_glob$age_y
T1_10_ad_value <- (pi*T1_10_ad$height_cm*(T1_10_ad$CBH_cm/2*pi)^2)/T1_10_ad$age_y

pairs_temp <- data.frame(NA, T1_10_glob_value, T1_10_ad_value, T1_10_ad_hybrid)
colnames(pairs_temp) <- c("Cam_Value", "Glob_Value", "Admix_Value", "Hybrid_Value")
pairs <- rbind(pairs, pairs_temp)

T4_2_cam <- T4_2[which(T4_2$h_posterior_mode==min(T4_2$h_posterior_mode)),]
T4_2_most_cam <- T4_2_cam[which(T4_2_cam$h_cred_int_upper==min(T4_2_cam$h_cred_int_upper)),]
T4_2_glob <- T4_2[which(T4_2$h_posterior_mode==max(T4_2$h_posterior_mode)),]
T4_2_most_glob <- T4_2_glob[which(T4_2_glob$h_cred_int_lower==max(T4_2_glob$h_cred_int_lower)),]
T4_2_ad <- T1_3[which(T4_2$hybrid_value2==max(T4_2$hybrid_value2)),]
T4_2_ad_hybrid <- T4_2_ad$hybrid_value2
T4_2_cam_value <- (pi*T4_2_most_cam$height_cm*(T4_2_most_cam$CBH_cm/2*pi)^2)/T4_2_most_cam$age_y
T4_2_glob_value <- (pi*T4_2_most_glob$height_cm*(T4_2_most_glob$CBH_cm/2*pi)^2)/T4_2_most_glob$age_y
T4_2_ad_value <- (pi*T4_2_ad$height_cm*(T4_2_ad$CBH_cm/2*pi)^2)/T4_2_ad$age_y

pairs_temp <- data.frame(T4_2_cam_value, T4_2_glob_value, T4_2_ad_value, T4_2_ad_hybrid)
colnames(pairs_temp) <- c("Cam_Value", "Glob_Value", "Admix_Value", "Hybrid_Value")
pairs <- rbind(pairs, pairs_temp)

T4_3_cam <- T4_3[which(T4_3$h_posterior_mode==min(T4_3$h_posterior_mode)),]
T4_3_most_cam <- T4_3_cam[which(T4_3_cam$h_cred_int_upper==min(T4_3_cam$h_cred_int_upper)),]
T4_3_ad <- T4_3[which(T4_3$hybrid_value2==max(T4_3$hybrid_value2)),]
T4_3_ad_hybrid <- T4_3_ad$hybrid_value2
T4_3_cam_value <- (pi*T4_3_most_cam$height_cm*(T4_3_most_cam$CBH_cm/2*pi)^2)/T4_3_most_cam$age_y
T4_3_ad_value <- (pi*T4_3_ad$height_cm*(T4_3_ad$CBH_cm/2*pi)^2)/T4_3_ad$age_y

pairs_temp <- data.frame(T4_3_cam_value, NA, T4_3_ad_value, T4_3_ad_hybrid)
colnames(pairs_temp) <- c("Cam_Value", "Glob_Value", "Admix_Value", "Hybrid_Value")
pairs <- rbind(pairs, pairs_temp)

T4_5_cam <- T4_5[which(T4_5$h_posterior_mode==min(T4_5$h_posterior_mode)),]
T4_5_most_cam <- T4_5_cam[which(T4_5_cam$h_cred_int_upper==min(T4_5_cam$h_cred_int_upper)),]
T4_5_ad <- T4_5[which(T4_5$hybrid_value2==max(T4_5$hybrid_value2)),]
T4_5_ad_hybrid <- T4_5_ad$hybrid_value2
T4_5_cam_value <- (pi*T4_5_most_cam$height_cm*(T4_5_most_cam$CBH_cm/2*pi)^2)/T4_5_most_cam$age_y
T4_5_ad_value <- (pi*T4_5_ad$height_cm*(T4_5_ad$CBH_cm/2*pi)^2)/T4_5_ad$age_y

pairs_temp <- data.frame(T4_5_cam_value, NA, T4_5_ad_value, T4_5_ad_hybrid)
colnames(pairs_temp) <- c("Cam_Value", "Glob_Value", "Admix_Value", "Hybrid_Value")
pairs <- rbind(pairs, pairs_temp)

T6_2_cam <- T6_2[which(T6_2$h_posterior_mode==min(T6_2$h_posterior_mode)),]
T6_2_most_cam <- T6_2_cam[which(T6_2_cam$h_cred_int_upper==min(T6_2_cam$h_cred_int_upper)),]
T6_2_ad <- T6_2[which(T6_2$hybrid_value2==max(T6_2$hybrid_value2)),]
T6_2_ad_hybrid <- T6_2_ad$hybrid_value2
T6_2_cam_value <- (pi*T6_2_most_cam$height_cm*(T6_2_most_cam$CBH_cm/2*pi)^2)/T6_2_most_cam$age_y
T6_2_ad_value <- (pi*T6_2_ad$height_cm*(T6_2_ad$CBH_cm/2*pi)^2)/T6_2_ad$age_y

pairs_temp <- data.frame(T6_2_cam_value, NA, T6_2_ad_value, T6_2_ad_hybrid)
colnames(pairs_temp) <- c("Cam_Value", "Glob_Value", "Admix_Value", "Hybrid_Value")
pairs <- rbind(pairs, pairs_temp)

T6_3_cam <- T6_3[which(T6_3$h_posterior_mode==min(T6_3$h_posterior_mode)),]
T6_3_most_cam <- T6_3_cam[which(T6_3_cam$h_cred_int_upper==min(T6_3_cam$h_cred_int_upper)),]
T6_3_ad <- T6_3[which(T6_3$hybrid_value2==max(T6_3$hybrid_value2)),]
T6_3_ad_hybrid <- T6_3_ad$hybrid_value2
T6_3_cam_value <- (pi*T6_3_most_cam$height_cm*(T6_3_most_cam$CBH_cm/2*pi)^2)/T6_3_most_cam$age_y
T6_3_ad_value <- (pi*T6_3_ad$height_cm*(T6_3_ad$CBH_cm/2*pi)^2)/T6_3_ad$age_y

pairs_temp <- data.frame(T6_3_cam_value, NA, T6_3_ad_value, T6_3_ad_hybrid)
colnames(pairs_temp) <- c("Cam_Value", "Glob_Value", "Admix_Value", "Hybrid_Value")
pairs <- rbind(pairs, pairs_temp)

T6_4_cam <- T6_4[which(T6_4$h_posterior_mode==min(T6_4$h_posterior_mode)),]
T6_4_most_cam <- T6_4_cam[which(T6_4_cam$h_cred_int_upper==min(T6_4_cam$h_cred_int_upper)),]
T6_4_ad <- T6_4[which(T6_4$hybrid_value2==max(T6_4$hybrid_value2)),]
T6_4_ad_hybrid <- T6_4_ad$hybrid_value2
T6_4_cam_value <- (pi*T6_4_most_cam$height_cm*(T6_4_most_cam$CBH_cm/2*pi)^2)/T6_4_most_cam$age_y
T6_4_ad_value <- (pi*T6_4_ad$height_cm*(T6_4_ad$CBH_cm/2*pi)^2)/T6_4_ad$age_y

pairs_temp <- data.frame(T6_4_cam_value, NA, T6_4_ad_value, T6_4_ad_hybrid)
colnames(pairs_temp) <- c("Cam_Value", "Glob_Value", "Admix_Value", "Hybrid_Value")
pairs <- rbind(pairs, pairs_temp)

T6_6_cam <- T6_6[which(T6_6$h_posterior_mode==min(T6_6$h_posterior_mode)),]
T6_6_most_cam <- T6_6_cam[which(T6_6_cam$h_cred_int_upper==min(T6_6_cam$h_cred_int_upper)),]
T6_6_ad <- T6_6[which(T6_6$hybrid_value2==max(T6_6$hybrid_value2)),]
T6_6_ad_hybrid <- T6_6_ad$hybrid_value2
T6_6_cam_value <- (pi*T6_6_most_cam$height_cm*(T6_6_most_cam$CBH_cm/2*pi)^2)/T6_6_most_cam$age_y
T6_6_ad_value <- (pi*T6_6_ad$height_cm*(T6_6_ad$CBH_cm/2*pi)^2)/T6_6_ad$age_y

pairs_temp <- data.frame(T6_6_cam_value, NA, T6_6_ad_value, T6_6_ad_hybrid)
colnames(pairs_temp) <- c("Cam_Value", "Glob_Value", "Admix_Value", "Hybrid_Value")
pairs <- rbind(pairs, pairs_temp)

T6_7_cam <- T6_7[which(T6_7$h_posterior_mode==min(T6_7$h_posterior_mode)),]
T6_7_most_cam <- T6_7_cam[which(T6_7_cam$h_cred_int_upper==min(T6_7_cam$h_cred_int_upper)),]
T6_7_ad <- T6_7[which(T6_7$hybrid_value2==max(T6_7$hybrid_value2)),]
T6_7_ad_hybrid <- T6_7_ad$hybrid_value2
T6_7_cam_value <- (pi*T6_7_most_cam$height_cm*(T6_7_most_cam$CBH_cm/2*pi)^2)/T6_7_most_cam$age_y
T6_7_ad_value <- (pi*T6_7_ad$height_cm*(T6_7_ad$CBH_cm/2*pi)^2)/T6_7_ad$age_y

pairs_temp <- data.frame(T6_7_cam_value, NA, T6_7_ad_value, T6_7_ad_hybrid)
colnames(pairs_temp) <- c("Cam_Value", "Glob_Value", "Admix_Value", "Hybrid_Value")
pairs <- rbind(pairs, pairs_temp)

T6_8_cam <- T6_8[which(T6_8$h_posterior_mode==min(T6_8$h_posterior_mode)),]
T6_8_most_cam <- T6_8_cam[which(T6_8_cam$h_cred_int_upper==min(T6_8_cam$h_cred_int_upper)),]
T6_8_ad <- T6_8[which(T6_8$hybrid_value2==max(T6_8$hybrid_value2)),]
T6_8_ad_hybrid <- T6_8_ad$hybrid_value2
T6_8_cam_value <- (pi*T6_8_most_cam$height_cm*(T6_8_most_cam$CBH_cm/2*pi)^2)/T6_8_most_cam$age_y
T6_8_ad_value <- (pi*T6_8_ad$height_cm*(T6_8_ad$CBH_cm/2*pi)^2)/T6_8_ad$age_y

pairs_temp <- data.frame(T6_8_cam_value, NA, T6_8_ad_value, T6_8_ad_hybrid)
colnames(pairs_temp) <- c("Cam_Value", "Glob_Value", "Admix_Value", "Hybrid_Value")
pairs <- rbind(pairs, pairs_temp)

T6_10_glob <- T6_10[which(T6_10$h_posterior_mode==max(T6_10$h_posterior_mode)),]
T6_10_most_glob <- T6_10_glob[which(T6_10_glob$h_cred_int_lower==max(T6_10_glob$h_cred_int_lower)),]
T6_10_ad <- T6_10[which(T6_10$hybrid_value2==max(T6_10$hybrid_value2)),]
T6_10_ad_hybrid <- T6_10_ad$hybrid_value2
T6_10_glob_value <- (pi*T6_10_most_glob$height_cm*(T6_10_most_glob$CBH_cm/2*pi)^2)/T6_10_most_glob$age_y
T6_10_ad_value <- (pi*T6_10_ad$height_cm*(T6_10_ad$CBH_cm/2*pi)^2)/T6_10_ad$age_y

pairs_temp <- data.frame(NA, T6_10_glob_value, T6_10_ad_value, T6_10_ad_hybrid)
colnames(pairs_temp) <- c("Cam_Value", "Glob_Value", "Admix_Value", "Hybrid_Value")
pairs <- rbind(pairs, pairs_temp)

####T-Tests, the actually ueful bits!####
pairs_cam <- pairs[which(!is.na(pairs$Cam_Value)),]
pairs_cam$cam_ad_diff <- pairs_cam$Cam_Value-pairs_cam$Admix_Value
t.test(pairs_cam$cam_ad_diff)

pairs_glob <- pairs[which(!is.na(pairs$Glob_Value)),]
pairs_glob$glob_ad_diff <- pairs_glob$Glob_Value-pairs_glob$Admix_Value
t.test(pairs_glob$glob_ad_diff)


#Housekeeping
rm(list=ls())
setwd('C:/Users/samra/Documents/My Documents/Uni/Imperial/Summer Project/R Scripts') #Standard setting of working directory

#Libraries
library(readxl)
library(ggplot2)

######
candidate_loci_an <- read_excel("Candidate_Adaptive_Loci.xlsx", sheet="Sheet1")
candidate_loci_annotated <- candidate_loci_an[which(candidate_loci_an$`NCBI Annotation`!="NA"),]
candidate_loci_exons <- candidate_loci_an[which(candidate_loci_an$`Exon?`=="Y"),]
candindate_loci_exons_scaled <- nrow(candidate_loci_exons)/nrow(candidate_loci_annotated)

random_loci_an_1 <- read_excel("Random_Loci_Samples.xlsx", sheet="Random_1")
random_loci_1_annotated <- random_loci_an_1[which(random_loci_an_1$`NCBI Annotation`!="NA"),]
random_loci_1_exons <- random_loci_an_1[which(random_loci_an_1$`Exon?`=="Y"),]
random_loci_1_exons_scaled <- nrow(random_loci_1_exons)/nrow(random_loci_1_annotated)

random_loci_an_2 <- read_excel("Random_Loci_Samples.xlsx", sheet="Random_2")
random_loci_2_annotated <- random_loci_an_2[which(random_loci_an_2$`NCBI Annotation`!="NA"),]
random_loci_2_exons <- random_loci_an_2[which(random_loci_an_2$`Exon?`=="Y"),]
random_loci_2_exons_scaled <- nrow(random_loci_2_exons)/nrow(random_loci_2_annotated)

random_loci_an_3 <- read_excel("Random_Loci_Samples.xlsx", sheet="Random_3")
random_loci_3_annotated <- random_loci_an_3[which(random_loci_an_3$`NCBI Annotation`!="NA"),]
random_loci_3_exons <- random_loci_an_3[which(random_loci_an_3$`Exon?`=="Y"),]
random_loci_3_exons_scaled <- nrow(random_loci_3_exons)/nrow(random_loci_3_annotated)

annotated_loci <- data.frame(c("candidate", "random_1", "random_2", "random_3"),
           c(nrow(candidate_loci_annotated), nrow(random_loci_1_annotated), nrow(random_loci_2_annotated), nrow(random_loci_3_annotated)),
           c(candindate_loci_exons_scaled, random_loci_1_exons_scaled, random_loci_2_exons_scaled, random_loci_3_exons_scaled))
colnames(annotated_loci) <- c("Subset", "Number", "Proportion_Exon")
barplot(annotated_loci$Number)
barplot(annotated_loci$Proportion_Exon)

######
candidate_loci_GO_class <- read_excel("Candidate_Introgressed_Loci.xlsx", sheet="Sheet6")
total <- nrow(candidate_loci_GO_class)
MF <- nrow(candidate_loci_GO_class[which(candidate_loci_GO_class$`Molecular Function`==1),])
BP <- nrow(candidate_loci_GO_class[which(candidate_loci_GO_class$`Biological Process`==1),])
CC <- nrow(candidate_loci_GO_class[which(candidate_loci_GO_class$`Cellular Component`==1),])
MF_freq <- MF/total
BP_freq <- BP/total 
CC_freq <- CC/total

candidate_loci_MF <- read_excel("Candidate_Introgressed_Loci.xlsx", sheet="Molecular_Function")
candidate_loci_MF[is.na(candidate_loci_MF)] <- 0
MF_annotations <- data.frame(matrix(NA,
                                    nrow=0,
                                    ncol=3))
colnames(MF_annotations) <- c("Annotation", "N", "Freq")
for (i in 4:ncol(candidate_loci_MF)) {
  temp <- data.frame(matrix(NA,
                                      nrow=1,
                                      ncol=3))
  colnames(temp) <- c("Annotation", "N", "Freq")
  temp$Annotation[1] <- colnames(candidate_loci_MF)[(i)]
  temp$N[1] <- sum(candidate_loci_MF[,(i)], na.rm=TRUE)
  if(temp$N>0){
    MF_annotations <- rbind(MF_annotations, temp)
  }
}
for(i in 1:nrow(MF_annotations)){
  MF_annotations$Freq[i] <- MF_annotations$N[i]/nrow(candidate_loci_MF)
}
barplot(MF_annotations$Freq[which(MF_annotations$Annotation!="molecular_function")],
        names.arg = MF_annotations$Annotation[which(MF_annotations$Annotation!="molecular_function")])

candidate_loci_BP <- read_excel("Candidate_Introgressed_Loci.xlsx", sheet="Biological_Process")
candidate_loci_BP[is.na(candidate_loci_BP)] <- 0
BP_annotations <- data.frame(matrix(NA,
                                    nrow=0,
                                    ncol=3))
colnames(BP_annotations) <- c("Annotation", "N", "Freq")
for (i in 4:ncol(candidate_loci_BP)) {
  temp <- data.frame(matrix(NA,
                            nrow=1,
                            ncol=3))
  colnames(temp) <- c("Annotation", "N", "Freq")
  temp$Annotation[1] <- colnames(candidate_loci_BP)[(i)]
  temp$N[1] <- sum(candidate_loci_BP[,(i)], na.rm=TRUE)
  if(temp$N>0){
    BP_annotations <- rbind(BP_annotations, temp)
  }
}
for(i in 1:nrow(BP_annotations)){
  BP_annotations$Freq[i] <- BP_annotations$N[i]/nrow(candidate_loci_BP)
}
barplot(BP_annotations$Freq[which(BP_annotations$Annotation!="biological_process")],
        names.arg = BP_annotations$Annotation[which(BP_annotations$Annotation!="biological_process")])

candidate_loci_CC <- read_excel("Candidate_Introgressed_Loci.xlsx", sheet="Cellular_Component")
candidate_loci_CC[is.na(candidate_loci_CC)] <- 0
CC_annotations <- data.frame(matrix(NA,
                                    nrow=0,
                                    ncol=3))
colnames(CC_annotations) <- c("Annotation", "N", "Freq")
for (i in 4:ncol(candidate_loci_CC)) {
  temp <- data.frame(matrix(NA,
                            nrow=1,
                            ncol=3))
  colnames(temp) <- c("Annotation", "N", "Freq")
  temp$Annotation[1] <- colnames(candidate_loci_CC)[(i)]
  temp$N[1] <- sum(candidate_loci_CC[,(i)], na.rm=TRUE)
  if(temp$N>0){
    CC_annotations <- rbind(CC_annotations, temp)
  }
}
for(i in 1:nrow(CC_annotations)){
  CC_annotations$Freq[i] <- CC_annotations$N[i]/nrow(candidate_loci_CC)
}
barplot(CC_annotations$Freq[which(CC_annotations$Annotation!="cellular_component")],
        names.arg = CC_annotations$Annotation[which(CC_annotations$Annotation!="cellular_component")])

####Graphs####
#Random 1
random_loci_1_GO_class <- read_excel("Random_Loci_Samples.xlsx", sheet="Random_1_GO")
random_1_total <- nrow(random_loci_1_GO_class)
random_1_MF <- nrow(random_loci_1_GO_class[which(random_loci_1_GO_class$`Molecular Function`==1),])
random_1_BP <- nrow(random_loci_1_GO_class[which(random_loci_1_GO_class$`Biological Process`==1),])
random_1_CC <- nrow(random_loci_1_GO_class[which(random_loci_1_GO_class$`Cellular Component`==1),])
random_1_MF_freq <- random_1_MF/random_1_total
random_1_BP_freq <- random_1_BP/random_1_total 
random_1_CC_freq <- random_1_CC/random_1_total

random_loci_1_MF <- read_excel("Random_Loci_Samples.xlsx", sheet="Random_1_Molecular_Function")
random_loci_1_MF[is.na(random_loci_1_MF)] <- 0
MF_annotations_random_1 <- data.frame(matrix(NA,
                                    nrow=0,
                                    ncol=3))
colnames(MF_annotations_random_1) <- c("Annotation", "N", "Freq")
for (i in 4:ncol(random_loci_1_MF)) {
  temp <- data.frame(matrix(NA,
                            nrow=1,
                            ncol=3))
  colnames(temp) <- c("Annotation", "N", "Freq")
  temp$Annotation[1] <- colnames(random_loci_1_MF)[(i)]
  temp$N[1] <- sum(random_loci_1_MF[,(i)], na.rm=TRUE)
  if(temp$N>0){
    MF_annotations_random_1 <- rbind(MF_annotations_random_1, temp)
  }
}
for(i in 1:nrow(MF_annotations_random_1)){
  MF_annotations_random_1$Freq[i] <- MF_annotations_random_1$N[i]/nrow(random_loci_1_MF)
}
barplot(MF_annotations_random_1$Freq[which(MF_annotations_random_1$Annotation!="molecular_function")],
        names.arg = MF_annotations_random_1$Annotation[which(MF_annotations_random_1$Annotation!="molecular_function")])

random_loci_1_BP <- read_excel("Random_Loci_Samples.xlsx", sheet="Random_1_Biological_Process")
random_loci_1_BP[is.na(random_loci_1_BP)] <- 0
BP_annotations_random_1 <- data.frame(matrix(NA,
                                    nrow=0,
                                    ncol=3))
colnames(BP_annotations_random_1) <- c("Annotation", "N", "Freq")
for (i in 4:ncol(random_loci_1_BP)) {
  temp <- data.frame(matrix(NA,
                            nrow=1,
                            ncol=3))
  colnames(temp) <- c("Annotation", "N", "Freq")
  temp$Annotation[1] <- colnames(random_loci_1_BP)[(i)]
  temp$N[1] <- sum(random_loci_1_BP[,(i)], na.rm=TRUE)
  if(temp$N>0){
    BP_annotations_random_1 <- rbind(BP_annotations_random_1, temp)
  }
}
for(i in 1:nrow(BP_annotations_random_1)){
  BP_annotations_random_1$Freq[i] <- BP_annotations_random_1$N[i]/nrow(random_loci_1_BP)
}
barplot(BP_annotations_random_1$Freq[which(BP_annotations_random_1$Annotation!="biological_process")],
        names.arg = BP_annotations_random_1$Annotation[which(BP_annotations_random_1$Annotation!="biological_process")])

random_loci_1_CC <- read_excel("Random_Loci_Samples.xlsx", sheet="Random_1_Cellular_Component")
random_loci_1_CC[is.na(random_loci_1_CC)] <- 0
CC_annotations_random_1 <- data.frame(matrix(NA,
                                    nrow=0,
                                    ncol=3))
colnames(CC_annotations_random_1) <- c("Annotation", "N", "Freq")
for (i in 4:ncol(random_loci_1_CC)) {
  temp <- data.frame(matrix(NA,
                            nrow=1,
                            ncol=3))
  colnames(temp) <- c("Annotation", "N", "Freq")
  temp$Annotation[1] <- colnames(random_loci_1_CC)[(i)]
  temp$N[1] <- sum(random_loci_1_CC[,(i)], na.rm=TRUE)
  if(temp$N>0){
    CC_annotations_random_1 <- rbind(CC_annotations_random_1, temp)
  }
}
for(i in 1:nrow(CC_annotations_random_1)){
  CC_annotations_random_1$Freq[i] <- CC_annotations_random_1$N[i]/nrow(random_loci_1_CC)
}
barplot(CC_annotations_random_1$Freq[which(CC_annotations_random_1$Annotation!="cellular_component")],
        names.arg = CC_annotations_random_1$Annotation[which(CC_annotations_random_1$Annotation!="cellular_component")])

#Random 2
random_loci_2_GO_class <- read_excel("Random_Loci_Samples.xlsx", sheet="Random_2_GO")
random_2_total <- nrow(random_loci_2_GO_class)
random_2_MF <- nrow(random_loci_2_GO_class[which(random_loci_2_GO_class$`Molecular Function`==1),])
random_2_BP <- nrow(random_loci_2_GO_class[which(random_loci_2_GO_class$`Biological Process`==1),])
random_2_CC <- nrow(random_loci_2_GO_class[which(random_loci_2_GO_class$`Cellular Component`==1),])
random_2_MF_freq <- random_2_MF/random_2_total
random_2_BP_freq <- random_2_BP/random_2_total 
random_2_CC_freq <- random_2_CC/random_2_total

random_loci_2_MF <- read_excel("Random_Loci_Samples.xlsx", sheet="Random_2_Molecular_Function")
random_loci_2_MF[is.na(random_loci_2_MF)] <- 0
MF_annotations_random_2 <- data.frame(matrix(NA,
                                             nrow=0,
                                             ncol=3))
colnames(MF_annotations_random_2) <- c("Annotation", "N", "Freq")
for (i in 4:ncol(random_loci_2_MF)) {
  temp <- data.frame(matrix(NA,
                            nrow=1,
                            ncol=3))
  colnames(temp) <- c("Annotation", "N", "Freq")
  temp$Annotation[1] <- colnames(random_loci_2_MF)[(i)]
  temp$N[1] <- sum(random_loci_2_MF[,(i)], na.rm=TRUE)
  if(temp$N>0){
    MF_annotations_random_2 <- rbind(MF_annotations_random_2, temp)
  }
}
for(i in 1:nrow(MF_annotations_random_2)){
  MF_annotations_random_2$Freq[i] <- MF_annotations_random_2$N[i]/nrow(random_loci_2_MF)
}
barplot(MF_annotations_random_2$Freq[which(MF_annotations_random_2$Annotation!="molecular_function")],
        names.arg = MF_annotations_random_2$Annotation[which(MF_annotations_random_2$Annotation!="molecular_function")])

random_loci_2_BP <- read_excel("Random_Loci_Samples.xlsx", sheet="Random_2_Biological_Process")
random_loci_2_BP[is.na(random_loci_2_BP)] <- 0
BP_annotations_random_2 <- data.frame(matrix(NA,
                                             nrow=0,
                                             ncol=3))
colnames(BP_annotations_random_2) <- c("Annotation", "N", "Freq")
for (i in 4:ncol(random_loci_2_BP)) {
  temp <- data.frame(matrix(NA,
                            nrow=1,
                            ncol=3))
  colnames(temp) <- c("Annotation", "N", "Freq")
  temp$Annotation[1] <- colnames(random_loci_2_BP)[(i)]
  temp$N[1] <- sum(random_loci_2_BP[,(i)], na.rm=TRUE)
  if(temp$N>0){
    BP_annotations_random_2 <- rbind(BP_annotations_random_2, temp)
  }
}
for(i in 1:nrow(BP_annotations_random_2)){
  BP_annotations_random_2$Freq[i] <- BP_annotations_random_2$N[i]/nrow(random_loci_2_BP)
}
barplot(BP_annotations_random_2$Freq[which(BP_annotations_random_2$Annotation!="biological_process")],
        names.arg = BP_annotations_random_2$Annotation[which(BP_annotations_random_2$Annotation!="biological_process")])

random_loci_2_CC <- read_excel("Random_Loci_Samples.xlsx", sheet="Random_2_Cellular_Component")
random_loci_2_CC[is.na(random_loci_2_CC)] <- 0
CC_annotations_random_2 <- data.frame(matrix(NA,
                                             nrow=0,
                                             ncol=3))
colnames(CC_annotations_random_2) <- c("Annotation", "N", "Freq")
for (i in 4:ncol(random_loci_2_CC)) {
  temp <- data.frame(matrix(NA,
                            nrow=1,
                            ncol=3))
  colnames(temp) <- c("Annotation", "N", "Freq")
  temp$Annotation[1] <- colnames(random_loci_2_CC)[(i)]
  temp$N[1] <- sum(random_loci_2_CC[,(i)], na.rm=TRUE)
  if(temp$N>0){
    CC_annotations_random_2 <- rbind(CC_annotations_random_2, temp)
  }
}
for(i in 1:nrow(CC_annotations_random_2)){
  CC_annotations_random_2$Freq[i] <- CC_annotations_random_2$N[i]/nrow(random_loci_2_CC)
}
barplot(CC_annotations_random_2$Freq[which(CC_annotations_random_2$Annotation!="cellular_component")],
        names.arg = CC_annotations_random_2$Annotation[which(CC_annotations_random_2$Annotation!="cellular_component")])

#Random 3
random_loci_3_GO_class <- read_excel("Random_Loci_Samples.xlsx", sheet="Random_3_GO")
random_3_total <- nrow(random_loci_3_GO_class)
random_3_MF <- nrow(random_loci_3_GO_class[which(random_loci_3_GO_class$`Molecular Function`==1),])
random_3_BP <- nrow(random_loci_3_GO_class[which(random_loci_3_GO_class$`Biological Process`==1),])
random_3_CC <- nrow(random_loci_3_GO_class[which(random_loci_3_GO_class$`Cellular Component`==1),])
random_3_MF_freq <- random_3_MF/random_3_total
random_3_BP_freq <- random_3_BP/random_3_total 
random_3_CC_freq <- random_3_CC/random_3_total

random_loci_3_MF <- read_excel("Random_Loci_Samples.xlsx", sheet="Random_3_Molecular_Function")
random_loci_3_MF[is.na(random_loci_3_MF)] <- 0
MF_annotations_random_3 <- data.frame(matrix(NA,
                                             nrow=0,
                                             ncol=3))
colnames(MF_annotations_random_3) <- c("Annotation", "N", "Freq")
for (i in 4:ncol(random_loci_3_MF)) {
  temp <- data.frame(matrix(NA,
                            nrow=1,
                            ncol=3))
  colnames(temp) <- c("Annotation", "N", "Freq")
  temp$Annotation[1] <- colnames(random_loci_3_MF)[(i)]
  temp$N[1] <- sum(random_loci_3_MF[,(i)], na.rm=TRUE)
  if(temp$N>0){
    MF_annotations_random_3 <- rbind(MF_annotations_random_3, temp)
  }
}
for(i in 1:nrow(MF_annotations_random_3)){
  MF_annotations_random_3$Freq[i] <- MF_annotations_random_3$N[i]/nrow(random_loci_3_MF)
}
barplot(MF_annotations_random_3$Freq[which(MF_annotations_random_3$Annotation!="molecular_function")],
        names.arg = MF_annotations_random_3$Annotation[which(MF_annotations_random_3$Annotation!="molecular_function")])

random_loci_3_BP <- read_excel("Random_Loci_Samples.xlsx", sheet="Random_3_Biological_Process")
random_loci_3_BP[is.na(random_loci_3_BP)] <- 0
BP_annotations_random_3 <- data.frame(matrix(NA,
                                             nrow=0,
                                             ncol=3))
colnames(BP_annotations_random_3) <- c("Annotation", "N", "Freq")
for (i in 4:ncol(random_loci_3_BP)) {
  temp <- data.frame(matrix(NA,
                            nrow=1,
                            ncol=3))
  colnames(temp) <- c("Annotation", "N", "Freq")
  temp$Annotation[1] <- colnames(random_loci_3_BP)[(i)]
  temp$N[1] <- sum(random_loci_3_BP[,(i)], na.rm=TRUE)
  if(temp$N>0){
    BP_annotations_random_3 <- rbind(BP_annotations_random_3, temp)
  }
}
for(i in 1:nrow(BP_annotations_random_3)){
  BP_annotations_random_3$Freq[i] <- BP_annotations_random_3$N[i]/nrow(random_loci_3_BP)
}
barplot(BP_annotations_random_3$Freq[which(BP_annotations_random_3$Annotation!="biological_process")],
        names.arg = BP_annotations_random_3$Annotation[which(BP_annotations_random_3$Annotation!="biological_process")])

random_loci_3_CC <- read_excel("Random_Loci_Samples.xlsx", sheet="Random_3_Cellular_Component")
random_loci_3_CC[is.na(random_loci_3_CC)] <- 0
CC_annotations_random_3 <- data.frame(matrix(NA,
                                             nrow=0,
                                             ncol=3))
colnames(CC_annotations_random_3) <- c("Annotation", "N", "Freq")
for (i in 4:ncol(random_loci_3_CC)) {
  temp <- data.frame(matrix(NA,
                            nrow=1,
                            ncol=3))
  colnames(temp) <- c("Annotation", "N", "Freq")
  temp$Annotation[1] <- colnames(random_loci_3_CC)[(i)]
  temp$N[1] <- sum(random_loci_3_CC[,(i)], na.rm=TRUE)
  if(temp$N>0){
    CC_annotations_random_3 <- rbind(CC_annotations_random_3, temp)
  }
}
for(i in 1:nrow(CC_annotations_random_3)){
  CC_annotations_random_3$Freq[i] <- CC_annotations_random_3$N[i]/nrow(random_loci_3_CC)
}
barplot(CC_annotations_random_3$Freq[which(CC_annotations_random_3$Annotation!="cellular_component")],
        names.arg = CC_annotations_random_3$Annotation[which(CC_annotations_random_3$Annotation!="cellular_component")])


######
candidate_loci_MF2 <- read_excel("Candidate_Introgressed_Loci.xlsx", sheet="Molecular_Function")
candidate_loci_MF2[is.na(candidate_loci_MF2)] <- 0
MF_annotations2 <- data.frame(matrix(NA,
                                    nrow=0,
                                    ncol=3))
colnames(MF_annotations2) <- c("Annotation", "N", "Freq")
for (i in 4:ncol(candidate_loci_MF2)) {
  temp <- data.frame(matrix(NA,
                            nrow=1,
                            ncol=3))
  colnames(temp) <- c("Annotation", "N", "Freq")
  temp$Annotation[1] <- colnames(candidate_loci_MF2)[(i)]
  temp$N[1] <- sum(candidate_loci_MF2[,(i)], na.rm=TRUE)
  MF_annotations2 <- rbind(MF_annotations2, temp)
}
for(i in 1:nrow(MF_annotations2)){
  MF_annotations2$Freq[i] <- MF_annotations2$N[i]/total
}
MF_annotations2$Data <- "Candidate"

random_loci_1_MF2 <- read_excel("Random_Loci_Samples.xlsx", sheet="Random_1_Molecular_Function")
random_loci_1_MF2[is.na(random_loci_1_MF2)] <- 0
MF_annotations_random_1_2 <- data.frame(matrix(NA,
                                           nrow=0,
                                           ncol=3))
colnames(MF_annotations_random_1_2) <- c("Annotation", "N", "Freq")
for (i in 4:ncol(random_loci_1_MF2)) {
  temp <- data.frame(matrix(NA,
                            nrow=1,
                            ncol=3))
  colnames(temp) <- c("Annotation", "N", "Freq")
  temp$Annotation[1] <- colnames(random_loci_1_MF2)[(i)]
  temp$N[1] <- sum(random_loci_1_MF2[,(i)], na.rm=TRUE)
    MF_annotations_random_1_2 <- rbind(MF_annotations_random_1_2, temp)
}
for(i in 1:nrow(MF_annotations_random_1_2)){
  MF_annotations_random_1_2$Freq[i] <- MF_annotations_random_1_2$N[i]/random_1_total
}
MF_annotations_random_1_2$Data <- "Random Subset 1" 

random_loci_2_MF2 <- read_excel("Random_Loci_Samples.xlsx", sheet="Random_2_Molecular_Function")
random_loci_2_MF2[is.na(random_loci_2_MF2)] <- 0
MF_annotations_random_2_2 <- data.frame(matrix(NA,
                                               nrow=0,
                                               ncol=3))
colnames(MF_annotations_random_2_2) <- c("Annotation", "N", "Freq")
for (i in 4:ncol(random_loci_2_MF2)) {
  temp <- data.frame(matrix(NA,
                            nrow=1,
                            ncol=3))
  colnames(temp) <- c("Annotation", "N", "Freq")
  temp$Annotation[1] <- colnames(random_loci_2_MF2)[(i)]
  temp$N[1] <- sum(random_loci_2_MF2[,(i)], na.rm=TRUE)
  MF_annotations_random_2_2 <- rbind(MF_annotations_random_2_2, temp)
}
for(i in 1:nrow(MF_annotations_random_2_2)){
  MF_annotations_random_2_2$Freq[i] <- MF_annotations_random_2_2$N[i]/random_2_total
}
MF_annotations_random_2_2$Data <- "Random Subset 2" 

random_loci_3_MF2 <- read_excel("Random_Loci_Samples.xlsx", sheet="Random_3_Molecular_Function")
random_loci_3_MF2[is.na(random_loci_3_MF2)] <- 0
MF_annotations_random_3_2 <- data.frame(matrix(NA,
                                               nrow=0,
                                               ncol=3))
colnames(MF_annotations_random_3_2) <- c("Annotation", "N", "Freq")
for (i in 4:ncol(random_loci_3_MF2)) {
  temp <- data.frame(matrix(NA,
                            nrow=1,
                            ncol=3))
  colnames(temp) <- c("Annotation", "N", "Freq")
  temp$Annotation[1] <- colnames(random_loci_3_MF2)[(i)]
  temp$N[1] <- sum(random_loci_3_MF2[,(i)], na.rm=TRUE)
  MF_annotations_random_3_2 <- rbind(MF_annotations_random_3_2, temp)
}
for(i in 1:nrow(MF_annotations_random_3_2)){
  MF_annotations_random_3_2$Freq[i] <- MF_annotations_random_3_2$N[i]/random_3_total
}
MF_annotations_random_3_2$Data <- "Random Subset 3" 

MF2 <- rbind(MF_annotations2, MF_annotations_random_1_2)
MF2 <- rbind(MF2, MF_annotations_random_2_2)
MF2 <- rbind(MF2, MF_annotations_random_3_2)

MF2$Group <- "Random"
for(i in 1:nrow(MF2)){
  if(MF2$Data[i]=="Candidate"){MF2$Group[i] <- "Candidate"}
}

Annotation <- as.data.frame(unique(MF2$Annotation))
for(i in 1:nrow(Annotation)){
  if(MF2$N[which(MF2$Annotation==Annotation$`unique(MF2$Annotation)`[i]&
               MF2$Data=="Candidate")]==0&
     MF2$N[which(MF2$Annotation==Annotation$`unique(MF2$Annotation)`[i]&
                 MF2$Data=="Random Subset 1")]==0&
     MF2$N[which(MF2$Annotation==Annotation$`unique(MF2$Annotation)`[i]&
                 MF2$Data=="Random Subset 2")]==0&
     MF2$N[which(MF2$Annotation==Annotation$`unique(MF2$Annotation)`[i]&
                 MF2$Data=="Random Subset 3")]==0){
    MF2 <- MF2[-c(which(MF2$Annotation==Annotation$`unique(MF2$Annotation)`[i])),]
  }
}
MF2 <- MF2[-c(which(MF2$Annotation=="molecular_function")),]
MF_Plot <- ggplot(data=MF2)+
  geom_col(aes(x=Annotation, y=Freq, fill=Data), position="dodge")+
  theme(axis.text.x=element_text(angle=70, vjust=1.01, hjust=1.1))+
  xlab("Gene Ontology Annotation")+
  ylab("Frequency")+
  ggtitle("A)")+
  ylim(0,1)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background= element_blank(),
        axis.line= element_line(colour="black"))
ggsave("MF_Plot.png", width = 15, height = 5)



candidate_loci_BP2 <- read_excel("Candidate_Introgressed_Loci.xlsx", sheet="Biological_Process")
candidate_loci_BP2[is.na(candidate_loci_BP2)] <- 0
BP_annotations2 <- data.frame(matrix(NA,
                                     nrow=0,
                                     ncol=3))
colnames(BP_annotations2) <- c("Annotation", "N", "Freq")
for (i in 4:ncol(candidate_loci_BP2)) {
  temp <- data.frame(matrix(NA,
                            nrow=1,
                            ncol=3))
  colnames(temp) <- c("Annotation", "N", "Freq")
  temp$Annotation[1] <- colnames(candidate_loci_BP2)[(i)]
  temp$N[1] <- sum(candidate_loci_BP2[,(i)], na.rm=TRUE)
  BP_annotations2 <- rbind(BP_annotations2, temp)
}
for(i in 1:nrow(BP_annotations2)){
  BP_annotations2$Freq[i] <- BP_annotations2$N[i]/total
}
BP_annotations2$Data <- "Candidate"

random_loci_1_BP2 <- read_excel("Random_Loci_Samples.xlsx", sheet="Random_1_Biological_Process")
random_loci_1_BP2[is.na(random_loci_1_BP2)] <- 0
BP_annotations_random_1_2 <- data.frame(matrix(NA,
                                               nrow=0,
                                               ncol=3))
colnames(BP_annotations_random_1_2) <- c("Annotation", "N", "Freq")
for (i in 4:ncol(random_loci_1_BP2)) {
  temp <- data.frame(matrix(NA,
                            nrow=1,
                            ncol=3))
  colnames(temp) <- c("Annotation", "N", "Freq")
  temp$Annotation[1] <- colnames(random_loci_1_BP2)[(i)]
  temp$N[1] <- sum(random_loci_1_BP2[,(i)], na.rm=TRUE)
  BP_annotations_random_1_2 <- rbind(BP_annotations_random_1_2, temp)
}
for(i in 1:nrow(BP_annotations_random_1_2)){
  BP_annotations_random_1_2$Freq[i] <- BP_annotations_random_1_2$N[i]/random_1_total
}
BP_annotations_random_1_2$Data <- "Random Subset 1" 

random_loci_2_BP2 <- read_excel("Random_Loci_Samples.xlsx", sheet="Random_2_Biological_Process")
random_loci_2_BP2[is.na(random_loci_2_BP2)] <- 0
BP_annotations_random_2_2 <- data.frame(matrix(NA,
                                               nrow=0,
                                               ncol=3))
colnames(BP_annotations_random_2_2) <- c("Annotation", "N", "Freq")
for (i in 4:ncol(random_loci_2_BP2)) {
  temp <- data.frame(matrix(NA,
                            nrow=1,
                            ncol=3))
  colnames(temp) <- c("Annotation", "N", "Freq")
  temp$Annotation[1] <- colnames(random_loci_2_BP2)[(i)]
  temp$N[1] <- sum(random_loci_2_BP2[,(i)], na.rm=TRUE)
  BP_annotations_random_2_2 <- rbind(BP_annotations_random_2_2, temp)
}
for(i in 1:nrow(BP_annotations_random_2_2)){
  BP_annotations_random_2_2$Freq[i] <- BP_annotations_random_2_2$N[i]/random_2_total
}
BP_annotations_random_2_2$Data <- "Random Subset 2" 

random_loci_3_BP2 <- read_excel("Random_Loci_Samples.xlsx", sheet="Random_3_Biological_Process")
random_loci_3_BP2[is.na(random_loci_3_BP2)] <- 0
BP_annotations_random_3_2 <- data.frame(matrix(NA,
                                               nrow=0,
                                               ncol=3))
colnames(BP_annotations_random_3_2) <- c("Annotation", "N", "Freq")
for (i in 4:ncol(random_loci_3_BP2)) {
  temp <- data.frame(matrix(NA,
                            nrow=1,
                            ncol=3))
  colnames(temp) <- c("Annotation", "N", "Freq")
  temp$Annotation[1] <- colnames(random_loci_3_BP2)[(i)]
  temp$N[1] <- sum(random_loci_3_BP2[,(i)], na.rm=TRUE)
  BP_annotations_random_3_2 <- rbind(BP_annotations_random_3_2, temp)
}
for(i in 1:nrow(BP_annotations_random_3_2)){
  BP_annotations_random_3_2$Freq[i] <- BP_annotations_random_3_2$N[i]/random_3_total
}
BP_annotations_random_3_2$Data <- "Random Subset 3" 

BP2 <- rbind(BP_annotations2, BP_annotations_random_1_2)
BP2 <- rbind(BP2, BP_annotations_random_2_2)
BP2 <- rbind(BP2, BP_annotations_random_3_2)

BP2$Group <- "Random"
for(i in 1:nrow(BP2)){
  if(BP2$Data[i]=="Candidate"){BP2$Group[i] <- "Candidate"}
}

Annotation <- as.data.frame(unique(BP2$Annotation))
for(i in 1:nrow(Annotation)){
  if(BP2$N[which(BP2$Annotation==Annotation$`unique(BP2$Annotation)`[i]&
                 BP2$Data=="Candidate")]==0&
     BP2$N[which(BP2$Annotation==Annotation$`unique(BP2$Annotation)`[i]&
                 BP2$Data=="Random Subset 1")]==0&
     BP2$N[which(BP2$Annotation==Annotation$`unique(BP2$Annotation)`[i]&
                 BP2$Data=="Random Subset 2")]==0&
     BP2$N[which(BP2$Annotation==Annotation$`unique(BP2$Annotation)`[i]&
                 BP2$Data=="Random Subset 3")]==0){
    BP2 <- BP2[-c(which(BP2$Annotation==Annotation$`unique(BP2$Annotation)`[i])),]
  }
}
BP2 <- BP2[-c(which(BP2$Annotation=="biological_process")),]
BP_Plot <- ggplot(data=BP2)+
  geom_col(aes(x=Annotation, y=Freq, fill=Data), position="dodge")+
  xlab("Gene Ontology Annotation")+
  ylab("Frequency")+
  ggtitle("B)")+
  ylim(0,1)+
  theme(axis.text.x=element_text(angle=70, vjust=1.01, hjust=1.1))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background= element_blank(),
        axis.line= element_line(colour="black"))
ggsave("BP_Plot.png", width = 15, height = 5)



candidate_loci_CC2 <- read_excel("Candidate_Introgressed_Loci.xlsx", sheet="Cellular_Component")
candidate_loci_CC2[is.na(candidate_loci_CC2)] <- 0
CC_annotations2 <- data.frame(matrix(NA,
                                     nrow=0,
                                     ncol=3))
colnames(CC_annotations2) <- c("Annotation", "N", "Freq")
for (i in 4:ncol(candidate_loci_CC2)) {
  temp <- data.frame(matrix(NA,
                            nrow=1,
                            ncol=3))
  colnames(temp) <- c("Annotation", "N", "Freq")
  temp$Annotation[1] <- colnames(candidate_loci_CC2)[(i)]
  temp$N[1] <- sum(candidate_loci_CC2[,(i)], na.rm=TRUE)
  CC_annotations2 <- rbind(CC_annotations2, temp)
}
for(i in 1:nrow(CC_annotations2)){
  CC_annotations2$Freq[i] <- CC_annotations2$N[i]/total
}
CC_annotations2$Data <- "Candidate"

random_loci_1_CC2 <- read_excel("Random_Loci_Samples.xlsx", sheet="Random_1_Cellular_Component")
random_loci_1_CC2[is.na(random_loci_1_CC2)] <- 0
CC_annotations_random_1_2 <- data.frame(matrix(NA,
                                               nrow=0,
                                               ncol=3))
colnames(CC_annotations_random_1_2) <- c("Annotation", "N", "Freq")
for (i in 4:ncol(random_loci_1_CC2)) {
  temp <- data.frame(matrix(NA,
                            nrow=1,
                            ncol=3))
  colnames(temp) <- c("Annotation", "N", "Freq")
  temp$Annotation[1] <- colnames(random_loci_1_CC2)[(i)]
  temp$N[1] <- sum(random_loci_1_CC2[,(i)], na.rm=TRUE)
  CC_annotations_random_1_2 <- rbind(CC_annotations_random_1_2, temp)
}
for(i in 1:nrow(CC_annotations_random_1_2)){
  CC_annotations_random_1_2$Freq[i] <- CC_annotations_random_1_2$N[i]/random_1_total
}
CC_annotations_random_1_2$Data <- "Random Subset 1" 

random_loci_2_CC2 <- read_excel("Random_Loci_Samples.xlsx", sheet="Random_2_Cellular_Component")
random_loci_2_CC2[is.na(random_loci_2_CC2)] <- 0
CC_annotations_random_2_2 <- data.frame(matrix(NA,
                                               nrow=0,
                                               ncol=3))
colnames(CC_annotations_random_2_2) <- c("Annotation", "N", "Freq")
for (i in 4:ncol(random_loci_2_CC2)) {
  temp <- data.frame(matrix(NA,
                            nrow=1,
                            ncol=3))
  colnames(temp) <- c("Annotation", "N", "Freq")
  temp$Annotation[1] <- colnames(random_loci_2_CC2)[(i)]
  temp$N[1] <- sum(random_loci_2_CC2[,(i)], na.rm=TRUE)
  CC_annotations_random_2_2 <- rbind(CC_annotations_random_2_2, temp)
}
for(i in 1:nrow(CC_annotations_random_2_2)){
  CC_annotations_random_2_2$Freq[i] <- CC_annotations_random_2_2$N[i]/random_2_total
}
CC_annotations_random_2_2$Data <- "Random Subset 2" 

random_loci_3_CC2 <- read_excel("Random_Loci_Samples.xlsx", sheet="Random_3_Cellular_Component")
random_loci_3_CC2[is.na(random_loci_3_CC2)] <- 0
CC_annotations_random_3_2 <- data.frame(matrix(NA,
                                               nrow=0,
                                               ncol=3))
colnames(CC_annotations_random_3_2) <- c("Annotation", "N", "Freq")
for (i in 4:ncol(random_loci_3_CC2)) {
  temp <- data.frame(matrix(NA,
                            nrow=1,
                            ncol=3))
  colnames(temp) <- c("Annotation", "N", "Freq")
  temp$Annotation[1] <- colnames(random_loci_3_CC2)[(i)]
  temp$N[1] <- sum(random_loci_3_CC2[,(i)], na.rm=TRUE)
  CC_annotations_random_3_2 <- rbind(CC_annotations_random_3_2, temp)
}
for(i in 1:nrow(CC_annotations_random_3_2)){
  CC_annotations_random_3_2$Freq[i] <- CC_annotations_random_3_2$N[i]/random_3_total
}
CC_annotations_random_3_2$Data <- "Random Subset 3" 

CC2 <- rbind(CC_annotations2, CC_annotations_random_1_2)
CC2 <- rbind(CC2, CC_annotations_random_2_2)
CC2 <- rbind(CC2, CC_annotations_random_3_2)

CC2$Group <- "Random"
for(i in 1:nrow(CC2)){
  if(CC2$Data[i]=="Candidate"){CC2$Group[i] <- "Candidate"}
}

Annotation <- as.data.frame(unique(CC2$Annotation))
for(i in 1:nrow(Annotation)){
  if(CC2$N[which(CC2$Annotation==Annotation$`unique(CC2$Annotation)`[i]&
                 CC2$Data=="Candidate")]==0&
     CC2$N[which(CC2$Annotation==Annotation$`unique(CC2$Annotation)`[i]&
                 CC2$Data=="Random Subset 1")]==0&
     CC2$N[which(CC2$Annotation==Annotation$`unique(CC2$Annotation)`[i]&
                 CC2$Data=="Random Subset 2")]==0&
     CC2$N[which(CC2$Annotation==Annotation$`unique(CC2$Annotation)`[i]&
                 CC2$Data=="Random Subset 3")]==0){
    CC2 <- CC2[-c(which(CC2$Annotation==Annotation$`unique(CC2$Annotation)`[i])),]
  }
}
CC2 <- CC2[-c(which(CC2$Annotation=="cellular_component")),]
CC_Plot <- ggplot(data=CC2)+
  geom_col(aes(x=Annotation, y=Freq, fill=Data), position="dodge")+
  xlab("Gene Ontology Annotation")+
  ylab("Frequency")+
  ylim(0,1)+
  ggtitle("C)")+
  theme(axis.text.x=element_text(angle=70, vjust=1.01, hjust=1.1))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background= element_blank(),
        axis.line= element_line(colour="black"))
ggsave("CC_Plot.png", width = 15, height = 5)

candidate_loci_an <- read_excel("Candidate_Adaptive_Loci.xlsx", sheet="Sheet1")
random_loci_an_1 <- read_excel("Random_Loci_Samples.xlsx", sheet="Random_1")
random_loci_an_2 <- read_excel("Random_Loci_Samples.xlsx", sheet="Random_2")
random_loci_an_3 <- read_excel("Random_Loci_Samples.xlsx", sheet="Random_3")

Genes <- data.frame(matrix(NA, ncol=3, nrow=4))
colnames(Genes) <- c("Dataset", "N_Genes", "N_Exons")
Genes$Dataset[1] <- "Candidate"
Genes$Dataset[2] <- "Random 1"
Genes$Dataset[3] <- "Random 2"
Genes$Dataset[4] <- "Random 3"

Genes$N_Genes[1] <- nrow(candidate_loci_an[which(candidate_loci_an$`NCBI Annotation`!="NA"),])
Genes$N_Genes[2] <- nrow(random_loci_an_1[which(random_loci_an_1$`NCBI Annotation`!="NA"),])
Genes$N_Genes[3] <- nrow(random_loci_an_2[which(random_loci_an_2$`NCBI Annotation`!="NA"),])
Genes$N_Genes[4] <- nrow(random_loci_an_3[which(random_loci_an_3$`NCBI Annotation`!="NA"),])

Genes$N_Exons[1] <- nrow(candidate_loci_an[which(candidate_loci_an$`Exon?`=="Y"),])
Genes$N_Exons[2] <- nrow(random_loci_an_1[which(random_loci_an_1$`Exon?`=="Y"),])
Genes$N_Exons[3] <- nrow(random_loci_an_2[which(random_loci_an_2$`Exon?`=="Y"),])
Genes$N_Exons[4] <- nrow(random_loci_an_3[which(random_loci_an_3$`Exon?`=="Y"),])

Gene_Plot <- ggplot(data=Genes)+
  geom_col(aes(y=N_Genes, x=Dataset), fill= "#0072B2")+
  geom_col(aes(y=N_Exons, x=Dataset), fill= "#D55E00")+
  xlab("Group")+
  ylab("Number of SNPs")+
  #ggtitle("A)")+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background= element_blank(),
        axis.line= element_line(colour="black"))
ggsave("Gene_Plot.png", width = 5, height = 5)

####Housekeeping####
rm(list=ls())
setwd('C:/Users/samra/Documents/My Documents/Uni/Imperial/Summer Project/R Scripts') #Usual setting of WD

####Libraries####
library(vcfR)
library(ggplot2)
library(tidyr)

####Encoding VCF as STRUCTURE input file####
vcf <- read.vcfR("eucal_grandisRef.postfilt.noIndels.vcf.gz", verbose=FALSE) #Loading in the vcf

vcf@gt <- subset(vcf@gt, select=-c(Euc_1620_1))

colnames(vcf@gt)[which(colnames(vcf@gt)=="Euc_1620_2")] <- "Euc_1620"

#Creating matching meta-data for vcfR data
meta.path="eucalyptus_sample_data_2022.csv" #Your metadata file, don't think this is strictly necessary, but makes it easer to change
meta <- read.csv(meta.path) #This just loads that metadata file in
head(meta) #Just to have a look at the meta, possibly not necessary

#This whole lump of stuff here (down to line 133) creates empty metadata rows for sequenced individuals missing matching meta-data records
meta <- subset(meta, !is.na(meta$Sample))
meta$ID <- NA
for(i in 1:nrow(meta)){
  meta$ID[i] <- paste0("Euc_", meta$Sample[i])
}
tree <- as.data.frame(colnames(vcf@gt))
meta3 <- data.frame()
for(i in 2:ncol(vcf@gt)){
  a <- subset(meta, meta$ID==colnames(vcf@gt)[i])
  meta3 <- rbind(meta3, a)
}
#112
new_row <- NA
meta3 <- rbind(meta3, new_row)
meta3$ID[629] <- "Euc_112"
#1186
new_row <- NA
meta3 <- rbind(meta3, new_row)
meta3$ID[630] <- "Euc_1186"
#1190
new_row <- NA
meta3 <- rbind(meta3, new_row)
meta3$ID[631] <- "Euc_1190"
#1209
new_row <- NA
meta3 <- rbind(meta3, new_row)
meta3$ID[632] <- "Euc_1209"
#1297
new_row <- NA
meta3 <- rbind(meta3, new_row)
meta3$ID[633] <- "Euc_1297"
#1574
new_row <- NA
meta3 <- rbind(meta3, new_row)
meta3$ID[634] <- "Euc_1574"
#1867
new_row <- NA
meta3 <- rbind(meta3, new_row)
meta3$ID[637] <- "Euc_1867"
#1891
new_row <- NA
meta3 <- rbind(meta3, new_row)
meta3$ID[638] <- "Euc_1891"
#1898
new_row <- NA
meta3 <- rbind(meta3, new_row)
meta3$ID[639] <- "Euc_1898"
#198
new_row <- NA
meta3 <- rbind(meta3, new_row)
meta3$ID[640] <- "Euc_198"
#2282
new_row <- NA
meta3 <- rbind(meta3, new_row)
meta3$ID[641] <- "Euc_2282"
#248
new_row <- NA
meta3 <- rbind(meta3, new_row)
meta3$ID[642] <- "Euc_248"
#2604
new_row <- NA
meta3 <- rbind(meta3, new_row)
meta3$ID[643] <- "Euc_2604"
#263
new_row <- NA
meta3 <- rbind(meta3, new_row)
meta3$ID[644] <- "Euc_263"
#278
new_row <- NA
meta3 <- rbind(meta3, new_row)
meta3$ID[645] <- "Euc_278"
#2817
new_row <- NA
meta3 <- rbind(meta3, new_row)
meta3$ID[646] <- "Euc_2817"
#2877
new_row <- NA
meta3 <- rbind(meta3, new_row)
meta3$ID[647] <- "Euc_2877"
#304
new_row <- NA
meta3 <- rbind(meta3, new_row)
meta3$ID[648] <- "Euc_304"
#373
new_row <- NA
meta3 <- rbind(meta3, new_row)
meta3$ID[649] <- "Euc_373"
#440
new_row <- NA
meta3 <- rbind(meta3, new_row)
meta3$ID[650] <- "Euc_440"
#447
new_row <- NA
meta3 <- rbind(meta3, new_row)
meta3$ID[651] <- "Euc_447"
#490
new_row <- NA
meta3 <- rbind(meta3, new_row)
meta3$ID[652] <- "Euc_490"

meta4 <- data.frame()
colname <- as.data.frame(colnames(vcf@gt))
for(i in 2:nrow(colname)){
  meta4 <- rbind(meta4, meta3[which(meta3$ID==colname[i,]),])
}
colnames(meta4)[1] <- "Transect"

#Getting right format for the diploid input into STRUCTURE
#This makes a file with a binary 0 for ancestral allele and 1 for derived allele for SNPs
ncol <- ncol(vcf@gt)
cols <- (ncol-1)*2 #Because the first column is the IDs
STR <- data.frame(matrix(0,
                         nrow=nrow(vcf@fix),
                         ncol=cols)) #Making a matrix of 0s of the right size
pb <- txtProgressBar(min= 0,
                     max= ncol(vcf@gt),
                     style= 3,
                     width= 100,
                     char= "=") #Making a progress bar so you can see how its progressing
for(m in 2:ncol(vcf@gt)){
  m1 <- m-1
  m2 <- m1*2
  m21 <- m2-1
  for(i in 1:nrow(vcf@fix)){
    a <- as.data.frame(strsplit(vcf@gt[i,m], ":"))
    if(a[1,]=="1/1"){geno1 <- "1"
    geno2 <- "1"
    STR[i, m21] <- geno1
    STR[i, m2] <- geno2}
    if(a[1,]=="1/0"){geno1 <- "0"
    geno2 <- "1"
    STR[i, m21] <- geno1
    STR[i, m2] <- geno2}
    if(a[1,]=="0/1"){geno1 <- "0"
    geno2 <- "1"
    STR[i, m21] <- geno1
    STR[i, m2] <- geno2}
    if(a[1,]=="./."){geno1 <- "-9"
    geno2 <- "-9"
    STR[i, m21] <- geno1
    STR[i, m2] <- geno2}}
  setTxtProgressBar(pb, m)
}
close(pb)

temp2 <- data.frame()
for(i in 1:nrow(vcf@fix)){
  x <- paste0((vcf@fix[i,1]), "_", (vcf@fix[i,2]))
  temp2 <- rbind(temp2, x)}
names <- temp2[,1]
row.names(STR) <- names #Adding the ID names (think I've pulled them out of the vcf)

vcf_gt <- as.data.frame(vcf@gt)
temp3 <- data.frame()
for(i in 2:ncol(vcf@gt)){
  x <- names(vcf_gt[i])
  temp3 <- rbind(temp3, x)
  temp3 <- rbind(temp3, x)
}
names2 <- temp3[,1]
names(STR) <- names2 #Adding the column names (Loci)

STR <- as.data.frame(t(STR)) #Turning the dataframe around
#So that loci are columns and individuals are rows

STR2 <- STR #Just caching away a spare here, incase I mess it up later

temp_pop <- data.frame()
for(i in 1:nrow(meta4)){
  if(!is.na(meta4$Transect[i]) & meta4$Transect[i]=="T1"){a <- 1}
  if(!is.na(meta4$Transect[i]) & meta4$Transect[i]=="T2"){a <- 2}
  if(!is.na(meta4$Transect[i]) & meta4$Transect[i]=="T3"){a <- 3}
  if(!is.na(meta4$Transect[i]) & meta4$Transect[i]=="T4"){a <- 4}
  if(!is.na(meta4$Transect[i]) & meta4$Transect[i]=="T5"){a <- 5}
  if(!is.na(meta4$Transect[i]) & meta4$Transect[i]=="T6"){a <- 6}
  if(!is.na(meta4$Transect[i]) & meta4$Transect[i]=="T7"){a <- 7}
  if(!is.na(meta4$Transect[i]) & meta4$Transect[i]=="T8"){a <- 8}
  if(is.na(meta4$Transect[i])){a <- 999}
  temp_pop <- rbind(temp_pop, a)
  temp_pop <- rbind(temp_pop, a)
}
STR <- cbind(temp_pop, STR) #This is adding a population column based on the transect my trees were on

STR_mat <- as.matrix(STR)
write.table(STR_mat, "Eucalyptus_STRUCTURE_input_1620.txt", sep="\t", quote=FALSE)
#This can then be opened in notepad
#And the first bit of test (the first column heading, probably called "temp" or "X" or something)
#Is the column heading for the "populations" column, so needs deleting


####Minor Allele Frequency Filters####
vcf <- read.vcfR("eucal_grandisRef.postfilt.noIndels.vcf.gz", verbose=FALSE) #Loading in the vcf

STR <- read.table("Eucalyptus_STRUCTURE_input_1620.txt", 
                  sep="\t", header=TRUE, row.names = 1) #This is reading in the text file I made in the other document

allele_freq <- maf(vcf, element=2)
allele_freq <- as.data.frame(allele_freq)

hist(allele_freq$Frequency, breaks=50)
abline(v=0.01)
abline(v=0.05)
abline(v=0.1)

cutoff_0.01 <- allele_freq[which(allele_freq$Frequency>=0.01),]
cutoff_0.05 <- allele_freq[which(allele_freq$Frequency>=0.05),]
cutoff_0.1 <- allele_freq[which(allele_freq$Frequency>=0.1),]

cut_0.01 <- row.names(cutoff_0.01)
cut_0.05 <- row.names(cutoff_0.05)
cut_0.1 <- row.names(cutoff_0.1)

subset_MAF_0.01 <- as.matrix(STR[,c("X.1", cut_0.01)])
subset_MAF_0.05 <- as.matrix(STR[,c("X.1", cut_0.05)])
subset_MAF_0.1 <- as.matrix(STR[,c("X.1", cut_0.1)])

write.table(subset_MAF_0.01, "Eucalyptus_STRUCTURE_input_MAF_0.01.txt", sep="\t", quote=FALSE)
write.table(subset_MAF_0.05, "Eucalyptus_STRUCTURE_input_MAF_0.05.txt", sep="\t", quote=FALSE)
write.table(subset_MAF_0.1, "Eucalyptus_STRUCTURE_input_MAF_0.1.txt", sep="\t", quote=FALSE)
rm(list = ls())
setwd('C:/Users/samra/Documents/My Documents/Uni/Imperial/Summer Project/R Scripts')

#Loading required libraries####
library(pcadapt)
library(vcfR)

#Loading datasets####
vcf.path="eucal_grandisRef.postfilt.noIndels.vcf/eucal_grandisRef.postfilt.noIndels.vcf"
meta.path="eucalyptus_sample_data_2022.csv"

vcf <- read.vcfR("eucal_grandisRef.postfilt.noIndels.vcf.gz", verbose=FALSE)

genos <- read.pcadapt(vcf.path, type=c("vcf"))
meta <- read.csv(meta.path)
head(meta)

#Subsetting meta to allow matching of VCF and meta-data####
meta2 <- subset(meta, !is.na(meta$Sample))
meta2$ID <- NA
for(i in 1:nrow(meta2)){
  meta2$ID[i] <- paste0("Euc_", meta2$Sample[i])
}

tree <- as.data.frame(colnames(vcf@gt))

meta3 <- data.frame()
for(i in 2:ncol(vcf@gt)){
  a <- subset(meta2, meta2$ID==colnames(vcf@gt)[i])
  meta3 <- rbind(meta3, a)
}

#Filling in missing IDs which lacked meta-data####
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

#1620_1
new_row <- NA
meta3 <- rbind(meta3, new_row)
meta3$ID[635] <- "Euc_1620_1"

#1620_2
new_row <- NA
meta3 <- rbind(meta3, new_row)
meta3$ID[636] <- "Euc_1620_2"

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

for(i in 1:nrow(meta4)){
  if(meta4$Landrace[i]=="Eucalyptus"|
     meta4$Landrace[i]=="Eucalyptus "|
     meta4$Landrace[i]=="Unknown"|
     is.na(meta4$Landrace[i])){
    meta4$Landrace[i]<-"Unknown"}
}

#PCA####
#Plotting the proportion of explained variance for each principal component
x <- pcadapt(input=genos,K=5)
plot(x, option="screeplot")+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background= element_blank(),
        axis.line= element_line(colour="black"))

#PC1 & PC2
plot(x, option="scores", pop=meta4$Landrace, i=1, j=2, plt.pkg="ggplot")+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background= element_blank(),
        axis.line= element_line(colour="black"))

#PC3 & PC4
plot(x, option="scores", pop=meta4$Landrace, i=3, j=4, plt.pkg="ggplot")


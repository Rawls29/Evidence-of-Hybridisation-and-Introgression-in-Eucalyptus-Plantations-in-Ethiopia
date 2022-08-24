####Housekeeping####
rm(list=ls())
setwd('C:/Users/samra/Documents/My Documents/Uni/Imperial/Summer Project/R Scripts') #Usual setting of WD

####Libraries####
#install.packages("devtools", repos = "http://cran.us.r-project.org")
#devtools::install_github("ribailey/gghybrid")
library(gghybrid)
library(ggplot2)
library(writexl)
library(readxl)
library(writexl)
library(gplots)

####Making gghybrid input files for complete set and each transect####
STR <- read.table("Eucalyptus_STRUCTURE_input_MAF_0.1.txt", header=TRUE)
a <- data.frame(matrix(NA,
                       nrow=1302,
                       ncol=1))
a[,1] <- row.names(STR)
names(a) <- c("INDLABEL")
STR <- cbind(a, STR)
row.names(STR) <- NULL
colnames(STR)[which(colnames(STR)=="X.1")] <- "Transect"
STR_out_2 <- read.table("combined_STR_admix.txt")
STR_out_2 <- subset(STR_out_2, select = -c(Transect))
colnames(STR_out_2) <- c("INDLABEL", "Q_Pop_1", "Q_Pop_2")
STR_out_2$Pop <- NA
for(i in 1:nrow(STR_out_2)){
  if(STR_out_2$Q_Pop_1[i]==1){STR_out_2$Pop[i]<-1}else{
    if(STR_out_2$Q_Pop_2[i]==1){STR_out_2$Pop[i]<-2}else{
      STR_out_2$Pop[i]<-3
    }
  }}
STR_out_2 <- STR_out_2[,-c(2,3)]
STR_out_3 <- data.frame(matrix(NA,
                               nrow=1302,
                               ncol=2))
colnames(STR_out_3) <- c("INDLABEL", "Pop")
for(i in 1:nrow(STR_out_2)){
  a2 <- i*2
  a1 <- a2-1
  STR_out_3$INDLABEL[a1] <- STR_out_2$INDLABEL[i]
  STR_out_3$INDLABEL[a2] <- paste0(STR_out_2$INDLABEL[i], ".1")
  STR_out_3$Pop[a1] <- STR_out_2$Pop[i]
  STR_out_3$Pop[a2] <- STR_out_2$Pop[i]
}
STR_out_3 <- merge(STR_out_3, STR, by=c("INDLABEL"))
colnames(STR_out_3)[which(colnames(STR_out_3)=="Pop")] <- "POPID"
STR_out_3$POPID <- as.integer(STR_out_3$POPID)
ID <- data.frame(matrix(NA,
                        nrow=1302,
                        ncol=1))
names(ID) <- c("ID")
nID <- nrow(STR_out_3)/2
for(i in 1:nID){
  a2 <- i*2
  a1 <- a2-1
  ID$ID[a1] <- STR_out_3$INDLABEL[a1]
  ID$ID[a2] <- STR_out_3$INDLABEL[a1]
}
STR_out_3 <- cbind(ID, STR_out_3)

for(i in 1:8){
  gghybrid_in <- STR_out_3[which(STR_out_3$Transect==i),]
  colnames <- as.data.frame(t(colnames(STR_out_3)))
  colnames(colnames) <- colnames(STR_out_3)
  gghybrid_in <- rbind(colnames, gghybrid_in)
  names(gghybrid_in)<- NULL
  rownames(gghybrid_in) <- NULL
  write.table(gghybrid_in, paste0("gghybrid_in_0.1_con_T", i, ".txt"), sep="\t", quote=FALSE)
}

colnames <- as.data.frame(t(colnames(STR_out_3)))
colnames(colnames) <- colnames(STR_out_3)
STR_out_3 <- rbind(colnames, STR_out_3)
names(STR_out_3)<- NULL
rownames(STR_out_3) <- NULL


write.table(STR_out_3, "gghybrid_0.1_df_con.txt", sep="\t", quote=FALSE)

###Investigating locus distribution####
colnames <- gghybrid_input[1,]
loci <- as.character(colnames[,-c(1:5)])
loci_split <- as.data.frame(strsplit(loci, "_"))
loci_split <- as.data.frame(t(loci_split))
colnames(loci_split) <- c("a", "Scaffold", "Location")
Loci_per_scaffold <- as.data.frame(table(loci_split$Scaffold))
barplot(Loci_per_scaffold$Freq)

####gghybrid Analysis for Whole Dataset####
gghybrid_data <- read.data(
  file = "gghybrid_0.1_df_con.txt",
  precol.headers = 1,
  EXTRACOL = 3,
  MISSINGVAL = -9,
  NUMINDS = 651,
  PLOIDY = 2,
  ONEROW = 0,
  nprecol = 5,
  INDLABEL = 1,
  POPID = 1,
  NUMLOCI= 6025
)
save(gghybrid_data, file="gghybrid_0.1_data.R")

prepped_data <- data.prep(gghybrid_data$data,
                          gghybrid_data$loci,
                          gghybrid_data$alleles,
                          marker.info = NULL, #Default
                          S0 = "1",
                          S1 = "2",
                          INDLABEL.name = "INDLABEL",
                          POPID.name = "POPID",
                          sourceAbsent = FALSE, #Default
                          precols=gghybrid_data$precols,
                          max.S.MAF = 0.5, #Default
                          min.diff = 0, #Default
                          min.allele.copies.S0 = 0, #Default
                          min.allele.copies.S1 = 0, #Default
                          AF.CIoverlap = TRUE, #Default
                          return.genotype.table = FALSE, #Default
                          return.locus.table = FALSE #Default
)
save(prepped_data, file="prepped_data_0.1.R")

esth_results <- esth(
  prepped_data$data.prep,
  gghybrid_data$precols,
  test.subject = "ID",
  nitt=3000,
  burnin=1000
)
save(esth_results, file="esth_results_b1000_n3000_0.1.R")

ggcline_results <- ggcline(
  prepped_data$data.prep,
  esth_results,
  esth.colname = "h_posterior_mode",
  test.subject = "locus",
  read.data.precols = gghybrid_data$precols,
  nitt=5000,
  burnin=2000
)
save(ggcline_results, file="ggcline_results_b2000_n5000_0.1.R")

####Histograms in ggplot####
load("esth_results_b1000_n3000_0.1.R")

hi <- as.data.frame(esth_results$hi)
hi_t1 <- subset(hi, hi$Transect==1)
hi_t2 <- subset(hi, hi$Transect==2)
hi_t3 <- subset(hi, hi$Transect==3)
hi_t4 <- subset(hi, hi$Transect==4)
hi_t5 <- subset(hi, hi$Transect==5)
hi_t6 <- subset(hi, hi$Transect==6)
hi_t7 <- subset(hi, hi$Transect==7)
hi_t8 <- subset(hi, hi$Transect==8)

all <- ggplot(data=hi, aes(x=h_posterior_mode ,fill=..x..))+
  geom_histogram(aes(x=h_posterior_mode), binwidth = 0.1)+
  scale_fill_gradient(low="#FF3366", high="#56B4E9")+
  xlab("Hybrid Index Posterior Mode")+
  ylab("Frequency")+
  ggtitle("A) All Transects")+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background= element_blank(),
        axis.line= element_line(colour="black"),
        legend.position="none")
t1 <- ggplot(data=hi_t1, aes(x=h_posterior_mode ,fill=..x..))+
  geom_histogram(aes(x=h_posterior_mode), binwidth = 0.1)+
  scale_fill_gradient(low="#FF3366", high="#56B4E9")+
  xlab("Hybrid Index Posterior Mode")+
  ylab("Frequency")+
  ggtitle("B) Transect 1")+
  ylim(0,90)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background= element_blank(),
        axis.line= element_line(colour="black"),
        legend.position="none")
t2 <- ggplot(data=hi_t2, aes(x=h_posterior_mode ,fill=..x..))+
  geom_histogram(aes(x=h_posterior_mode), binwidth = 0.1)+
  scale_fill_gradient(low="#FF3366", high="#56B4E9")+
  xlab("Hybrid Index Posterior Mode")+
  ylab("Frequency")+
  ggtitle("C) Transect 2")+
  ylim(0,90)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background= element_blank(),
        axis.line= element_line(colour="black"),
        legend.position="none")
t3 <- ggplot(data=hi_t3, aes(x=h_posterior_mode ,fill=..x..))+
  geom_histogram(aes(x=h_posterior_mode), binwidth = 0.1)+
  scale_fill_gradient(low="#FF3366", high="#56B4E9")+
  xlab("Hybrid Index Posterior Mode")+
  ylab("Frequency")+
  ggtitle("D) Transect 3")+
  ylim(0,90)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background= element_blank(),
        axis.line= element_line(colour="black"),
        legend.position="none")
t4 <- ggplot(data=hi_t4, aes(x=h_posterior_mode ,fill=..x..))+
  geom_histogram(aes(x=h_posterior_mode), binwidth = 0.1)+
  scale_fill_gradient(low="#FF3366", high="#56B4E9")+
  xlab("Hybrid Index Posterior Mode")+
  ylab("Frequency")+
  ggtitle("E) Transect 4")+
  ylim(0,90)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background= element_blank(),
        axis.line= element_line(colour="black"),
        legend.position="none")
t5 <- ggplot(data=hi_t5, aes(x=h_posterior_mode ,fill=..x..))+
  geom_histogram(aes(x=h_posterior_mode), binwidth = 0.1)+
  scale_fill_gradient(low="#FF3366", high="#56B4E9")+
  xlab("Hybrid Index Posterior Mode")+
  ylab("Frequency")+
  ggtitle("F) Transect 5")+
  ylim(0,90)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background= element_blank(),
        axis.line= element_line(colour="black"),
        legend.position="none")
t6 <- ggplot(data=hi_t6, aes(x=h_posterior_mode ,fill=..x..))+
  geom_histogram(aes(x=h_posterior_mode), binwidth = 0.1)+
  scale_fill_gradient(low="#FF3366", high="#56B4E9")+
  xlab("Hybrid Index Posterior Mode")+
  ylab("Frequency")+
  ggtitle("G) Transect 6")+
  ylim(0,90)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background= element_blank(),
        axis.line= element_line(colour="black"),
        legend.position="none")
t7 <- ggplot(data=hi_t7, aes(x=h_posterior_mode ,fill=..x..))+
  geom_histogram(aes(x=h_posterior_mode), binwidth = 0.1)+
  scale_fill_gradient(low="#FF3366", high="#56B4E9")+
  xlab("Hybrid Index Posterior Mode")+
  ylab("Frequency")+
  ggtitle("H) Transect 7")+
  ylim(0,90)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background= element_blank(),
        axis.line= element_line(colour="black"),
        legend.position="none")
t8 <- ggplot(data=hi_t8, aes(x=h_posterior_mode ,fill=..x..))+
  geom_histogram(aes(x=h_posterior_mode), binwidth = 0.1)+
  scale_fill_gradient(low="#FF3366", high="#56B4E9")+
  xlab("Hybrid Index Posterior Mode")+
  ylab("Frequency")+
  ggtitle("I) Transect 8")+
  ylim(0,90)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background= element_blank(),
        axis.line= element_line(colour="black"),
        legend.position="none")

# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

multiplot(all, t3, t6, t1, t4, t7, t2, t5, t8, cols=3)

####gghybrid Analysis for Transect 1 and Identification of Putative Introgressed Loci####
gghybrid_data <- read.data(
  file = "gghybrid_in_0.1_con_T1.txt",
  precol.headers = 1,
  EXTRACOL = 3,
  MISSINGVAL = -9,
  NUMINDS = 125,
  PLOIDY = 2,
  ONEROW = 0,
  nprecol = 5,
  INDLABEL = 1,
  POPID = 1,
  NUMLOCI= 6025
)
save(gghybrid_data, file="gghybrid_0.1_T1_data.R")

prepped_data <- data.prep(gghybrid_data$data,
                          gghybrid_data$loci,
                          gghybrid_data$alleles,
                          marker.info = NULL, #Default
                          S0 = "1",
                          S1 = "2",
                          INDLABEL.name = "INDLABEL",
                          POPID.name = "POPID",
                          sourceAbsent = FALSE, #Default
                          precols=gghybrid_data$precols,
                          max.S.MAF = 0.5, #Default
                          min.diff = 0, #Default
                          min.allele.copies.S0 = 0, #Default
                          min.allele.copies.S1 = 0, #Default
                          AF.CIoverlap = TRUE, #Default
                          return.genotype.table = FALSE, #Default
                          return.locus.table = FALSE #Default
)
save(prepped_data, file="prepped_data_0.1_T1.R")

esth_results <- esth(
  prepped_data$data.prep,
  gghybrid_data$precols,
  test.subject = "ID",
  nitt=3000,
  burnin=1000
)
save(esth_results, file="esth_results_b1000_n3000_0.1_T1.R")

ggcline_results <- ggcline(
  prepped_data$data.prep,
  esth_results,
  esth.colname = "h_posterior_mode",
  test.subject = "locus",
  read.data.precols = gghybrid_data$precols,
  nitt=5000,
  burnin=2000
)
save(ggcline_results, file="ggcline_results_b2000_n5000_0.1_T1.R")

load("ggcline_results_b2000_n5000_0.1_T1.R")
big_jumps <- ggcline_results$gc[which(
  (ggcline_results$gc$S0.prop_1-ggcline_results$gc$S1.prop_1)>0.7|
    (ggcline_results$gc$S0.prop_1-ggcline_results$gc$S1.prop_1)<(-0.7)),]
big_jumps_sig_clines <- big_jumps[which((big_jumps$v_lower_95>1&
                                           big_jumps$v_upper_95>1)|
                                          (big_jumps$v_lower_95<1&
                                             big_jumps$v_upper_95<1)),]
big_jumps_sig_clines_centres <- big_jumps_sig_clines[which((big_jumps_sig_clines$centre_lower_95>0.5&
                                                              big_jumps_sig_clines$centre_upper_95>0.5)|
                                                             (big_jumps_sig_clines$centre_lower_95<0.5&
                                                                big_jumps_sig_clines$centre_upper_95<0.5)),]
big_jumps_sig_clines_centres$slope <- NA
for(i in 1:nrow(big_jumps_sig_clines_centres)){
  if(big_jumps_sig_clines_centres$exp_mean_log_v[i]>1){big_jumps_sig_clines_centres$slope[i] <- "steep"}
  if(big_jumps_sig_clines_centres$exp_mean_log_v[i]<1){big_jumps_sig_clines_centres$slope[i] <- "shallow"}
}

big_jumps_sig_clines_centres$fav_pop <- NA
for(i in 1:nrow(big_jumps_sig_clines_centres)){
  if(big_jumps_sig_clines_centres$invlogit_mean_logit_centre[i]>0.5){big_jumps_sig_clines_centres$fav_pop[i] <- "Pop_1"}
  if(big_jumps_sig_clines_centres$invlogit_mean_logit_centre[i]<0.5){big_jumps_sig_clines_centres$fav_pop[i] <- "Pop_0"}
}

table(big_jumps_sig_clines_centres$slope, big_jumps_sig_clines_centres$fav_pop)

write.table(big_jumps_sig_clines_centres, file="sig_SNPs_T1.txt", sep="\t", quote=FALSE)

####gghybrid Analysis for Transect 4 and Identification of Putative Introgressed Loci####
gghybrid_data <- read.data(
  file = "gghybrid_in_0.1_con_T4.txt",
  precol.headers = 1,
  EXTRACOL = 3,
  MISSINGVAL = -9,
  NUMINDS = 90,
  PLOIDY = 2,
  ONEROW = 0,
  nprecol = 5,
  INDLABEL = 1,
  POPID = 1,
  NUMLOCI= 6025
)
save(gghybrid_data, file="gghybrid_0.1_T4_data.R")

prepped_data <- data.prep(gghybrid_data$data,
                          gghybrid_data$loci,
                          gghybrid_data$alleles,
                          marker.info = NULL, #Default
                          S0 = "1",
                          S1 = "2",
                          INDLABEL.name = "INDLABEL",
                          POPID.name = "POPID",
                          sourceAbsent = FALSE, #Default
                          precols=gghybrid_data$precols,
                          max.S.MAF = 0.5, #Default
                          min.diff = 0, #Default
                          min.allele.copies.S0 = 0, #Default
                          min.allele.copies.S1 = 0, #Default
                          AF.CIoverlap = TRUE, #Default
                          return.genotype.table = FALSE, #Default
                          return.locus.table = FALSE #Default
)
save(prepped_data, file="prepped_data_0.1_T4.R")

esth_results <- esth(
  prepped_data$data.prep,
  gghybrid_data$precols,
  test.subject = "ID",
  nitt=3000,
  burnin=1000
)
save(esth_results, file="esth_results_b1000_n3000_0.1_T4.R")

ggcline_results <- ggcline(
  prepped_data$data.prep,
  esth_results,
  esth.colname = "h_posterior_mode",
  test.subject = "locus",
  read.data.precols = gghybrid_data$precols,
  nitt=5000,
  burnin=2000
)
save(ggcline_results, file="ggcline_results_b2000_n5000_0.1_T4.R")

load("ggcline_results_b2000_n5000_0.1_T4.R")

big_jumps <- ggcline_results$gc[which(
  (ggcline_results$gc$S0.prop_1-ggcline_results$gc$S1.prop_1)>0.7|
    (ggcline_results$gc$S0.prop_1-ggcline_results$gc$S1.prop_1)<(-0.7)),]
big_jumps_sig_clines <- big_jumps[which((big_jumps$v_lower_95>1&
                                           big_jumps$v_upper_95>1)|
                                          (big_jumps$v_lower_95<1&
                                             big_jumps$v_upper_95<1)),]
big_jumps_sig_clines_centres <- big_jumps_sig_clines[which((big_jumps_sig_clines$centre_lower_95>0.5&
                                                              big_jumps_sig_clines$centre_upper_95>0.5)|
                                                             (big_jumps_sig_clines$centre_lower_95<0.5&
                                                                big_jumps_sig_clines$centre_upper_95<0.5)),]
big_jumps_sig_clines_centres$slope <- NA
for(i in 1:nrow(big_jumps_sig_clines_centres)){
  if(big_jumps_sig_clines_centres$exp_mean_log_v[i]>1){big_jumps_sig_clines_centres$slope[i] <- "steep"}
  if(big_jumps_sig_clines_centres$exp_mean_log_v[i]<1){big_jumps_sig_clines_centres$slope[i] <- "shallow"}
}

big_jumps_sig_clines_centres$fav_pop <- NA
for(i in 1:nrow(big_jumps_sig_clines_centres)){
  if(big_jumps_sig_clines_centres$invlogit_mean_logit_centre[i]>0.5){big_jumps_sig_clines_centres$fav_pop[i] <- "Pop_1"}
  if(big_jumps_sig_clines_centres$invlogit_mean_logit_centre[i]<0.5){big_jumps_sig_clines_centres$fav_pop[i] <- "Pop_0"}
}

table(big_jumps_sig_clines_centres$slope, big_jumps_sig_clines_centres$fav_pop)

write.table(big_jumps_sig_clines_centres, file="sig_SNPs_T4.txt", sep="\t", quote=FALSE)

####gghybrid Analysis for Transect 6 and Identification of Putative Introgressed Loci####
gghybrid_data <- read.data(
  file = "gghybrid_in_0.1_con_T6.txt",
  precol.headers = 1,
  EXTRACOL = 3,
  MISSINGVAL = -9,
  NUMINDS = 101,
  PLOIDY = 2,
  ONEROW = 0,
  nprecol = 5,
  INDLABEL = 1,
  POPID = 1,
  NUMLOCI= 6025
)
save(gghybrid_data, file="gghybrid_0.1_T6_data.R")

prepped_data <- data.prep(gghybrid_data$data,
                          gghybrid_data$loci,
                          gghybrid_data$alleles,
                          marker.info = NULL, #Default
                          S0 = "1",
                          S1 = "2",
                          INDLABEL.name = "INDLABEL",
                          POPID.name = "POPID",
                          sourceAbsent = FALSE, #Default
                          precols=gghybrid_data$precols,
                          max.S.MAF = 0.5, #Default
                          min.diff = 0, #Default
                          min.allele.copies.S0 = 0, #Default
                          min.allele.copies.S1 = 0, #Default
                          AF.CIoverlap = TRUE, #Default
                          return.genotype.table = FALSE, #Default
                          return.locus.table = FALSE #Default
)
save(prepped_data, file="prepped_data_0.1_T6.R")

esth_results <- esth(
  prepped_data$data.prep,
  gghybrid_data$precols,
  test.subject = "ID",
  nitt=3000,
  burnin=1000
)
save(esth_results, file="esth_results_b1000_n3000_0.1_T6.R")

ggcline_results <- ggcline(
  prepped_data$data.prep,
  esth_results,
  esth.colname = "h_posterior_mode",
  test.subject = "locus",
  read.data.precols = gghybrid_data$precols,
  nitt=5000,
  burnin=2000
)
save(ggcline_results, file="ggcline_results_b2000_n5000_0.1_T6.R")

load("ggcline_results_b2000_n5000_0.1_T6.R")
big_jumps <- ggcline_results$gc[which(
  (ggcline_results$gc$S0.prop_1-ggcline_results$gc$S1.prop_1)>0.7|
    (ggcline_results$gc$S0.prop_1-ggcline_results$gc$S1.prop_1)<(-0.7)),]

big_jumps_sig_clines <- big_jumps[which((big_jumps$v_lower_95>1&
                                           big_jumps$v_upper_95>1)|
                                          (big_jumps$v_lower_95<1&
                                             big_jumps$v_upper_95<1)),]


big_jumps_sig_clines_centres <- big_jumps_sig_clines[which((big_jumps_sig_clines$centre_lower_95>0.5&
                                                              big_jumps_sig_clines$centre_upper_95>0.5)|
                                                             (big_jumps_sig_clines$centre_lower_95<0.5&
                                                                big_jumps_sig_clines$centre_upper_95<0.5)),]
plot_clinecurve(
  big_jumps_sig_clines_centres,
  cline.locus = big_jumps_sig_clines_centres$locus[1],
  locus.column = "locus"
)
plot_clinecurve(
  big_jumps_sig_clines_centres,
  cline.locus = big_jumps_sig_clines_centres$locus[31],
  locus.column = "locus"
)

big_jumps_sig_clines_centres$slope <- NA
for(i in 1:nrow(big_jumps_sig_clines_centres)){
  if(big_jumps_sig_clines_centres$exp_mean_log_v[i]>1){big_jumps_sig_clines_centres$slope[i] <- "steep"}
  if(big_jumps_sig_clines_centres$exp_mean_log_v[i]<1){big_jumps_sig_clines_centres$slope[i] <- "shallow"}
}

big_jumps_sig_clines_centres$fav_pop <- NA
for(i in 1:nrow(big_jumps_sig_clines_centres)){
  if(big_jumps_sig_clines_centres$invlogit_mean_logit_centre[i]>0.5){big_jumps_sig_clines_centres$fav_pop[i] <- "Pop_1"}
  if(big_jumps_sig_clines_centres$invlogit_mean_logit_centre[i]<0.5){big_jumps_sig_clines_centres$fav_pop[i] <- "Pop_0"}
}

table(big_jumps_sig_clines_centres$slope, big_jumps_sig_clines_centres$fav_pop)

write.table(big_jumps_sig_clines_centres, file="sig_SNPs_T6.txt", sep="\t", quote=FALSE)

####Idenifying SNPs that are significant in more than 1 Transect####
T1_snps <- read.table("sig_SNPs_T1.txt", sep="\t")
T4_snps <- read.table("sig_SNPs_T4.txt", sep="\t")
T6_snps <- read.table("sig_SNPs_T6.txt", sep="\t")

T1_4 <- merge(T1_snps, T4_snps, by="locus")
T1_6 <- merge(T1_snps, T6_snps, by="locus")
T4_6 <- merge(T4_snps, T6_snps, by="locus")

full <- rbind(T1_4, T1_6)
full <- rbind(full, T4_6)
candidate_loci <- unique(full$locus)
candidate_loci <- as.data.frame(strsplit(candidate_loci, "_"))
candidate_loci <- as.data.frame(t(candidate_loci))
write_xlsx(candidate_loci, "candidate_loci.xlsx")

####Analyses of Candidate SNPs####
candidate_loci <- read_excel("candidate_loci.xlsx")
candidate_loci$V4 <- paste0(candidate_loci$V1, "_",
                            candidate_loci$V2, "_",
                            candidate_loci$V3)

T1_candidates <- subset(T1_snps, T1_snps$locus %in% candidate_loci$V4)
T4_candidates <- subset(T4_snps, T4_snps$locus %in% candidate_loci$V4)
T6_candidates <- subset(T6_snps, T6_snps$locus %in% candidate_loci$V4)

nrow(T1_candidates[which(T1_candidates$fav_pop=="Pop_1"),])

T1_candidates_reduced <- T1_candidates[c("locus", "slope", "fav_pop")]
colnames(T1_candidates_reduced) <- c("locus", "T1_slope", "T1_fav_pop")

T4_candidates_reduced <- T4_candidates[c("locus", "slope", "fav_pop")]
colnames(T4_candidates_reduced) <- c("locus", "T4_slope", "T4_fav_pop")

T6_candidates_reduced <- T6_candidates[c("locus", "slope", "fav_pop")]
colnames(T6_candidates_reduced) <- c("locus", "T6_slope", "T6_fav_pop")

colnames(candidate_loci) <- c("NW", "Scaffold", "Location", "locus")

candidate_snps1 <- merge(candidate_loci, T1_candidates_reduced, by="locus")
candidate_snps1_4 <- merge(candidate_snps1, T4_candidates_reduced, by="locus")
disagree1_4 <- candidate_snps1_4[which(candidate_snps1_4$T1_slope!=candidate_snps1_4$T4_slope|candidate_snps1_4$T1_fav_pop!=candidate_snps1_4$T4_fav_pop),]
candidate_snps1_6 <- merge(candidate_snps1, T6_candidates_reduced, by="locus")
disagree1_6 <- candidate_snps1_6[which(candidate_snps1_6$T1_slope!=candidate_snps1_6$T6_slope|candidate_snps1_6$T1_fav_pop!=candidate_snps1_6$T6_fav_pop),]
candidate_snps4 <- merge(candidate_loci, T4_candidates_reduced, by="locus")
candidate_snps4_6 <- merge(candidate_snps4, T6_candidates_reduced, by="locus")
disagree4_6 <- candidate_snps4_6[which(candidate_snps4_6$T1_slope!=candidate_snps4_6$T6_slope|candidate_snps4_6$T1_fav_pop!=candidate_snps4_6$T6_fav_pop),]

candidate_snps <- merge(candidate_loci, T1_candidates_reduced, by="locus", all=TRUE)
candidate_snps <- merge(candidate_snps, T4_candidates_reduced, by="locus", all=TRUE)
candidate_snps <- merge(candidate_snps, T6_candidates_reduced, by="locus", all=TRUE)

candidate_snps$slope <- NA
candidate_snps$fav_pop <- NA
for(i in 1:nrow(candidate_snps)){
  if(!is.na(candidate_snps$T1_slope[i])&candidate_snps$T1_slope[i]=="steep"|
     !is.na(candidate_snps$T4_slope[i])&candidate_snps$T4_slope[i]=="steep"|
     !is.na(candidate_snps$T6_slope[i])&candidate_snps$T6_slope[i]=="steep"){candidate_snps$slope[i]<-"steep"}
  if(!is.na(candidate_snps$T1_slope[i])&candidate_snps$T1_slope[i]=="shallow"|
     !is.na(candidate_snps$T4_slope[i])&candidate_snps$T4_slope[i]=="shallow"|
     !is.na(candidate_snps$T6_slope[i])&candidate_snps$T6_slope[i]=="shallow"){candidate_snps$slope[i]<-"shallow"}
  if(!is.na(candidate_snps$T1_fav_pop[i])&candidate_snps$T1_fav_pop[i]=="Pop_0"|
     !is.na(candidate_snps$T4_fav_pop[i])&candidate_snps$T4_fav_pop[i]=="Pop_0"|
     !is.na(candidate_snps$T6_fav_pop[i])&candidate_snps$T6_fav_pop[i]=="Pop_0"){candidate_snps$fav_pop[i]<-"Pop_0"}
  if(!is.na(candidate_snps$T1_fav_pop[i])&candidate_snps$T1_fav_pop[i]=="Pop_1"|
     !is.na(candidate_snps$T4_fav_pop[i])&candidate_snps$T4_fav_pop[i]=="Pop_1"|
     !is.na(candidate_snps$T6_fav_pop[i])&candidate_snps$T6_fav_pop[i]=="Pop_1"){candidate_snps$fav_pop[i]<-"Pop_1"}
}

candidate_snps$adaptive <- "N"
for(i in 1:nrow(candidate_snps)){
  if(candidate_snps$slope[i]=="steep"){candidate_snps$adaptive[i]<-"Y"}
}

adaptive_snps <- candidate_snps[which(candidate_snps$adaptive=="Y"),]

write_xlsx(adaptive_snps, "candidate_adaptive_snps.xlsx")

####Candidates Venn####
candidates <- read_excel("candidate_adaptive_snps.xlsx")

x <- list('Transect_1'=c(candidates$locus[which(!is.na(candidates$T1_slope))]), 
          'Transect_4'=c(candidates$locus[which(!is.na(candidates$T4_slope))]), 
          'Transect_6'=c(candidates$locus[which(!is.na(candidates$T6_slope))]))
venn(x)

####Analysis of All Differential Allele SNPs####
load("ggcline_results_b2000_n5000_0.1_T1.R")
big_jumps_1 <- ggcline_results$gc[which(
  (ggcline_results$gc$S0.prop_1-ggcline_results$gc$S1.prop_1)>0.7|
    (ggcline_results$gc$S0.prop_1-ggcline_results$gc$S1.prop_1)<(-0.7)),]
write.table(big_jumps_1, file="diff_SNPs_T1.txt", sep="\t", quote=FALSE)

load("ggcline_results_b2000_n5000_0.1_T4.R")
big_jumps_4 <- ggcline_results$gc[which(
  (ggcline_results$gc$S0.prop_1-ggcline_results$gc$S1.prop_1)>0.7|
    (ggcline_results$gc$S0.prop_1-ggcline_results$gc$S1.prop_1)<(-0.7)),]
write.table(big_jumps_4, file="diff_SNPs_T4.txt", sep="\t", quote=FALSE)

load("ggcline_results_b2000_n5000_0.1_T6.R")
big_jumps_6 <- ggcline_results$gc[which(
  (ggcline_results$gc$S0.prop_1-ggcline_results$gc$S1.prop_1)>0.7|
    (ggcline_results$gc$S0.prop_1-ggcline_results$gc$S1.prop_1)<(-0.7)),]
write.table(big_jumps_6, file="diff_SNPs_T6.txt", sep="\t", quote=FALSE)

T1_snps <- read.table("diff_SNPs_T1.txt", sep="\t")
T4_snps <- read.table("diff_SNPs_T4.txt", sep="\t")
T6_snps <- read.table("diff_SNPs_T6.txt", sep="\t")

T1_4 <- merge(T1_snps, T4_snps, by="locus")
T1_6 <- merge(T1_snps, T6_snps, by="locus")
T4_6 <- merge(T4_snps, T6_snps, by="locus")

no_T2 <- rbind(T1_4, T1_6)
no_T2 <- rbind(no_T2, T4_6)
diff_loci_no_T2 <- as.data.frame(unique(no_T2$locus))
colnames(diff_loci_no_T2) <- c("locus")

T1_diff <- subset(T1_snps, T1_snps$locus %in% diff_loci_no_T2$locus)
T1_diff$slope <- NA
for(i in 1:nrow(T1_diff)){
  if(T1_diff$v_lower_95[i]>1&T1_diff$v_upper_95[i]>1){T1_diff$slope[i]<-"steep"}else{
    if(T1_diff$v_lower_95[i]<1&T1_diff$v_upper_95[i]<1){T1_diff$slope[i]<-"shallow"}else{
      T1_diff$slope[i]<-"null"
    }
  }
}
T1_diff$fav_pop <- NA
for(i in 1:nrow(T1_diff)){
  if(T1_diff$centre_lower_95[i]>0.5&T1_diff$centre_upper_95[i]>0.5){T1_diff$fav_pop[i]<-"pop_1"}else{
    if(T1_diff$centre_lower_95[i]<0.5&T1_diff$centre_upper_95[i]<0.5){T1_diff$fav_pop[i]<-"pop_0"}else{
      T1_diff$fav_pop[i]<-"null"
    }
  }
}

T4_diff <- subset(T4_snps, T4_snps$locus %in% diff_loci_no_T2$locus)
T4_diff$slope <- NA
for(i in 1:nrow(T4_diff)){
  if(T4_diff$v_lower_95[i]>1&T4_diff$v_upper_95[i]>1){T4_diff$slope[i]<-"steep"}else{
    if(T4_diff$v_lower_95[i]<1&T4_diff$v_upper_95[i]<1){T4_diff$slope[i]<-"shallow"}else{
      T4_diff$slope[i]<-"null"
    }
  }
}
T4_diff$fav_pop <- NA
for(i in 1:nrow(T4_diff)){
  if(T4_diff$centre_lower_95[i]>0.5&T4_diff$centre_upper_95[i]>0.5){T4_diff$fav_pop[i]<-"pop_1"}else{
    if(T4_diff$centre_lower_95[i]<0.5&T4_diff$centre_upper_95[i]<0.5){T4_diff$fav_pop[i]<-"pop_0"}else{
      T4_diff$fav_pop[i]<-"null"
    }
  }
}

T6_diff <- subset(T6_snps, T6_snps$locus %in% diff_loci_no_T2$locus)
T6_diff$slope <- NA
for(i in 1:nrow(T6_diff)){
  if(T6_diff$v_lower_95[i]>1&T6_diff$v_upper_95[i]>1){T6_diff$slope[i]<-"steep"}else{
    if(T6_diff$v_lower_95[i]<1&T6_diff$v_upper_95[i]<1){T6_diff$slope[i]<-"shallow"}else{
      T6_diff$slope[i]<-"null"
    }
  }
}
T6_diff$fav_pop <- NA
for(i in 1:nrow(T6_diff)){
  if(T6_diff$centre_lower_95[i]>0.5&T6_diff$centre_upper_95[i]>0.5){T6_diff$fav_pop[i]<-"pop_1"}else{
    if(T6_diff$centre_lower_95[i]<0.5&T6_diff$centre_upper_95[i]<0.5){T6_diff$fav_pop[i]<-"pop_0"}else{
      T6_diff$fav_pop[i]<-"null"
    }
  }
}

T1_diff_reduced <- T1_diff[c("locus", "slope", "fav_pop")]
colnames(T1_diff_reduced) <- c("locus", "T1_slope", "T1_fav_pop")

T4_diff_reduced <- T4_diff[c("locus", "slope", "fav_pop")]
colnames(T4_diff_reduced) <- c("locus", "T4_slope", "T4_fav_pop")

T6_diff_reduced <- T6_diff[c("locus", "slope", "fav_pop")]
colnames(T6_diff_reduced) <- c("locus", "T6_slope", "T6_fav_pop")

diff_snps1 <- merge(diff_loci_no_T2, T1_diff_reduced, by="locus")
diff_snps1_4 <- merge(diff_snps1, T4_diff_reduced, by="locus")
disagree1_4 <- diff_snps1_4[which(diff_snps1_4$T1_slope!=diff_snps1_4$T4_slope|diff_snps1_4$T1_fav_pop!=diff_snps1_4$T4_fav_pop),]
diff_snps1_6 <- merge(diff_snps1, T6_diff_reduced, by="locus")
disagree1_6 <- diff_snps1_6[which(diff_snps1_6$T1_slope!=diff_snps1_6$T6_slope|diff_snps1_6$T1_fav_pop!=diff_snps1_6$T6_fav_pop),]
diff_snps4 <- merge(diff_loci_no_T2, T4_diff_reduced, by="locus")
diff_snps4_6 <- merge(diff_snps4, T6_diff_reduced, by="locus")
disagree4_6 <- diff_snps4_6[which(diff_snps4_6$T1_slope!=diff_snps4_6$T6_slope|diff_snps4_6$T1_fav_pop!=diff_snps4_6$T6_fav_pop),]

diff_snps <- merge(diff_loci_no_T2, T1_diff_reduced, by="locus", all=TRUE)
diff_snps <- merge(diff_snps, T4_diff_reduced, by="locus", all=TRUE)
diff_snps <- merge(diff_snps, T6_diff_reduced, by="locus", all=TRUE)

diff_snps$slope <- NA
diff_snps$fav_pop <- NA
for(i in 1:nrow(diff_snps)){
  if(!is.na(diff_snps$T1_slope[i])&diff_snps$T1_slope[i]=="steep"|
     !is.na(diff_snps$T4_slope[i])&diff_snps$T4_slope[i]=="steep"|
     !is.na(diff_snps$T6_slope[i])&diff_snps$T6_slope[i]=="steep"){diff_snps$slope[i]<-"steep"}else{
  if(!is.na(diff_snps$T1_slope[i])&diff_snps$T1_slope[i]=="shallow"|
     !is.na(diff_snps$T4_slope[i])&diff_snps$T4_slope[i]=="shallow"|
     !is.na(diff_snps$T6_slope[i])&diff_snps$T6_slope[i]=="shallow"){diff_snps$slope[i]<-"shallow"}else{
       diff_snps$slope<-"null"
     }}
  if(!is.na(diff_snps$T1_fav_pop[i])&diff_snps$T1_fav_pop[i]=="Pop_0"|
     !is.na(diff_snps$T4_fav_pop[i])&diff_snps$T4_fav_pop[i]=="Pop_0"|
     !is.na(diff_snps$T6_fav_pop[i])&diff_snps$T6_fav_pop[i]=="Pop_0"){diff_snps$fav_pop[i]<-"Pop_0"}else{
  if(!is.na(diff_snps$T1_fav_pop[i])&diff_snps$T1_fav_pop[i]=="Pop_1"|
     !is.na(diff_snps$T4_fav_pop[i])&diff_snps$T4_fav_pop[i]=="Pop_1"|
     !is.na(diff_snps$T6_fav_pop[i])&diff_snps$T6_fav_pop[i]=="Pop_1"){diff_snps$fav_pop[i]<-"Pop_1"}else{
       diff_snps$fav_pop <- "null"
     }}
}

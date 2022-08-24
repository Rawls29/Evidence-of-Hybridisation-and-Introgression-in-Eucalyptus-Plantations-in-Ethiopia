rm(list=ls())
setwd('C:/Users/samra/Documents/My Documents/Uni/Imperial/Summer Project/R Scripts') #Standard setting of working directory
###Libraries####
library("writexl")

###Input####
gghybrid_input <- read.table("gghybrid_0.1_df_con.txt")

###Randomly sampling 91 loci####
colnames <- gghybrid_input[1,]
loci <- colnames[,-c(1:5)]

random_numbers <- sample(1:ncol(loci), 91)
random_sample <- as.character(loci[,c(random_numbers)])
random_sample <- as.data.frame(strsplit(random_sample, "_"))
random_sample <- as.data.frame(t(random_sample))

write_xlsx(random_sample, "gghybrid_input_random_sample.xlsx")

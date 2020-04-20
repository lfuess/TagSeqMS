##This script will be used to prepare data for DESEQ2##

library(dplyr)
library(janitor)
library(data.table)

initialdata = read.csv("Masterdata_TagSeq_SDH_23July2018.csv", check.names= FALSE)
dim(initialdata);
names(initialdata);

##Now to select the data we want##

data = as.data.frame(initialdata[c(2,4,25:28,30:31,34,54)]);
names(data)

dim(data)
head(data)

##input the list of samples you are using##

##we use a file hand currated which lists all the samples that passed sequencing and assembly##
samples = read.csv("Good_Samples.csv", check.names = FALSE)

names(samples)
dim(samples)

##and combine the files to reduce down to just the good samples##

merged=merge(samples, data, by.x = "Sample", by.y = "sample_ID")
dim(merged) ##should have 408 lines remaining##

##Select everything but F1s
GBC=merged[merged[,5] == "GBC",]
RBC=merged[merged[,5] == "RBC",]
F2=merged[merged[,5] == "F2",]

final=rbind(GBC,RBC,F2)

dim(final) ##should be 393##

##order##
final1=final[order(final$Sample),]

##change NAs in sex to Unkown##
levels <- levels(final1$Sex)
levels[length(levels) + 1] <- "U"
final1$Sex <- factor(final1$Sex, levels = levels)
final1[c("Sex")][is.na(final1[c("Sex")])] <- "U"

##remove rows with NAs##
final = final1[complete.cases(final1), ]

##write out the experimental design file##
write.csv(final, file = "ExpDesign.csv", row.names=FALSE)

##last step is to refine our read count matrix##
##get a list of new good samples (samples which didn't have NAs for any of our model factors)##
gs = final[,c(1:2)]

##merge that with your read count matrix, to get a read count matrix with only good samples and no F1s##
reads = as.data.frame(t(read.csv("allcounts.csv", check.names = FALSE)))
reads[1:10,1:10]
  ##fix names##
  reads = reads %>%
    row_to_names(row_number = 1)
  reads = setDT(reads, keep.rownames = TRUE)[]
  ##merge##
  merged=merge(reads, gs, by.x = "rn", by.y = "Sample")
  merged=merged[,c(1:25846)]
  ##write it out##
  finalcounts=as.data.frame(t(merged))
  finalcounts[1:10,1:10]
  finalcounts = finalcounts %>%
    row_to_names(row_number = 1)
  write.csv(finalcounts, "allcounts_final.csv")  
  
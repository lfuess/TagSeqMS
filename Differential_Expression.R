##this script will take our processed read count file and run differential expression analysis

#BiocManager::install("DESeq2")
library("DESeq2")

##input our data##

countData <- read.csv("allcounts_noF1noNA.csv", row.names="Gene", check.names= FALSE)
colData <- read.csv("ExpDesign_NoNA.csv") ##note I set all NAs here for Sex to U for unknown (in Excel.. need to add that to code)##

##read normalization- if needed for other applications (WGCNA, heatmaps)##

#dds <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~ Fused.Organs) ##doesn't matter##
#dds <- estimateSizeFactors(dds) 
#dds <- estimateDispersions(dds)

#vst <- getVarianceStabilizedData(dds)

#write.csv(vst, file = "normalizedreads_test.csv")


##Alright let's try our first model##

names(colData)
#colnames(countData)=NULL

##build a model matrix with 0 rows removed##
##set the model##
dds <- model.matrix( ~ library_name + WormSource + CrossDir + worm_present + Fused.Organs + BG_Gadj_MFI + Sex + CrossDir*worm_present + CrossDir*Fused.Organs, colData)
##use code from the manual vignette remove rows with 0s##
colnames(dds)
unname(dds)
all.zero <- apply(dds, 2, function(x) all (x==0))
all.zero
idx <- which(all.zero)
dds <- dds[,-idx]
unname(dds)

##now model it, following suggestions for large models with rows which do not converge in beta##
object = DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~ library_name + WormSource + CrossDir + worm_present + Fused.Organs + BG_Gadj_MFI + Sex + CrossDir*worm_present + CrossDir*Fused.Organs, ignoreRank = TRUE)
  ##filter out rows where gene is expressed in less than half of the samples##
  object <- estimateSizeFactors(object)
  nc <- counts(object, normalized=TRUE)
  object <- object[ rowSums(nc > 0) >= 195 ]

##actually run the analysis; do it in two steps to increase power and get more rows to converge##
object <- estimateDispersions(object, modelMatrix = dds)
#save(colData, countData, dds, object, file = "Stickeblack-Frozed.RData")
#load("Stickeblack-Frozed.RData")
dds1 <- nbinomWaldTest(object, maxit=500, modelMatrix = dds)

#save(colData, countData, dds, dds1, object, file = "Stickeblack-DEseq-Full-StringParam.RData")
#load("Stickeblack-DEseq-Full-StringParam_nosexworm.RData")

##converte the results and remove any rows which still didn't converge##
res <- results(dds1)
res = res[mcols(dds1)$fullBetaConv,]

##Look at Infection##
resinfect <- results(dds1,contrast=list(c("worm_presentTRUE")))
resorderd <-resinfect[order(resinfect$padj),] 
head(resorderd, 310)
write.csv(resorderd, file = "DESeq_Infection_Full_StringParam_noNAs.csv") ##71 affected

##Next GBC v F2##
resGvF2 <- results(dds1,contrast=list(c("CrossDirGBC")))
resorderd <-resGvF2[order(resGvF2$padj),] 
head(resorderd,6000)
write.csv(resorderd, file = "DESeq_GBCvF2_Full_StringParam_noNAs.csv") ##4132 affected

##Next RBC v F2##
resRvF2 <- results(dds1,contrast=list(c("CrossDirRBC")))
resorderd <-resRvF2[order(resRvF2$padj),] 
head(resorderd, 8975) #8974 affected##
write.csv(resorderd, file = "DESeq_RBCvF2_Full_StringParam_noNAs.csv") ## write CSV ##

##Next GBC v RBC##
resGvR <- results(dds1,contrast=list(c("CrossDirGBC"),  c("CrossDirRBC")))
resorderd <-resGvR[order(resGvR$padj),] 
head(resorderd, 846) #845 affected##
write.csv(resorderd, file = "DESeq_RBCvGBC_Full_StringParam_noNAs.csv") ## write CSV ##

##Next Fibrosis##
resFib <- results(dds1,contrast=list(c("Fused.OrgansTRUE")))
resorderd <-resFib[order(resFib$padj),] 
head(resorderd, 10) #0 affected##
write.csv(resorderd, file = "DESeq_Fibrosis_Full_StringParam_noNAs.csv") ## write CSV ##

##ROS##
resROS <- results(dds1,contrast=list(c("BG_Gadj_MFI")))
resorderd <-resROS[order(resROS$padj),] 
head(resorderd, 10) #0 affected##
write.csv(resorderd, file = "DESeq_ROS_Full_StringParam_noNAs.csv")

##Sex##
resSex <- results(dds1,contrast=list(c("SexM")))
resorderd <-resSex[order(resSex$padj),] 
head(resorderd, 500) #498 affected##
write.csv(resorderd, file = "DESeq_Sex_Full_StringParam_noNAs.csv")

##Worm##
resWorm <- results(dds1,contrast=list(c("WormSourceGos")))
resorderd <-resWorm[order(resWorm$padj),] 
head(resorderd, 9405) #9404 affected##
write.csv(resorderd, file = "DESeq_Worm_Full_StringParam_noNAs.csv")

##Moving on to Interactive Effects, starting with Cross by Infection##
resCrossInfect <- results(dds1,contrast=list(c("CrossDirGBC.worm_presentTRUE")))
resorderd <-resCrossInfect[order(resCrossInfect$padj),] 
head(resorderd, 10) #0 affected##
write.csv(resorderd, file = "DESeq_FvGbyInfect_Full_StringParam_noNAs.csv")

resCrossInfect1 <- results(dds1,contrast=list(c("CrossDirRBC.worm_presentTRUE")))
resorderd <-resCrossInfect1[order(resCrossInfect1$padj),] 
head(resorderd, 105) #2 affected##
write.csv(resorderd, file = "DESeq_FvRbyInfect_Full_StringParam_noNAs.csv")## write CSV ##

resROSCross3 <- results(dds1,contrast=list(c("CrossDirGBC.worm_presentTRUE"),  c("CrossDirRBC.worm_presentTRUE")))
resorderd <-resROSCross3[order(resROSCross3$padj),] 
head(resorderd, 10) #2 affected##
write.csv(resorderd, file = "DESeq_GvRbyInfect_Full_StringParam_noNAs.csv")

##Cross by Fibrosis##
resCrossFibrosis <- results(dds1,contrast=list(c("CrossDirRBC.Fused.OrgansTRUE")))
resorderd <-resCrossFibrosis[order(resCrossFibrosis$padj),] 
head(resorderd, 10) #0 affected##
write.csv(resorderd, file = "DESeq_CrossbyFibrosis_Full_StringParam_noNAs.csv")

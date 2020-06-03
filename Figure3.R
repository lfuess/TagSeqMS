##lets make figure 3 for the MS##
library(tidyverse)
library(data.table)

##first we have to process the data##
  ##input reads##
  reads = read.csv("normalizedreads.csv", check.names = FALSE)
  colnames(reads)[1] <- "Gene"
  row.names(reads) <- reads$Gene
  reads[1] <- NULL
  reads[1:10,1:10]
  ##select sig genes##
  sig = reads[c("ENSGACT00000015170.1", "ENSGACT00000018024.1","ENSGACT00000024555.1","ENSGACT00000003040.1"),]
  sig = as.data.frame(t(sig))
  sig=setDT(sig, keep.rownames = TRUE)[]
  colnames(sig)[1] <- "Sample"
  ##input metadata##
  meta = read.csv("ExpDesign.csv")
  meta = meta[,c(1,3,6)]
  ##select RBC##
  RBC= meta[meta$CrossDir == "RBC", ]
  RBC_Infect = RBC[RBC$worm_present == "TRUE", ]  
  RBC_Infect <- as.vector(RBC_Infect['Sample'])
  RBC_UnInfect = RBC[RBC$worm_present == "FALSE", ]  
  RBC_UnInfect <- as.vector(RBC_UnInfect['Sample'])
  ##Select GBC##
  GBC= meta[meta$CrossDir == "GBC", ]
  GBC_Infect = GBC[GBC$worm_present == "TRUE", ]  
  GBC_Infect <- as.vector(GBC_Infect['Sample'])
  GBC_UnInfect = GBC[GBC$worm_present == "FALSE", ]  
  GBC_UnInfect <- as.vector(GBC_UnInfect['Sample'])
  
##now we can subset the data to create averages and stderr##
  RBC_Infect = merge(sig,RBC_Infect, by="Sample")
  RBC_Infect = RBC_Infect %>% remove_rownames %>% column_to_rownames(var="Sample")
  RBC_Infect[c("average"),]=colMeans(RBC_Infect, na.rm = FALSE, dims = 1)  
  RBC_Infect[c("se"),]=sapply(RBC_Infect,function(x)sd(x)/sqrt(length(x)))
  RBC_Infect[c("Cross"),]="RBC"  
  RBC_Infect[c("Infect"),]="TRUE"
  RBC_Infect[nrow(RBC_Infect) + 1,] = c("CCN3","Unannotated","SH2D1A","LSM10")
  rownames(RBC_Infect)[29]<-"name"
  RBC_Infect=RBC_Infect[c(25:29),]  
  RBC_Infect=as.data.frame(t(RBC_Infect))
  
  
  GBC_Infect = merge(sig,GBC_Infect, by="Sample")
  GBC_Infect = GBC_Infect %>% remove_rownames %>% column_to_rownames(var="Sample")
  GBC_Infect[c("average"),]=colMeans(GBC_Infect, na.rm = FALSE, dims = 1)  
  GBC_Infect[c("se"),]=sapply(GBC_Infect,function(x)sd(x)/sqrt(length(x)))
  GBC_Infect[c("Cross"),]="GBC"  
  GBC_Infect[c("Infect"),]="TRUE"
  GBC_Infect[nrow(GBC_Infect) + 1,] = c("CCN3","Unannotated","SH2D1A","LSM10")
  rownames(GBC_Infect)[25]<-"name"
  GBC_Infect=GBC_Infect[c(21:25),]  
  GBC_Infect=as.data.frame(t(GBC_Infect))
  
  
  
  RBC_UnInfect = merge(sig,RBC_UnInfect, by="Sample")
  RBC_UnInfect = RBC_UnInfect %>% remove_rownames %>% column_to_rownames(var="Sample")
  RBC_UnInfect[c("average"),]=colMeans(RBC_UnInfect, na.rm = FALSE, dims = 1)  
  RBC_UnInfect[c("se"),]=sapply(RBC_UnInfect,function(x)sd(x)/sqrt(length(x)))
  RBC_UnInfect[c("Cross"),]="RBC"  
  RBC_UnInfect[c("Infect"),]="FALSE"
  RBC_UnInfect[nrow(RBC_UnInfect) + 1,] = c("CCN3","Unannotated","SH2D1A","LSM10")
  rownames(RBC_UnInfect)[73]<-"name"
  RBC_UnInfect=RBC_UnInfect[c(69:73),]  
  RBC_UnInfect=as.data.frame(t(RBC_UnInfect))
  
  
  GBC_UnInfect = merge(sig,GBC_UnInfect, by="Sample")
  GBC_UnInfect = GBC_UnInfect %>% remove_rownames %>% column_to_rownames(var="Sample")
  GBC_UnInfect[c("average"),]=colMeans(GBC_UnInfect, na.rm = FALSE, dims = 1)  
  GBC_UnInfect[c("se"),]=sapply(GBC_UnInfect,function(x)sd(x)/sqrt(length(x)))
  GBC_UnInfect[c("Cross"),]="GBC"  
  GBC_UnInfect[c("Infect"),]="FALSE"
  GBC_UnInfect[nrow(GBC_UnInfect) + 1,] = c("CCN3","Unannotated","SH2D1A","LSM10")
  rownames(GBC_UnInfect)[70]<-"name"
  GBC_UnInfect=GBC_UnInfect[c(66:70),]  
  GBC_UnInfect=as.data.frame(t(GBC_UnInfect))

  ##merge them all together to get your final data frame for graphing!!
  data = rbind(RBC_Infect, RBC_UnInfect, GBC_Infect, GBC_UnInfect)

  
##Actual Graphing Part##
  levels(data$Infect)[levels(data$Infect)=="FALSE"] <- "Uninfected"
  levels(data$Infect)[levels(data$Infect)=="TRUE"] <- "Infected"
  data$se=as.numeric(as.character(data$se))
  data$average=as.numeric(as.character(data$average))

  data$Infect <- factor(data$Infect, levels = rev(levels(factor(data$Infect))))
  
  gp1 <- ggplot(data, aes(x=Infect, y=average, colour=Cross, group=Cross))
  gp1 + geom_line(aes(linetype=Cross), size=.6) + scale_colour_manual(values=c("#f3c558","#6dae90")) +
    geom_point(aes(shape=Cross), size=3) + scale_linetype_manual(values=c("solid", "solid")) + 
    scale_shape_manual(values=c(17, 18)) + geom_errorbar(aes(ymax=average+se, ymin=average-se), width=.1) + 
    facet_wrap(~name, ncol=2, scales="free") + theme_bw() +
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
    ylab("Normalized Expression") + xlab("Infection") 

##This script will parse results from DESeq2 to make files of just significant genes, and files for IPA##
##input is the files output by the previous script: Differential_Expression.R##

##load anotations##
annos = read.csv("sticklemasterannos_219.csv")
names(annos)
annos = annos[,c(2,1,7)]

##go factor by factor##
  
  ##read in data##
  infect = read.csv("DESeq_Infection_Full_StringParam_new.csv")
  infect_annotated = merge(infect, annos, by.x = "X", by.y = "gene", all.x = TRUE)
  ##select significant genes##
  infect_annotated_sig = infect_annotated[(infect_annotated[,7]<.1),]
  infect_annotated_sig = infect_annotated_sig[rowSums(is.na(infect_annotated_sig)) != ncol(infect_annotated_sig), ]
  ##write it out##
  write.csv(infect_annotated_sig, "Infection_Sig_Genes.csv", row.names = FALSE)
  
  ##read in data##
  FvG = read.csv("DESeq_GBCvF2_Full_StringParam_new.csv")
  FvG_annotated = merge(FvG, annos, by.x = "X", by.y = "gene", all.x = TRUE)
  ##select significant genes##
  FvG_annotated_sig = FvG_annotated[(FvG_annotated[,7]<.1),]
  FvG_annotated_sig = FvG_annotated_sig[rowSums(is.na(FvG_annotated_sig)) != ncol(FvG_annotated_sig), ]
  ##write it out##
  write.csv(FvG_annotated_sig, "FvG_Sig_Genes.csv", row.names = FALSE)
  
  ##read in data##
  FvR = read.csv("DESeq_F2vsRBC_Full_StringParam_new.csv")
  FvR_annotated = merge(FvR, annos, by.x = "X", by.y = "gene", all.x = TRUE)
  ##select significant genes##
  FvR_annotated_sig = FvR_annotated[(FvR_annotated[,7]<.1),]
  FvR_annotated_sig = FvR_annotated_sig[rowSums(is.na(FvR_annotated_sig)) != ncol(FvR_annotated_sig), ]
  ##write it out##
  write.csv(FvR_annotated_sig, "FvR_Sig_Genes.csv", row.names = FALSE)
  
  ##read in data##
  RvG = read.csv("DESeq_RBCvGBC_Full_StringParam_new.csv")
  RvG_annotated = merge(RvG, annos, by.x = "X", by.y = "gene", all.x = TRUE)
  ##select significant genes##
  RvG_annotated_sig = RvG_annotated[(RvG_annotated[,7]<.1),]
  RvG_annotated_sig = RvG_annotated_sig[rowSums(is.na(RvG_annotated_sig)) != ncol(RvG_annotated_sig), ]
  ##write it out##
  write.csv(RvG_annotated_sig, "RvG_Sig_Genes.csv", row.names = FALSE)

  ##read in data##
  Fib = read.csv("DESeq_Fibrosis_Full_StringParam_new.csv")
  Fib_annotated = merge(Fib, annos, by.x = "X", by.y = "gene", all.x = TRUE)
  ##select significant genes##
  Fib_annotated_sig = Fib_annotated[(Fib_annotated[,7]<.1),]
  Fib_annotated_sig = Fib_annotated_sig[rowSums(is.na(Fib_annotated_sig)) != ncol(Fib_annotated_sig), ]
  ##write it out##
  write.csv(Fib_annotated_sig, "Fibrosis_Sig_Genes.csv", row.names = FALSE)
  
  ##read in data##
  FGI = read.csv("DESeq_FvGbyInfect_Full_StringParam_new.csv")
  FGI_annotated = merge(FGI, annos, by.x = "X", by.y = "gene", all.x = TRUE)
  ##select significant genes##
  FGI_annotated_sig = FGI_annotated[(FGI_annotated[,7]<.1),]
  FGI_annotated_sig = FGI_annotated_sig[rowSums(is.na(FGI_annotated_sig)) != ncol(FGI_annotated_sig), ]
  ##write it out##
  write.csv(FGI_annotated_sig, "F2vGBCxInfect_Sig_Genes.csv", row.names = FALSE)
  
  ##read in data##
  FRI = read.csv("DESeq_FvRbyInfect_Full_StringParam_new.csv")
  FRI_annotated = merge(FRI, annos, by.x = "X", by.y = "gene", all.x = TRUE)
  ##select significant genes##
  FRI_annotated_sig = FRI_annotated[(FRI_annotated[,7]<.1),]
  FRI_annotated_sig = FRI_annotated_sig[rowSums(is.na(FRI_annotated_sig)) != ncol(FRI_annotated_sig), ]
  ##write it out##
  write.csv(FRI_annotated_sig, "F2vRBCxInfect_Sig_Genes.csv", row.names = FALSE)
  
  ##read in data##
  RGI = read.csv("DESeq_GvRbyInfect_Full_StringParam_new.csv")
  RGI_annotated = merge(RGI, annos, by.x = "X", by.y = "gene", all.x = TRUE)
  ##select significant genes##
  RGI_annotated_sig = RGI_annotated[(RGI_annotated[,7]<.1),]
  RGI_annotated_sig = RGI_annotated_sig[rowSums(is.na(RGI_annotated_sig)) != ncol(RGI_annotated_sig), ]
  ##write it out##
  write.csv(RGI_annotated_sig, "RBCvGBCxInfect_Sig_Genes.csv", row.names = FALSE)
  
  ##read in data##
  FRF = read.csv("DESeq_CrossbyFibrosis_Full_StringParam_new.csv")
  FRF_annotated = merge(FRF, annos, by.x = "X", by.y = "gene", all.x = TRUE)
  ##select significant genes##
  FRF_annotated_sig = FRF_annotated[(FGI_annotated[,7]<.1),]
  FRF_annotated_sig = FRF_annotated_sig[rowSums(is.na(FRF_annotated_sig)) != ncol(FRF_annotated_sig), ]
  ##write it out##
  write.csv(FRF_annotated_sig, "CrossXFibrosis_Sig_Genes.csv", row.names = FALSE)

  
##one last step to make the IPA input- we have to make a single file with all the data##
  ##use infection as a base and remove what you want##
  names(infect_annotated)
  infect_IPA = infect_annotated[c(1,8,3,6,7)]
  ##rename the columns##
  names(infect_IPA) = c("Gene", "spID", "Infect_LFC", "Infect_pval", "Infect_padj")
  ##add in the rest of the data##
  names(FvG_annotated)
  FvG_IPA = FvG_annotated[c(1,3,6,7)]
  names(FvG_IPA) = c("Gene", "FvG_LFC", "FvG_pval", "FvG_padj")
  FvR_IPA = FvR_annotated[c(1,3,6,7)]
  names(FvR_IPA) = c("Gene", "FvR_LFC", "FvR_pval", "FvR_padj")
  RvG_IPA = RvG_annotated[c(1,3,6,7)]
  names(RvG_IPA) = c("Gene", "RvG_LFC", "RvG_pval", "RvG_padj")
  Fib_IPA = Fib_annotated[c(1,3,6,7)]
  names(Fib_IPA) = c("Gene", "Fib_LFC", "Fib_pval", "Fib_padj")
  ##merge them all together##
  step1 = merge(infect_IPA, FvG_IPA, by = "Gene")
  step2 = merge(step1, FvR_IPA, by = "Gene")  
  step3 = merge(step2, RvG_IPA, by = "Gene")  
  final = merge(step3, Fib_IPA, by = "Gene")  

  ##and write it out##
  write.csv(final, "TagSeqMS_IPA.csv", row.names = FALSE)
  
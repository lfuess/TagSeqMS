##This script helps us sort things out for our Venn Diagram and other comparasion stats##
##It looks for overlap in signficant DEGs across different factors from our original model##

library(tidyr)

##infection genes##
infect=read.csv("Infection_Sig_Genes.csv")
infect = infect %>% drop_na(X)
infect = infect[,c(1,3,7,8,9)]
names(infect) <- c("Gene","Infect_LFC","Infect_padj","spID","Protein")

##fibrosis genes##
fibrosis=read.csv("Fibrosis_Sig_Genes.csv")
fibrosis = fibrosis %>% drop_na(X)
fibrosis = fibrosis[,c(1,3,7,8,9)]
names(fibrosis) <- c("Gene","Fibrosis_LFC","Fibrosis_padj","spID","Protein")

##infect vs. fibrosis##
infectfib=merge(infect,fibrosis, by = "Gene", all.x = TRUE, all.y = TRUE)
write.csv(infectfib,"infectvfib.csv", row.names = FALSE)

##FvR genes##
FvR=read.csv("FvR_Sig_Genes.csv")
FvR = FvR %>% drop_na(X)
FvR = FvR[,c(1,3,7)]
names(FvR) <- c("Gene","FvR_LFC","FvR_padj")

##FvG genes##
FvG=read.csv("FvG_Sig_Genes.csv")
FvG = FvG %>% drop_na(X)
FvG = FvG[,c(1,3,7)]
names(FvG) <- c("Gene","FvG_LFC","FvG_padj")

##RvG genes##
RvG=read.csv("RvG_Sig_Genes.csv")
RvG = RvG %>% drop_na(X)
RvG = RvG[,c(1,3,7)]
names(RvG) <- c("Gene","RvG_LFC","RvG_padj")

##infect vs. crosses##
infectcross1=merge(infect,FvG, by = "Gene", all.x = TRUE)
infectcross2=merge(infectcross1,FvR, by = "Gene", all.x = TRUE)
infectcross3=merge(infectcross2,RvG, by = "Gene", all.x = TRUE)
write.csv(infectcross3,"infectvcross.csv", row.names = FALSE)

##all crosses##
cross1=merge(FvR,FvG, by = "Gene", all.x = TRUE, all.y = TRUE)
cross2=merge(cross1,RvG, by = "Gene", all.x = TRUE, all.y = TRUE)
write.csv(cross2,"crosscomp.csv", row.names = FALSE)


##create something for figuring out venn diagram comparing many factors##
  ##create a list of DEGs for cross type##
  FvG_genes = FvG[,c(1:2)]
  names(FvG_genes) <- c("Gene","LFC")
  FvR_genes = FvR[,c(1:2)]
  names(FvR_genes) <- c("Gene","LFC")
  RvG_genes = RvG[,c(1:2)]
  names(RvG_genes) <- c("Gene","LFC")
  ##combine
  all_cross = rbind(FvG_genes, FvR_genes, RvG_genes)
  ##add column##
  all_cross['Cross_Sig']='Y'
  ##remove suprious columns##
  all_cross=all_cross[c(1,3)]
  ##remove duplicate rows  
  all_cross=unique(all_cross)

  ##do the same for interaction DEGs##
  FGI=read.csv("F2vGBCxInfect_Sig_Genes.csv")
  FGI=FGI[,c(1,2)]  
  names(FGI) <- c("X","baseMean_FGI")
  FRI=read.csv("F2vRBCxInfect_Sig_Genes.csv")
  FRI=FRI[,c(1,2)] 
  names(FRI) <- c("X","baseMean_FRI")
  RGI=read.csv("RBCvGBCxInfect_Sig_Genes.csv")
  RGI=RGI[,c(1,2)]  
  names(RGI) <- c("X","baseMean_RGI")
  ##make something for a standalone venn
  cross1=merge(FGI,RGI, by = "X", all.x = TRUE, all.y = TRUE)
  cross2=merge(cross1,FRI, by = "X", all.x = TRUE, all.y = TRUE)
  write.csv(cross2,"crossinfectcomp.csv", row.names = FALSE)
  ##combine
  all_interact = rbind(FGI, FRI, RGI)
  ##add column##
  all_interact['Interact_Sig']='Y'
  ##remove suprious columns##
  all_interact=all_interact[c(1,3)]
  ##remove duplicate rows  
  all_interact=unique(all_interact)
  names(all_interact) = c("Gene", "Interact_Sig")

  ##do the same for infection##
  infect_genes = infect[,c(1,2)]
  ##add column##
  infect_genes['Infect_Sig']='Y'
  ##remove suprious columns##
  infect_genes=infect_genes[c(1,3)]
  
  ##do the same for fibrosis##
  fib_genes = fibrosis[,c(1,2)]
  ##add column##
  fib_genes['Fibrosis_Sig']='Y'
  ##remove suprious columns##
  fib_genes=fib_genes[c(1,3)]

  ##and then we merge them all together## 
  merge = merge(all_cross, all_interact, by = "Gene", all.x = TRUE, all.y = TRUE)
  merge1 = merge(merge,infect_genes, by = "Gene", all.x = TRUE, all.y = TRUE)  
  final = merge(merge1, fib_genes, by = "Gene", all.x = TRUE, all.y = TRUE)  

  ##write it out##
  write.csv(final,"VennDiagram_Info.csv", row.names = FALSE)
  
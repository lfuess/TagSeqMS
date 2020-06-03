#install.packages('VennDiagram')
library(VennDiagram)


draw.triple.venn(area1=240, area2=313, area3=94, n12=231, n13=74, n23=85, n123=74, category = c("F2vGBC", "F2vRBC", "RBCvGBC"), lty = "blank", fill=c("#6dae90", "#f3c558", "#30b4cc"), alpha=rep(.7, 3), cex=rep(1.7, 7), fontfamily=rep("Arial Black", 7), cat.cex=rep(1, 3), cat.fontface=rep("bold.italic", 3), cat.fontfamily=rep("Arial",3))


1: FvG
2: FvR
3: RvG
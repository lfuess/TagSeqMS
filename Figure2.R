#install.packages('VennDiagram')
library(VennDiagram)

draw.triple.venn(area1=7745, area2=10601, area3=1202, n12=7175, n13=399, n23=856, n123=203, category = c("F2vGBC", "F2vRBC", "RBCvGBC"), lty = "blank", fill=c("#6dae90", "#f3c558", "#30b4cc"), alpha=rep(.7, 3), cex=rep(1.7, 7), fontfamily=rep("Arial Black", 7), cat.cex=rep(1, 3), cat.fontface=rep("bold.italic", 3), cat.fontfamily=rep("Arial",3))


1: FvG
2: FvR
3: RvG
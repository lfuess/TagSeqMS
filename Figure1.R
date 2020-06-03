##This code makes Figure 1 from the MS##
##It uses a hand currated file from the IPA output##

library(ggplot2)

infect = read.csv("Infect_bar_all.csv")
ggplot(infect, aes(y=z.score, x=reorder(Pathway, z.score))) + 
  geom_bar(position="dodge", stat="identity", fill = "#cc3d24") + coord_flip()+
  theme_classic() + geom_abline(slope=0, intercept=0,  col = "black")+
  xlab("Pathway") + ylab("z-score")
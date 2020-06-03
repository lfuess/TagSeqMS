##This code makes Figure 4 from the MS##
##It uses a hand currated file from the IPA output##

library(ggplot2)

fibrosis = read.csv("Fib_bar_all.csv")
ggplot(fibrosis, aes(y=z.score, x=reorder(Pathway, z.score))) + 
  geom_bar(position="dodge", stat="identity", fill = "#004f7a" ) + coord_flip()+
  theme_classic() + geom_abline(slope=0, intercept=0,  col = "black")+
  xlab("Pathway") + ylab("z-score")
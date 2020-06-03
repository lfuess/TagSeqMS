##Scatter plot figure##
##uses handcurrated output from IPA##

library(ggplot2)

scat = read.csv("Infect_Fib_path_scat.csv")

ggplot(scat, aes(y=Fibrosis, x=Infection, color=Color)) + 
  geom_point() + geom_abline(slope=1, intercept=0, linetype="dashed") +
  theme_classic() +
  theme(legend.background = element_rect(fill="white",
                                         size=0.5, linetype="solid", 
                                         colour ="black")) +
  scale_color_manual(values=c("#004f7a", "#cc3d24")) +
  xlab("Infection z-score") + ylab("Fibrosis z-score") +
  theme(legend.position="none")

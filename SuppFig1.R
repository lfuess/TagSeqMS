##let's make a regulator comparison bar graph

library(ggplot2)

scat = read.csv("Infect_Fib_reg_scat.csv")

p1 = ggplot(scat, aes(y=Fibrosis, x=Infection, color=Color)) + 
  geom_point() + geom_abline(slope=1, intercept=0, linetype="dashed") +
  theme_classic() +
  theme(legend.background = element_rect(fill="white",
                                         size=0.5, linetype="solid", 
                                         colour ="black")) +
  scale_color_manual(values=c("#004f7a", "#cc3d24")) +
  xlab("Infection z-score") + ylab("Fibrosis z-score") +
  theme(legend.position="none")

##we can also make unique plots for each fibrosis and infection##

infect = read.csv("Infect_reg_bar.csv")
p2 = ggplot(infect, aes(y=z.score, x=reorder(Regulator, z.score))) + 
  geom_bar(position="dodge", stat="identity", fill = "#cc3d24") + coord_flip()+
  theme_classic() + geom_abline(slope=0, intercept=0,  col = "black")+
  xlab("Regulator") + ylab("z-score")


fibrosis = read.csv("Fib_reg_bar.csv")
p3 = ggplot(fibrosis, aes(y=z.score, x=reorder(Regulator, z.score))) + 
  geom_bar(position="dodge", stat="identity", fill = "#004f7a" ) + coord_flip()+
  theme_classic() + geom_abline(slope=0, intercept=0,  col = "black")+
  xlab("Regulator") + ylab("z-score")

##make them into a multiplot##
library("cowplot")
p4 = plot_grid(p2, p3, ncol = 1, align = "v")
ggdraw() +
  draw_plot(p4, x = 0, y = 0, width = .5, height = 1) +
  draw_plot(p1, x = .5, y = .5, width = .5, height = .5) +
  draw_plot_label(label = c("A", "C", "B"), size = 15,
                  x = c(0, 0.5, 0), y = c(1, 1, 0.5))

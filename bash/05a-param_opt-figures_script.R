# Parameter optimization----

## subset by species

## import data

## make sure factor variables are in correct order for plots
parameter_opt$Change <- factor(parameter_opt$Change, 
                               levels=c("2/3", "3/4", "4/5","5/6", "6/7", "7/8",
                                        "8/9", "9/10", "10/11","11/12"))

CD_opt <- subset(parameter_opt, species == "CD")
DS_opt <- subset(parameter_opt, species == "DS")
EC_opt <- subset(parameter_opt, species == "EC")
LS_opt <- subset(parameter_opt, species == "LS")
PA_opt <- subset(parameter_opt, species == "PA")
TO_opt <- subset(parameter_opt, species == "TO")

CD_plot <- ggplot(aes(Change, R60.loci, group=1), data=CD_opt)  +
  geom_line(color="blue", size =0.8) +  
  geom_point(color="black", size =1.5) +
  theme_bw() +
  theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank(),
        axis.title = element_text(size=7),
        axis.text = element_text(size=6),
        plot.title = element_text(size = 7)) +
  ylab("Change in R60 loci") + xlab("Change in Ustacks M &\n Cstacks n") +
  ggtitle("Parameter optimization: \nCynodon dactylon") +
  ylim(-20, 80) + geom_hline(yintercept = 0, color ="grey", linetype="dashed")

DS_plot <- ggplot(aes(Change, R60.loci, group=1), data=DS_opt)  +
  geom_line(color="blue", size =0.8) +  
  geom_point(color="black", size =1.5) +
  theme_bw() +
  theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank(),
        axis.title = element_text(size=7),
        axis.text = element_text(size=6),
        plot.title = element_text(size = 7)) +
  ylab("Change in R60 loci") + xlab("Change in Ustacks M &\n Cstacks n") +
  ggtitle("Parameter optimization: \nDigitaria sanguinalis") +
  ylim(-20, 80) + geom_hline(yintercept = 0, color ="grey", linetype="dashed")

EC_plot <- ggplot(aes(Change, R60.loci, group=1), data=EC_opt)  +
  geom_line(color="blue", size =0.8) +  
  geom_point(color="black", size =1.5) +
  theme_bw() +
  theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank(),
        axis.title = element_text(size=7),
        axis.text = element_text(size=6),
        plot.title = element_text(size = 7)) +
  ylab("Change in R60 loci") + xlab("Change in Ustacks M &\n Cstacks n") +
  ggtitle("Parameter optimization: \nErigeron canadensis") +
  ylim(-20, 80) + geom_hline(yintercept = 0, color ="grey", linetype="dashed")

LS_plot <- ggplot(aes(Change, R60.loci, group=1), data=LS_opt)  +
  geom_line(color="blue", size =0.8) +  
  geom_point(color="black", size =1.5) +
  theme_bw() +
  theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank(),
        axis.title = element_text(size=7),
        axis.text = element_text(size=6),
        plot.title = element_text(size = 7)) +
  ylab("Change in R60 loci") + xlab("Change in Ustacks M &\n Cstacks n") +
  ggtitle("Parameter optimization: \nLactuca serriola") +
  ylim(-20, 80) + geom_hline(yintercept = 0, color ="grey", linetype="dashed")

PA_plot <- ggplot(aes(Change, R60.loci, group=1), data=PA_opt)  +
  geom_line(color="blue", size =0.8) +  
  geom_point(color="black", size =1.5) +
  theme_bw() +
  theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank(),
        axis.title = element_text(size=7),
        axis.text = element_text(size=6),
        plot.title = element_text(size = 7)) +
  ylab("Change in R60 loci") + xlab("Change in Ustacks M &\n Cstacks n") +
  ggtitle("Parameter optimization: \nPoa annua") +
  ylim(-20, 80) + geom_hline(yintercept = 0, color ="grey", linetype="dashed")

TO_plot <- ggplot(aes(Change, R60.loci, group=1), data=TO_opt)  +
  geom_line(color="blue", size =0.8) +  
  geom_point(color="black", size =1.5) +
  theme_bw() +
  theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank(),
        axis.title = element_text(size=7),
        axis.text = element_text(size=6),
        plot.title = element_text(size = 7)) +
  ylab("Change in R60 loci") + xlab("Change in Ustacks M &\n Cstacks n") +
  ggtitle("Parameter optimization: \nTaraxacum officinale") +
  ylim(-20, 80) + geom_hline(yintercept = 0, color ="grey", linetype="dashed")


#plot using Rmisc multiplot function

multiplot(CD_plot, DS_plot, EC_plot, LS_plot, PA_plot, TO_plot, cols = 2)
ppi <- 300
pdf("plot.pdf", width=9, height=8)
multiplot(CD_plot, DS_plot, EC_plot, LS_plot, PA_plot, TO_plot, cols = 3)
dev.off()

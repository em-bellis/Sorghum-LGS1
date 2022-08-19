#updated 8/19 to include gene models
library(ggplot2)

td.lgs <- read.table('/Users/emilywork/Desktop/LGS1/Figs_12.3.18/data/LGS1.1Mb.miss.0.07.bi.maf0.01.recode.Tajima.D', header=T)
td.p14 <- read.table('/Users/emilywork/Desktop/LGS1/Figs_12.3.18/data/P14.1Mb.miss.0.07.bi.maf0.01.recode.Tajima.D', header=T)
td.smax1 <- read.table('/Users/emilywork/Desktop/LGS1/LGS1_revision_8.19/Chr7.20Mb.miss0.07.bi.maf0.01.TajimaD', header=T)

td.l<-ggplot(td.lgs, aes(x=BIN_START/1000000, y=TajimaD)) +geom_point(size=1)  + geom_line(lwd=0.2)+geom_vline(xintercept=69986146/1000000, col="grey") + theme_classic() + xlab("") + ylab("Tajima's D")+ xlim(c((69986146-100000)/1000000, (69986146+100000)/1000000)) + ylim(c(-2.5,6))+ geom_hline(yintercept=ml, lty=1, col="black", lwd=0.5)#+ annotate("text", x=69.775, y=5, label="LGS1", fontface="italic", size=5)+ annotate("text", x=(69986146+100000)/1000000,y=5, label="S5_69986146", col="grey")
td.p<-ggplot(td.p14, aes(x=BIN_START/1000000, y=TajimaD)) + geom_point(size=1) + geom_line(lwd=0.2)+ ylab("")+ geom_vline(xintercept=21521798/1000000, col="grey") + xlim(c((21521798-100000)/1000000, (21521798+100000)/1000000))+ theme_classic()+ xlab("")+ ylim(c(-2.5,6))+ geom_hline(yintercept=0.71, lty=1, lwd=0.5)+ annotate("text", x=21.37, y=5, label="pectinesterase",  size=5)#+ annotate("text", x=(21521798+100000)/1000000,y=5, label="S2_21521798", col="grey")
td.s<-ggplot(td.smax1, aes(x=BIN_START/1000000, y=TajimaD)) + geom_point(size=1) + geom_line(lwd=0.2)+ ylab("")+ geom_vline(xintercept=14459084/1000000, col="grey") + xlim(c((14458102-100000)/1000000, (14459084+100000)/1000000))+ theme_classic()+ xlab("Position (Mb)")+ ylim(c(-2.5,6)) + geom_hline(yintercept=0.71, lty=1, lwd=0.5)#+ annotate("text", x=14.25, y=5, label="SMAX1", fontface="italic", size=5)+ annotate("text", x=(14458102+100000)/1000000,y=5, label="S7_14459084", col="grey")

pdf('/Users/emilybellis/Desktop/LGS1_figs/Fig5.pdf', width=3.41, height=4)
multiplot(td.p, td.l, td.s,cols=1)
dev.off()

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }

 if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}
####TD confidence intervals
td.ran1k <- read.table('/Users/emilywork/Desktop/LGS1/LGS1_revision_8.19/randomTD.txt')
ll <- quantile(td.ran1k$V1, c(0.025,0.5, 0.975))[1]
ml <- quantile(td.ran1k$V1, c(0.025,0.5, 0.975))[2]
ul <- quantile(td.ran1k$V1, c(0.025,0.5, 0.975))[3]


#####include gene models in fig?
gm.lgs <- read.table('/Users/emilywork/Desktop/LGS1/LGS1_revision_8.19/lgs1.genemodels.txt', header=F)

td.l2 <- td.l + geom_segment(aes(x=gm.lgs$V4/1000000, xend=gm.lgs$V5/1000000, y=-2.5, yend=-2.5), data=gm.lgs, size=3)+annotate("rect", xmin=-Inf, xmax=Inf, ymin=ll, ymax=ul, alpha=0.2)#+ geom_hline(yintercept=ul, lty=2, col="black", lwd=0.5)+ geom_hline(yintercept=ll, lty=2, col="black", lwd=0.5)

gm.smx <- read.table('/Users/emilywork/Desktop/LGS1/LGS1_revision_8.19/smx.genemodels.txt', header=F)
td.s2 <- td.s + geom_segment(aes(x=gm.smx$V4/1000000, xend=gm.smx$V5/1000000, y=-2.5, yend=-2.5), data=gm.smx, size=3)+annotate("rect", xmin=-Inf, xmax=Inf, ymin=ll, ymax=ul, alpha=0.2)#+ geom_hline(yintercept=ul, lty=2, col="black", lwd=0.5)+ geom_hline(yintercept=ll, lty=2, col="black", lwd=0.5)

gm.pec <- read.table('/Users/emilywork/Desktop/LGS1/LGS1_revision_8.19/pec.genemodels.txt', header=F)
td.p2 <- td.p + geom_segment(aes(x=gm.pec$V4/1000000, xend=gm.pec$V5/1000000, y=-2.5, yend=-2.5), data=gm.pec, size=3)+annotate("rect", xmin=-Inf, xmax=Inf, ymin=ll, ymax=ul, alpha=0.2)

pdf('/Users/emilywork/Desktop/LGS1/LGS1_revision_8.19/Fig4.pdf', width=3.41, height=4)
multiplot(td.p2, td.l2, td.s2,cols=1)
dev.off()

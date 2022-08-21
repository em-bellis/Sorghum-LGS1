library('tidyr')
library('ggplot2')
library('dplyr')

gwas_gemma <- read.table('../Fig3/global_200_0.01.assoc.txt', header=T) #kinship matrix maf 0.01; snps 0.01; this file was made for global set, setting prob of occurrence outside of 200 km to 0; 10.30.18

##combine gemma output with allele counts
counts <- read.table('../Fig3/lines.geo.frq.count', header=T, row.names=NULL)
counts$rs <- paste0("S",counts$row.names,"_",counts$CHROM)
comb <- inner_join(gwas_gemma, counts)
comb <- comb %>% tidyr::separate(ALLELE.COUNT, c("Allele","Count"),sep=":") %>% filter(Count >= 10)
comb$Count <- as.numeric(comb$Count)
comb <- subset(comb, Count > 12)
gwas <- comb %>% tidyr::separate(rs,c("Chr","Pos"),sep="_")
colnames(gwas)[1] <- "Chromosome"
colnames(gwas)[3] <- "Position"
gwas <- gwas %>% dplyr::select(Chromosome, Position, p_score)
gwas$Position<-as.numeric(gwas$Position)
gwas$Chromosome<-factor(gwas$Chromosome, levels=c(1:10))

##prep data for manhattan plot
don <- gwas %>% 
  
  # Compute chromosome size
  group_by(Chromosome) %>% 
  summarise(chr_len=max(Position)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  dplyr::select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(gwas, ., by=c("Chromosome"="Chromosome")) %>%
  
  # Add a cumulative position of each SNP
  arrange(Chromosome, Position) %>%
  mutate( BPcum=Position+tot)

colnames(don)[3] <- "P.value"

#color sorghum qtl from OZ atlas
qtl <- read.csv('../Fig3/OZ_Sorghum.csv', head=T)
qtl <- qtl %>% tidyr::separate(LG.Start.End..v3.0., c("Chromosome","Position"), sep=":") %>% separate(Position, c("Start","Stop"), sep="-")
qtl$Start <- as.numeric(qtl$Start)
qtl$Stop <- as.numeric(qtl$Stop)

is9830 <- subset(qtl, qtl$Population =="IS9830/E36-1") #low germination
n13 <- subset(qtl, qtl$Population =='N13/E36-1')		#mechanical resistance

don.n13 <- NULL
for (i in 1:nrow(n13)){
	qtl.curr <- subset(don, Chromosome==n13$Chromosome[i] & Position >= n13$Start[i] & Position <= n13$Stop[i] )
	don.n13 <- rbind.data.frame(don.n13, qtl.curr)
}

don.is9830 <- NULL
for (i in 1:nrow(is9830)){
	qtl.curr <- subset(don, Chromosome==is9830$Chromosome[i] & Position >= is9830$Start[i] & Position <= is9830$Stop[i] )
	don.is9830 <- rbind.data.frame(don.is9830, qtl.curr)
}

#color SL biosynthesis gene models
don.sl <- NULL
sl <- read.csv('../Fig3/SL_genes.csv', header=T)
sl <- subset(sl, Chromosome != "NA")
sl <- sl %>%tidyr::separate(Start, c("StartPos","NA"),"-")
sl$StartPos <- as.numeric(gsub(",","",sl$StartPos))
sl$Stop <- as.numeric(gsub(",","",sl$Stop))
colnames(sl)[5] <- "Start"
for (i in 1:nrow(sl)){
	qtl.curr <- subset(don, Chromosome==sl$Chr[i] & Position >= sl$Start[i] & Position <= sl$Stop[i] )
	don.sl <- rbind.data.frame(don.sl, qtl.curr)
}

axisdf = don %>% group_by(Chromosome) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

##determine p-value cut-off
#gemma p.adjust
gwas_gemma$adj <- p.adjust(gwas_gemma$p_score, method="fdr", n=length(gwas_gemma$p_score))
tmp <- subset(gwas_gemma, adj <0.05) #note that these include SNPs nearby to each other removed in the final table
tmp[order(tmp$adj),]
thresh <- -log10(max(tmp$p_score))

p <- ggplot(don, aes(x=BPcum, y=-log10(P.value), col=Chromosome)) + geom_point(alpha=0.4, size=1) + theme_classic(base_size=20) + xlab("Chromosome") + guides(col=FALSE) + ylab("-log10(P)")+scale_x_continuous(label = axisdf$Chromosome, breaks=axisdf$center)+scale_colour_manual(values=rep(c("grey35","grey55"),5))+ geom_hline(yintercept=thresh, lty=2) + ylab(expression(paste("-log"["10"],"(",italic("p"),")")))

##add qtl
qtl.cols <- c('#56B4E9','#E69F00','#D55E00')
pdf('/Users/emilywork/Desktop/LGS1/Figs_12.3.18/Fig2.gwas.pdf', width=7, height=2)
p + geom_point(data=don.n13,col=qtl.cols[1], alpha=0.7, size=1) + geom_point(data=don.is9830, col=qtl.cols[2], alpha=0.7, size=1) + geom_point(data=don.sl, col=qtl.cols[3], size=1) + ylim(c(0,9)) + annotate("rect", xmin=525000000, xmax=545000000, ymin = 6.3, ymax=6.8, fill=qtl.cols[1])+annotate("rect", xmin=525000000, xmax=545000000, ymin = 7.0, ymax=7.5, fill=qtl.cols[2]) + annotate("rect", xmin=525000000, xmax=545000000, ymin = 7.7, ymax=8.2, fill=qtl.cols[3])
#+ annotate("text", x=620000000, y=8, label="SL biosynthesis", size=5) + annotate("text", x=620000000, y=7.3, label="QTL(low germ.)", size=5) + annotate("text", x=629000000, y=6.6, label="QTL(mechanical)", size=5)
dev.off()

##loci over threshold in each qtl
subset(don.n13, -log10(P.value) > thresh)
subset(don.is9830, -log10(P.value) > thresh)
subset(don.sl, -log10(P.value) > thresh)

##in case a lower res version is needed...
install.packages('devtools')
devtools::install_github('VPetukhov/ggrastr')
library(ggrastr)
library(Cairo)
p <- ggplot(don, aes(x=BPcum, y=-log10(P.value), col=Chromosome)) + geom_point_rast(alpha=0.4, size=1) + theme_classic(base_size=20) + xlab("Chromosome") + guides(col=FALSE) + ylab("-log10(P)")+scale_x_continuous(label = axisdf$Chromosome, breaks=axisdf$center)+scale_colour_manual(values=rep(c("grey35","grey55"),5))+ geom_hline(yintercept=thresh, lty=2) + ylab(expression(paste("-log"["10"],"(",italic("p"),")")))
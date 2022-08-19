##visualize 'niche envelopes' for sorghum-only and all occurrence models
##2.13.19

library(dismo)
library(raster)
library(gdalUtils)
library(rgdal)
library(ggplot2)

all <- raster("/Users/emilybellis/Desktop/LGS1_revision_2019/enm_figs/preds.all.tif")
sorg <- raster("/Users/emilybellis/Desktop/LGS1_revision_2019/enm_figs/preds.sorg.tif")
occ <- read.csv('/Users/emilybellis/Desktop/striga_specificity/sorg_mill_map/curated_occ_ESB_7.26.18.csv', head=T)

##what is the distribution of hss for actual occurrences?
occ$hss <- extract(all, cbind(occ$lon, occ$lat))
occ$hss_sorg <- extract(sorg, cbind(occ$lon, occ$lat))
my.lab <- expression(paste)

pdf('/Users/emilybellis/Desktop/LGS1_revision_2019/enm_figs/hss_dist.pdf', width = 2, height=2)
ggplot(occ, aes(x=hss)) + geom_histogram(bins=20, fill="grey90", col="black") + theme_classic() + xlab("Parasite HS") + ylab("Parasite\noccurrences")
dev.off()

#similar plot for sorghum
striga_score <- read.table('/Users/emilybellis/Desktop/LGS1_revision_2019/enm_figs/global_200.fam', header=F)

pdf('/Users/emilybellis/Desktop/LGS1_revision_2019/enm_figs/hss_sorg_dist.pdf', width = 2, height=2)
ggplot(striga_score, aes(x=V6)) + geom_histogram(bins=20, fill="grey90", col="black") + theme_classic() + xlab("Parasite HS") + ylab("Sorghum\naccessions")
dev.off()


nrow(occ)
quantile(occ$hss, na.rm=T, p=0.1)
##1369
ggplot(occ, aes(x=hss_sorg)) + geom_histogram(bins=20)

##visualize core, suitable, and overall SSA distribution
#create masks 
core <- calc(all, fun=function(x){ x[x < 0.5] <- NA; return(x)} )
core.25 <- calc(all, fun=function(x){ x[x < 0.1] <- NA; return(x)} )

score <- calc(sorg, fun=function(x){ x[x < 0.5] <- NA; return(x)} )
score.25 <- calc(sorg, fun=function(x){ x[x < 0.1] <- NA; return(x)} )

#1-(cellStats(core, stat="countNA")/86369001) #gets you the proportion of SSA that is core or suitable

#get data for each environment, mask based on HSS
env <- raster("/Users/emilybellis/Desktop/LGS1_revision_2019/enm_figs/env12.tif")
env.core <- mask(env, core)
#cellStats(env.core, stat="mean")
env.core25 <- mask(env, core.25)
env.all <- mask(env, all)

quantile(env.core, p=c(0.1,0.9,0.5))
quantile(env.all, p=c(0.1,0.9,0.5))


#change to dataframe
env.core.df <- as.data.frame(env.core, na.rm =T)
colnames(env.core.df) <- "env"
env.core25.df <- as.data.frame(env.core25, na.rm =T)
colnames(env.core25.df) <- "env"
env.all.df <- as.data.frame(env.all, na.rm =T)
colnames(env.all.df) <- "env"

ks.test(env.all.df$env,env.core25.df$env)
ks.test(env.all.df$env, env.core.df$env)

pdf("/Users/emilybellis/Desktop/LGS1_revision_2019/enm_figs/pet.pdf", width = 3, height =3, pointsize=16)
ggplot() +geom_histogram(data=env.all.df,aes(x=env), fill="grey50",col="black", size=0.1)+ geom_histogram(data=env.core25.df,aes(x=env), fill="#EDB48EFF",col="black", size=0.5)+ geom_histogram(data=env.core.df,aes(x=env), fill="#00A600FF", col="black", size=0.5)+ theme_bw() + xlab("Phosphorus (mg/kg)") + ylab("Number grid cells") +xlim(c(0,3000))+geom_vline(xintercept=quantile(env.core, p=0.5), lty=2)
dev.off()

############plot for sorghum model
env.score <- mask(env, score)
#cellStats(env12.score, stat="mean")
env.score25 <- mask(env, score.25)
env.s <- mask(env, sorg)

env.core.df <- as.data.frame(env.score, na.rm =T)
colnames(env.core.df) <- "env"
env.core25.df <- as.data.frame(env.score25, na.rm =T)
colnames(env.core25.df) <- "env"
env.all.df <- as.data.frame(env.s, na.rm =T)
colnames(env.all.df) <- "env"

pdf("/Users/emilybellis/Desktop/LGS1_revision_2019/enm_figs/cly.sorg.pdf", width = 3, height =3, pointsize=16)
ggplot() +geom_histogram(data=env.all.df,aes(x=env), fill="grey50",col="black", size=0.1)+ geom_histogram(data=env.core25.df,aes(x=env), fill="#EDB48EFF",col="black", size=0.5)+ geom_histogram(data=env.core.df,aes(x=env), fill="#00A600FF", col="black", size=0.5)+ theme_bw() + xlab("Annual precipitation (mm/yr)") + ylab("Number grid cells")+geom_vline(xintercept=quantile(env.score, p=0.5), lty=2)
dev.off()

quantile(env.score, p=c(0.1,0.9,0.5))

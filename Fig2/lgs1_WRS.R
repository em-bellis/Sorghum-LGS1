###updated 8.7.19
library(raster)
library(rgdal)
library(ggplot2)
library(maptools)
library(maps)
data(wrld_simpl)
library(gridExtra)
library(plyr)
library(dismo)
library(viridis)
library(RColorBrewer)
data(wrld_simpl)
library(rgeos)

sorg <- raster('/Users/emilywork/Desktop/ENM/distModel7.30/preds.sorg.tif', header=T)
all<- raster('/Users/emilywork/Desktop/ENM/distModel7.30/preds.all.tif')

del9067 <- read.table('/Users/emilywork/Desktop/LGS1/MSfigs/createDistmaps/9506.coords.v2.txt', header=F)
del9067$V4 <- "lgs1-3"
del9519 <- read.table('/Users/emilywork/Desktop/LGS1/MSfigs/createDistmaps/9519.coords.v2.txt', header=F)
del9519$V4 <- "lgs1-2"
fs.a <- read.table('/Users/emilywork/Desktop/LGS1/MSfigs/createDistmaps/FS.pres.coords.txt', header=F)
fs.a$V4 <- "Frameshift"
fs.t <- read.table('/Users/emilywork/Desktop/LGS1/MSfigs/createDistmaps/FS.abs.coords.txt', header=F)
fs.t$V4 <- "Intact"

lgs <- rbind.data.frame(del9067, del9519, fs.a, fs.t)
lgs <- subset(lgs, V2 != "NA")
lgs$V4 <- as.factor(lgs$V4)

##extract probability of striga occurrence at location of each landrace
lgs.points <- cbind.data.frame(lgs$V3, lgs$V2)
colnames(lgs.points) <- c("lon","lat")
coordinates(lgs.points) <- ~lon + lat
proj4string(lgs.points)=proj4string(all)
lgs$all.prob <- extract(all, lgs.points)

#apply 200km mask
shgeo <- read.csv('/Users/emilywork/Desktop/ENM/Final_occurrence_data/curated_occ_ESB_7.26.18.csv', header=T,stringsAsFactors=F) ##1116 occurrences
coordinates(shgeo) <- ~lon+lat
projection(shgeo) <- CRS('+proj=longlat +datum=WGS84')
projection(lgs.points) <- CRS('+proj=longlat +datum=WGS84')
x <- circles(shgeo, d=200000, lonlat=T) #sample within radius of 200km
pol <- polygons(x)
ovr <- over(lgs.points, pol)
i <- which(is.na(ovr))

lgs.df <- lgs
lgs.df$within200 <- "N"
for (j in 1:nrow(lgs.df)) {
	lgs.df$within200[j] <- j
}

within200 <- lgs.df[ ! lgs.df$within200 %in% i, ] #619
notwithin200 <- lgs.df[lgs.df$within200 %in% i, ] #356
notwithin200$all.prob<- 0
within200 <- as.data.frame(within200)
notwithin200 <- as.data.frame(notwithin200)
replaced <- rbind.data.frame(within200, notwithin200)
replaced$V4 <- as.factor(replaced$V4)

l3.wrs <- subset(replaced, V4 == "lgs1-3" | V4=="Intact")
l2.wrs <- subset(replaced, V4 == "lgs1-2" | V4=="Intact")
fs.wrs <- subset(replaced, V4 == "Frameshift" | V4=="Intact")

wilcox.test(l3.wrs$all.prob~l3.wrs$V4)
wilcox.test(l2.wrs$all.prob~l2.wrs$V4)
wilcox.test(fs.wrs$all.prob~fs.wrs$V4)

######### what about if subset to within Africa?
map(database="world", xlim=c(-20,130),ylim=c(-50,60))
points(replaced$V3, replaced$V2, col="blue", cex=0.5, pch=19) #975 total
afr_200 <- subset(replaced, V2 < 30 & V3 < 47)
afr_200 <- subset(afr_200, V3 < 43 | V2 < 12)
points(afr_200$V3, afr_200$V2, col="red", cex=0.5, pch=19)

l3.wrs <- subset(afr_200, V4 == "lgs1-3" | V4=="Intact")
l2.wrs <- subset(afr_200, V4 == "lgs1-2" | V4=="Intact")
fs.wrs <- subset(afr_200, V4 == "Frameshift" | V4=="Intact")

wilcox.test(l3.wrs$all.prob~l3.wrs$V4)
wilcox.test(l2.wrs$all.prob~l2.wrs$V4)
wilcox.test(fs.wrs$all.prob~fs.wrs$V4)

######### what about if subset to within west Africa?
wafr_200 <- subset(afr_200, V3 < 19 & V2 > 0)
points(wafr_200$V3, wafr_200$V2, col="turquoise", cex=0.5, pch=19)

l3.wrs <- subset(wafr_200, V4 == "lgs1-3" | V4=="Intact")
l2.wrs <- subset(wafr_200, V4 == "lgs1-2" | V4=="Intact")
fs.wrs <- subset(wafr_200, V4 == "Frameshift" | V4=="Intact")

wilcox.test(l3.wrs$all.prob~l3.wrs$V4)
wilcox.test(l2.wrs$all.prob~l2.wrs$V4)
wilcox.test(fs.wrs$all.prob~fs.wrs$V4)

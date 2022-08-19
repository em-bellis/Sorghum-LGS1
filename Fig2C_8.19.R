###distribution model fig
library(raster)
library(rgdal)
library(ggplot2)
library(maps)
data(worldMapEnv)
library(gridExtra)
library(dismo)

all<- raster('/Users/emilywork/Desktop/ENM/distModel7.30/preds.all.tif')

del9067 <- read.table('/Users/emilywork/Desktop/LGS1/MSfigs/createDistmaps/9506.coords.v2.txt', header=F)
del9067$V4 <- "lgs1-3"
del9519 <- read.table('/Users/emilywork/Desktop/LGS1/MSfigs/createDistmaps/9519.coords.v2.txt', header=F)
del9519$V4 <- "lgs1-2"

#position 69985710
fs.a <- read.table('/Users/emilywork/Desktop/LGS1/MSfigs/createDistmaps/FS.pres.coords.txt', header=F)
fs.a$V4 <- "Frameshift"
fs.t <- read.table('/Users/emilywork/Desktop/LGS1/MSfigs/createDistmaps/FS.abs.coords.txt', header=F)
fs.t$V4 <- "Intact"

lgs <- rbind.data.frame(del9067, del9519)
lgs <- subset(lgs, V2 != "NA")
lgs$V4 <- as.factor(lgs$V4)
fs <- rbind.data.frame(fs.a, fs.t)
fs <- subset(fs, V2 != "NA")

#map
pdf('/Users/emilywork/Desktop/LGS1/LGS1_revision_8.19/Fig2C.pdf')
map(database="world", xlim=c(-20,90),ylim=c(-35,40))
box()
points(fs.t$V3, fs.t$V2, col=rgb(190, 190, 190, max = 255, alpha = 140), pch=19, cex=0.5)
points(fs.a$V3, fs.a$V2, col=rgb(117,112,179, max=255, alpha=140), pch=19, cex=0.75) #purple; fs
points(del9067$V3, del9067$V2,col=rgb(217,95,2, max=255, alpha=140), pch=19, cex=0.75) #red; lgs1-3
points(del9519$V3,del9519$V2,col=rgb(27,158,119, max=255, alpha=140), pch=19, cex=0.75) #teal; lgs1-2
scalebar(2000, xy = c(-15,-30), type = "bar", divs = 2, lonlat = TRUE)
dev.off()

##extract probability of striga occurrence at location of each landrace
lgs.points <- cbind.data.frame(lgs$V3, lgs$V2)
fs.points <- cbind.data.frame(fs$V3, fs$V2)
colnames(fs.points) <- c("lon","lat")
colnames(lgs.points) <- c("lon","lat")
lgs$all.prob <- raster::extract(all, lgs.points)
fs$all.prob <- raster::extract(all, fs.points)
lgs$all.prob[is.na(lgs$all.prob)] <- 0
fs$all.prob[is.na(fs$all.prob)] <- 0

#apply 200km mask
shgeo <- read.csv('/Users/emilywork/Desktop/ENM/Final_occurrence_data/curated_occ_ESB_7.26.18.csv', header=T,stringsAsFactors=F) ##1116 occurrences
coordinates(shgeo) <- ~lon+lat
coordinates(lgs.points) <- ~lon+lat
coordinates(fs.points) <- ~lon+lat
projection(shgeo) <- CRS('+proj=longlat +datum=WGS84')
projection(fs.points) <- CRS('+proj=longlat +datum=WGS84')
projection(lgs.points) <- CRS('+proj=longlat +datum=WGS84')
x <- circles(shgeo, d=200000, lonlat=T) #sample within radius of 200km
pol <- polygons(x)
ovr <- over(fs.points, pol)
ovr.l <- over(lgs.points, pol)
i <- which(is.na(ovr))
i.l <- which(is.na(ovr.l))

fs$within200 <- "N"
lgs$within200 <- "N"
for (j in 1:nrow(lgs)) {
	lgs$within200[j] <- j
}

for (j in 1:nrow(fs)) {
	fs$within200[j] <- j
}

within200 <- fs[ ! fs$within200 %in% i, ] #1246
notwithin200 <- fs[fs$within200 %in% i, ] #363
notwithin200$all.prob<- 0
within200 <- as.data.frame(within200)
notwithin200 <- as.data.frame(notwithin200)
fs.replaced <- rbind.data.frame(within200, notwithin200)

within200 <- lgs[ ! lgs$within200 %in% i.l, ] #1246
notwithin200 <- lgs[lgs$within200 %in% i.l, ] #363
notwithin200$all.prob<- 0
within200 <- as.data.frame(within200)
notwithin200 <- as.data.frame(notwithin200)
lgs.replaced <- rbind.data.frame(within200, notwithin200)

#####probability density functions
s3.points <- rbind(fs.replaced, lgs.replaced)
s3.points <- subset(s3.points, V2 != "NA")
coordinates(s3.points) <- ~V3 + V2
proj4string(s3.points)=proj4string(all)
s3.points <- as.data.frame(s3.points)
p.a <- ggplot(s3.points, aes(x=all.prob, lty=V4, group=V4, col=V4)) + geom_density(size=1.25) + theme_classic(base_size=16) + scale_color_manual(values=c("#7570b3", rgb(190,190,190, max=255, alpha=220), "#1b9e77", "#d95f02"), name="")+ scale_linetype_manual(values=c(1,2,1,1), name="") +xlab("Parasite occurrence") + ylab("Density")+theme(legend.position=c(0.8,0.8))

pdf('/Users/emilywork/Desktop/LGS1/LGS1_revision_8.19/Fig2D.pdf', width=3.5, height=3,pointsize=16)
p.a
dev.off()
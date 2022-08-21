library(maptools)
data(wrld_simpl)
library(raster)
library(dplyr)
library(rgdal)
library(dismo)
library(maps)
library(viridis)
library(rgeos)
library(plotrix)
library(tidyverse)

shgeo <- read.csv('curated_occ_ESB_7.26.18.csv', header=T,stringsAsFactors=F) ##1369 S. hermonthica occurrences
all<- raster('/Users/emilywork/Desktop/STRIGA/ENM/distModel7.30/preds.all.tif') ## too big to host on GitHub; users need to download locally from https://scholarsphere.psu.edu/resources/484b3bb5-22bf-4fe1-a93f-b73569b31151 
coords <- read.table('hmp_ne.lat.lon.txt', header=F)

data(worldMapEnv)

#####################
## Fig 1A
#####################
# no masking
maps::map(database="world", xlim=c(-20,60),ylim=c(-40,45),axes=TRUE)
plot(all, add=T, legend=F, col=viridis(n=100)) #joyce votes for viridis
maps::map(database="world", xlim=c(-20,60),ylim=c(-40,45),add=T, col="white", lwd=0.75)
#points(shgeo$lon, shgeo$lat, col='black',cex=0.5, pch=21)
scalebar(1000, xy = c(35,-33), type = "bar", divs = 2, lonlat = TRUE)
plot(all, legend.only=T, add=T, smallplot=c(0.2,0.25, 0.2,0.4), breaks=c(0,1), col=viridis(n=20))
plot(all, legend.only=T, add=T, smallplot=c(0.2,0.25, 0.2,0.4),col=viridis(n=100))
text(-7,-6.5, "Habitat suitability")
text(31,-34, "km")

# replace outside of 200 km w/value of 0
coordinates(shgeo) <- ~lon+lat
projection(shgeo) <- CRS('+proj=longlat +datum=WGS84')
x <- circles(shgeo, d=200000, lonlat=T) #sample within radius of 200km
pol <- polygons(x)

# make transparent
mygrey <- rgb(190, 190, 190, max = 255, alpha = 130, names = "mygrey50")

# crop to polygon
all_within <- mask(all, pol)
all_out <- mask(all, pol, inverse=T)

# fig w/masking
pdf('/Users/emilywork/Desktop/LGS1/LGS1_revision_8.19/Fig1a.pdf')
maps::map(database="world", xlim=c(-20,60),ylim=c(-40,45), lwd=0.5)
plot(all_within, add=T, legend=F, col=viridis(n=100)) #joyce votes for viridis
maps::map(database="world", xlim=c(-20,60),ylim=c(-40,45),add=T, col="black", lwd=0.5)
plot(all_out, legend=F, col=viridis(n=100), alpha=0.5, add=T)
scalebar(1000, xy = c(35,-33), type = "bar", divs = 2, lonlat = TRUE)
color.legend(-10, -30, -5, -10, rect.col=viridis(n=100), c("0.0","0.2","0.4","0.6","0.8","1.0"), gradient="y", align="rb")
#text(-7,-6.5, "Habitat suitability")
#text(31,-34, "km")
dev.off()

#####################
## Fig 1B
#####################
striga_score <- read.table('global_200.fam', header=F) # striga score in V6
bap <- read.table('referenced_bap.txt', header=T)
bap <- subset(bap, Latitude <30)

# histogram
pdf('/Users/emilywork/Desktop/LGS1/LGS1_revision_8.19/hss_sorg_dist.pdf', width = 2, height=2)
ggplot(striga_score, aes(x=V6)) + geom_histogram(bins=20, fill="grey90", col="black") + theme_classic() + xlab("") + ylab("") + xlim(c(-0.05,1))
dev.off()

# map
pdf('/Users/emilywork/Desktop/LGS1/LGS1_revision_8.19/Fig1.D.pdf')
par(fig = c(0,1,0,1))
maps::map(database="world", xlim=c(-20,100),ylim=c(-40,45), col="grey70", lwd=0.5)
points(coords$V3, coords$V2, col='grey20',cex=0.1, pch=20)
points(bap$Longitude, bap$Latitude, col='green3',cex=0.3, pch=20)
#scalebar(5000, xy = c(0,-39), type = "line", divs = 2, lonlat = TRUE)
dev.off()

#####################
## Fig 1C
#####################
# map
pdf('/Users/emilywork/Desktop/LGS1/LGS1_revision_8.19/Fig1.E.pdf')
par(fig = c(0,1,0,1))
maps::map(database="world", xlim=c(-20,100),ylim=c(-40,45), col="grey70", lwd=0.5)
points(shgeo$lon, shgeo$lat, col='grey20',cex=0.1, pch=20)
scalebar(5000, xy = c(0,-39), type = "line", divs = 2, lonlat = TRUE)
dev.off()

# histogram
pdf('/Users/emilybellis/Desktop/LGS1_revision_2019/enm_figs/hss_dist.pdf', width = 2, height=2)
ggplot(occ, aes(x=hss)) + geom_histogram(bins=20, fill="grey90", col="black") + theme_classic() + xlab("") + ylab("")
dev.off()

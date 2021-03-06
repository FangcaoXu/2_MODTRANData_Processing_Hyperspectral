#
#
#  "`-''-/").___..--''"`-._
# (`6_ 6  )   `-.  (     ).`-.__.`)   WE ARE ...
# (_Y_.)'  ._   )  `._ `. ``-..-'    PENN STATE!
#   _ ..`--'_..-_/  /--'_.' ,'
# (il),-''  (li),'  ((!.-'
#
#
# Author: Guido Cervone (cervone@psu.edu) and Fangcao Xu (xfangcao@psu.edu)
#         Geoinformatics and Earth Observation Laboratory (http://geolab.psu.edu)
#         Department of Geography and Institute for CyberScience
#         The Pennsylvania State University
#

library(colorRamps)
library(data.table)
library(dplyr)
library(gtable)
library(gtools)
library(stringr)
library(scales)
library(tidyverse)
library(viridis) 
library(RColorBrewer)

### This script is for the predicted matrix csv files with the local time

source("Z_global_variables.R")
source("Z_MODTRAN_functions.R")
source("Z_plot_functions.R")
source("Z_general_functions.R")


### -----------------------------------------------------Test Test Test -----------------------------------------------------------------------------
### -----------------------------------------------------Test Test Test -----------------------------------------------------------------------------
### ------------------------------------ Understand the difference of files with suffix *_scan.csv and .csv -----------------------------------------
# result1<- head(read.csv("/home/graduate/fbx5002/MODTRAN6/Radianceout_scan.csv", skip=6, header = FALSE, stringsAsFactors = FALSE),-1)
# result2<- head(read.csv("/home/graduate/fbx5002/MODTRAN6/Radianceout.csv", skip=6, header = FALSE, stringsAsFactors = FALSE),-1)
# colnames(result1) <- c("freq", "path_trans", "path_emiss","path_thermal_scat","surface_emiss", "path_multi_scat",
#                        "sing_scat","grnd_rflt","drct_rflt","total_rad","irrad_ref","irrad_obs","-natlog_path_trans","direct_emiss",
#                        "ToA_irrad", "bbody_temp")
# colnames(result2) <- c("freq", "path_trans", "path_emiss","path_thermal_scat","surface_emiss", "path_multi_scat",
#                        "sing_scat","grnd_rflt","drct_rflt","total_rad","irrad_ref","irrad_obs","-natlog_path_trans","direct_emiss",
#                        "ToA_irrad", "bbody_temp")
# result2$wavelength <- 1e4/result2$freq
# result2 <- result2[which(result2$wavelength<=13 & result2$wavelength>=7),-1]
# result2 <- result2[order(result2$wavelength),]
# result2 <- cbind(result2[,16], result2[,-16])
# colnames(result2)[1] <- "wavelength"
# rownames(result2) <- 1: nrow(result2)


### -----------------------------------------------------Codes start from here -----------------------------------------------------------------------------
### -----------------------------------------------------Codes start from here -----------------------------------------------------------------------------
### -----------------------------------------------------Codes start from here -----------------------------------------------------------------------------

# Section 1: Load all files -------------------------------------------------------------------
# load all csv results generated by each set of geometry parameters
# the order of listed files follows the sequence of 1,10,11,12,...2,20,21,22....
### Codes below did the same as mixedsort
# re <- regexpr("([0-9])+\\.csv$",temp.digit) # \\ means escape the dot with a backslapsh (which is escaped by another backslash)
# rm <- regmatches(temp.digit,re)
# dig <- as.numeric(gsub("\\.csv","",rm))
# o <- order(dig)
# temp.digit <- temp.digit[o]
temp.digit <- list.files('~/geolab_storage_V3/data/DARPA/7.5_12um/FixHeight',pattern="*[[:digit:]].csv",full.names=TRUE)
temp.digit <- mixedsort(temp.digit)
### load all scan.csv files --------------------------------------------------
temp.scan <- list.files('~/geolab_storage_V3/data/DARPA/7.5_12um/FixHeight',pattern="*scan.csv",full.names=TRUE)
temp.scan <- mixedsort(temp.scan)

# Section 2: Combine geometric information with radiance information ----------------------------
files.digit <- MODTRAN.readfiles.digit(temp.digit)
files.scan <- MODTRAN.readfiles.scan(temp.scan)
table.digit <- MODTRAN.combineGeo.digit(files.digit, locs[,c(1,2)])[[1]]
table.scan <- MODTRAN.combineGeo.scan(files.scan, locs[,c(1,2)])[[1]]
# reformat the data can split into different groups: split(), group_by(), gather(), spread(), slice()
# Fix phi
dataplot.theta.digit <- split(table.digit, table.digit$phi)
dataplot.theta.scan <- split(table.scan, table.scan$phi)
# Fix theta
dataplot.phi.digit <- split(table.digit, table.digit$theta)
dataplot.phi.scan <- split(table.scan, table.scan$theta)

# Section 3: Have fun with different plots -------------------------------------------------------
### margin (top,right,bottom,left)
# legend.justification is the center point of the legend
# legend.key.size: the area size of keys = legend.key.height & legend.key.width
label <- c(`total_rad`="Total radiance at sensor",`path_emiss` = "Path thermal emission radiance", `surface_emiss`="Surface thermal emission radiance",
           `grnd_rflt`="Ground reflected radiance", `path_thermal_scat`="Path thermal scattering radiance", `sing_scat`="Single scattering solar radiance", 
           `path_multi_scat`="Path multiple scattering solar radiance")
layout <- matrix(c(1,1,2,5,3,4,6,7),nrow=4,ncol=2, byrow=TRUE)
layout <- matrix(c(1,2,3,4,5,6),nrow=2,ncol=3, byrow=TRUE)
layout <- matrix(c(2,5,4,3,6),nrow=1,ncol=5, byrow=TRUE)
# Plot using the facet
# Fix phi or theta for plotting
# Only the data.plot.theta.digit[[1]] has theta 0 since we didn't read files with theta 0, azimuths > 0
boxplotfunctheta.digit(dataplot.theta.digit[[1]], layout, 3500,5000,
             title=as.expression(bquote(phi ~ "=" ~.(phi.list[1])^degree ~ theta ~ "= [30-90"^degree*"]"~"2:00 pm")), filename="FixedPhiDigitBoxplot2PM.jpg",legend.com=F)
boxplotfunctheta.scan(dataplot.theta.scan[[1]], layout, 4000,2400, "gray", 0.02,
             title=as.expression(bquote(phi ~ "=" ~.(phi.list[1])^degree ~ theta ~ "= [30-90"^degree*"]"~"2:00 pm")), filename="FixedPhiScanBoxplot2PM.jpg",legend.com=F)
boxplotfuncphi.digit(dataplot.phi.digit[[1]], layout, 3500,5000, 
             title=as.expression(bquote(theta ~ "=" ~.(rev(theta.list)[1])^degree ~ phi ~ "= [0-360"^degree*"]"~"2:00 pm")), filename="FixedThetaDigitBoxplot2PM.jpg")
boxplotfuncphi.scan(dataplot.phi.scan[[1]], layout, 4000,2400, "gray", 0.02,
             title=as.expression(bquote(theta ~ "=" ~.(rev(theta.list)[1])^degree ~ phi ~ "= [0-360"^degree*"]"~"2:00 pm")), filename="FixedThetaScanBoxplot2PM.jpg")

# two rows layout
boxplotfunctheta.scan(dataplot.theta.scan[[1]], layout, 3600,2200, "gray", 0.02, title=NULL, filename="FixedPhiBoxplot.jpg",legend.com=F)
boxplotfuncphi.scan(dataplot.phi.scan[[1]], layout, 3600,2200, "gray", 0.02, title=NULL, filename="FixedThetaBoxplot.jpg")
# one row layout
boxplotfunctheta.scan(dataplot.theta.scan[[1]], layout, 5600,1000, "gray", 0.02, filename="FixedHeightPhiScanBoxplot2PM.jpg",legend.com=F)
boxplotfuncphi.scan(dataplot.phi.scan[[1]], layout, 5600,1000, "gray", 0.02, filename="FixedHeightThetaScanBoxplot2PM.jpg",legend.com=F)


# Section 4: Single Scattering Solar Radiance ----------------------------------------------------------------------------------
# Check the sun azimuth influence on the single solar scattering radiance given the location of the sensor 
singscat.ex<- dataplot.phi.scan[[1]][,c(1,7,18)]
singscat.ex<-singscat.ex[which(singscat.ex$wavelength == singscat.ex$wavelength[which.max(singscat.ex$sing_scat)]),]
singscat.ex<-singscat.ex[-nrow(singscat.ex),]  # Remove 360 because it's the same value as 0
n   <- round(normalize1nmax( singscat.ex[,2], nmax=nrow(singscat.ex) ) )
#col <- rev(colorRampPalette(brewer.pal(9, "RdGy"))(length(n)))
col <- colorRampPalette(rev(brewer.pal(11,"Spectral")))(length(n))
angle = 203
p203 <-xu.pie.xy(angle,r=1.2)
p23 <- xu.pie.xy(angle-180,r=1.2)
t203 <- xu.pie.xy(angle,r=.6)
legend.text<- seq(min(singscat.ex[,2]), max(singscat.ex[,2]), length=5)
legend.col <- seq(1,length(col), length=5)
# Plot pie chart
jpeg("piechart_SingleSolarScattering.jpg", width= 800, height=1000, units = "px", quality = 100, res=300)
par(mai = c(0.1, 0, 0, 0)) #(bottom, left, top, right)
xu.pie(rep(1,length(n)), col=col[n], clockwise = TRUE, labels=seq(0,355,5), border = FALSE)
lines(c(p203[,1], p23[,1]), c(p203[,2], p23[,2]), lty=2,col="grey")
points(p203,col="orange",pch=19,cex=2)
text(t203, paste("Sun azimuth = 203"),pos=4,srt=64,col="grey",offset=0.5, cex=0.4)
legend('bottom',horiz=T,fill=col[legend.col], legend=round(legend.text*1E9,2), bty="n", yjust = 0, cex=0.5)
mtext(expression(paste("(W cm-2 sr-1 ", mu, "m-1) x ",10^-09)), 1, line=-0.5, cex=0.5)
title("2019-04-18 14:00",line = -1,cex.main=0.5)
dev.off()

### ---------------------------------------------------------New Section Starts from here------------------------------------------------------------------------------
### ---------------------------------------------------------New Section Starts from here------------------------------------------------------------------------------
### ---------------------------------------------------------New Section Starts from here------------------------------------------------------------------------------
label <- c(`total_rad`="Total radiance at sensor",`path_emiss` = "Path thermal emission radiance", `surface_emiss`="Surface thermal emission radiance",
           `grnd_rflt`="Ground reflected radiance", `path_thermal_scat`="Path thermal scattering radiance", `sing_scat`="Single scattering solar radiance", 
           `path_multi_scat`="Path multiple scattering solar radiance")
layout <- matrix(c(1,1,2,5,3,4,6,7),nrow=4,ncol=2, byrow=TRUE)

# Section 5: Plots for the fixed range  --------------------------------------------------
# Daytime 2:00 pm
temp.digit1 <- list.files('~/geolab_storage_V3/data/DARPA/7.5_12um/FixRange',pattern="*[[:digit:]].csv",full.names=TRUE)
temp.digit1 <- mixedsort(temp.digit1)
temp.scan1 <- list.files('~/geolab_storage_V3/data/DARPA/7.5_12um/FixRange',pattern="*scan.csv",full.names=TRUE)
temp.scan1 <- mixedsort(temp.scan1)
files.digit1 <- MODTRAN.readfiles.digit(temp.digit1)
files.scan1 <- MODTRAN.readfiles.scan(temp.scan1)
table.digit1 <- MODTRAN.combineGeo.digit(files.digit1, eles)[[1]]
table.scan1 <- MODTRAN.combineGeo.scan(files.scan1, eles)[[1]]
# For Error in as.vector(x, mode) : cannot coerce type 'closure' to vector of type 'any', check whether you have set the layout
# function(plotdata, layout, imgwidth, imgheight, boxcolor="black",boxlwd=0.2, title=NULL, filename="", 
#          fixScatteringSolar=F, legend.pos=c(0.005,0.99) ,legend.com = T)
boxplotfunctheta.digit(table.digit1, layout, title=as.expression(bquote(phi ~ "= 0"^degree ~ theta ~ "= [30-90"^degree*"]"~"2:00 pm")), legend.com=F)
boxplotfunctheta.scan(table.scan1, layout, boxcolor="gray", boxlwd= 0.02, title=as.expression(bquote(phi ~ "= 0"^degree ~ theta ~ "= [30-90"^degree*"]"~"2:00 pm")), legend.com=F)
boxplotfunctrans.comp.digit(table.digit1, dataplot.theta.digit[[1]], boxcolor="gray", boxlwd= 0.1,title=as.expression(bquote(phi ~ "=" ~.(phi.list[1])^degree ~ theta ~ "= [30-90"^degree*"]"~"2:00 pm")))
boxplotfunctrans.comp.scan(table.scan1, dataplot.theta.scan[[1]], boxcolor= "gray", boxlwd= 0.02, title=as.expression(bquote(phi ~ "=" ~.(phi.list[1])^degree ~ theta ~ "= [30-90"^degree*"]"~"2:00 pm")))
# export the plots to jpegs
boxplotfunctheta.digit(table.digit1, layout, 3500,5000,
                  title=as.expression(bquote(phi ~ "= 0"^degree ~ theta ~ "= [30-90"^degree*"]"~"2:00 pm")), filename="FixedRangeDigitBoxplot2PM.jpg",legend.com=F)
boxplotfunctheta.scan(table.scan1, layout, 3500,5000, boxcolor="gray", boxlwd= 0.02,
                  title=as.expression(bquote(phi ~ "= 0"^degree ~ theta ~ "= [30-90"^degree*"]"~"2:00 pm")), filename="FixedRangeScanBoxplot2PM.jpg",legend.com=F)
# transmission comparison between fixed range and fixed height
boxplotfunctrans.comp.digit(table.digit1, dataplot.theta.digit[[1]], 5000,2000, boxcolor="gray", boxlwd= 0.1,
                      title=as.expression(bquote(phi ~ "=" ~.(phi.list[1])^degree ~ theta ~ "= [30-90"^degree*"]"~"2:00 pm")), filename="TransmissionDigitComp.jpg")
boxplotfunctrans.comp.scan(table.scan1, dataplot.theta.scan[[1]], 5000,2000, boxcolor= "gray", boxlwd= 0.02,
                      title=as.expression(bquote(phi ~ "=" ~.(phi.list[1])^degree ~ theta ~ "= [30-90"^degree*"]"~"2:00 pm")), filename="TransmissionScanComp.jpg")

### ---------------------------------------------------- Day & Nighttime
temp.digit2<- list.files('~/disk10TB/DARPA/MODTRANSimulated/MidLatitude/springsummer',pattern="*107.*10.csv",full.names=TRUE)
temp.digit2 <- mixedsort(temp.digit2)
temp.scan2 <- list.files('~/disk10TB/DARPA/MODTRANSimulated/MidLatitude/springsummer',pattern="*107.*10.csv",full.names=TRUE)
temp.scan2 <- mixedsort(temp.scan2)
# Radiance11_2.csv: 11 means 11th geometric setting and 2 means 2:00 GMT
for(i in 1:6){
  time <- Localtime[i]
  temp.scan2.bytime <- temp.scan2[seq(i, 78, 6)]
  files.scan2 <- MODTRAN.readfiles.scan(temp.scan2.bytime)
  table.scan2 <-  MODTRAN.combineGeo.scan(files.scan2, eles)[[1]]
  boxplotfunctheta.scan(table.scan2, layout, 3500,5000, "gray", 0.02,
                    title=as.expression(bquote(phi ~ "= 0"^degree ~ theta ~ "= [30-90"^degree*"]"~.(time)*":00")), filename=paste("FixedRangeScanBoxplot",time,".jpg",sep = ""), fixScatteringSolar =TRUE, legend.com=F)
}

### 
temp.scan3 <- list.files('~/disk10TB/DARPA/MODTRANSimulated/MidLatitude_VNIR_SWIR/OneYear_2hours', pattern="*107_16_100_scan.csv",full.names=TRUE)
filename.pattern=c("Geometry", "azimuth","Day","Time","Reflectivity")
temp.scan3 <- mixedsort(temp.scan3)
vars <- return.index(temp.scan3, range=1:length(filename.pattern), names = filename.pattern)
try <- MODTRAN.calculate.downwelling(MODTRAN.combineGeo(MODTRAN.readfiles.scan(temp.scan3), eles))
plot(x=try$wavelength, y=try$path_trans, col=try$theta)

# Section 6: Research Sun locations for given time ----------------------------------------------------------------------------------
# https://www.usno.navy.mil/USNO/ry <- atory Website is out of service. It is undergoing modernization and the expected completion of work is estimated as 30 April 2020 
day <- 18
month <- 4
time <- "14:00"
state <- "PA"
city <- "State+College"
url <- "https://aa.usno.navy.mil/cgi-bin/aa_altazw.pl?form=1&body=10&year=2019&month=#MONTH#&day=#DAY#&intv_mag=10&state=#STATE#&place=#CITY#"
myurl <- url
myurl <-gsub("#MONTH#",month,myurl)
myurl <-gsub("#DAY#",day,myurl)
myurl <-gsub("#STATE#",state,myurl)
myurl <-gsub("#CITY#",city,myurl)
sun.pos <- scan(myurl, what="character")
sun.pos <-sun.pos[-(1:85)]
sun.pos <-sun.pos[-( (length(sun.pos)-6):length(sun.pos)) ]
sun.pos <- matrix(sun.pos, ncol=3, byrow=T)
sun.df  <- data.frame(hour=as.POSIXct(paste("2019-",month,"-",day," ",sun.pos[,1],sep=""),origin="19070-01-01",tz="EST5DT"),
                                     alt=as.numeric(sun.pos[,2]), azim=as.numeric(sun.pos[,3]))
target.time       <- as.POSIXct(paste("2019-",month,"-",day," ",time,sep=""),origin="19070-01-01",tz="EST5DT")
target.time.index <- which.min( abs(sun.df$hour - target.time) )
target.azimuth <- sun.df$azim[target.time.index]
#plot(sun.df$hour, sun.df$alt)
plot(sun.df$hour, sun.df$azim,pch=16,col="grey",ylab="Azimuth",xlab="Time")
points(target.time, target.azimuth, col="red2",cex=2)
lines(c(target.time,target.time),c(0,360),lty=2,col="red2")
lines(c(min(sun.df$hour), max(sun.df$hour)), c(target.azimuth, target.azimuth), col="red2",lty=2)
text(min(sun.df$hour)+2000, target.azimuth, target.azimuth,pos=3, col="red2")
text(target.time, min(sun.df$azim)+10,time, pos=2,srt=90,col="red2")


# Section 7: Check the Distribution of the Simulated Data over Year --------
temp.scan.differentdays <- list.files('~/disk10TB/DARPA/MODTRANSimulated/MidLatitude/springsummer',pattern="*[A-z]0_[0-9]+_14_10_scan.csv",full.names=TRUE)
temp.scan.differentdays <- mixedsort(temp.scan.differentdays)
files.scan.differentdays <- MODTRAN.readfiles.scan(temp.scan.differentdays)
table.scan.differentdays <-  as.data.frame(MODTRAN.combineGeo.scan(files.scan.differentdays, data.frame("days"=matrix(79:265)))[[1]])
boxplotfunc.scan(table.scan.differentdays, "days", layout, 3500,5000, "white", 0.02, legend.com=T, legend.continous = TRUE, 
                 title="Distribution of Simulated Data over Days", filename="Distribution_ovdays.jpg")


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
library(gtools)
library(RColorBrewer)
library(scales)
library(viridis) 

source("Z_global_variables.R")
source("Z_MODTRAN_functions.R")
source("Z_plot_functions.R")
source("Z_general_functions.R")

# Valid output in our MODTRAN simulations
#[1] "wavelength"         "path_trans"         "path_emiss"         "path_thermal_scat"  "surface_emiss"      "path_multi_scat"   
#[7] "sing_scat"          "grnd_rflt"          "drct_rflt"          "total_rad"          "irrad_ref"          "irrad_obs"         
#[13] "-natlog_path_trans" "direct_emiss"       "ToA_irrad"          "bbody_temp"         "theta"              "elevation"  

variable <- "direct_emiss"
layout <- matrix(c(1,1,2,5,3,4,6,7),nrow=4,ncol=2, byrow=TRUE)

temp.scan <- list.files('~/disk10TB/DARPA/MODTRANSimulated/MidLatitude/springsummer',pattern="*107.*scan.csv",full.names=TRUE)
temp.scan <- mixedsort(temp.scan) # geometric elevation angle, day, time, reflectivity (i.e., 13*365*6*7)

# degree 90, 2 pm 
files90 <- temp.scan[which(sapply(temp.scan, function(x) grepl("Radiance1",x)))]
files90.2pm <- files90[which(sapply(files90, function(x) grepl("Radiance1_107_18",x)))]
files90.scan <- MODTRAN.readfiles.scan(files90.2pm)
dataplot90.2pm <- as.data.frame(do.call(rbind, files90.scan))

boxplotfunc.scan(dataplot90.2pm, variable, layout, boxcolor="white", boxlwd=0.02, title=as.expression(bquote("day:107"~phi ~ "= 0"^degree ~ theta ~ "= 90"^degree~"14:00")),legend.com=F)
boxplotfunc.scan(dataplot90.2pm, variable, layout, 3500,5000, "white", 0.02,
                      title=as.expression(bquote("day:107"~phi ~ "= 0"^degree ~ theta ~ "= 90"^degree~"14:00")), filename=paste("EmissivityBoxplot1",".jpg",sep = ""),legend.com=F)

# degree 40, 2 pm
files40 <- temp.scan[which(sapply(temp.scan, function(x) grepl("Radiance10",x)))]
files40.2pm <- files40[which(sapply(files40, function(x) grepl("Radiance10_107_18",x)))]
files40.scan <- MODTRAN.readfiles.scan(files40.2pm)
dataplot40.2pm <- as.data.frame(do.call(rbind, files40.scan))

boxplotfunc.scan(dataplot40.2pm, variable, layout, boxcolor="white", boxlwd=0.02, title=as.expression(bquote("day:107"~phi ~ "= 0"^degree ~ theta ~ "= 40"^degree~"14:00")),legend.com=F)
boxplotfunc.scan(dataplot40.2pm, variable, layout, 3500,5000, boxcolor="white", boxlwd=0.02,
                 title=as.expression(bquote("day:107"~phi ~ "= 0"^degree ~ theta ~ "= 40"^degree~"14:00")), filename=paste("EmissivityBoxplot2",".jpg",sep = ""),legend.com=F)


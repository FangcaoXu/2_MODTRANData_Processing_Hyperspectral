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


################# polyethylene reflectivity retrieved from https://webbook.nist.gov/cgi/cbook.cgi?ID=C9002884&Units=SI&Type=IR-SPEC&Index=1#
source("Z_general_functions.R")

library(dplyr)  
library(hyperSpec)

# Grass
# reflec <- read.delim('~/disk10TB/DARPA/BH_analysis/MaterialsReflectivity/Full/vegetation.grass.unknown.unknown.all.grass.jhu.becknic.spectrum.txt', header = FALSE, stringsAsFactors = FALSE)
# reflec <- reflec[-which(reflec==""),] # vector
# reflec <- data.frame(cbind(as.numeric(reflec[seq(21, length(reflec), 2)]), 1-as.numeric(reflec[seq(22, length(reflec), 2)])/100))
reflec <- read.csv('~/disk10TB/DARPA/MatrixCSV/Stage2_BHData/MODTRAN/DifferentMaterials/Reflectivity/LAMB_GRASS_reflec.csv', 
                   header = T, stringsAsFactors = FALSE, row.names = 1)
colnames(reflec) <- c("wavelength", "reflectivity")
load('~/disk10TB/DARPA/BH_analysis/SpectralGeoDataTable/wavelength_side.RData')
# find nearest wavelengths
idx1<- sapply(wavelength_side[1, -1], function(x) findnearest(reflec$wavelength, x))
idx2 <- sapply(wavelength_side[2, -1], function(x) findnearest(reflec$wavelength, x))
reflec.side1<- reflec[idx1, 2]
reflec.side2<- reflec[idx2, 2]
grass.reflec <- rbind(reflec.side1, reflec.side2)
write.csv(grass.reflec, "~/disk10TB/DARPA/BH_analysis/TrainingData/Grass_reflec.csv", row.names = F)


# Polyethylene,  CAS Registry Number: C9002884
jdx.file <- read.jdx('~/disk10TB/DARPA/jdxfiles_NISTBook/polyethylene.jdx')
wavenumber <- jdx.file@wavelength # length: 7153, cm-1
wavelength <- 10000/wavenumber    # from 1.333684 to 16.6689 um
# jdx.file@data: (2,2) dataframe. colnames: "spc", "filename"
spectral.info <- jdx.file@data$spc # (2, 7153) matrix
# weird column names of spectral.info
reflec.1 <- unname(spectral.info[1, ]) # length: 7153
reflec.2 <- unname(spectral.info[2,])  # length: 7153
# fg="transparent", color for the foreground (i.e., axes, boxes around plots)
# bg="red", background color
seg.index <- which(wavelength>=7.5&wavelength<12)
wavelength.subset <- wavelength
wavelength.subset[seg.index] <- NA

jpeg("polyethylene.jpg", width= 2000, height=1250, units = "px", quality = 100, res=300)
#par(mfrow=c(2,1), mar=c(2.5,3,1,1), oma=c(0,0,0,0), mgp=c(1.5, 0.5, 0)) 
#plot(wavenumber, reflec.1, type="l", xlab="Wavenumber (cm-1)", ylab="Reflectivity [0,1]")
par(mar=c(2.5,3,1,1), oma=c(0,0,0,0), mgp=c(1.5, 0.5, 0)) 
plot(wavelength.subset, reflec.1, type="l", xlab=expression(paste("Wavelength (", mu, "m)")), ylab="Reflectivity [0,1]")
lines(wavelength[seg.index], reflec.1[seg.index], type="l", col="red")
dev.off()
# # write the reflectivity to the csv
reflec.poly <- reflec.1[findnearest(wv, wavelength)] # 256
plot(wv, reflec.poly, type="l", ylim=c(0,1),xlab=expression(paste("Wavelength (", mu, "m)")), ylab="reflectivity")
write.csv(cbind(wv, reflec.poly), '~/disk10TB/DARPA/MatrixCSV/DifferentMaterials/Polyethylene/polyethylene_info/polyethtlene_reflectivity.csv')

# transmissivity + absorptance + reflectivity = 1
# most solid objects exhibit very low transmission of infrared energy
jdx.file <- read.jdx('~/disk10TB/DARPA/jdxfiles_NISTBook/ammonia.jdx')
wavenumber <- jdx.file@wavelength   # length(wavenumber) 
wavelength <- 10000/wavenumber    
# jdx.file@data: dataframe. colnames: "spc", "filename"
spectral.info <- jdx.file@data$spc 
# weird column names of spectral.info
absorp <- spectral.info[1, ]
jpeg("ammonia.jpg", width= 2000, height=2500, units = "px", quality = 100, res=300)
par(mfrow=c(2,1), mar=c(2.5,3,1,1), oma=c(0,0,0,0), mgp=c(1.5, 0.5, 0)) 
plot(wavenumber, absorp, type="l", xlab="Wavenumber (cm-1)", ylab="Absorptance")
grid()
plot(wavelength, absorp, type="l", xlab=expression(paste("Wavelength (", mu, "m)")), ylab="Absorptance")
grid()
dev.off()
# write the absorptance to the csv
absorp <- absorp[findnearest(wv, wavelength)] # 256
plot(wv, absorp, type="l", ylim=c(0,1),xlab=expression(paste("Wavelength (", mu, "m)")), ylab="Absorptance")
write.csv(cbind(wv, absorp), '~/disk10TB/DARPA/MatrixCSV/DifferentMaterials/Ammonia/ammonia_info/ammonia_absorptance.csv')

# acetone
jdx.file <- read.jdx('~/disk10TB/DARPA/jdxfiles_NISTBook/acetone.jdx')
wavenumber <- jdx.file@wavelength   # length(wavenumber) 
wavelength <- 10000/wavenumber    
# jdx.file@data: dataframe. colnames: "spc", "filename"
spectral.info <- jdx.file@data$spc 
# weird column names of spectral.info
absorp <- spectral.info[1, ]
jpeg("acetone.jpg", width= 2000, height=2500, units = "px", quality = 100, res=300)
par(mfrow=c(2,1), mar=c(2.5,3,1,1), oma=c(0,0,0,0), mgp=c(1.5, 0.5, 0)) 
plot(wavenumber, absorp, type="l", xlab="Wavenumber (cm-1)", ylab="Absorptance")
grid()
plot(wavelength, absorp, type="l", xlab=expression(paste("Wavelength (", mu, "m)")), ylab="Absorptance")
grid()
dev.off()
# write the absorptance to the csv
absorp <- absorp[findnearest(wv, wavelength)] # 256
plot(wv, absorp, type="l", ylim=c(0,1),xlab=expression(paste("Wavelength (", mu, "m)")), ylab="Absorptance")
write.csv(cbind(wv, absorp), '~/disk10TB/DARPA/MatrixCSV/DifferentMaterials/Acetone/acetone_info/acetone_absorptance.csv')


#####################################################################################  ExoScan
# pipe operator
exoscan <- list.files("~/disk10TB/DARPA/jdxfiles_NISTBook/exoscan", full.names = TRUE) %>% .[-length(.)]
# * at least 0 time, * at least one time
# [^_] ^: except
# ([^_]*_){2}: /home/graduate/fbx5002/disk10TB/DARPA/jdxfiles_; NISTBook/exoscan/2019-04-18T14-17-56_
# {}, ? only works for () which means match at most one time or 
index <- data.frame("site"=gsub("^([^_]*_){2}([^_]*)?_(.*)\\.csv", "\\2", exoscan),
                   "material" = gsub("^([^_]*_){2}([^_]*)?_(.*)\\.csv", "\\3", exoscan)) 

# read the file; change the column names; convert the wavenumber to wavelength; set the reflectivity between [0,1]
exoscan.read <- lapply(exoscan, function(x) read.csv(x, skip = 1, header = TRUE)[, 1:2]) %>% 
  lapply(., setNames, c("wavenumber", "reflectivity")) %>% 
  lapply(., function(x){x$wavelength=10000/x$wavenumber; return(x)}) %>%
  lapply(., function(x){x$reflectivity = ifelse(x$reflectivity>100, 1, x$reflectivity/100); return(x)})

exoscan.read[[1]][which(exoscan.read[[1]]$wavelength >= 7.5 & exoscan.read[[1]]$wavelength <=12),]

# write to csv
for(i in 1:length(exoscan.read)){
  temp <- exoscan.read[[i]]
  # use the site and material as the filename
  basename <-paste(index[i,1], index[i,2],sep='_')
  # plot the reflectivity against the wavelength and wavenumber
  jpeg(paste('~/disk10TB/DARPA/MatrixCSV/DifferentMaterials/exoscan/exoscan_info/', basename, '.jpg', sep=''), width= 4700, height=3200,
       units = "px", quality = 100, res=300)
  par(mar=c(3,3,2,1), oma = c(2,2,0,0), mgp=c(2, 0.5, 0))
  # par(mfrow=c(2,1), mar=c(2.5,3,1,1), oma=c(0,0,0,0), mgp=c(1.5, 0.5, 0)) 
  # plot(temp$wavenumber,temp$reflectivity, type="l", xlab="Wavenumber (cm-1)", ylab="reflectivity")
  # grid(lwd=0.7)
  plot(temp$wavelength, temp$reflectivity, type="l", xlab=expression(paste("Wavelength (", mu, "m)")), ylab="Reflectivity [0,1]")
  grid(lwd=0.1)
  dev.off()
  # # plot only the longwave hyperspectral reflectivity
  # reflec <- temp$reflectivity[findnearest(wv, temp$wavelength)]
  # # plot(wv, reflec, type="l", ylim=c(0,1),xlab=expression(paste("Wavelength (", mu, "m)")), ylab="reflectivity")
  # write.csv(cbind(wv, reflec), paste('~/disk10TB/DARPA/MatrixCSV/DifferentMaterials/exoscan/exoscan_info/', basename, '.csv', sep = ""))
}

# 
# #full spectral error added
# # plot retrieved reflectivity for different elevation angles in one panel as the boxplot
# plot_randomerror_reflectivity_boxplot <- function(true_reflectivity, filename){
#   # expand the true reflectivity to a matrix
#   poly.mat <- matrix(rep(true_reflectivity[,2],8), ncol=8)
#   wv.error <- runif(nrow(true_reflectivity), min = 0, max = 0.01)
#   angle.error <- lapply(1:nrow(true_reflectivity), function(x) rnorm(8, mean=true_reflectivity[,2][x], sd = wv.error[x]))
#   true_reflectivity.randomerror <- matrix(unlist(angle.error), nrow=nrow(true_reflectivity), byrow=TRUE)
#   true_reflectivity.randomerror[which(true_reflectivity.randomerror < 0)] <- 0
#   # plot the retrieved reflectivity versus true reflectivity
#   jpeg(filename, width= 4700, height=3200, units = "px", quality = 100, res=300)
#   par(mar=c(2,3,2,1), oma = c(2,2,0,0), mgp=c(2, 0.5, 0))
#   boxplot(t(true_reflectivity.randomerror), xlab=expression(paste("Wavelength (", mu, "m)")), ylab="Reflectivity [0,1]",
#           main="Boxplot for Reflectivity at Different Angles with random error added")
#   lines(as.factor(true_reflectivity[,1]), true_reflectivity[,2], lwd=1, col="blue")
#   legend('top', bty='n', legend="True", col="blue", lty=1, seg.len=1, cex=1, x.intersp=0.4)
#   dev.off()
# }
# 
# 
# for(i in 1:length(exoscan.read)){
#   reflectivity <- exoscan.read[[i]][, c(3,2)]
#   plot_randomerror_reflectivity_boxplot(reflectivity, paste('~/disk10TB/DARPA/MatrixCSV/DifferentMaterials/exoscan/exoscan_info/', 
#                                                             sub('\\.csv$', '', basename(reflectivities.exoscan))[i], '_randomError.jpg', sep=""))
# }

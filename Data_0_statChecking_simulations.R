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
###################### 
library(abind)
library(gtools)
library(plyr)
source("Z_general_functions.R")
source("Z_global_variables.R")
source("Z_plot_functions.R")
#str(vars.total)

wv=seq(0.4, 3.0075, 0.0175); an=an.full
reflec = c(5,10,30,50,80,100)
# time: 6, 8, 10, 12, 14, 16, 18, 20
# reflectivity: 5, 10, 30, 50, 80, 100
# 365*8*6 = 17520


###############################################################################################################################################################
################################################## Way 1: Read all simulated files, do statistics #############################################################
###############################################################################################################################################################
total.rad <- mixedsort(list.files("~/disk10TB/DARPA/MatrixCSV/MidLatitude_VNIR_SWIR/OneYear_2hours/TotalRadianceCSV",full.names=TRUE))
vars.total <-return.index(total.rad)[,2:4]
total.rad.read <- read.originalmatrix(total.rad, wv, an)

surf.rad <-  mixedsort(list.files("~/disk10TB/DARPA/MatrixCSV/MidLatitude_VNIR_SWIR/OneYear_2hours/SurfaceEmissionCSV",full.names=TRUE))
vars.surf <-return.index(surf.rad)[,2:4]
surf.rad.read <- read.originalmatrix(surf.rad, wv, an)

trans.rad <- mixedsort(list.files("~/disk10TB/DARPA/MatrixCSV/MidLatitude_VNIR_SWIR/OneYear_2hours/TransmissionCSV",full.names=TRUE))
vars.trans <-return.index(trans.rad)[,2:4]
trans.rad.read <- read.originalmatrix(trans.rad, wv, an)

down.rad <- mixedsort(list.files("~/disk10TB/DARPA/MatrixCSV/MidLatitude_VNIR_SWIR/OneYear_2hours/DownwellingCSV",full.names=TRUE))
vars.down <-return.index(down.rad)[,2:4]
down.rad.read <- read.originalmatrix(down.rad, wv, an)

solar1.rad <- mixedsort(list.files("~/disk10TB/DARPA/MatrixCSV/MidLatitude_VNIR_SWIR/OneYear_2hours/SolarSingleScatteringCSV",full.names=TRUE))
vars.solar1 <-return.index(solar1.rad)[,2:4]
solar1.rad.read <- read.originalmatrix(solar1.rad, wv, an)

solar2.rad <- mixedsort(list.files("~/disk10TB/DARPA/MatrixCSV/MidLatitude_VNIR_SWIR/OneYear_2hours/SolarMultiScatteringCSV",full.names=TRUE))
vars.solar2 <-return.index(solar2.rad)[,2:4]
solar2.rad.read <- read.originalmatrix(solar2.rad, wv, an)

up1.rad <- mixedsort(list.files("~/disk10TB/DARPA/MatrixCSV/MidLatitude_VNIR_SWIR/OneYear_2hours/PathThermalEmissionCSV",full.names=TRUE))
vars.up1 <-return.index(up1.rad)[,2:4]
up1.rad.read <- read.originalmatrix(up1.rad, wv, an)

up2.rad <- mixedsort(list.files("~/disk10TB/DARPA/MatrixCSV/MidLatitude_VNIR_SWIR/OneYear_2hours/PathThermalScatteringCSV",full.names=TRUE))
vars.up2 <-return.index(up2.rad)[,2:4]
up2.rad.read <- read.originalmatrix(up2.rad, wv, an)

all.equal(vars.total, vars.down) # True, 17520
all.equal(vars.total, vars.up1)  # True

# return 3 dimensional array [x, y, z]; (x,y): radiance matrix size; z: filelength 17520
total.marray <- statistic.prep(total.rad.read);     dim(total.marray)
trans.marray <- statistic.prep(trans.rad.read);     surf.marray <- statistic.prep(surf.rad.read)     
down.marray <- statistic.prep(down.rad.read);       
solar1.marray <- statistic.prep(solar1.rad.read);   solar2.marray <- statistic.prep(solar2.rad.read);  
up1.marray <- statistic.prep(up1.rad.read);         up2.marray <- statistic.prep(up2.rad.read);         

# # average at certain dimensions
# dim_aver = c(2,3) # 13, 17520 (vars.down: day, time, reflec)
# total.average <- apply(total.marray, dim_aver, mean) 
# trans.average <- apply(trans.marray, dim_aver, mean)
# surf.average <- apply(surf.marray, dim_aver, mean) 
# down.average <- apply(down.marray, dim_aver, mean)
# solar1.average <- apply(solar1.marray, dim_aver, mean) 
# solar2.average <- apply(solar2.marray, dim_aver, mean)
# up1.average <- apply(up1.marray,dim_aver, mean) 
# up2.average <- apply(up2.marray, dim_aver, mean)

# seq
check <- solar1.marray; vars <- vars.solar1
seq <- split(vars, vars['Time'])
seq <- lapply(seq, function(x) as.numeric(row.names(x)))
as.numeric(names(seq))
dim(check[,,seq[[1]][1:600]]) # 150, 13, 2190 (365* 6)
dim(abind(check[,,seq[[1]]], check[,,seq[[2]]], along=4))
# winter
tmp <- do.call(abind, list(lapply(seq, function(x) check[,,x[1:600]]), along=4))
tmp.average <- apply(tmp, c(2,4), mean)
tmp.average <- tmp.average[,c(7,8, 1:6)]
# summer
tmp1 <- do.call(abind, list(lapply(seq, function(x) check[,,x[1100:1700]]), along=4))
tmp1.average <- apply(tmp1, c(2,4), mean)
tmp1.average <- tmp1.average[,c(7,8, 1:6)]
dim(t(tmp.average))

col = rev(brewer.pal(11,"Spectral"))
zlim <- range(tmp.average)
zlim1 <- range(tmp1.average)
zlim <- c(min(zlim[1], zlim1[1]), max(zlim[2], zlim1[2]))
ztitle <- bquote(atop(NA,atop("Radiance","(W cm-2 sr-1"~mu*"m-1)")))
par(mar = c(2, 3, 2, 3))
xu.filled.contour(seq(6,20,2), an, t(tmp.average), color.palette=colorRampPalette(col), nlevels = 200, zlim=zlim,
                  plot.title=title(xlab="Time", ylab="Angle (degrees)",line=1.3,cex.lab=0.7), 
                  plot.axes={axis(1, cex.axis=0.7, padj =-1.5); axis(2, cex.axis=0.7, hadj=0.4)}, 
                  key.title = title(main= ztitle, line=1.2, cex.main = 0.8, outer = FALSE, font.main = 1),
                  key.axes=axis(4, at=seq(zlim[1],zlim[2],length.out =5),labels=format(seq(zlim[1],zlim[2],length.out =5),scientific = TRUE,digits=2), cex.axis=0.7))

xu.filled.contour(seq(6,20,2), an, t(tmp1.average), color.palette=colorRampPalette(col), nlevels = 200, zlim=zlim,
                  plot.title=title(xlab="Time", ylab="Angle (degrees)",line=1.3,cex.lab=0.7), 
                  plot.axes={axis(1, cex.axis=0.7, padj =-1.5); axis(2, cex.axis=0.7, hadj=0.4)}, 
                  key.title = title(main= ztitle, line=1.2, cex.main = 0.8, outer = FALSE, font.main = 1),
                  key.axes=axis(4, at=seq(zlim[1],zlim[2],length.out =5),labels=format(seq(zlim[1],zlim[2],length.out =5),scientific = TRUE,digits=2), cex.axis=0.7))

  
# Create a dataframe to store all average values over files
average.profile <- cbind(vars.total, total.average, trans.average, surf.average, down.average, solar1.average, solar2.average, up1.average, up2.average)

# # When downwelling and other components have different length, Add the downwelling to average profile with NA values for reflectivity 0
# average.profile <- cbind(vars.total, total.average, surf.average, trans.average, solar1.average, solar2.average, up1.average, up2.average)
# down.index.inprofile <- sapply(1:nrow(vars.down), function(i) which(average.profile$Day == vars.down$Day[i]
#                                                                     & average.profile$Time== vars.down$Time[i]
#                                                                     & average.profile$Reflectivity == vars.down$Reflectivity[i]))
# length(down.average); length(down.index.inprofile)
# # Day, Time, Reflectivity, Total, Surf, Trans, Solar1, Solar2, Up1, Up2, Down
# average.profile$down.average[down.index.inprofile] <- down.average
# # Day, Time, Reflectivity, Total, Surf, Trans, Down, Solar1, Solar2, Up1, Up2
# average.profile<- average.profile[, c(1:6, 11, 7:10)]


# select one time and one reflectivity for checking 
# change at days 79, 265
average.subset <- average.profile[which(average.profile$Reflectivity==30 & average.profile$Time==14),]
par(mfrow=c(2,4))
for(i in 4:ncol(average.subset)){ # first three columns are day, time, reflectivity
  name <- colnames(average.subset)[i]
  plot(average.subset[, c(1,i)],pch=16, cex=0.5, main=paste(name,"GMT 18:00, State College (2:00 PM)", sep="\n"), 
       xlab="Day", ylab="Average Radiance over Shortwave")
  abline(v=79, col="red")
  text(x=74, y= min(average.subset[, i])+ 0.5*(max(average.subset[, i])-min(average.subset[, i])), "Day=79", col="red", srt=90)
  abline(v=265, col="red")
  text(x=270, y= min(average.subset[, i])+ 0.5*(max(average.subset[, i])-min(average.subset[, i])), "Day=265", col="red", srt=-90)
}


###############################################################################################################################################################
############################################ Way 2: Checking MODTRAN built-in materials with training dataset #################################################
################################################################ Day:108, Time:14:00 ##########################################################################
day = 108; time=14
# 42 folders, each folder is for one material, containing 9 files
originalfolders.test <- list.files('~/disk10TB/DARPA/MatrixCSV/MidLatitude_VNIR_SWIR/DifferentMaterials/MODTRAN', full.names = TRUE)
# 9 folders, each folder is for each component
originalfolders.training <- '~/disk10TB/DARPA/MatrixCSV/MidLatitude_VNIR_SWIR/OneYear_2hours'

training.trans.files <- find.originalfiles(paste("Radiance_{day}_{time}_", reflectivity.default[-1], "_trans.csv", sep='') %>% gsub("\\{day\\}", day, .)  
                                           %>% gsub("\\{time\\}", time, .),  originalfolders.training)
training.down.files <- find.originalfiles(paste("Radiance_{day}_{time}_", reflectivity.default[-1], "_down.csv", sep='') %>% gsub("\\{day\\}", day, .) 
                                          %>% gsub("\\{time\\}", time, .),  originalfolders.training)
training.up1.files <- find.originalfiles(paste("Radiance_{day}_{time}_", reflectivity.default[-1], "_up1.csv", sep='') %>% gsub("\\{day\\}", day, .) 
                                         %>% gsub("\\{time\\}", time, .),  originalfolders.training)
training.up2.files <- find.originalfiles(paste("Radiance_{day}_{time}_", reflectivity.default[-1], "_up2.csv", sep='') %>% gsub("\\{day\\}", day, .) 
                                         %>% gsub("\\{time\\}", time, .),  originalfolders.training)
training.solar1.files <- find.originalfiles(paste("Radiance_{day}_{time}_", reflectivity.default[-1], "_solar1.csv", sep='') %>% gsub("\\{day\\}", day, .) 
                                            %>% gsub("\\{time\\}", time, .),  originalfolders.training)
training.solar2.files <- find.originalfiles(paste("Radiance_{day}_{time}_", reflectivity.default[-1], "_solar2.csv", sep='') %>% gsub("\\{day\\}", day, .) 
                                            %>% gsub("\\{time\\}", time, .),  originalfolders.training)
training.total.files <- find.originalfiles(paste("Radiance_{day}_{time}_", reflectivity.default[-1], "_total.csv", sep='') %>% gsub("\\{day\\}", day, .) 
                                           %>% gsub("\\{time\\}", time, .),  originalfolders.training)
training.trans <- read.originalmatrix(training.trans.files, wv, an)
training.down <- read.originalmatrix(training.down.files, wv, an)
training.up1 <- read.originalmatrix(training.up1.files, wv, an)
training.up2 <- read.originalmatrix(training.up2.files, wv, an)
training.solar1 <- read.originalmatrix(training.solar1.files, wv, an)
training.solar2 <- read.originalmatrix(training.solar2.files, wv, an)
training.total <- read.originalmatrix(training.total.files, wv, an)
# 1 dimension: wavelength ; 2 dimension: angle; 3 dimension: file
# dim(statistic.prep(training.trans)) # dimensions: 150, 13, 6
# calculate the maximum value at each wavelength
training.trans.max <- apply(statistic.prep(training.trans),1, max)
training.down.max <- apply(statistic.prep(training.down),1, max)
training.up1.max <-apply(statistic.prep(training.up1),1, max)
training.up2.max <-apply(statistic.prep(training.up2),1, max)
training.solar1.max <-apply(statistic.prep(training.solar1),1, max)
training.solar2.max <-apply(statistic.prep(training.solar2),1, max)
training.total.max <-apply(statistic.prep(training.total),1, max)

#############################################################################################################################################################
outliers <- list()
data.read <- list()
materials <- basename(originalfolders.test)
for(i in 1:length(originalfolders.test)){
  originalfolder.test <- originalfolders.test[i]
  test.trans.file <- paste(originalfolder.test, "/Radiance_108_14_trans.csv", sep="")
  if(!file.exists(test.trans.file)){
    stop("transmission file not found")
  }else{
    test.trans<- read.originalmatrix(test.trans.file, wv, an)[[1]] #unnamed dataframe
    test.trans.max <- apply(test.trans, 1, max)
    index <- which((training.trans.max-test.trans.max)<0)
    trans.outlier <- cbind(training=training.trans.max[index], test=test.trans.max[index])
  }
  ### downwelling component
  test.down.file <- paste(originalfolder.test, "/Radiance_108_14_down.csv", sep="")
  if(!file.exists(test.down.file)){
    print("downwelling file not found (downwelling NA?)")
    print(paste(i, ":", originalfolders.test[i]))
    down.outlier <- NA
  }else{
    test.down<- read.originalmatrix(test.down.file, wv, an)[[1]]
    test.down.max <- apply(test.down, 1, max)
    index <- which((training.down.max -test.down.max)<0)
    down.outlier <- cbind(training=training.down.max[index], test=test.down.max[index])
  }
  ### up1 component
  test.up1.file <- paste(originalfolder.test, "/Radiance_108_14_up1.csv", sep="")
  if(!file.exists(test.up1.file)){
    stop("path thermal emission file (up1) not found")
  }else{
    test.up1<- read.originalmatrix(test.up1.file, wv, an)[[1]]
    test.up1.max <- apply(test.up1, 1, max)
    index <- which((training.up1.max -test.up1.max)<0)
    up1.outlier <- cbind(training=training.up1.max[index], test=test.up1.max[index])
  }
  ### up2 component
  test.up2.file <- paste(originalfolder.test, "/Radiance_108_14_up2.csv", sep="")
  if(!file.exists(test.up2.file)){
    stop("path thermal scattering file (up2) file not found")
  }else{
    test.up2<- read.originalmatrix(test.up2.file, wv, an)[[1]]
    test.up2.max <- apply(test.up2, 1, max)
    index <- which((training.up2.max - test.up2.max)<0)
    up2.outlier <- cbind(training=training.up2.max[index], test=test.up2.max[index])
  }
  ### solar1 component
  test.solar1.file <- paste(originalfolder.test, "/Radiance_108_14_solar1.csv", sep="")
  if(!file.exists(test.solar1.file)){
    stop("solar single scattering (solar1) file not found")
  }else{
    test.solar1<- read.originalmatrix(test.solar1.file, wv, an)[[1]]
    test.solar1.max <- apply(test.solar1, 1, max)
    index <- which((training.solar1.max -test.solar1.max)<0)
    solar1.outlier <- cbind(training=training.solar1.max[index], test=test.solar1.max[index])
  }
  ### solar2 component
  test.solar2.file <- paste(originalfolder.test, "/Radiance_108_14_solar2.csv", sep="")
  if(!file.exists(test.solar2.file)){
    stop("solar multi scattering (solar2) file file not found")
  }else{
    test.solar2<- read.originalmatrix(test.solar2.file, wv, an)[[1]]
    test.solar2.max <- apply(test.solar2, 1, max)
    index <- which((training.solar2.max -test.solar2.max)<0)
    solar2.outlier <- cbind(training=training.solar2.max[index], test=test.solar2.max[index])
  }
  ### total component
  test.total.file <- paste(originalfolder.test, "/Radiance_108_14_total.csv", sep="")
  if(!file.exists(test.total.file)){
    stop("total radiance file not found")
  }else{
    test.total <- read.originalmatrix(test.total.file, wv, an)[[1]]
    test.total.max <- apply(test.total, 1, max)
    index <- which((training.total.max -test.total.max)<0)
    total.outlier <- cbind(training=training.total.max[index], test=test.total.max[index])
  }
  
  data.read[[i]] <- list(trans=test.trans, down=test.down, up1=test.up1, up2=test.up2, solar1=test.solar1, solar2=test.solar2, total=test.total)
  outliers[[i]] <- list(trans=trans.outlier, down=down.outlier, up1=up1.outlier, up2=up2.outlier, solar1=solar1.outlier, solar2=solar2.outlier, total=total.outlier)
}
names(data.read)<-materials
names(outliers)<-materials


#############################################################################################################################################################
# let's calculate the number of outliers for each material and each component
stat.outlier.number <- matrix(nrow=length(outliers), ncol=7)
rownames(stat.outlier.number) <- materials
colnames(stat.outlier.number) <- c("trans", "down", "up1", "up2", "solar1", "solar2", "total")
for(i in 1:length(outliers)){
  stat.outlier.number[i,] <- unlist(lapply(outliers[[i]], function(x) if(is.null(dim(x))) NA else dim(x)[[1]]))
}

#############################################################################################################################################################
# let's calculate the difference by their magnitude compared to original signals
outliers.magnitude <- list()
for(i in 1:length(outliers)){
  temp <- lapply(outliers[[i]], function(x) if(is.null(dim(x))) NA else if(dim(x)[[1]]==0) NA else (x[,'test']-x[,'training'])/x[,'training']) # list
  outliers.magnitude[[i]] <- temp[which(!is.na(temp))]        # list
}
names(outliers.magnitude) <- materials
# given a confidence level, filter out the outliers
confidence.level <- 0.1
outliers.magnitude.conf <- lapply(outliers.magnitude, function(x) lapply(x, function(y) if(length(which(y >=confidence.level))==0) NA else y[which(y >=confidence.level)]) )
outliers.magnitude.conf <- lapply(outliers.magnitude.conf, function(x) x[!is.na(x)])
outliers.magnitude.conf <- outliers.magnitude.conf[sapply(outliers.magnitude.conf, function(x) length(x)!=0)]
# find all outliers material and component
temp <- unlist(outliers.magnitude.conf, recursive = FALSE )
material.component <- do.call(rbind, strsplit(names(temp), "\\."))
colnames(material.component) <- c("material", "component")
# find all outlier wavelengths
cols.wv <- sort(unique(as.numeric(unlist(sapply(temp, function(x) names(x))))))
value.wv <- matrix(nrow=nrow(material.component), ncol=length(cols.wv))
colnames(value.wv) <- cols.wv
for(i in 1:length(temp)){
  temp.wv <-as.numeric(names(temp[[i]]))
  temp.values <- unname(temp[[i]])
  value.wv[i,which(cols.wv %in% temp.wv)] <- temp.values
}
number.nna <- apply(value.wv, 1, function(x) length(which(!is.na(x))))
stat.outlier.material.comp.wv <- cbind(material.component, number=number.nna, value.wv)  


#############################################################################################################################################################
# check result
outliers; stat.outlier.number; stat.outlier.material.comp.wv;

data.read$LAMB_CLOUD_DECK$down['2.5525',]
training.down$Radiance_108_14_100_down['2.5525',]

data.read$LAMB_CLOUD_DECK$down['2.605',]
training.down$Radiance_108_14_100_down['2.605',]

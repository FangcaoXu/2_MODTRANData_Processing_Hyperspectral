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
###################### output file name represents for the local time
library(colorRamps)
library(data.table)
library(RColorBrewer)
library(scales)


source("Z_general_functions.R")
source("Z_MODTRAN_functions.R")
# Valid output in our MODTRAN simulations
#[1] "wavelength"         "path_trans"         "path_emiss"         "path_thermal_scat"  "surface_emiss"      "path_multi_scat"   
#[7] "sing_scat"          "grnd_rflt"          "drct_rflt"          "total_rad"          "irrad_ref"          "irrad_obs"         
#[13] "-natlog_path_trans" "direct_emiss"       "ToA_irrad"          "bbody_temp"         "theta"              "elevation"  

# writeMODTRAN2Matrix.singleset('~/MODTRAN6/try', "~/MODTRAN6/try/try1", writeinsame=TRUE, day=191, time=22, reflec=0, solar=TRUE)

################################################### Generate matrix output ###################################################
######################################################################################## Mid-latitude Spring and Summer Model
writeMODTRAN2Matrix('~/disk10TB/DARPA/MODTRANSimulated/MidLatitude/springsummer', "~/disk10TB/DARPA/MatrixCSV/MidLatitude")
######################################################################################## Mid-latitude Autumn and Winter Model
writeMODTRAN2Matrix('~/disk10TB/DARPA/MODTRANSimulated/MidLatitude/autumnwinter', "~/disk10TB/DARPA/MatrixCSV/MidLatitude")

# ######################################################################################## Tropical atmospheric model
# writeMODTRAN2Matrix('~/disk10TB/DARPA/MODTRANSimulated/Tropical', "~/disk10TB/DARPA/MatrixCSV/Tropical")
# ######################################################################################## US_Standard_1976 atmospheric model
# writeMODTRAN2Matrix('~/disk10TB/DARPA/MODTRANSimulated/US_Standard_1976', "~/disk10TB/DARPA/MatrixCSV/US_Standard_1976")
# ######################################################################################## Sub-Arctic Summer Model
# writeMODTRAN2Matrix('~/disk10TB/DARPA/MODTRANSimulated/Subarc/summer', "~/disk10TB/DARPA/MatrixCSV/SubArctic")
# ######################################################################################## Sub-Arctic Winter Model
# writeMODTRAN2Matrix('~/disk10TB/DARPA/MODTRANSimulated/Subarc/winter', "~/disk10TB/DARPA/MatrixCSV/SubArctic")

######################################################################################## Simulations for BH data 
writeMODTRAN2Matrix('~/disk10TB/DARPA/MODTRANSimulated/Stage2_BHData', '~/disk10TB/DARPA/MatrixCSV/Stage2_BHData/BHSimulatedMatrix', 
                    eles=seq(60,30,-5), filename.pattern=c("Geometry", "azimuth","Day","Time","Reflectivity"))

######################################################################################## Mid-latitude One-year Model VNIR/SWIR
writeMODTRAN2Matrix('~/disk10TB/DARPA/MODTRANSimulated/MidLatitude_VNIR_SWIR/OneYear_2hours', '~/disk10TB/DARPA/MatrixCSV/MidLatitude_VNIR_SWIR/OneYear_2hours', 
                    solar=TRUE, filename.pattern=c("Geometry", "azimuth","Day","Time","Reflectivity"))



#################################################################################################################################
############################################### Task 1: SWIR MODTRAN built-in materials ########################################
#################################################################################################################################
######################################################################################## One Year VNIR/SWIR MODTRAN materials
folders.modtran <- list.files('~/disk10TB/DARPA/MODTRANSimulated/MidLatitude_VNIR_SWIR/DifferentMaterials/OneYear_2hours_MODTRANBuiltin_Materials', full.names = TRUE)
materials.MODTRAN <- basename(folders.modtran)
# 365*8*13 37960
# extract reflectivity for different materials 
for(i in length(materials.MODTRAN)){
  print(materials.MODTRAN[i])
  temp <- list.files(folders.modtran[i], pattern='Radiance0_108_10_scan.csv', full.names=TRUE)
  temp <- head(read.csv(temp, skip=6, header = FALSE, stringsAsFactors = FALSE),-1)
  temp <- setNames(temp, c("wavelength", "path_trans", "path_emiss","path_thermal_scat","surface_emiss", "path_multi_scat",
              "sing_scat","grnd_rflt","drct_rflt","total_rad","irrad_ref","irrad_obs","-natlog_path_trans","direct_emiss",
              "ToA_irrad", "bbody_temp"))
  reflec <- cbind(as.numeric(temp$wavelength), 1-as.numeric(temp$direct_emiss))
  write.csv(reflec, paste('~/disk10TB/DARPA/MatrixCSV/MidLatitude_VNIR_SWIR/DifferentMaterials/MODTRAN_OneYear/', materials.MODTRAN[i], "_Reflectivity.csv", sep=""))
}
# LWIRRadiance2Radiance(files.ithfolder) # change folder name from LWIRRadiance to Radiance
############################################################
for(i in 1:length(folders.modtran)){
  if(length(list.files(folders.modtran[i], pattern='scan.csv', full.names=TRUE))==37960){
    writeMODTRAN2Matrix(folders.modtran[i],
                        paste('~/disk10TB/DARPA/MatrixCSV/MidLatitude_VNIR_SWIR/DifferentMaterials/MODTRAN_OneYear/', materials.MODTRAN[i], sep=""), 
                        solar=TRUE, filename.pattern=c("Geometry", "Day","Time"), downNA = TRUE)
  }
}

#################################################################################################################################
####################################### Task 2: BH LWIR MODTRAN Built-in Materials ##############################################
#################################################################################################################################
# read elevation angle
geofile <- read.csv('/amethyst/s0/fbx5002/PythonWorkingDir/sensorLocation/sensorAngleEle.csv', row.names = 1)
# change zenith to elevation angles
geofile$theta <- 90 - geofile$theta 

folders.modtran <- list.files('~/disk10TB/DARPA/MODTRANSimulated/Stage2_BHData/DifferentMaterials', full.names = TRUE)
materials.MODTRAN <- basename(folders.modtran)
# 1. 26784 scan.csv for each material: 6 temperature* 4 times of a day* 1116 geometry (31 theta, 36 ranges)
done <- NULL
for(i in 1:length(materials.MODTRAN)){
  print(i)
  scan.len <- length(list.files(folders.modtran[i], pattern='scan.csv', full.names=TRUE))
  print(paste(materials.MODTRAN[i], "has", scan.len, "files"))
  if(scan.len == 26784) done<-c(done, i)
}
done; materials.MODTRAN[done]
# 2. write all data into matrix for each material
for(i in 1:length(materials.MODTRAN)){
  writeMODTRAN2Matrix(folders.modtran[i], paste('~/disk10TB/DARPA/MatrixCSV/Stage2_BHData/MODTRAN/DifferentMaterials/', materials.MODTRAN[i], sep=""), 
                      filename.pattern=c("Geometry", "Day","Time","Temperature"), 
                      elevation.angles=geofile, downNA = TRUE, bh=TRUE)
  #free up RAM
  gc() 
}
# 3. check whether all materials are writen into matrix csv
donefolders <- list.files('~/disk10TB/DARPA/MatrixCSV/Stage2_BHData/MODTRAN/DifferentMaterials', full.names = T)
for (i in 1:length(materials.MODTRAN)){
  comp <- list.files(donefolders[i], full.names = T)
  print(paste(i, ", All simulations are done: ", all((sapply(comp, function(x) length(list.files(x))) == 24)[-c(5,6)]), 
              " for folder ", basename(donefolders[i]), sep = "")) # solar are not included 
}
# 4. generate reflectivity of each material
for(i in 1:length(materials.MODTRAN)){
  print(materials.MODTRAN[i])
  temp <- list.files(folders.modtran[i], pattern='Radiance0_107_14_295_scan.csv', full.names=TRUE)
  temp <- head(read.csv(temp, skip=6, header = FALSE, stringsAsFactors = FALSE),-1)
  temp <- setNames(temp, c("wavelength", "path_trans", "path_emiss","path_thermal_scat","surface_emiss", "path_multi_scat",
                           "sing_scat","grnd_rflt","drct_rflt","total_rad","irrad_ref","irrad_obs","-natlog_path_trans","direct_emiss",
                           "ToA_irrad", "bbody_temp"))
  reflec <- cbind(as.numeric(temp$wavelength), 1-as.numeric(temp$direct_emiss))
  write.csv(reflec, paste('~/disk10TB/DARPA/MatrixCSV/Stage2_BHData/MODTRAN/DifferentMaterials/Reflectivity/', materials.MODTRAN[i], "_reflec.csv", sep=""))
}

# all 29 materials at different temperature from 270 to 320 by every 1K at one geometries
geofile <- geofile[which(geofile$theta == 45 & geofile$range==5000),]
folders.modtran <- list.files('~/disk10TB/DARPA/MODTRANSimulated/Stage2_BHData/DifferentMaterials/TemperatureChecking', full.names = TRUE)
materials.MODTRAN <- basename(folders.modtran)
for(i in 1:length(materials.MODTRAN)){
  writeMODTRAN2Matrix(folders.modtran[i], paste('~/disk10TB/DARPA/MatrixCSV/Stage2_BHData/MODTRAN/DifferentMaterials/TemperatureChecking/', materials.MODTRAN[i], sep=""), 
                      filename.pattern=c("Geometry", "Day","Time","Temperature"), 
                      elevation.angles=geofile, downNA = T, bh=T)
  #free up RAM
  gc() 
}



##################################################################################################################################
##################################################################################################################################
# calculate solar blackbody radiance for normalization
wv <- seq(0.4, 3.0075, 0.0175) #um
solar_bbradiance <- blackbody_radiance_calculation(wv, temp = 5778)
solar_bbradiance <- as.data.frame(cbind(wv, radiance=solar_bbradiance*10^-4)) # [W cm^-2 sr^-1 um^-1]
write.csv(solar_bbradiance, "shortwaveSolarbbRadiance.csv")

for(i in c(350, 400, 450, 500)){
  bbradiance <- blackbody_radiance_calculation(wv, temp = i)
  bbradiance <- as.data.frame(cbind(wv, radiance=bbradiance*10^-4)) # [W cm^-2 sr^-1 um^-1]
  write.csv(bbradiance, paste("bbradiance_",i,".csv", sep=""))
}





# ######################################### generate matrix output for a specific day and time ####################################
# ######################################################################################## Singleset Ammonia
# writeMODTRAN2Matrix.singleset('~/disk10TB/DARPA/MODTRANSimulated/DifferentMaterials/Ammonia', "~/disk10TB/DARPA/MatrixCSV/DifferentMaterials/Ammonia",writeinsame=TRUE)
# ######################################################################################## Singleset Acetone
# writeMODTRAN2Matrix.singleset('~/disk10TB/DARPA/MODTRANSimulated/DifferentMaterials/Acetone', "~/disk10TB/DARPA/MatrixCSV/DifferentMaterials/Acetone",writeinsame=TRUE)
# ######################################################################################## Singleset exoscan
# folders.exoscan <- list.files('~/disk10TB/DARPA/MODTRANSimulated/Stage1_7.5_12um/DifferentMaterials/exoscan', full.names = TRUE)
# for(folder in folders.exoscan){
#   writeMODTRAN2Matrix.singleset(folder, paste('~/disk10TB/DARPA/MatrixCSV/Stage1_7.5_12um/DifferentMaterials/exoscan/originalCSV/', basename(folder), sep=""),writeinsame=TRUE)
# }
# 
# ######################################################################################## Singleset MODTRAN materials
# folders.modtran <- list.files('~/disk10TB/DARPA/MODTRANSimulated/MidLatitude_VNIR_SWIR/DifferentMaterials/MODTRAN', full.names = TRUE)
# folders.modtran.output <- paste('~/disk10TB/DARPA/MatrixCSV/MidLatitude_VNIR_SWIR/DifferentMaterials/MODTRAN/', basename(folders.modtran), sep='')
# sapply(folders.modtran.output, function(x) if(!dir.exists(x)) dir.create(x))
# for(i in 1:length(folders.modtran)){
#   writeMODTRAN2Matrix.singleset(folders.modtran[i], folders.modtran.output[i], writeinsame=TRUE, day=108, time=14, solar=TRUE)
# }

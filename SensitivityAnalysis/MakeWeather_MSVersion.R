#Create weather for Model Experiments w/ weather generator###
#packages
library(ggplot2)
library(geosphere)
setwd("") #set to location of terrestrial-aquatic-LCT folder

###create gradients####
pptGrads<-c(25, 50, 75, 100, 125, 150, 175, 200, 225, 250) #10 items
tmaxGrads<-c(12, 15, 18, 21, 24, 27) # 6 items

gradsDF<-data.frame(matrix(ncol= 6, nrow=60))
colnames(gradsDF)<-c("runNum", "Site", "Lat", "CorrT", "pptGrads", "tmaxGrads")
gradsDF[,1]<-1:nrow(gradsDF)
gradsDF[,2]<-rep("HighLat", times=60)
gradsDF[,3]<-rep("35", times=60)
gradsDF[,4]<-rep(1, times=60)
gradsDF[,5]<-rep(pptGrads, times=6)
gradsDF[,6]<-c(rep(tmaxGrads[1], times=10), rep(tmaxGrads[2], times=10), rep(tmaxGrads[3], times=10), rep(tmaxGrads[4], times=10), rep(tmaxGrads[5], times=10), rep(tmaxGrads[6], times=10))

####Make fake weather
#get the weather generator
source("Functions/genWeather_MSVersion.R")
source("Functions/EvaporationFunc_MSVersion.R")
dfOut<-data.frame(matrix(ncol=14, nrow=0))
colnames(dfOut)<-c("Month", "Day", "Precip", "AnnPpt", "MeanTemp", "DOY", "PAR", "VPD", "Evap", "SoilTemp", "runNum", "Site", "AnnPpt", "Tmax")
#run weather generator
for(i in 1:nrow(gradsDF)){
  dfOut1<-genWeather2(annPpt=gradsDF$pptGrads[i], fracWet=0.5, pptSeason=c(0.24166, -0.166),  tempSeason=c(0.2573, -0.2727), corrT=gradsDF$CorrT[i], Tmax=gradsDF$tmaxGrads[i], lat=gradsDF$Lat[i])
  
  dfOut1$runNum<-i
  dfOut1$Site<-gradsDF$Site[i]
  dfOut1$AnnPpt<-gradsDF$pptGrads[i]
  dfOut1$Tmax<-gradsDF$tmaxGrads[i]
  dfOut<-rbind(dfOut, dfOut1)
  
}

dfOut[1:(ncol(dfOut)-2)] <- lapply(dfOut[1:(ncol(dfOut)-2)], as.character) #make data character
dfOut[1:(ncol(dfOut)-2)] <- lapply(dfOut[1:(ncol(dfOut)-2)], as.numeric) #this will have a warning message, it's okay. It is telling you that the blanks were replaced with "NA"
dfOut[1:(ncol(dfOut)-2)]<-lapply(dfOut[1:(ncol(dfOut)-2)], round, digits=3) #round to the number of digits in the original data

#write.csv(dfOut, "SensitivityAnalysis/dfOut.csv, row.names=F)

#phenology (Don and Doff parameters)
source("Functions/gddSum_MSVersion.R")
phenolDF<-data.frame(matrix(ncol=3, nrow=nrow(gradsDF)))
colnames(phenolDF)<-c("runNum", "LeafOn", "LeafOff")

for(i in 1:nrow(gradsDF)){
  dfSub<-dfOut[dfOut$runNum==i,]
  phenolDF$runNum[i]<-i
  p<-gddSum(Tmean=as.numeric(as.character(dfSub$MeanTemp)), baseT_on=5, latit=as.numeric(as.character(gradsDF$Lat[i])))
  phenolDF$LeafOn[i]<-p[1]
  phenolDF$LeafOff[i]<-p[2]
}

#write.csv(phenolDF, "SensitivityAnalysis/PhenolDF.csv, row.names=F)
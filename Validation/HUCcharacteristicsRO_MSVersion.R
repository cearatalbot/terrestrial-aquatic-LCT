#####Get area, PFT, latitude, mean monthly temperature, and total annual precip for each HUC###
setwd("") #set to location of terrestrial-aquatic-LCT folder
source("Functions/genWeatherMOD_MSVersion.R") #source modified weather generator 
source("Functions/gddSum_MSVersion.R") #gddsum calculator

library(sirad)
library(geosphere)
#NOTES FOR METHODS: I removed catchments that extend beyond the US (only affected Region04; MI), I combined HAY and CROP land covers, will assume that mixed forest is 50% decid and 50% eg
#list of regions
ListRegions<-c("01", "02", "03W", "03S", "03N","04", "05", "06", 
               "07", "08", "09", "10L", "10U", "11", "12", "13", 
               "14", "15", "16", "17", "18")
ListMonths<-c("01", "02", "03", "04", "05", "06",
              "07", "08", "09", "10", "11", "12")
#loop to read in data files (Q: How many different files per hydroregion do I need? may want to read in and collate everything here)
df_PFT1<-data.frame(matrix(ncol=7, nrow=0)) #to store PFT for all
colnames(df_PFT1)<-c("Region", "con_weighted", "dec_weighted", "gr_weighted", "sh_weighted", "mx_weighted", "cr_weighted")
dfOut<-data.frame(matrix(ncol=14, nrow=0)) #df to store weather for all
colnames(dfOut)<-c("Month", "Day", "Precip", "AnnPpt", "MeanTemp", "DOY", "PAR", "VPD", "Evap", "SoilTemp", "runNum", "Site", "AnnPpt", "Tmax")
phenolDF<-data.frame(matrix(ncol=3, nrow=length(ListRegions))) #DF to store phenology
colnames(phenolDF)<-c("Region", "LeafOn", "LeafOff") #name columns
m2Region<-data.frame(matrix(ncol=2, nrow=21))
colnames(m2Region)<-c("Region", "m2area")
runoff<-data.frame(matrix(ncol=4, nrow=21))
colnames(runoff)<-c("Region", "RunoffMA", "Precip","RunoffRatio")#cm

for(i in 1:length(ListRegions)){
  df_NLCD<-read.csv(paste("streamCatDataTheory/NLCD2016_Region", ListRegions[i], ".csv", sep=""), stringsAsFactors=F) #read in NLCD
  df_NLCDx<-df_NLCD[order(df_NLCD$COMID),] #sort by comid
  df_PRECIP<-read.csv(paste("NHDdata/NHDPlusHydroRegion/NHDPlus", ListRegions[i], "/VPUAttributeExtension/IncrPrecipMA.txt", sep=""), stringsAsFactors=F) #precip=mm*100
  df_PRECIPx<-df_PRECIP[order(df_PRECIP$FeatureID),] #sort by feature id
  df_LAT<-read.csv(paste("NHDdata/NHDPlusHydroRegion/NHDPlus", ListRegions[i], "/VPUAttributeExtension/IncrLat.txt", sep=""), stringsAsFactors=F) #read in latitude
  df_LATx<-df_LAT[order(df_LAT$FeatureID),]
  df_TEMPx<-data.frame(matrix(ncol=13, nrow=nrow(df_LAT)))
  df_ROMA<-read.csv(paste("NHDdata/NHDPlusHydroRegion/NHDPlus", ListRegions[i], "/VPUAttributeExtension/ROMA.txt", sep=""), stringsAsFactors=F) #mean annual runoff=mm yr-1
  df_ROMAx<-df_ROMA[order(df_ROMA$FeatureID),] #sort by feature id
  colnames(df_TEMPx)<-c("Temp01", "Temp02", "Temp03", "Temp04", "Temp05", "Temp06", "Temp07", "Temp08", "Temp09", "Temp10", "Temp11", "Temp12", "FeatureID")
  for(z in 1:12){
    df_TEMP<-read.csv(paste("NHDdata/NHDPlusHydroRegion/NHDPlus", ListRegions[i], "/VPUAttributeExtension/IncrTempMM", ListMonths[z], ".txt", sep=""), stringsAsFactors=F) #read in PRISM temp for each month
    df_TEMP<-df_TEMP[order(df_TEMP$FeatureID),]
    df_TEMPx[,z]<-df_TEMP$TempV
  }
  df_TEMPx$FeatureID<-df_TEMP$FeatureID
  #combine these
  if(i==6){
    df_matching<-data.frame(cbind(df_PRECIPx$FeatureID, df_PRECIPx$PrecipV, df_TEMPx[,1:12], df_LATx$LatV, df_ROMAx$RunOffV))
    colnames(df_matching)<-c("FeatureID", "PrecipV", colnames(df_TEMPx[,1:12]), "LatV", "RunoffV")
    df<-merge(df_NLCDx, df_matching, by.x="COMID", by.y="FeatureID", all.x=T)
  } else {
    df<-cbind(df_NLCDx, df_PRECIPx$PrecipV, df_TEMPx[,1:12], df_LATx$LatV, df_ROMAx$RunOffV) #only if cols match
    colnames(df)[38:52]<-c("PrecipV", colnames(df_TEMPx[,1:12]), "LatV", "RunoffV") #change last two colnames
  }     
  fiName<-paste("Region", ListRegions[i], sep="") #make table name
  assign(fiName, df) #name table for region
  write.csv(df, paste("combined3/", fiName, ".csv", sep=""), row.names=F)
  #get weighted PFT land use
  df$conCat<-df$CatAreaSqKm*df$PctConif2016Cat #sqkm of each PFT in each sub catchment
  df$decCat<-df$CatAreaSqKm*df$PctDecid2016Cat
  df$grCat<-df$CatAreaSqKm*df$PctGrs2016Cat
  df$shCat<-df$CatAreaSqKm*df$PctShrb2016Cat
  df$mxCat<-df$CatAreaSqKm*df$PctMxFst2016Cat
  df$crCat<-df$CatAreaSqKm*(df$PctCrop2016Cat+df$PctHay2016Cat)
  df_PFT<-data.frame(matrix(ncol=7, nrow=1))
  colnames(df_PFT)<-c("Region", "con_weighted", "dec_weighted", "gr_weighted", "sh_weighted", "mx_weighted", "cr_weighted")
  df_PFT$Region<-paste("Region", ListRegions[i])
  df_PFT$con_weighted<-sum(df$conCat, na.rm=T)/sum(df$CatAreaSqKm, na.rm=T) #weighted average pct PFT
  df_PFT$dec_weighted<-sum(df$decCat, na.rm=T)/sum(df$CatAreaSqKm, na.rm=T)
  df_PFT$gr_weighted<-sum(df$grCat, na.rm=T)/sum(df$CatAreaSqKm, na.rm=T)
  df_PFT$sh_weighted<-sum(df$shCat, na.rm=T)/sum(df$CatAreaSqKm, na.rm=T)
  df_PFT$mx_weighted<-sum(df$mxCat, na.rm=T)/sum(df$CatAreaSqKm, na.rm=T)
  df_PFT$cr_weighted<-sum(df$crCat, na.rm=T)/sum(df$CatAreaSqKm, na.rm=T)
  df_PFT1<-rbind(df_PFT1, df_PFT) #combine for all regions
  #Generate weather
  #run weather generator
  Tmaxdf<-data.frame(matrix(ncol=12, nrow=1)) #for weather generator
  colnames(Tmaxdf)<-c("Tmax01", "Tmax02", "Tmax03", "Tmax04", "Tmax05", "Tmax06", 
                      "Tmax07", "Tmax08", "Tmax09", "Tmax10", "Tmax11", "Tmax12")
  for(z in 1:12){
    df[,(52+z)]<-df$CatAreaSqKm*(df[,(38+z)]/100) #area weighted temp
    Tmaxdf[1,z]<-sum(df[,(51+z)])/sum(df$CatAreaSqKm)
  }
  colnames(df)[53:64]<-c("Tmax01", "Tmax02", "Tmax03", "Tmax04", "Tmax05", "Tmax06", 
                         "Tmax07", "Tmax08", "Tmax09", "Tmax10", "Tmax11", "Tmax12")
  df$CatPrecip<-df$CatAreaSqKm*(df$PrecipV/1000)#area weighted precip w/ conversion from mm to cm
  df$CatLat<-df$CatAreaSqKm*(df$LatV) #area weighted latitude
  
  df$CatRunoff<-df$CatAreaSqKm*(df$RunoffV/10)#convert from mm to cm
  runoff$Region[i]<-paste("Region", ListRegions[i], sep="")
  runoff$RunoffMA[i]<-sum(df$CatRunoff)/sum(df$CatAreaSqKm)#in cm
  runoff$Precip[i]<-sum(df$CatPrecip)/sum(df$CatAreaSqKm) #in cm
  runoff$RunoffRatio[i]<-runoff$RunoffMA[i]/runoff$Precip[i]
    
  dfOut1<-genWeather2(annPpt=(sum(df$CatPrecip)/sum(df$CatAreaSqKm)), fracWet=0.5, pptSeason=c(0.24166, -0.166),  tempSeason=c(0.2573, -0.2727), corrT=1, Tmaxdf=Tmaxdf, lat=(sum(df$CatLat)/sum(df$CatAreaSqKm)))
  dfOut1$runNum<-i
  dfOut1$Site<-paste("Region", ListRegions[i], sep="")
  dfOut1$AnnPpt<-sum(df$CatPrecip)/sum(df$CatAreaSqKm)
  dfOut1$Tmax<-Tmaxdf$Tmax07[1]
  dfOut<-rbind(dfOut, dfOut1)
  
  dfOut[1:(ncol(dfOut)-2)] <- lapply(dfOut[1:(ncol(dfOut)-2)], as.character) #make data character
  dfOut[1:(ncol(dfOut)-2)] <- lapply(dfOut[1:(ncol(dfOut)-2)], as.numeric) #this will have a warning message, it's okay. It is telling you that the blanks were replaced with "NA"
  dfOut[1:(ncol(dfOut)-2)]<-lapply(dfOut[1:(ncol(dfOut)-2)], round, digits=3) #round to the number of digits in the original data
  #calculate phenology
  
  p<-gddSum(Tmean=as.numeric(as.character(dfOut1$MeanTemp)), baseT_on=5, latit=sum(df$CatLat)/sum(df$CatAreaSqKm)) #Tthres=100
  phenolDF$Region[i]<-paste("Region", ListRegions[i], sep="")
  phenolDF$LeafOn[i]<-p[1]
  phenolDF$LeafOff[i]<-p[2]
  
  #create a table with total area of each catchment
  m2Region$m2area[i]<-sum(df$CatAreaSqKm)*1000000
  m2Region$Region[i]<-paste("Region", ListRegions[i], sep="")
  
} #end of i loop; now we have 1 year of weather for each region :) 

#write.csv(dfOut, "Validation/regionsWeather.csv", row.names=F)#write out weather to keep
#write.csv(phenolDF, "Validation/regionsPhenology.csv", row.names=F)#write out phenology to keep
#write.csv(df_PFT1, "Validation/regionsPFT.csv", row.names=F)#write out PFT to keep
#write.csv(m2Region, "Validation/regionsArea.csv", row.names=F)#write out areas to keep
#write.csv(runoff, "Validation/regionsRunoff.csv", row.names=F)#write out areas to keep

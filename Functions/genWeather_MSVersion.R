######Function for generating weather data####
#
#load packages
library(sirad)
#daily precip, mean monthly temperature, total daily PAR, mean VPD, aquatic evap

#usage: genWeather(lon, lat, annPpt, fracWet, pptSeason(b,c), corrT, Tmax, tempSeason(b,c)   ..  ))
genWeather2<-function(annPpt, fracWet, pptSeason, corrT, Tmax, lat, tempSeason, corrV){
  tempDF<-data.frame(matrix(ncol=4, nrow = 12)) #new df for monthly info
  colnames(tempDF)<-c("Month", "ndays", "monthMeanT", "monthMeanTScaled")
  tempDF[,2]<-c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
  #generate monthly T
  YearTemp<-data.frame(matrix(ncol=4, nrow=0)) #make 4 for loop
  colnames(YearTemp)<-c("Month", "Day", "MeanTemp", "meanAppt") #add "fracwet" for loop
  
  for(t in 1:12){
    for(z in 1:12){
      tempDF[z,1]<- paste(z)
      monthMeanT<-(1*sin(tempSeason[1]*z + tempSeason[2])-0.2)*Tmax #seasonality
      tempDF[z,3]<-monthMeanT
    }
    tempDF[t,4]<-(tempDF[t,3]/max(tempDF[,3]))*Tmax 
    #monthDays<-tempDF[t,2]   # days in a month
    newTemp<-data.frame(matrix(ncol=4, nrow=tempDF$ndays[t])) #make 4 col for loop
    colnames(newTemp)<-c("Month", "Day", "MeanTemp", "meanAppt") #add "fracwet" for loop
    newTemp[,1]<-rep(t)
    newTemp[,2]<-seq(1,tempDF$ndays[t], by=1)
    newTemp[,3]<-rep(tempDF[t,4], tempDF$ndays[t])
    newTemp[,4]<-paste(annPpt)
    YearTemp<-rbind(YearTemp, newTemp)
  }
  
  ###soil surface temperature###
  #generate monthly T
  soilDF<-data.frame(matrix(ncol=3, nrow = 12)) #new df for monthly info
  colnames(soilDF)<-c("Month", "ndays", "monthMeanT")
  soilDF[,2]<-c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
  
  SoilTemp<-data.frame(matrix(ncol=4, nrow=0)) #make 4 for loop
  colnames(SoilTemp)<-c("Month", "Day", "SoilTemp", "meanAppt") #add "fracwet" for loop
  
  for (i in 1:12){
    soilDF[i,1]<- paste(i)
    soilDF[i,3]<-ifelse(tempDF[i,4] < 0, 0, tempDF[i,4])
    soilDF[i,4]<-paste(annPpt)
    SoilTemp1<-data.frame(matrix(ncol=4, nrow=soilDF$ndays[i])) #make 4 col for loop
    colnames(SoilTemp1)<-c("Month", "Day", "SoilTemp", "meanAppt") #add "fracwet" for loop
    SoilTemp1[,1]<-rep(i)
    SoilTemp1[,2]<-seq(1,soilDF$ndays[i], by=1)
    SoilTemp1[,3]<-rep(soilDF[i,3], soilDF$ndays[i])
    SoilTemp1[,4]<-paste(annPpt)
    SoilTemp<-rbind(SoilTemp, SoilTemp1)
  }
  
  #####precipitation#####
  monthDF<-data.frame(matrix(ncol=4, nrow = 12)) #new df for monthly info
  colnames(monthDF)<-c("Month", "ndays", "monthTotal", "monthTotScaled")
  monthDF[,2]<-c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
  
  YearPrecip<-data.frame(matrix(ncol=4, nrow=0)) #make 4 for loop
  colnames(YearPrecip)<-c("Month", "Day", "Precip", "meanAppt") #add "fracwet" for loop
  
  #loop to make monthly proportions from quadratic fit
  for (i in 1:12){
    monthDF[,1]<- paste(i)
    for(b in 1:12){
      monthTotal<-if(corrT < 0){
        1+(-1*sin(pptSeason[1]*b+pptSeason[2])-(corrT*sin(pptSeason[1]*b+pptSeason[2])))
      } else if (corrT > 0){
        1+corrT*sin(pptSeason[1]*b+pptSeason[2])
      } else if (corrT==1){
        corrT*sin(pptSeason[1]*b+pptSeason[2])
      }
      monthDF[b,3]<-monthTotal
    }
    monthDF[i,4]<-(monthDF[i,3]/sum(monthDF[,3]))*annPpt
    monthDays<-monthDF[i,2]   # days in a month
    wetDays<-floor(monthDays*fracWet)  # number of days with rain in the month
    monthMean<-monthDF[i,4]/wetDays  # mean amount of rain on days that it rained in the month
    meanlog<-log(monthMean)-(1^2)/2 # transforming mean for parameter rlnorm() needs
    events<-rlnorm(wetDays,meanlog,1) # random draws from lognormal with specified parameters
    events<-(events/sum(events))*monthDF[i,4] # scale random draws to month total
    precip<-sample(x=c(events,rep(0,(monthDays-wetDays))),size=monthDays,replace=FALSE)  # add in no rain days and randomly order them
    newPrecip<-data.frame(matrix(ncol=4, nrow=monthDays)) #make 4 col for loop
    colnames(newPrecip)<-c("Month", "Day", "Precip", "meanAppt") #add "fracwet" for loop
    newPrecip[,1]<-rep(i)
    newPrecip[,2]<-seq(1,monthDays, by=1)
    newPrecip[,3]<-precip
    newPrecip[,4]<-paste(annPpt)
    YearPrecip<-rbind(YearPrecip, newPrecip)
  }
  
  ######PAR######
  #convert decimal degrees to radians
  pi<-3.14159265 #pi
  lat<-as.numeric(lat) #latitude
  lat_radian<-lat * pi/180   #decimal deg to radians
  
  dfPar<-data.frame(matrix(ncol=3, nrow=365)) #empty df
  colnames(dfPar)<-c("lat", "DOY", "ETrad_MJ")
  dfPar$DOY<-seq(1, 365, by=1)
  dfPar$lat<-lat
  #get ET radiation based on latitude
  n<-extrat(dfPar$DOY, lat_radian)
  dfPar$ETrad_MJ<-n$ExtraTerrestrialSolarRadiationDaily #units are MJm-2
  #convert units from MJm-2 to Em-2
  dfPar$ETrad_Em2<-dfPar$ETrad_MJ*4.6 #full conversion written out elsewhere
  #convert ET radiation to PAR
  dfPar$PARest<-dfPar$ETrad_Em2*0.45 #source is Britton & Dodd 1993
  #for evap.. convert ETrad to W/m2
  dfPar$ETrad_W<-dfPar$ETrad_MJ/0.0864 
  #Introduce variation from cloudcover? 
  for(p in 1:365){
    dfPar$clouds[p]<-ifelse(YearPrecip$Precip[p] > 0, 0.6, 1) #could also use random numbers from normal distribution
    dfPar$PAR_clouds<-dfPar$PARest*dfPar$clouds
  }
  
  YearPAR<-data.frame(cbind(lat=dfPar$lat, DOY=dfPar$DOY, PAR=dfPar$PAR_clouds, clouds=dfPar$clouds, ETrad_W=dfPar$ETrad_W))

  #source for VPD equations... http://www.fao.org/3/X0490E/x0490e07.htm
  #vp = RH/100 * ((SVPmax_kPa + SVPmin_kPa)/2)
  dfVPD<-data.frame(matrix(ncol=9, nrow=0))
  colnames(dfVPD)<-c("Month", "Day", "DOY", "VPDmin_kPa", "VPDmax_kPa", "VPD_meanKpa", "SVP_meanKpa" , "VP_meanKpa", "aqEvap")
  
  #monthDF
  VPDmonth<-data.frame(matrix(ncol=5, nrow = 12)) #new df for monthly info
  colnames(VPDmonth)<-c("Month", "ndays", "VPDmin_kPa", "VPDmax_kPa", "VPDmean_kPa")
  VPDmonth[,2]<-c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
  
  VPDmin_kPa<-(4.6683+0.595*mean(tempDF[,4])+0.03974*annPpt-0.30026*max(tempDF[,4])-0.00355*mean(tempDF[,4])*annPpt+0.007778*mean(tempDF[,4])*max(tempDF[,4])-0.002819*annPpt*max(tempDF[,4])+0.00015*mean(tempDF[,4])*annPpt*max(tempDF[,4]))/10
  VPDmax_kPa<-(-0.1826+0.28864*mean(tempDF[,4])+0.15539*annPpt+1.4749*max(tempDF[,4])+0.00703*mean(tempDF[,4])*annPpt-0.0077*mean(tempDF[,4])*max(tempDF[,4])-0.0172*annPpt*max(tempDF[,4])+0.000111*mean(tempDF[,4])*annPpt*max(tempDF[,4]))/10
  
  for(x in 1:12){
    for(i in 1:12){
      VPDmonth[i,1]<- paste(i)
      monthminVPD<-(1*sin(tempSeason[1]*i + tempSeason[2]))*VPDmin_kPa#seasonality
      VPDmonth[i,3]<-monthminVPD
      monthmaxVPD<-(1*sin(tempSeason[1]*i + tempSeason[2]))*VPDmax_kPa #seasonality
      VPDmonth[i,4]<-monthmaxVPD
      VPDmonth[i,5]<-max(monthminVPD+sin(tempSeason[1]*i + tempSeason[2]),0) #begin with monthminVPD+
    }
    dfVPD1<-data.frame(matrix(ncol=9, nrow=VPDmonth$ndays[x]))
    colnames(dfVPD1)<-c("Month", "Day", "DOY", "VPDmin_kPa", "VPDmax_kPa", "VPD_meanKpa", "SVP_meanKpa" , "VP_meanKpa", "aqEvap")
    dfVPD1[,1]<-rep(x)
    dfVPD1[,2]<-seq(1,VPDmonth$ndays[x], by=1)
    dfVPD1[,4]<-rep(VPDmonth[x,3], VPDmonth$ndays[x])
    dfVPD1[,5]<-rep(VPDmonth[x,4], VPDmonth$ndays[x])
    dfVPD1[,6]<-rep(VPDmonth[x,5], VPDmonth$ndays[x])
    dfVPD1[,7]<-rep(0.6108*exp((17.27*tempDF[x,4])/(tempDF[x,4]+273.3)), VPDmonth$ndays[x]) 
    dfVPD<-rbind(dfVPD, dfVPD1)
  }
  dfVPD$DOY<-seq(1, 365, by=1)
  dfVPD$VP_meanKpa<-ifelse(dfVPD$SVP_meanKpa-dfVPD$VPD_meanKpa >= 0, dfVPD$SVP_meanKpa-dfVPD$VPD_meanKpa, 0)
  #loop evap calculation over each DOY
  #evap in mm day^-1
  for(b in 1:365){
    n<-evap_calc1(tmean=YearTemp$MeanTemp[b],svp=dfVPD$SVP_meanKpa[b], 
                  vp=dfVPD$VP_meanKpa[b], vpd=dfVPD$VPD_meanKpa[b], ETrad=YearPAR$ETrad_W[b])
    dfVPD$aqEvap[b]<-n}
  dfVPD<-as.data.frame(dfVPD)
  
  ########Combine weather#####
  df<-as.data.frame(cbind(YearPrecip$Month, YearPrecip$Day, YearPrecip$Precip, YearPrecip$meanAppt , YearTemp$MeanTemp, YearPAR$DOY, YearPAR$PAR, dfVPD$VPD_meanKpa, dfVPD$aqEvap, SoilTemp$SoilTemp))
  colnames(df)<-c("Month", "Day", "Precip", "AnnPpt", "MeanTemp", "DOY", "PAR", "VPD", "Evap", "SoilTemp")
  
  return(df)}
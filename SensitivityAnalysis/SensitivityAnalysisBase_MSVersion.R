###Run simulations across Temp/Precip gradients and base case for sensitivity analysis#####
#load packages
library(ggplot2)
library(nord)
library(geosphere)
library(deSolve)
#set working directory/source required functions
setwd("") #set to location of terrestrial-aquatic-LCT folder
source("Functions/TAM_MSVersion.R") #model code
source("Functions/summaryTAM_MSVersion.R") #summary function
source("Functions/findInitials_MSVersion.R") #initial conditions finder
dfOut<-read.csv("SensitivityAnalysis/dfOut.csv", stringsAsFactors = F) #used to be "WeatherData.csv"
phenolDF<-read.csv("SensitivityAnalysis/phenolDF.csv", stringsAsFactors = F)
trait_i0<-read.csv("SensitivityAnalysis/SASoil_trait_i0.csv", stringsAsFactors=F)
trait_i0_cp<-read.csv("SensitivityAnalysis/SASoil_trait_i0_cp.csv", stringsAsFactors=F)

###create gradients####
pptGrads<-c(25, 50, 75, 100, 125, 150, 175, 200, 225, 250) #10 items
tmaxGrads<-c(12, 15, 18, 21, 24, 27) # 6 items
gradsDF<-data.frame(matrix(ncol= 6, nrow=60))
colnames(gradsDF)<-c("runNum", "Site", "Lat", "CorrT", "pptGrads", "tmaxGrads")
gradsDF[,1]<-1:nrow(gradsDF)
gradsDF[,2]<-rep("HighLat", times=60)
gradsDF[,3]<-rep("35", times=60) #static latitude
gradsDF[,4]<-rep(1, times=60)
gradsDF[,5]<-rep(pptGrads, times=6)
gradsDF[,6]<-c(rep(tmaxGrads[1], times=10), rep(tmaxGrads[2], times=10), rep(tmaxGrads[3], times=10), rep(tmaxGrads[4], times=10), rep(tmaxGrads[5], times=10), rep(tmaxGrads[6], times=10))
    #####Set up for initial conditions####
    #trait_i0<-data.frame(matrix(nrow=60, ncol=1))
    #colnames(trait_i0)<-c("Val1")
    #Init<-c(4859.803, 2930.978, 7331.52, 10601.826, 14698.041, 17090.070,
            #16997.888, 19900.974 , 19641.915, 21889.722, -6844.758, 2145.273,
           # 7669.535, 11953.119, 15112.136, 18668.272, 22155.090, 23326.660,
            #24538.827, 25534.967, -9508.868, 1215.833, 6253.837, 11422.513,
            #15222.179, 19993.454, 21434.238, 23670.975, 27172.354, 29712.514,
            #-12453.754, -272.375, 4759.756, 8003.553, 13917.583, 17901.010,
            #22620.222, 25220.809, 26901.464, 28353.966, -16455.070, -4020.001,
           # 3278.430, 8062.284, 12205.615, 14034.055, 19171.822, 21048.932,
           # 25766.577, 29487.398, -21774.411, -9051.321, 1180.587, 5136.328,
            #7693.662, 11560.745, 16034.825, 17903.226, 21150.429, 23826.311)
    #trait_i0[,1]<-Init 
    #intital conditions for cp
      #trait_i0_cp<-data.frame(matrix(nrow=60, ncol=1))
      #colnames(trait_i0_cp)<-c("Val1")
      #Init<-5000
      #trait_i0_cp[,1]<-rep(Init, times=60)
    stored_PFT<-data.frame(matrix(ncol=67, nrow=(365*length(gradsDF)*2)))#DF to store results from looped sims
    colnames(stored_PFT)<-c("time", "Cw", "Cl", "Cs1", "Cs2",
                            "Cs3", "Cs4", "Cdoc1", "Cdoc2", "W1",
                            "W2", "Ca", "Cr", "Ccwd", "GPP", "Q1",         
                            "Q2", "Rf", "Ra", "NPP", "LCT1", "Rs3",         
                            "Dwater", "Wa", "TpotAreal" , "i0", "T",           
                            "Dvpd", "Dtemp", "Dlightbar", "GPPmax", "LAIareal",   
                            "GPPpotAreal", "WUE", "V", "Lw", "Ll", "L", "Rs1",
                            "Qout", "RfAreal", "Bdoc", "Q12", "LCT2", "Lr",          
                            "alCr", "lout", "Lg", "Bs1", "Ds1", "Ds3", "Ds4",        
                            "Ls3", "Ls1", "Rhdoc", "Rhdoc2", "Bs3", "Rs2",         
                            "Ls2", "Ds2", "Lf1", "Rm", "Bs2", "runNum",      
                            "pptGrads", "tmaxGrads", "trait")
      ##to make a loop for params
      newsim<-data.frame(matrix(ncol=66, nrow=nrow(gradsDF)))
      colnames(newsim)<-c("time", "Cw", "Cl", "Cs1", "Cs2",
                          "Cs3", "Cs4", "Cdoc1", "Cdoc2", "W1",
                          "W2", "Ca", "Cr", "Ccwd", "GPP", "Q1",         
                          "Q2", "Rf", "Ra", "NPP", "LCT1", "Rs3",         
                          "Dwater", "Wa", "TpotAreal" , "i0", "T",           
                          "Dvpd", "Dtemp", "Dlightbar", "GPPmax", "LAIareal",   
                          "GPPpotAreal", "WUE", "V", "Lw", "Ll", "L", "Rs1",
                          "Qout", "RfAreal", "Bdoc", "Q12", "LCT2", "Lr",          
                          "alCr", "lout", "Lg", "Bs1", "Ds1", "Ds3", "Ds4",        
                          "Ls3", "Ls1", "Rhdoc", "Rhdoc2", "Bs3", "Rs2",         
                          "Ls2", "Ds2", "Lf1", "Rm", "Bs2", "runNum",      
                          "pptGrads", "tmaxGrads")
      
      egON="NO" #this is for seasonal or EG leaves
      for(i in 1:60){
        repeat{ S0<-c(Cw=trait_i0[i,1],Cl=0,Cs1=5,Cs2=6,Cs3=300, Cs4=trait_i0_cp[i,1], Cdoc1=10, Cdoc2=7, W1=8, W2=5, Ca=500000, Cr=400, Ccwd=20)
       
         params=c(Amax=112, Ad=0.75, Kf=0.1, Tmin=4, Tmax=40, Topt=24,
                 Kvpd=0.05, PARhalf=17, k=0.58, Don=phenolDF$LeafOn[i],
                 Doff=phenolDF$LeafOff[i], Lmax=3, Ka=0.006, Q10v=2,
                 Kh=0.03, Q10s=2, f=0.04, Kwue=10.9, Wmax1=10, SLW=70,
                 Cfrac=0.45, Kw=0.005, Aa= 1e5, Ac= 2e5, zbar= 3, 
                 Cprecip = 1, deltaA = 0.01, deltaS1 = 0.37, Ks=1.2,
                 lambdaS1 = 0.40, lambdaS3= 0.0001, bi = 0.3, Tstar = 50, W20 = 5, 
                 Kr= 0.55, aw= 0.1, l = 1/365, lt = 0.5, r=1, Bp= 0.16, 
                 Wmax2=50, Kcwd=0.01, ag=0.1, rhoS1=0.4, rhoS2=0.5, rhoS3=0.55, 
                 deltaS3=0.0014, deltaS4=0.0001, deltaS2=0.0090, lambdaS2=0.04, fS1=0.6)
   
        subDF<-dfOut[dfOut$runNum==i, ]
        weatherDF<-subDF[rep(seq_len(nrow(subDF)), 200), ] #repeat one year of weather for 200 yrs
        weatherDF$runDay<-1:nrow(weatherDF)
        ### setting boundary conditions, parameters, and forcings
        # years to run the model for
        years2run=200
        times=seq(1,365*years2run)
        
        # define forcing approx functions
        PARapprox=approxfun(x=weatherDF$runDay, y = weatherDF$PAR) 
        Papprox=approxfun(x =weatherDF$runDay, y = weatherDF$Precip)  
        VPDapprox=approxfun(x=weatherDF$runDay, y = weatherDF$VPD)  
        Tair_approx=approxfun(x=weatherDF$runDay, y = weatherDF$MeanTemp)
        Tsoil_approx=approxfun(x=weatherDF$runDay, y = weatherDF$SoilTemp)
        Evap_approx=approxfun(x=weatherDF$runDay, y= weatherDF$Evap)
        Lapprox=approxfun(times,rep(c(rep(0,params[10]),(params[12]*params[20]*params[21]),rep(0,365-params[10]-1)),years2run))   
        Ll_approx=approxfun(times,rep(c(rep(0,params[11]),(params[12]*params[20]*params[21]),rep(0,365-params[11]-1)),years2run))  
        
        # runs simulation
        sim2=ode(y=S0,times=times,func=tamStep,parms=params,method="euler")
        sim2=data.frame(sim2)
        #label run with run number and gradient info
        sim2$runNum<-weatherDF$runNum[1]
        sim2$pptGrads<-weatherDF$AnnPpt[1]
        sim2$tmaxGrads<-weatherDF$Tmax[1]
            if(i==1){
              newsim[1:365,]<-sim2[(nrow(sim2)-364):nrow(sim2),]
            } else {
              newsim[(i*365-364):(i*365),]<-sim2[(nrow(sim2)-364):nrow(sim2),] 
            }
        iold_Cw=trait_i0[i,1]  #initial conditions finder for Cw
        iold_Cp=trait_i0_cp[i,1]  #initial conditions finder for Cs4
        inew_Cw<-findInitials(i0=iold_Cw, i1=sim2$Cw[nrow(sim2)])
        inew_Cp<-findInitials(i0=iold_Cp, i1=sim2$Cs4[nrow(sim2)])
        trait_i0[i,1]<-inew_Cw
        trait_i0_cp[i,1]<-inew_Cp
        
        if((inew_Cw < 0 | inew_Cw==iold_Cw)&(inew_Cp < 0 | inew_Cp==iold_Cp)){
          break
        } 
        }
        print(i)
        }
      stored_PFT[1:21900,]<-newsim
      stored_PFT$trait<-"Val1"
      
  #summarize
  subDF<-stored_PFT
  sumDF1<-data.frame(matrix(ncol=6, nrow=(60*37))) #store summary (37 vars) for each gradient (60) and each trait value
  colnames(sumDF1)<-c("runNum", "Variable", "Value", "Tmax", "Precip", "trait")
  for(x in 1:60){
    gradSum<-data.frame(matrix(ncol=6, nrow=37))
    colnames(gradSum)<-c("runNum", "Variable", "Value", "Tmax", "Precip")
    subDF1<-subDF[subDF$runNum==x,]
    gradSum[,1]<-x
    gradSum[,2:3]<-summaryTAM(df=subDF1) #should be subRegion
    gradSum[,4]<-subDF1$tmaxGrads[1]
    gradSum[,5]<-subDF1$pptGrads[1]
    gradSum[,6]<-"Val1" #label base
    if(x==1){
      sumDF1[1:37,]<-gradSum
    } else{
      sumDF1[((x*37)-36):(x*37),]<-gradSum
    }
  }#end x loop
  
  #make a two column df for initial conditions 
  #trait_i0_cp$Val3<-NA
  #trait_i0_cp$Val3<-trait_i0_cp$Val2
  #colnames(trait_i0_cp)<-c("Val2", "Val3")
  #trait_i0$Val3<-NA
  #trait_i0$Val3<-trait_i0$Val2
  #colnames(trait_i0)<-c("Val2", "Val3")
#write.csv(trait_i0, "SensitivityAnalysis/SASoil_trait_i0.csv", row.names=F) #write out initial conditions for Cw
#write.csv(trait_i0_cp, "SensitivityAnalysis/SASoil_trait_i0_cp.csv", row.names=F) #write out initial conditions for Cs4
#write.csv(sumDF1, "SensitivityAnalysis/SASoil_base.csv", row.names=F) #write file

############################################################
####RUNNING MODEL FOR PFTS and CLIMATES FOR EACH REGION#####
############################################################
library(deSolve)
setwd("") #set to location of terrestrial-aquatic-LCT folder
#Source functions
source("Functions/TAM_MSVersion.R") #model code
source("Functions/summaryTAM_MSVersion.R") #summary function
source("Functions/findInitials_MSVersion.R") #initial conditions finder
#read in files
dfWeather<-read.csv("Validation/regionsWeather.csv")#weather for each HUC
dfPhenol<-read.csv("Validation/regionsPhenology.csv")#phenology each HUC (based on weather)
dfPFT<-read.csv("Validation/regionsPFT.csv")#PFT converage each HUC
PFT_i0<-read.csv("Validation/PFT_i0.csv", stringsAsFactors=F)#wood initials
PFT_i0_cp<-read.csv("Validation/PFT_i0_cp.csv", stringsAsFactors=F)#CS4 initials
PFT_i0_cr<-read.csv("Validation/PFT_i0_cr.csv", stringsAsFactors=F)#Cr initials
PFT_i0_cm<-read.csv("Validation/PFT_i0_cm.csv", stringsAsFactors=F)#CS3 initials
PFT_i0_cf2<-read.csv("Validation/PFT_i0_cf2.csv", stringsAsFactors=F)#CS2 initials
PFT_i0_co<-read.csv("Validation/PFT_i0_co.csv", stringsAsFactors=F)#CS1 initials

ListRegions<-c("01", "02", "03W", "03S", "03N","04", "05", "06", 
               "07", "08", "09", "10L", "10U", "11", "12", "13", 
               "14", "15", "16", "17", "18")
pftNam<-c("EG", "DE", "GR", "SH") 

#PFT_i0_co<-data.frame(matrix(nrow=21, ncol=4))
#colnames(PFT_i0_co)<-c("EGi0", "DEi0", "GRi0", "SHi0")
#PFT_i0_co$EGi0<-15
#PFT_i0_co$DEi0<-15
#PFT_i0_co$GRi0<-15
#PFT_i0_co$SHi0<-15

stored_PFT<-data.frame(matrix(ncol=67, nrow=0))#DF to store results from looped sims
colnames(stored_PFT)<-c("time", "Cw", "Cl", "Cs1",
                        "Cs2", "Cs3", "Cs4", "Cdoc1",      
                        "Cdoc2", "W1", "W2", "Ca",         
                        "Cr", "Ccwd", "GPP", "Q1",         
                        "Q2", "Rf", "Ra", "NPP",        
                        "LCT1", "Rs3", "Dwater", "Wa",         
                        "TpotAreal", "i0", "T", "Dvpd",       
                        "Dtemp", "Dlightbar", "GPPmax", "LAIareal",   
                        "GPPpotAreal", "WUE", "V", "Lw",         
                        "Ll", "L", "Rs1", "Qout",       
                        "RfAreal", "Bdoc", "Q12", "LCT2",       
                        "Lr", "alCr", "lout", "Lg",          
                        "Bs1", "Ds1", "Ds3", "Ds4",         
                        "Ls3", "Ls1", "Rhdoc", "Rhdoc2",     
                        "Bs3", "Rs2",  "Ls2", "Ds2",        
                        "Lf1", "Rm", "Bs2", "Region",     
                        "pptGrads", "tmaxGrads", "PFT")
#21=number of HUCs
for(i in 1:21){
  #set up parameter sets for each PFT
  PFT_EG=c(Amax=60,Ad=0.75,Kf=0.1,Tmin=0,Tmax=40, Topt=20,
            Kvpd=0.21,PARhalf=17,k=0.5,Don=0,Doff=0,Lmax=2,
            Ka=0.006,Q10v=2,Kh=0.03,Q10s=2,f=0.04,Kwue=10.9,
            Wmax1=10,SLW=90,Cfrac=0.45,Kw=0.003,
            Aa= 1e5, Ac= 2e5,zbar= 3,Cprecip = 1,deltaA = 0.01, 
            deltaS1 = 0.37, lambdaS1 = 0.40, lambdaS3=0.0001, bi = 0.3,Tstar = 50, 
            W20 = 5, Kr= 0.55, aw= 0.1,l = 1/912.5, lt = 0.2, 
            r=1, Bp= 0.16, Wmax2=50, Kcwd=0.01, ag=0.1, Ks=1.2,
            rhoS1=0.4, rhoS2=0.5, rhoS3=0.55, deltaS3=0.0015, deltaS4=0.0001, 
            deltaS2=0.0090, lambdaS2=0.04, fS1=0.4) 
  PFT_DE=c(Amax=112, Ad=0.75, Kf=0.1, Tmin=4, Tmax=40, Topt=24,
           Kvpd=0.05, PARhalf=17, k=0.58, Don=dfPhenol$LeafOn[i],
           Doff=dfPhenol$LeafOff[i], Lmax=3, Ka=0.006, Q10v=2,
           Kh=0.03, Q10s=2, f=0.04, Kwue=10.9, Wmax1=10, SLW=70,
           Cfrac=0.45, Kw=0.005,  Aa= 1e5, 
           Ac= 2e5, zbar= 3, Cprecip = 1, deltaA = 0.01, deltaS1 = 0.37, Ks=1.2,
           lambdaS1 = 0.40, lambdaS3= 0.0001, bi = 0.3, Tstar = 50, W20 = 5, 
           Kr= 0.55, aw= 0.1, l = 1/365, lt = 0.5, r=1, 
           Bp= 0.16, Wmax2=50, Kcwd=0.01, ag=0.1, 
           rhoS1=0.4, rhoS2=0.5, rhoS3=0.55, deltaS3=0.0015, deltaS4=0.0001, 
           deltaS2=0.0090, lambdaS2=0.04, fS1=0.4)
  PFT_SH=c(Amax=66.5, Ad=0.75, Kf=0.1, Tmin=4, Tmax=40, Topt=24,
           Kvpd=0.2, PARhalf=17, k=0.58, Don=dfPhenol$LeafOn[i], 
           Doff=dfPhenol$LeafOff[i], Lmax=1, Ka=0.006, Q10v=2,
           Kh=0.03, Q10s=2, f=0.04, Kwue=6, Wmax1=10, SLW=100,
           Cfrac=0.45, Kw=0.03,  Aa= 1e5, 
           Ac=2e5, zbar= 3, Cprecip=1, deltaA=0.01, deltaS1=0.37, Ks=1.2,
           lambdaS1 = 0.40, lambdaS3= 0.0001, bi=0.3,Tstar=50, W20=5, Kr=0.35,
           aw=0.1, l = 1/365, lt=0.5, r=1, Bp= 0.16, 
           Wmax2=50, Kcwd=0.01, ag=0.1, rhoS1=0.4, rhoS2=0.5,rhoS3=0.55, deltaS3=0.0015, 
           deltaS4=0.0001, deltaS2=0.0090, lambdaS2=0.04, fS1=0.4) 
  
  PFT_GR=c(Amax=190, Ad=0.8, Kf=0.1, Tmin=4, Tmax=40, Topt=28,
           Kvpd=0.05, PARhalf=26, k=0.68, Don=dfPhenol$LeafOn[i],
           Doff=dfPhenol$LeafOff[i], Lmax=1, Ka=0, Q10v=2, 
           Kh=0.03, Q10s=2, f=0.04, Kwue=10.9, Wmax1=10, SLW=70,
           Cfrac=0.45, Kw=0.03, Aa= 1e5, 
           Ac=2e5, zbar=3, Cprecip=1, deltaA=0.01, deltaS1=0.37,
           lambdaS1 = 0.40, lambdaS3= 0.0001, bi=0.3, Tstar=50, W20=5, 
           Kr=0.55, aw=0, l=1/365, lt=0.5, r=1, Ks=1.2,
           Bp= 0.16, Wmax2=50, Kcwd=0, ag=0.1, bo2=0.00001, 
           rhoS1=0.4, rhoS2=0.5, rhoS3=0.55, deltaS3=0.0015, deltaS4=0.0001, 
           deltaS2=0.0090, lambdaS2=0.04, fS1=0.4) 
  
  subDF<-dfWeather[dfWeather$Site==paste("Region", ListRegions[i],sep=""), ]
  weatherDF<-subDF[rep(seq_len(nrow(subDF)), 200), ]
  weatherDF$runDay<-1:nrow(weatherDF)
  years2run=200 # years to run the model for
  times=seq(1,365*years2run)
  
  newsim<-data.frame(matrix(ncol=67, nrow=1460))
  colnames(newsim)<-c("time", "Cw", "Cl", "Cs1",
                      "Cs2", "Cs3", "Cs4", "Cdoc1",      
                      "Cdoc2", "W1", "W2", "Ca",         
                      "Cr", "Ccwd", "GPP", "Q1",         
                      "Q2", "Rf", "Ra", "NPP",        
                      "LCT1", "Rs3", "Dwater", "Wa",         
                      "TpotAreal", "i0", "T", "Dvpd",       
                      "Dtemp", "Dlightbar", "GPPmax", "LAIareal",   
                      "GPPpotAreal", "WUE", "V", "Lw",         
                      "Ll", "L", "Rs1", "Qout",       
                      "RfAreal", "Bdoc", "Q12", "LCT2",       
                      "Lr", "alCr", "lout", "Lg",          
                      "Bs1", "Ds1", "Ds3", "Ds4",         
                      "Ls3", "Ls1", "Rhdoc", "Rhdoc2",     
                      "Bs3", "Rs2",  "Ls2", "Ds2",        
                      "Lf1", "Rm", "Bs2", "Region",     
                      "pptGrads", "tmaxGrads", "PFT") 
  
  for(b in 1:4){
    repeat{  
      S0_EG<-c(Cw=PFT_i0$EGi0[i],Cl=100,Cs1=PFT_i0_co$EGi0[i],Cs2=PFT_i0_cf2$EGi0[i], Cs3=PFT_i0_cm$EGi0[i], Cs4=PFT_i0_cp$EGi0[i], Cdoc1=10, Cdoc2=7, W1=8, W2=5, Ca=500000, Cr=PFT_i0_cr$EGi0[i], Ccwd=20)
      S0_DE<-c(Cw=PFT_i0$DEi0[i],Cl=0,Cs1=PFT_i0_co$DEi0[i], Cs2=PFT_i0_cf2$DEi0[i],Cs3=PFT_i0_cm$DEi0[i], Cs4=PFT_i0_cp$DEi0[i], Cdoc1=10, Cdoc2=7, W1=8, W2=5, Ca=500000, Cr=PFT_i0_cr$DEi0[i], Ccwd=20)
      S0_GR<-c(Cw=PFT_i0$GRi0[i],Cl=0,Cs1=PFT_i0_co$GRi0[i], Cs2=PFT_i0_cf2$GRi0[i] ,Cs3=PFT_i0_cm$GRi0[i], Cs4=PFT_i0_cp$GRi0[i], Cdoc1=10, Cdoc2=7, W1=8, W2=5,  Ca=500000, Cr=PFT_i0_cr$GRi0[i], Ccwd=0)
      S0_SH<-c(Cw=PFT_i0$SHi0[i],Cl=0,Cs1=PFT_i0_co$SHi0[i], Cs2=PFT_i0_cf2$SHi0[i],Cs3=PFT_i0_cm$SHi0[i], Cs4=PFT_i0_cp$SHi0[i], Cdoc1=10, Cdoc2=7, W1=8, W2=5, Ca=500000, Cr=PFT_i0_cr$SHi0[i], Ccwd=10)
      
      PFT<-pftNam[[b]] 
      
      if(PFT=="EG"){
        S0=S0_EG
        params=PFT_EG
        egON="EG"
      } else if(PFT=="DE"){
        S0=S0_DE
        params=PFT_DE
        egON="NO"
      } else if(PFT=="GR"){
        S0=S0_GR
        params=PFT_GR
        egON="NO"
      } else if(PFT=="SH"){
        S0=S0_SH
        params=PFT_SH
        egON="NO"
      }
      
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
      
      sim2$Region<-weatherDF$Site[1]
      sim2$pptGrads<-subDF$AnnPpt[1]
      sim2$tmaxGrads<-subDF$Tmax[1]
      sim2$PFT<-pftNam[b]
      if(b==1){
        newsim[1:365,]<-sim2[72636:73000,]
      } else {
        newsim[(b*365-364):(b*365),]<-sim2[72636:73000,] 
      }
      iold_Cw=PFT_i0[i,b]  #initial conditions finder for Cw
      iold_Cp=PFT_i0_cp[i,b]  #initial conditions finder for Cs4
      iold_Cr=PFT_i0_cr[i,b]
      iold_Cm=PFT_i0_cm[i,b]
      inew_Cw<-findInitials(i0=iold_Cw, i1=sim2$Cw[nrow(sim2)])
      inew_Cp<-findInitials(i0=iold_Cp, i1=sim2$Cs4[nrow(sim2)])
      inew_Cr<-findInitials(i0=iold_Cr, i1=sim2$Cr[nrow(sim2)])
      inew_Cm<-findInitials(i0=iold_Cm, i1=sim2$Cs3[nrow(sim2)])
      PFT_i0[i,b]<-inew_Cw
      PFT_i0_cp[i,b]<-inew_Cp
      PFT_i0_cr[i,b]<-inew_Cr
      PFT_i0_cm[i,b]<-inew_Cm

      if((inew_Cw < 0 | inew_Cw==iold_Cw)&(inew_Cp < 0 | inew_Cp==iold_Cp)&(inew_Cm < 0 | inew_Cm==iold_Cm)&(inew_Cr < 0 | inew_Cr==iold_Cr)){
        break
      } 
    } #end repeat
  } #end b loop through PFTs
  stored_PFT<-rbind(stored_PFT,newsim)
  print(i)
} #end i loop through regions

#write out files
#write.csv(stored_PFT, "Validation/storedRegions.csv", row.names=F)
#write.csv(PFT_i0, "Validation/PFT_i0.csv", row.names=F)
#write.csv(PFT_i0_cp, "Validation/PFT_i0_cp.csv", row.names=F) # initial conditions CS4
#write.csv(PFT_i0_cr, "Validation/PFT_i0_cr.csv", row.names=F) #initial conditions Cr
#write.csv(PFT_i0_cm, "Validation/PFT_i0_cm.csv", row.names=F) #initial conditions CS3
#write.csv(PFT_i0_cf2, "Validation/PFT_i0_cf2.csv", row.names=F) #initial conditions CS2
#write.csv(PFT_i0_co, "Validation/PFT_i0_co.csv", row.names=F) #initial conditions CS1

##########Analyze output and comparison plots#############
#####################read in files made earlier in script#####
library(ggplot2)
library(reshape2)
RegionsOut<-read.csv("Validation/storedRegions.csv", stringsAsFactors=F)
dfPFT<-read.csv("Validation/regionsPFT.csv", stringsAsFactors=F)#write out PFT to keep
m2Region<-read.csv("Validation/regionsArea.csv", stringsAsFactors=F)#read in areas
runoff<-read.csv("Validation/regionsRunoff.csv", stringsAsFactors=F)#read in areas

ListRegions<-c("01", "02", "03W", "03S", "03N","04", "05", "06", 
               "07", "08", "09", "10L", "10U", "11", "12", "13", 
               "14", "15", "16", "17", "18")
pftNam<-c("EG", "DE", "GR", "SH") 

#Annual Summary calculation
sumRegions<-data.frame(matrix(ncol=6, nrow=0))
colnames(sumRegions)<-c("pft", "Variable", "Value", "Tmax", "Precip", "Region")
for(x in 1:21){
  subRegion<-RegionsOut[RegionsOut$Region==paste("Region",ListRegions[x],sep=""),]  
  
  pftSum1<-data.frame(matrix(ncol=6, nrow=0))
  colnames(pftSum1)<-c("pft", "Variable", "Value", "Tmax", "Precip")
  for(i in 1:4){
    pftSum<-data.frame(matrix(ncol=5, nrow=37))
    colnames(pftSum)<-c("pft", "Variable", "Value", "Tmax", "Precip")
    subPFT<-subRegion[subRegion$PFT==pftNam[i],]  
    pftSum[,1]<-pftNam[i]
    pftSum[,2:3]<-summaryTAM(df=subPFT) 
    pftSum[,4]<-subRegion$tmaxGrads[1] 
    pftSum[,5]<-subRegion$pptGrads[1]
    pftSum1<-rbind(pftSum1, pftSum)
  }#end i loop
  pftSum1$Region<-paste("Region", ListRegions[x], sep="")
  sumRegions<-rbind(sumRegions, pftSum1)
}#end x loop
sumRegions$Value<-as.numeric(sumRegions$Value)
data_wide <- dcast(sumRegions, pft + Region + Tmax + Precip ~ Variable, value.var = "Value", fun.aggregate = mean)

#make LCT 0 if unrealistic neg NPP/Cw/Cr is neg
for(i in 1:nrow(data_wide)){
  if(data_wide$Cwmin[i] < 0){
    data_wide$SumLdoc[i]<-0
  } else if (data_wide$NPP < 0){
    data_wide$SumLdoc[i]<-0
} else if (data_wide$Cr < 0){
  data_wide$SumLdoc[i]<-0
}}

totalRegions<-data.frame(matrix(ncol=5, nrow=21))
colnames(totalRegions)<-c("Region", "EG_gC", "DE_gC", "SH_gC", "GR_gC")
for(x in 1:21){
  regional<-data_wide[data_wide$Region==paste("Region", ListRegions[x], sep=""),]
  totalRegions$Region[x]<-paste("Region", ListRegions[x], sep="")
  
  totalRegions$DE_gC[x]<-(((dfPFT$dec_weighted[x]+(0.5*dfPFT$mx_weighted[x]))/100)*m2Region$m2area[x])*as.numeric(regional$SumLdoc[1])
  totalRegions$EG_gC[x]<-(((dfPFT$con_weighted[x]+(0.5*dfPFT$mx_weighted[x]))/100)*m2Region$m2area[x])*as.numeric(regional$SumLdoc[2]) 
  totalRegions$GR_gC[x]<-(((dfPFT$gr_weighted[x]+dfPFT$cr_weighted[x])/100)*m2Region$m2area[x])*as.numeric(regional$SumLdoc[3])
  totalRegions$SH_gC[x]<-((dfPFT$sh_weighted[x]/100)*m2Region$m2area[x])*as.numeric(regional$SumLdoc[4])
}
totalRegions$total<-(totalRegions$EG_gC+totalRegions$DE_gC+totalRegions$SH_gC+totalRegions$GR_gC)/1000000000000 #to TgC
totalR03<-((totalRegions$total[3]*m2Region$m2area[3])+(totalRegions$total[4]*m2Region$m2area[4])+(totalRegions$total[5]*m2Region$m2area[5]))/(m2Region$m2area[3]+m2Region$m2area[4]+m2Region$m2area[5])
totalR10<-((totalRegions$total[13]*m2Region$m2area[13])+(totalRegions$total[14]*m2Region$m2area[14]))/(m2Region$m2area[13]+m2Region$m2area[14])
subTotals<-data.frame(cbind(totalRegions$Region, totalRegions$total))
colnames(subTotals)<-c("Region", "TgC")
subTotals$Region[3]<-"Region03"
subTotals$TgC[3]<-totalR03
subTotals<-subTotals[-4,]
subTotals<-subTotals[-4,]
subTotals$Region[10]<-"Region10"
subTotals$TgC[10]<-totalR10
subTotals<-subTotals[-11,]
#write.csv(subTotals, "subTotalsHUC2.csv", row.names=F)

ButmanTable<-read.csv("Validation/ButmanTable.csv", stringsAsFactors=F)
ButmanTable$LCT<-ButmanTable$Lake.efflux..Tg.C.yr.1.*.30+ButmanTable$Steam.efflux..Tg.C.yr.1.*.30+(ButmanTable$Total.lateral.flux..TgC.yr.1.)*(ButmanTable$PCT_OC/100)+ButmanTable$Lake.burial..Tg.C.yr.1.
#Hotchkiss et al 2015: ~28% of stream/river emissions are from terrestrial DOC
AllLCT<-data.frame(cbind(subTotals, ButmanTable$LCT))
colnames(AllLCT)<-c("Region", "model", "Butman")
AllLCT$model<-as.numeric(AllLCT$model)
AllLCT$RegNum<-1:nrow(AllLCT)

borderTheme0.5<-theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
                      panel.background = element_rect(colour = "black", size=1, fill=NA), 
                      strip.background = element_blank(), 
                      panel.spacing = unit(0.3, "lines"), axis.ticks.length = unit(0.2, "lines"), 
                      text = element_text(size = 12),legend.background=element_blank(), 
                      legend.key = element_rect(fill = NA), aspect.ratio = 1)

library(ggrepel)
plot1<-ggplot(data=AllLCT, aes(x=Butman, y=model))+geom_point(size=0.7)+
  #geom_text(label=AllLCT$Region)+
  scale_x_continuous(limits=c(0,12), breaks=seq(0, 12, 3), expand=c(0,0))+
  scale_y_continuous(limits=c(0,12), breaks=seq(0, 12, 3), expand=c(0,0))+
  geom_abline(slope=1, intercept = 0, linetype="dashed", color="grey")+
  labs(y=expression('Modeled LCT'~(Tg~C~yr^-1)),x=expression('Literature LCT'~(Tg~C~yr^-1)))+
  geom_text(label=AllLCT$RegNum, hjust = -0.1, nudge_x = 0.2,size=2.5)+
  annotate("text", x=1, y=11, label= "A.", size=4.5)+
  borderTheme0.5
plot1
#ggsave(filename = "Validation/ButmanPointsLCT.png", plot=plot1, width = 3.5, height = 2.5, units= "in", device='png', dpi=320)
colnames(AllLCT)<-c("Region", "ModelLCT", "ButmanLCT", "RegNum")
#write.csv(AllLCT, "Validation/CompareLCT.csv", row.names = F)

###NPP###
totalRegionsNPP<-data.frame(matrix(ncol=6, nrow=21)) #empty to store calculations
colnames(totalRegionsNPP)<-c("Region", "EG_gC", "DE_gC", "SH_gC", "GR_gC", "total")
for(x in 1:21){
  regional<-data_wide[data_wide$Region==paste("Region", ListRegions[x], sep=""),]
  totalRegionsNPP$Region[x]<-paste("Region", ListRegions[x], sep="")
  
  #total NPP in g C yr-1 for each HUC
  totalRegionsNPP$DE_gC[x]<-(((dfPFT$dec_weighted[x]+(0.5*dfPFT$mx_weighted[x]))/100)*m2Region$m2area[x])*as.numeric(regional$NPP[1])
  totalRegionsNPP$EG_gC[x]<-(((dfPFT$con_weighted[x]+(0.5*dfPFT$mx_weighted[x]))/100)*m2Region$m2area[x])*as.numeric(regional$NPP[2]) 
  totalRegionsNPP$GR_gC[x]<-(((dfPFT$gr_weighted[x]+dfPFT$cr_weighted[x])/100)*m2Region$m2area[x])*as.numeric(regional$NPP[3])
  totalRegionsNPP$SH_gC[x]<-((dfPFT$sh_weighted[x]/100)*m2Region$m2area[x])*as.numeric(regional$NPP[4])
  #to g C m-2 yr-1
  totalRegionsNPP$total[x]<-(totalRegionsNPP$EG_gC[x]+totalRegionsNPP$DE_gC[x]+totalRegionsNPP$SH_gC[x]+totalRegionsNPP$GR_gC[x])/((((dfPFT$dec_weighted[x]+(0.5*dfPFT$mx_weighted[x]))/100)*m2Region$m2area[x])+((((dfPFT$con_weighted[x]+(0.5*dfPFT$mx_weighted[x]))/100)*m2Region$m2area[x])+(((dfPFT$gr_weighted[x]+dfPFT$cr_weighted[x])/100)*m2Region$m2area[x])+((dfPFT$sh_weighted[x]/100)*m2Region$m2area[x]))) 
}
#calculate NPP for modified HUCs (to match Butman et al. 2016)
totalR03<-(((totalRegionsNPP$total[3]*m2Region$m2area[3])+(totalRegionsNPP$total[4]*m2Region$m2area[4])+(totalRegionsNPP$total[5]*m2Region$m2area[5]))/(m2Region$m2area[3]+m2Region$m2area[4]+m2Region$m2area[5]))
totalR10<-(((totalRegionsNPP$total[13]*m2Region$m2area[13])+(totalRegionsNPP$total[14]*m2Region$m2area[14]))/(m2Region$m2area[13]+m2Region$m2area[14]))
subTotalsNPP<-data.frame(cbind(totalRegionsNPP$Region, totalRegionsNPP$total))
colnames(subTotalsNPP)<-c("Region", "TgC")
subTotalsNPP$Region[3]<-"Region03"
subTotalsNPP$TgC[3]<-totalR03
subTotalsNPP<-subTotalsNPP[-4,]
subTotalsNPP<-subTotalsNPP[-4,]
subTotalsNPP$Region[10]<-"Region10"
subTotalsNPP$TgC[10]<-totalR10
subTotalsNPP<-subTotalsNPP[-11,]

AllLCT<-cbind(subTotalsNPP, ButmanTable$NPP, ButmanTable$NPPsd)
colnames(AllLCT)<-c("Region", "model", "Butman", "Butmansd")
AllLCT$model<-as.numeric(AllLCT$model)
AllLCT$RegNum<-1:nrow(AllLCT)
plot2<-ggplot(data=AllLCT, aes(x=Butman, y=model))+
  scale_x_continuous(limits=c(0,1300), breaks=seq(0, 1200, 300), expand=c(0,0))+
  scale_y_continuous(limits=c(0,1300), breaks=seq(0, 1200, 300), expand=c(0,0))+geom_abline(slope=1, intercept = 0, linetype="dashed", color="grey")+
  labs(y=expression('Modeled NPP'~(g~C~m^-2~yr^-1)),x=expression('Literature NPP'~(g~C~m^-2~yr^-1)))+
  geom_pointrange(aes(xmin=Butman-Butmansd, xmax=Butman+Butmansd), size=0.1, color="black")+
  geom_point(size=0.7, color="black")+
  geom_text(label=AllLCT$RegNum, hjust = -0.5, nudge_x = 0.3,size=2.5)+
  annotate("text", x=100, y=1200, label= "B.", size=4.5)+
  borderTheme0.5
plot2
#ggsave(filename = "Validation/ButmanPointsNPP.png", plot=plot2, width = 3.6, height = 2.6, units= "in", device='png', dpi=320)
colnames(AllLCT)<-c("Region", "ModelNPP", "ButmanNPP", "ButmanSD", "RegNum")
#write.csv(AllLCT, "Validation/CompareNPP.csv", row.names = F)

totalRegionsLCTareal<-data.frame(matrix(ncol=5, nrow=21))
colnames(totalRegionsLCTareal)<-c("Region", "EG_gC", "DE_gC", "SH_gC", "GR_gC")
for(x in 1:21){
  regional<-data_wide[data_wide$Region==paste("Region", ListRegions[x], sep=""),]
  totalRegionsLCTareal$Region[x]<-paste("Region", ListRegions[x], sep="")
  #area weighted LCT g C m-2 yr-1 
  totalRegionsLCTareal$DE_gC[x]<-(((dfPFT$dec_weighted[x]+(0.5*dfPFT$mx_weighted[x]))/100)*m2Region$m2area[x])*as.numeric(regional$SumLdoc[1])/m2Region$m2area[x]
  totalRegionsLCTareal$EG_gC[x]<-(((dfPFT$con_weighted[x]+(0.5*dfPFT$mx_weighted[x]))/100)*m2Region$m2area[x])*as.numeric(regional$SumLdoc[2])/m2Region$m2area[x]
  totalRegionsLCTareal$GR_gC[x]<-(((dfPFT$gr_weighted[x]+dfPFT$cr_weighted[x])/100)*m2Region$m2area[x])*as.numeric(regional$SumLdoc[3])/m2Region$m2area[x]
  totalRegionsLCTareal$SH_gC[x]<-((dfPFT$sh_weighted[x]/100)*m2Region$m2area[x])*as.numeric(regional$SumLdoc[4])/m2Region$m2area[x]
}
totalRegionsLCTareal$total<-(totalRegionsLCTareal$EG_gC+totalRegionsLCTareal$DE_gC+totalRegionsLCTareal$SH_gC+totalRegionsLCTareal$GR_gC)

for(i in 1:nrow(totalRegionsLCTareal)){
  subdwide<-data_wide[data_wide$Region==totalRegionsLCTareal$Region[i], ]
  totalRegionsLCTareal$Precip[i]<-subdwide$Precip[1]
  totalRegionsLCTareal$Tmax[i]<-subdwide$Tmax[1]
}
arealplusclim<-data.frame(cbind(totalRegionsLCTareal$Region, totalRegionsLCTareal$total, totalRegionsLCTareal$Precip, totalRegionsLCTareal$Tmax))
colnames(arealplusclim)<-c("Region", "LCT", "Precip", "Tmax")
for(i in 1:nrow(arealplusclim)){
  arealplusclim$RegNum[i]<-unlist(strsplit(arealplusclim$Region[i], split="Region"))[2]
}
arealplusclim$Tmax<-as.numeric(arealplusclim$Tmax)
arealplusclim$Precip<-as.numeric(arealplusclim$Precip)
arealplusclim$LCT<-as.numeric(arealplusclim$LCT)
plot3<-ggplot(data=arealplusclim, aes(x=Precip, y=LCT, group=Tmax, colour=Tmax))+
  scale_x_continuous(limits=c(25,162), breaks=seq(25, 150, 25), expand=c(0,0))+
  scale_y_continuous(limits=c(0,25), breaks=seq(0, 25, 5), expand=c(0,0))+
  labs(y=expression('LCT'~(g~C~m^-2~yr^-1)),x=expression('Precipitaiton'~(cm~yr^-1)))+
  scale_color_viridis_c("Temp (C)", option="B", alpha = 0.8, begin = 0, end = 0.9)+
  geom_point(size=1)+
  geom_text_repel(label=arealplusclim$RegNum,size=2, color="black")+
  borderTheme0.5
plot3
#ggsave(filename = "Validation/ArealLCT.png", plot=plot3, width = 3.6, height = 2.6, units= "in", device='png', dpi=320)
#write.csv(arealplusclim, "Validation/ModeledArealLCT.csv", row.names = F)
#######RUnoff ratio
totalRegions<-data.frame(matrix(ncol=6, nrow=21))
colnames(totalRegions)<-c("Region", "EG_cm", "DE_cm", "SH_cm", "GR_cm", "total")
for(x in 1:21){
  regional<-data_wide[data_wide$Region==paste("Region", ListRegions[x], sep=""),]
  totalRegions$Region[x]<-paste("Region", ListRegions[x], sep="")
  
  totalRegions$DE_cm[x]<-(((dfPFT$dec_weighted[x]+(0.5*dfPFT$mx_weighted[x]))/100)*m2Region$m2area[x])*as.numeric(regional$SumD[1])/m2Region$m2area[x]
  totalRegions$EG_cm[x]<-(((dfPFT$con_weighted[x]+(0.5*dfPFT$mx_weighted[x]))/100)*m2Region$m2area[x])*as.numeric(regional$SumD[2])/m2Region$m2area[x] 
  totalRegions$GR_cm[x]<-(((dfPFT$gr_weighted[x]+dfPFT$cr_weighted[x])/100)*m2Region$m2area[x])*as.numeric(regional$SumD[3])/m2Region$m2area[x]
  totalRegions$SH_cm[x]<-((dfPFT$sh_weighted[x]/100)*m2Region$m2area[x])*as.numeric(regional$SumD[4])/m2Region$m2area[x]
  #totalRegions$total[x]<-(totalRegions$EG_cm[x]+totalRegions$DE_cm[x]+totalRegions$SH_cm[x]+totalRegions$GR_cm[x])/((((dfPFT$dec_weighted[x]+(0.5*dfPFT$mx_weighted[x]))/100)*m2Region$m2area[x])+((((dfPFT$con_weighted[x]+(0.5*dfPFT$mx_weighted[x]))/100)*m2Region$m2area[x])+(((dfPFT$gr_weighted[x]+dfPFT$cr_weighted[x])/100)*m2Region$m2area[x])+((dfPFT$sh_weighted[x]/100)*m2Region$m2area[x]))) 
}

totalRegions$total<-(totalRegions$EG_cm+totalRegions$DE_cm+totalRegions$SH_cm+totalRegions$GR_cm)

totalR03<-((totalRegions$total[3]*m2Region$m2area[3])+(totalRegions$total[4]*m2Region$m2area[4])+(totalRegions$total[5]*m2Region$m2area[5]))/(m2Region$m2area[3]+m2Region$m2area[4]+m2Region$m2area[5])
totalR10<-((totalRegions$total[13]*m2Region$m2area[13])+(totalRegions$total[14]*m2Region$m2area[14]))/(m2Region$m2area[13]+m2Region$m2area[14])
subTotals<-data.frame(cbind(totalRegions$Region, totalRegions$total))
colnames(subTotals)<-c("Region", "cm_m2")
subTotals$Region[3]<-"Region03"
subTotals$cm_m2[3]<-totalR03
subTotals<-subTotals[-4,]
subTotals<-subTotals[-4,]
subTotals$Region[10]<-"Region10"
subTotals$cm_m2[10]<-totalR10
subTotals<-subTotals[-11,]

##runoff
totalPrecipR03<-((runoff$Precip[3]*m2Region$m2area[3])+(runoff$Precip[4]*m2Region$m2area[4])+(runoff$Precip[5]*m2Region$m2area[5]))/(m2Region$m2area[3]+m2Region$m2area[4]+m2Region$m2area[5])
totalPrecipR10<-((runoff$Precip[13]*m2Region$m2area[13])+(runoff$Precip[14]*m2Region$m2area[14]))/(m2Region$m2area[13]+m2Region$m2area[14])
totalRatioR03<-((runoff$RunoffRatio[3]*m2Region$m2area[3])+(runoff$RunoffRatio[4]*m2Region$m2area[4])+(runoff$RunoffRatio[5]*m2Region$m2area[5]))/(m2Region$m2area[3]+m2Region$m2area[4]+m2Region$m2area[5])
totalRatioR10<-((runoff$RunoffRatio[13]*m2Region$m2area[13])+(runoff$RunoffRatio[14]*m2Region$m2area[14]))/(m2Region$m2area[13]+m2Region$m2area[14])
runoff$Region[3]<-"Region03"
runoff$Precip[3]<-totalPrecipR03
runoff$RunoffRatio[3]<-totalRatioR03
runoff<-runoff[-4,]
runoff<-runoff[-4,]
runoff$Region[10]<-"Region10"
runoff$Precip[10]<-totalPrecipR10
runoff$RunoffRatio[10]<-totalRatioR10
runoff<-runoff[-11,]


AllRunoff<-cbind(runoff, subTotals$cm_m2)
colnames(AllRunoff)[5]<-c("ModelRunoff")
AllRunoff$ModelRunoff<-as.numeric(AllRunoff$ModelRunoff)
AllRunoff$ModelRatio<-AllRunoff$ModelRunoff/AllRunoff$Precip
AllRunoff$RegNum<-1:nrow(AllRunoff)

library(ggrepel)
plot4<-ggplot(data=AllRunoff, aes(x=RunoffRatio, y=ModelRatio))+geom_point(size=0.7)+
  #geom_text(label=AllLCT$Region)+
  scale_x_continuous(limits=c(0,0.8), breaks=seq(0, 0.8, 0.2), expand=c(0,0))+
  scale_y_continuous(limits=c(0,0.8), breaks=seq(0, 0.8, 0.2), expand=c(0,0))+
  geom_abline(slope=1, intercept = 0, linetype="dashed", color="grey")+
  labs(y=expression('Modeled runoff ratio'),x=expression('Runoff ratio'))+
  geom_text(label=AllLCT$RegNum, hjust = -0.1, nudge_x = 0.007,size=2.5)+
  #annotate("text", x=1, y=11, label= "A.", size=4.5)+
  borderTheme0.5
plot4
#ggsave(filename = "Validation/RunoffRatios.png", plot=plot4, width = 3.5, height = 2.5, units= "in", device='png', dpi=320)

#write.csv(AllRunoff, "Validation/RunoffRatios.csv", row.names = F)

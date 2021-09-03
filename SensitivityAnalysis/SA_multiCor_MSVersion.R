######Run sensitivity analysis in parallel######
#load packages
library(reshape2)
library(parallel)
library(deSolve)
#set up working directory, which contains all input files & output folder
setwd("") #set to location of terrestrial-aquatic-LCT folder
#source functions
source("Functions/findInitials_MSVersion.R") #find initial conditions
source("Functions/TAM_MSVersion.R") #model code
#read in files
dfOut<-read.csv("SensitivityAnalysis/dfOut.csv", stringsAsFactors=F) #weather forcings
phenolDF<-read.csv("SensitivityAnalysis/phenolDF.csv", stringsAsFactors=F) #plant phenology forcings
trait_i0<-read.csv("SensitivityAnalysis/SASoil_trait_i0.csv", stringsAsFactors=F) #Cw initial conditions
trait_i0_cp<-read.csv("SensitivityAnalysis/SASoil_trait_i0_cp.csv", stringsAsFactors=F) #Cs4 initial conditions
traitList<-c("Val2", "Val3") #Val2 = -10%, Val3 = +10%, Val1 is base case.. 

#Loop to make a list for each parameter, where it is decreased by 10% (1:51) and
#increased by 10% (52:102)
paramAll<-list() #to store list for each parameter
for(p in 1:51){
  params=c(Amax=112, Ad=0.75, Kf=0.1, Tmin=4, Tmax=40, Topt=24,
          Kvpd=0.05, PARhalf=17, k=0.58, Don=100, Doff=280, Lmax=3, 
          Ka=0.006, Q10v=2,Kh=0.03, Q10s=2, f=0.04, Kwue=10.9, Wmax1=10, SLW=70,
          Cfrac=0.45, Kw=0.005, Aa= 1e5, Ac= 2e5, zbar= 3, 
          Cprecip = 1, deltaA = 0.01, deltaS1 = 0.37, Ks=1.2,
          lambdaS1 = 0.40, lambdaS3= 0.0001, bi = 0.3, Tstar = 50, W20 = 5, 
          Kr= 0.55, aw= 0.1, l = 1/365, lt = 0.5, r=1, Bp= 0.16, 
          Wmax2=50, Kcwd=0.01, ag=0.1, rhoS1=0.4, rhoS2=0.5, rhoS3=0.55, 
          deltaS3=0.0014, deltaS4=0.0001, deltaS2=0.0090, lambdaS2=0.04, fS1=0.6)
  paramList<-names(params)
  paramVals<-c(rep(params, times=2))
  paramVals[[p]]<-params[[p]]*.90 #-10%
  paramVals[[p+51]]<-params[[p]]+params[[p]]*.10 #+10%
  
  paramAll<-c(paramAll, list(paramVals))
  #names(paramAll)[p]<-paste(names(params[p]))
  #mclapply(FUN=SA_fun(paramsIn=paramVals, ptype=names(paramVals[p])))
  
} #end p loop

#A function to run simulations for sensitivity analysis 
#over climate gradients for use with mclapply
SA_fun<-function(paramsIn){
  
  stored_PFT<-data.frame(matrix(ncol=67, nrow=43800))#DF to store results from looped sims
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
  
  V1<-paramsIn[1:51] #first param list with one param decreased by 10%
  V2<-paramsIn[52:102] #second param list with one param increased by 10%
  pDiff=c() 
  for(z in 1:51){
    pDiff[z]<-V1[[z]]-V2[[z]]
  }
  currParam<-names(paramsIn[which(pDiff!=0)]) #identifies which param the mclapply is on 
  for(b in 1:2){
    #to make a loop for params
    newsim<-data.frame(matrix(ncol=66, nrow=21900)) #60*365=21900; number of climate gradients * number of days in last year of simulations
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
    phenolDFsub<-phenolDF[phenolDF$trait=="Val1",]
    egON="NO" #if "NO", then deciduous phenology is on
    for(i in 1:60){
      #select params depending on b
      #this will replace the phenology parameters to match the climate 
      if(currParam=="Don"){
        paramsIn[10]<-as.numeric(phenolDF$LeafOn[i])-7
        paramsIn[61]<-as.numeric(phenolDF$LeafOn[i+60])+7
      } else if(currParam=="Doff"){
        paramsIn[11]<-as.numeric(phenolDF$LeafOff[i])-7
        paramsIn[62]<-as.numeric(phenolDF$LeafOff[i+60])+7
      } else{
        paramsIn[[10]]<-as.numeric(phenolDF$LeafOn[i])
        paramsIn[[61]]<-as.numeric(phenolDF$LeafOn[i])
        paramsIn[[11]]<-as.numeric(phenolDF$LeafOff[i])
        paramsIn[[62]]<-as.numeric(phenolDF$LeafOff[i])
      }
      #select for b 
      if(b==1){
        params=paramsIn[1:51]
      } else{
        params=paramsIn[52:104]}
      
      repeat{ 
        S0<-c(Cw=trait_i0[i,b],Cl=0,Cs1=5,Cs2=6,Cs3=300, Cs4=trait_i0_cp[i,b], Cdoc1=10, Cdoc2=7, W1=8, W2=5, Ca=500000, Cr=400, Ccwd=20)
        
        subDF<-dfOut[dfOut$runNum==i, ] #subset weather
        weatherDF<-subDF[rep(seq_len(nrow(subDF)), 200), ] #repeat data for 200 years
        weatherDF$runDay<-1:nrow(weatherDF) #assign run day
        # years to run the model for
        years2run=200 
        times=seq(1,365*years2run)
        
        # define forcing approx functions
        PARapprox<<-approxfun(x=weatherDF$runDay, y = weatherDF$PAR) 
        Papprox<<-approxfun(x =weatherDF$runDay, y = weatherDF$Precip) 
        VPDapprox<<-approxfun(x=weatherDF$runDay, y = weatherDF$VPD)  
        Tair_approx<<-approxfun(x=weatherDF$runDay, y = weatherDF$MeanTemp)
        Tsoil_approx<<-approxfun(x=weatherDF$runDay, y = weatherDF$SoilTemp)
        Evap_approx<<-approxfun(x=weatherDF$runDay, y= weatherDF$Evap)
        Lapprox<<-approxfun(times,rep(c(rep(0,params[10]),(params[12]*params[20]*params[21]),rep(0,365-params[10]-1)),years2run))   
        Ll_approx<<-approxfun(times,rep(c(rep(0,params[11]),(params[12]*params[20]*params[21]),rep(0,365-params[11]-1)),years2run))  
        
        # runs simulation
        sim2=ode(y=S0,times=times,func=tamStep,parms=params,method="euler")
        sim2=data.frame(sim2)
        
        #label simulations with climate info
        sim2$runNum<-weatherDF$runNum[1] 
        sim2$pptGrads<-weatherDF$AnnPpt[1]
        sim2$tmaxGrads<-weatherDF$Tmax[1]
        #store last year of simulations only
        if(i==1){
          newsim[1:365,]<-sim2[72636:73000,] 
        } else {
          newsim[(i*365-364):(i*365),]<-sim2[72636:73000,]
        }
        
        #find initial conditions 
        iold_Cw=trait_i0[i,b]  # Cw
        iold_Cp=trait_i0_cp[i,b]  # Cs4
        inew_Cw<-findInitials(i0=iold_Cw, i1=sim2$Cw[nrow(sim2)])
        inew_Cp<-findInitials(i0=iold_Cp, i1=sim2$Cs4[nrow(sim2)])
        trait_i0[i,b]<-inew_Cw
        trait_i0_cp[i,b]<-inew_Cp
        
        if((inew_Cw < 0 | inew_Cw==iold_Cw)&(inew_Cp < 0 | inew_Cp==iold_Cp)){
          break
        } 
      }
      
      #print(i)
    } # end i loop
    newsim$trait<-traitList[b] #label the trait for the simulation
    if(b==1){
      stored_PFT[1:21900,]<-newsim
    }else {
      stored_PFT[(b*21900-21899):(b*21900),]<-newsim
    }
  } #end b loop     
  
  fiName<-paste("SASoil",currParam, sep="_")
  fiName1<-paste("Out/", fiName, ".csv", sep="")
  write.csv(stored_PFT, fiName1, row.names=F)
  fiName2<-paste("Out/", fiName,"_i0_cp", ".csv", sep="")
  write.csv(trait_i0_cp, fiName2, row.names=F)
  fiName3<-paste("Out/", fiName, "_i0_cw",".csv", sep="")
  write.csv(trait_i0, fiName3, row.names=F)
  
} #end SA_fun

#parallelize to run SA_fun for each param in paramAll(n=51)
mclapply(X=paramAll, FUN=SA_fun,mc.preschedule = T, mc.cores=48)


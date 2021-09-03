###running SIPNETPLUSLAKES with traits###
#packages
library(ggplot2)
library(viridis)
library(geosphere)
library(deSolve)
#set working directory
setwd("") #set to location of terrestrial-aquatic-LCT folder
source("Functions/TAM_MSVersion.R") #model code
source("Functions/summaryTAM_MSVersion.R") #summary function
source("Functions/findInitials_MSVersion.R") #initial conditions finder
#read in files
dfOut<-read.csv("SensitivityAnalysis/dfOut.csv") #data from weather generator
phenolDF<-read.csv("SensitivityAnalysis/phenolDF.csv") #leaf on/off phenology
#initial conditions for wood pool
TOPT_i0<-read.csv("ModelExperiments/Topt_i0_cw_NewS.csv") #for Topt param
AMAX_i0<-read.csv("ModelExperiments/Amax_i0_cw_NewS.csv") #for Amax param
WUE_i0<-read.csv("ModelExperiments/Kwue_i0_cw_NewS.csv") #for Kwue param
AW_i0<-read.csv("ModelExperiments/aw_i0_cw_NewS.csv") #for aw param
#initial conditions for S4 (passive soil C pool)
TOPT_i0_cp<-read.csv("ModelExperiments/Topt_i0_cp_NewS.csv") #for Topt param
AW_i0_cp<-read.csv("ModelExperiments/aw_i0_cp_NewS.csv") #for aw param
AMAX_i0_cp<-read.csv("ModelExperiments/Amax_i0_cp_NewS.csv") #for Amax param
WUE_i0_cp<-read.csv("ModelExperiments/Kwue_i0_cw_NewS.csv") #for Kwue param


###create gradients####
#We later remove the 25 cm yr-1 precip gradient because it tends to create
#unrealistic negative NPP at most temperatures.
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

  #list of param/trait values 
  AwList<-c(0, 0.1, 0.2, 0.3, 0.4) 
  kwueList<-c(6, 8, 10, 12, 14)
  AmaxList<-c(60, 90, 120, 150, 180)
  ToptList<-c(20, 22, 24, 26, 28)

######Run simulations
listP<-c("Amax", "Topt", "Kwue", "aw") #parameters in analysis

#loops over each param in listP
for(m in 1:length(listP)){
  if(m==1){
    traitList<-AmaxList
    cw_init<-AMAX_i0
    cp_init<-AMAX_i0_cp
  } else if(m==2){
    traitList<-ToptList
    cw_init<-TOPT_i0
    cp_init<-TOPT_i0_cp
  } else if(m==3){
    traitList<-kwueList
    cw_init<-WUE_i0
    cp_init<-WUE_i0_cp
  } else if(m==4){
    traitList<-AwList
    cw_init<-AW_i0
    cp_init<-AW_i0_cp
  }
  
stored_PFT<-data.frame(matrix(ncol=67, nrow=109500))#DF to store results from looped sims
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

for(b in 1:length(traitList)){
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
 repeat{ S0<-c(Cw=cw_init[i,b],Cl=0,Cs1=5,Cs2=6,Cs3=300, Cs4=cp_init[i,b], Cdoc1=10, Cdoc2=7, W1=8, W2=5, Ca=500000, Cr=400, Ccwd=20)
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
    
    if(m==1){
      params[[1]]<-traitList[b] #amax
    }else if(m==2){
      params[[6]]<-traitList[b] #topt
    }else if(m==3){
      params[[18]]<-traitList[b] #kwue
    }else if (m==4){
      params[[36]]<-traitList[b] #aw
    }
    subDF<-dfOut[dfOut$runNum==i, ]
    weatherDF<-subDF[rep(seq_len(nrow(subDF)), 100), ]
    weatherDF$runDay<-1:nrow(weatherDF)
    ### setting boundary conditions, parameters, and forcings
    trait<-traitList[[b]]
    # years to run the model for
    years2run=100
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
    
    sim2$runNum<-weatherDF$runNum[1] #label run number
    sim2$pptGrads<-weatherDF$AnnPpt[1] #label precip
    sim2$tmaxGrads<-weatherDF$Tmax[1] #label tmax
    #store last year of simulations (at equilibrium)
    if(i==1){
      newsim[1:365,]<-sim2[36136:36500,] #89426:91250
    } else {
      newsim[(i*365-364):(i*365),]<-sim2[36136:36500,] 
    }
    
    iold_Cw=cw_init[i,b]  #initial conditions finder for Cw
    iold_Cp=cp_init[i,b]  #initial conditions finder for Cp
    inew_Cw<-findInitials(i0=iold_Cw, i1=sim2$Cw[nrow(sim2)])
    inew_Cp<-findInitials(i0=iold_Cp, i1=sim2$Cs4[nrow(sim2)])
    if(inew_Cw!=cw_init[i,b] | inew_Cp!=cp_init[i,b]) {
      cw_init[i,b]<-inew_Cw
      cp_init[i,b]<-inew_Cp}
    if((inew_Cw < 0 | inew_Cw==iold_Cw)&(inew_Cp < 0 | inew_Cp==iold_Cp)){
      break
    } 
 }  
   print(i)
   }
  newsim$trait<-traitList[b] #label the trait value for the simulation
  if(b==1){
    stored_PFT[1:21900,]<-newsim
  }else {
    stored_PFT[(b*21900-21899):(b*21900),]<-newsim
  }
}
#name and write out last year of simulations
finame_df<-paste("ModelExperiments/",listP[m], "Out_NewS.csv", sep="")
finame_init_cp<-paste("ModelExperiments/",listP[m], "_i0_cp_NewS.csv", sep="")
finame_init_cw<-paste("ModelExperiments/",listP[m], "_i0_cw_NewS.csv", sep="")

#write.csv(stored_PFT, finame_df, row.names=F)
#write.csv(cp_init, finame_init_cp, row.names=F)
#write.csv(cw_init, finame_init_cw, row.names=F)

#Annual summary 
sumDF<-data.frame(matrix(ncol=6, nrow=0))#store summary (37 vars) for each gradient (60) and each trait value
colnames(sumDF)<-c("runNum", "Variable", "Value", "Tmax", "Precip", "trait")

for(b in 1:length(traitList)){
subDF<-stored_PFT[stored_PFT$trait==traitList[b],]
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
    gradSum[,6]<-traitList[b]
    if(x==1){
      sumDF1[1:37,]<-gradSum
    } else{
      sumDF1[((x*37)-36):(x*37),]<-gradSum
    }
  }#end x loop
    sumDF<-rbind(sumDF, sumDF1)
}#end b loop

sumDF$Value<-as.numeric(sumDF$Value)
data_wide <- dcast(sumDF, trait + runNum + Tmax + Precip ~ Variable, value.var = "Value", fun.aggregate = mean)
data_wide <-data_wide[(data_wide$Cw >= 0) & (data_wide$NPP > 0) &(data_wide$Cr > 0), ] #remove unrealistic negative NPP and sims w/neg wood or roots

finame_wide<-paste("ModelExperiments/", listP[m], "_wide.csv", sep="") 
#write.csv(data_wide, finame_wide, row.names=F) #write out processed data
assign(finame_wide, data_wide)
} #end m loop

###density plots
borderTheme0.5<-theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
                      panel.background = element_rect(colour = "black", size=1.5, fill=NA), 
                      strip.background = element_blank(), 
                      panel.spacing = unit(0.3, "lines"), axis.ticks.length = unit(0.2, "lines"), 
                      text = element_text(size = 14),legend.background=element_blank(), 
                      legend.key = element_rect(fill = NA),)

mycols<-viridis(6, alpha = 0.6, end = 1, option = "D")
#for kwue
Kwue_wide$trait_f<-as.factor(Kwue_wide$trait)
Kwue_wide<-Kwue_wide[Kwue_wide$Precip!=25,]
bp<-ggplot(Kwue_wide, aes(x=SumLdoc, fill=trait_f, group=trait_f))+
  scale_fill_manual(values=mycols, expression("K"[wue]))+
  geom_density()+borderTheme0.5+
  labs(x=expression('LCT'~(g~C~m^-2~yr^-1)), y=expression('Density'))+
  scale_x_continuous(limits=c(0,125), breaks=seq(0, 125, 25), expand=c(0.004,0))+
  scale_y_continuous(limits=c(0, 0.025), breaks=seq(0, 0.025, 0.005), expand=c(0.008, 0))+
  annotate("text", x=8, y=0.0235, label="A.", size=4.5)
bp  
#ggsave(filename = "KwueDensity_NewS.png", plot=bp, width = 3.5, height = 2.5, units= "in", device='png', dpi=320)

#for Amax
Amax_wide$trait_f<-as.factor(Amax_wide$trait)
Amax_wide<-Amax_wide[Amax_wide$Precip!=25,]
bp<-ggplot(Amax_wide, aes(x=SumLdoc, fill=trait_f, group=trait_f))+
  scale_fill_manual(values=mycols, expression("A"[max]))+
  geom_density()+borderTheme0.5+
  labs(x=expression('LCT'~(g~C~m^-2~yr^-1)), y=expression('Density'))+
  scale_x_continuous(limits=c(0,125), breaks=seq(0, 125, 25), expand=c(0.004,0))+
  scale_y_continuous(limits=c(0, 0.025), breaks=seq(0, 0.025, 0.005), expand=c(0.008, 0))+
  annotate("text", x=8, y=0.0235, label="B.", size=4.5)
bp  
#ggsave(filename = "AmaxDensity_NewS.png", plot=bp, width = 3.5, height = 2.5, units= "in", device='png', dpi=320)

#for Topt
Topt_wide$trait_f<-as.factor(Topt_wide$trait)
Topt_wide<-Topt_wide[Topt_wide$Precip!=25,]
bp<-ggplot(Topt_wide, aes(x=SumLdoc, fill=trait_f, group=trait_f))+
  scale_fill_manual(values=mycols, expression("T"[opt]))+
  geom_density()+borderTheme0.5+
  labs(x=expression('LCT'~(g~C~m^-2~yr^-1)), y=expression('Density'))+
  scale_x_continuous(limits=c(0,125), breaks=seq(0, 125, 25), expand=c(0.004,0))+
  scale_y_continuous(limits=c(0, 0.025), breaks=seq(0, 0.025, 0.005), expand=c(0.008, 0))+
  annotate("text", x=8, y=0.0235, label="C.", size=4.5)
bp  
#ggsave(filename = "ToptDensity_NewS.png", plot=bp, width = 3.5, height = 2.5, units= "in", device='png', dpi=320)

#for aw
aw_wide$trait_f<-as.factor(aw_wide$trait)
aw_wide<-aw_wide[Aw_wide$Precip!=25,]
bp<-ggplot(Aw_wide, aes(x=SumLdoc, fill=trait_f, group=trait_f))+
  scale_fill_manual(values=mycols, expression("a"[w]))+
  geom_density()+borderTheme0.5+
  labs(x=expression('LCT'~(g~C~m^-2~yr^-1)), y=expression('Density'))+
  scale_x_continuous(limits=c(0,125), breaks=seq(0, 125, 25), expand=c(0.004,0))+
  scale_y_continuous(limits=c(0, 0.025), breaks=seq(0, 0.025, 0.005), expand=c(0.008, 0))+
  annotate("text", x=8, y=0.0235, label="D.", size=4.5)
bp  
#ggsave(filename = "AwDensity_NewS.png", plot=bp, width = 3.5, height = 2.5, units= "in", device='png', dpi=320)

###statistics
ScoreVals<-function(Vals){
  my_sd<- sd(Vals)*sqrt((length(Vals)-1)/(length(Vals)))
  my_mean <- mean(Vals)
  my_scored<-(Vals-my_mean)/my_sd
  return(my_scored)}

Kwue_wide$Precip_scored<-ScoreVals(Kwue_wide$Precip)
Kwue_wide$Tmax_scored<-ScoreVals(Kwue_wide$Tmax)
Kwue_wide$trait_scored<-ScoreVals(Kwue_wide$trait)
myfit_kwue<-lm(SumLdoc~trait_scored+Precip_scored+Tmax_scored, data=Kwue_wide)
#               Estimate Std. Error t value  Pr(>|t|) 
#trait_scored   12.7347     0.4711   27.03   <2e-16 ***
#Precip_scored  26.5822     0.4742   56.06   <2e-16 ***
#Tmax_scored    -4.8505     0.4699  -10.32   <2e-16 ***

Amax_wide$Precip_scored<-ScoreVals(Amax_wide$Precip)
Amax_wide$Tmax_scored<-ScoreVals(Amax_wide$Tmax)
Amax_wide$trait_scored<-ScoreVals(Amax_wide$trait)
myfit_amax<-lm(SumLdoc~Precip_scored+trait_scored+Tmax_scored, data=Amax_wide)
#               Estimate Std. Error t value  Pr(>|t|) 
#Precip_scored  27.6443     0.4868  56.788  < 2e-16 ***
#trait_scored   -1.9387     0.4842  -4.004 8.19e-05 ***
#Tmax_scored    -4.3805     0.4822  -9.085  < 2e-16 ***
aw_wide$Precip_scored<-ScoreVals(aw_wide$Precip)
aw_wide$Tmax_scored<-ScoreVals(aw_wide$Tmax)
aw_wide$trait_scored<-ScoreVals(aw_wide$trait)
myfit_aw<-lm(SumLdoc~Precip_scored*trait_scored*Tmax_scored, data=aw_wide)
#               Estimate Std. Error t value  Pr(>|t|) 
#Precip_scored  24.1789     0.5684  42.536  < 2e-16 ***
#trait_scored  -10.8015     0.5660 -19.085  < 2e-16 ***
#Tmax_scored    -4.2562     0.5684  -7.487 1.03e-12 ***

Topt_wide$Precip_scored<-ScoreVals(Topt_wide$Precip)
Topt_wide$Tmax_scored<-ScoreVals(Topt_wide$Tmax)
Topt_wide$trait_scored<-ScoreVals(Topt_wide$trait)
myfit_topt<-lm(SumLdoc~Precip_scored*trait_scored*Tmax_scored, data=Topt_wide)
#               Estimate Std. Error t value  Pr(>|t|) 
#Precip_scored  28.2004     0.3261  86.473   <2e-16 ***
#trait_scored    2.8745     0.3247   8.854   <2e-16 ***
#Tmax_scored    -4.7599     0.3255 -14.624   <2e-16 ***




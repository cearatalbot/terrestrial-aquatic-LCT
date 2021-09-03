#packages
library(ggplot2) 
library(viridis)
library(reshape2)
library(gridExtra)
setwd("")
#source summary function
source("Functions/summaryTAM_MSVersion.R") #current version has 37 rows
#list traits
traitList<-c("Val1","Val2", "Val3") #var2: -10%, var3: +10%

#parameter list/names
params=c(Amax=112, Ad=0.75, Kf=0.1, Tmin=4, Tmax=40, Topt=24,
         Kvpd=0.05, PARhalf=17, k=0.58, Don=100,
         Doff=200, Lmax=3, Ka=0.006, Q10v=2,
         Kh=0.03, Q10s=2, f=0.04, Kwue=10.9, Wmax1=10, SLW=70,
         Cfrac=0.45, Kw=0.005, Aa= 1e5, Ac= 2e5, zbar= 3, 
         Cprecip = 1, deltaA = 0.01, deltaS1 = 0.37, Ks=1.2,
         lambdaS1 = 0.40, lambdaS3= 0.0001, bi = 0.3, Tstar = 50, W20 = 5, 
         Kr= 0.55, aw= 0.1, l = 1/365, lt = 0.5, r=1, Bp= 0.16, 
         Wmax2=50, Kcwd=0.01, ag=0.1, rhoS1=0.4, rhoS2=0.5, rhoS3=0.55, 
         deltaS3=0.0014, deltaS4=0.0001, deltaS2=0.0090, lambdaS2=0.04, fS1=0.6)
paramList<-names(params)
SA_base<-read.csv("SensitivityAnalysis/SASoil_base.csv", stringsAsFactors = F)

allParams<-data.frame(matrix(ncol=42, nrow=9540)) 
mycols<-c("runNum", "Tmax", "Precip", "trait", 
          "Param", "Cw", "Cr", "Cl", "Ccwd", "Cs1","Cs3", "Cdoc1", "Cdoc2", "W1", "W2",
          "Ca", "Ra", "Rr", "Rs1", "Rs3", "LCT1", "LCT2", "SumLdoc",
          "Rhdoc", "Rhdoc2", "SumRSoil", "Q2", "Q1", "GPP", "T", 
          "SumD", "NPP",  "NEE", "LCE/NEE", "Cwmin", "Cs4", "Ds1", "Ds3", 
          "Ds4", "Ls1", "Ls3", "Cs2")
colnames(allParams)<-mycols
dfChange1<-data.frame(matrix(ncol=43, nrow=0))
Climruns<-c(1:60)
#read in data
for(i in 1:51){
  fiName<-paste("SASoil",paramList[i], sep="_")
  fiName1<-paste("Out/",fiName, ".csv", sep="")
  df1<-read.csv(fiName1, stringsAsFactors = F) #assign file to generated name
  #summarize 
  sumDF<-data.frame(matrix(ncol=6, nrow=0))#store summary (36 vars) for each gradient (60) and each trait value
  colnames(sumDF)<-c("runNum", "Variable", "Value", "Tmax", "Precip", "trait")
  
  for(b in 2:length(traitList)){
    subDF<-df1[df1$trait==traitList[b],]
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
  } #end b loop
  sumDF$Value<-as.numeric(sumDF$Value)
  SA_base$Value<-as.numeric(SA_base$Value)
  df<-rbind(SA_base, sumDF)
  df$Param<-paramList[i]
  data_wide <- dcast(df, runNum + Tmax + Precip + trait + Param~ Variable, value.var = "Value", fun.aggregate = mean)
  data_wide$`LCE/NEE`<-NULL
  data_wide$SumD<-NULL
 
  assign(fiName, value=data_wide)
  if(i==1){
      allParams[1:180,]<-data_wide
    } else{
    allParams[(i*180-179):(i*180),]<-data_wide
    }
  
  #pct change
  dfChange<-data.frame(matrix(ncol=(ncol(allParams)+1), nrow=120))
  colnames(dfChange)<-c(names(allParams), "GoodCw")
  for(x in 1:60){
    subwide<-data_wide[data_wide$runNum==Climruns[x],]
    if(x==1){
        for(b in 6:42){
            dfChange[1,b]<-(((subwide[2,b]-subwide[1,b])/subwide[1,b])*100)
            dfChange[2,b]<-(((subwide[3,b]-subwide[1,b])/subwide[1,b])*100)
            dfChange[1:2,1]<-subwide[1,1]
            dfChange[1:2,2]<-subwide[1,2]
            dfChange[1:2,3]<-subwide[1,3]
            dfChange[1:2,4]<-c("Val2", "Val3")
            dfChange[1:2,43]<-ifelse(min(subwide$Cw) < 0, -999, 0)
            }}
      else {
          for(b in 6:42){
               dfChange[(x*2-1),b]<-(((subwide[2,b]-subwide[1,b])/subwide[1,b])*100) #calculate differences
               dfChange[(x*2),b]<-(((subwide[3,b]-subwide[1,b])/subwide[1,b])*100)
                 dfChange[(x*2-1):(x*2),1]<-subwide[1,1] #assign run num
                 dfChange[(x*2-1):(x*2),2]<-subwide[1,2] #assign temp
                 dfChange[(x*2-1):(x*2),3]<-subwide[1,3] #assign precip
                 dfChange[(x*2-1):(x*2),4]<-c("Val2", "Val3") #assign param change
                 dfChange[(x*2-1):(x*2),43]<-ifelse(min(subwide$Cw) < 0, -999, 0) #flag neg Cw
               }}
    
    dfChange[1:120,5]<-subwide[1,5]
  } #end x loop
  colnames(dfChange1)<-names(dfChange)
  dfChange1<-rbind(dfChange1, dfChange)
} #end i loop
dfChange2<-dfChange1[dfChange1[,41]==0,] #remove neg Cw
dfChange2$GoodCw<-NULL #can now remove flag column
change_long<-melt(dfChange2, id.vars=c("runNum", "Tmax", "Precip","trait", "Param"))
change_long$value<-round(change_long$value, digits=3)
change_long<-change_long[change_long$Precip !="25", ] #we no longer use 25 as part of our precipitation gradient

#rm aquatic parameters as they don't act on LCT
removeList<-c("Don", "Doff", "zbar", "Aa", "Ac", "Cprecip", "deltaA")
for(i in 1:length(removeList)){
change_long<-change_long[change_long$Param != removeList[i], ]
}

###########SUPPLEMENTARY FIGS#############
#set wd to store plots
setwd("")
#colors for group plots
GroupCols<-viridis(3, alpha = 1, begin = 0, end=0.8, option = "D", direction=-1)

#theme
borderTheme0.5<-theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
                      panel.background = element_rect(colour = "black", size=1, fill=NA), 
                      strip.background = element_blank(), 
                      panel.spacing = unit(0.3, "lines"), axis.ticks.length = unit(0.2, "lines"), 
                      text = element_text(size = 14),legend.background=element_blank(), 
                      legend.key = element_rect(fill = NA), axis.text.x=element_text(angle = -45, hjust = 0))


#assign params to a group
ParamGroup<-data.frame(matrix(ncol=3, nrow=44))
colnames(ParamGroup)<-c("Param", "reName","Group")
ParamGroup$Param<-c("Amax", "Ad", "Kf", "Tmin", "Tmax", "Topt",
                    "Kvpd", "PARhalf", "k", "Lmax", "Ka", "Q10v",
                    "Kh", "Q10s", "f", "Kwue"," Wmax1", "SLW",
                    "Cfrac", "Kw", "deltaS1", "Ks",
                    "lambdaS1", "lambdaS3", "bi", "Tstar", "W20", 
                    "Kr", "aw", "l", "lt", "r", "Bp", 
                    "Wmax2", "Kcwd", "ag", "rhoS1", "rhoS2", "rhoS3", 
                    "deltaS3", "deltaS4", "deltaS2"," lambdaS2", "fS1") 
ParamGroup$reName<-c("A[max]", "A[d]", "K[f]", "T[min]", "T[max]", "T[opt]", 
                     "K[vpd]", "PAR[1/2]", "k", "L[max]",
                     "K[a]", "Q[10,v]", "K[h]", "Q[10,s]", "f",  "K[wue]",   
                     "W[max,1]", "SLW", "C[frac]", "K[w]",
                     "delta[S1]",   
                     "K[s]", "lambda[S1]", "lambda[S3]", "b[i]", "T*", "W[2]^o", "K[r]",
                     "a[w]", "l", "l[t]", "r", "B[p]", "W[max,2]", "K[cwd]",  "a[g]", 
                     "rho[S1]","rho[S2]","rho[S3]",    
                     "deltaS3", "deltaS4", "deltaS2", "lambdaS2", "fS1")
ParamGroup$Group<-c("Plant", "Plant", "Plant", "Plant", "Plant", "Plant", 
                    "Plant", "Plant", "Plant", "Plant", "Plant", "Plant", 
                    "Soil", "Soil", "Soil", "Plant",
                    "Water", "Plant", "Plant", "Plant", "Soil", 
                    "Water", "Soil", "Soil", "Water", "Water", "Soil",
                    "Plant", "Plant", "Plant", "Plant", "Water", "Water", "Water",
                    "Soil", "Plant", "Soil", "Soil", "Soil", "Soil", "Soil", 
                    "Soil", "Soil", "Soil")

varList<-c("Ddoc", "Qbdoc", "NPP", "Cd1", "Cd2", "Co", "Cf2","Cm", "Cp")  
cornerLab<-c("A.", "B.", "C.",
             "D.", "E.", "F.",
             "G.", "H.", "I.")
df1<-data.frame(matrix(ncol=ncol(change_long), nrow=0))

for(i in 1:length(varList)){
df<-change_long[change_long$variable==varList[i], ]
df1<-rbind(df1, df)
} #df1 is the change data frame for only the variables i need
df1$value_abs<-abs(df1$value)
dfTop<-data.frame(matrix(ncol=2, nrow=0)) #df to store top 10 params for each variable 
colnames(dfTop)<-c("variable", "TopParam")
for(i in 1:length(varList)){
  VarInt<-df1[df1$variable==varList[i],]
  for(x in 1:nrow(VarInt)){
    #VarInt$Param_f<-factor(VarInt$Param, levels=VarOrder)
    VarInt$Group[x]<-ParamGroup[ParamGroup$Param==VarInt$Param[x], 3]
    VarInt$reName[x]<-ParamGroup[ParamGroup$Param==VarInt$Param[x], 2]
  }
  VarInt2<-aggregate(VarInt$value_abs, by=list(VarInt$reName), FUN=median)
  colnames(VarInt2)<-c("Param", "median")
  VarOrder<-VarInt2[order(VarInt2$median, decreasing=T),]
  VarOrder<-as.character(unique(VarOrder$Param)) 
    VarInta<-data.frame(matrix(ncol=6, nrow=0))
    for(y in 1:10){
      VarIntb<-VarInt[VarInt$reName==VarOrder[y], ]
      VarInta<-rbind(VarInta, VarIntb)
    }

VarInta$ParamNew_f<-factor(VarInta$reName, levels=VarOrder[1:10])
VarInta$Group_f<-factor(VarInta$Group, levels=c("Plant", "Water", "Soil"))

#loop to make plots
if (i==1){
gb<-ggplot(VarInta, aes(x=ParamNew_f, y=value_abs))+geom_boxplot(aes(colour=Group_f), position=position_dodge(0.4), width=0.6)+
  scale_y_continuous(limits=c(0, 20), breaks=c(seq(0, 20, 5)))+
  labs(y=" Mean absolute \n change (%)", x="Parameter")+
  scale_x_discrete(labels=c(expression("K"[wue]),expression(lambda[S1]), expression("f"[S1]), expression("W"[MAX1]), expression("Q"[10]["s"]), expression("K"[h]),expression("T"[opt]),expression("K"[f]), expression("f"), expression("W"[2]^o)))+
  scale_color_manual(values=c(GroupCols), " Parameter \n type")+
  annotate("text", y=19, x=0.8, label=cornerLab[i], size=4.5)
  gb<-gb+borderTheme0.5
print(gb)
fiName<-paste(varList[i], "Box",".png", sep="")
#ggsave(filename = fiName, plot=gb, width = 6, height = 2.9, units= "in", device='png', dpi=320)
} else if (i==2){
  gb<-ggplot(VarInta, aes(x=ParamNew_f, y=value_abs))+geom_boxplot(aes(colour=Group_f), position=position_dodge(0.4), width=0.6)+
    scale_y_continuous(limits=c(0, 20), breaks=c(seq(0, 20, 5)))+
    labs(y=" Mean absolute \n change (%)", x="Parameter")+
    scale_x_discrete(labels=c(expression("K"[h]),expression("K"[wue]), expression(lambda[S1]), expression("f"[S1]), expression("Q"[10]["s"]), expression("B"[p]),expression("b"[i]),expression("K"[s]), expression("T"[opt]), expression("T"^"*")))+
    scale_color_manual(values=c(GroupCols), " Parameter \n type")+
    annotate("text", y=19, x=0.8, label=cornerLab[i], size=4.5)
  gb<-gb+borderTheme0.5
  print(gb)
  fiName<-paste(varList[i], "Box",".png", sep="")
  #ggsave(filename = fiName, plot=gb, width = 6, height = 2.9, units= "in", device='png', dpi=320)
} else if (i==3){
  gb<-ggplot(VarInta, aes(x=ParamNew_f, y=value_abs))+geom_boxplot(aes(colour=Group_f), position=position_dodge(0.4), width=0.6)+
    scale_y_continuous(limits=c(0, 20), breaks=c(seq(0, 20, 5)))+
    labs(y=" Mean absolute \n change (%)", x="Parameter")+
    scale_x_discrete(labels=c(expression("K"[wue]),expression("T"[opt]), expression("K"[f]), expression("f"), expression("W"[2]^o), expression("B"[p]),expression("L"[max]),expression("W"[MAX1]), expression("b"[i]), expression("Q"[10]["v"])))+
    scale_color_manual(values=c(GroupCols), " Parameter \n type")+
    annotate("text", y=19, x=0.8, label=cornerLab[i], size=4.5)
  gb<-gb+borderTheme0.5
  print(gb)
  fiName<-paste(varList[i], "Box",".png", sep="")
  #ggsave(filename = fiName, plot=gb, width = 6, height = 2.9, units= "in", device='png', dpi=320)
  
} else if (i==4){
  gb<-ggplot(VarInta, aes(x=ParamNew_f, y=value_abs))+geom_boxplot(aes(colour=Group_f), position=position_dodge(0.4), width=0.6)+
    scale_y_continuous(limits=c(0, 20), breaks=c(seq(0, 20, 5)))+
    labs(y=" Mean absolute \n change (%)", x="Parameter")+
    scale_x_discrete(labels=c(expression("K"[wue]),expression(lambda[S1]), expression("f"[S1]), expression("W"[MAX1]), expression("K"[h]), expression("T"[opt]),expression("Q"[10]["s"]),expression("K"[f]), expression("f"), expression("W"[2]^o)))+
    scale_color_manual(values=c(GroupCols), " Parameter \n type")+
    annotate("text", y=19, x=0.8, label=cornerLab[i], size=4.5)
  gb<-gb+borderTheme0.5
  print(gb)
  fiName<-paste(varList[i], "Box",".png", sep="")
  #ggsave(filename = fiName, plot=gb, width = 6, height = 2.9, units= "in", device='png', dpi=320)
  
} else if (i==5){
  gb<-ggplot(VarInta, aes(x=ParamNew_f, y=value_abs))+geom_boxplot(aes(colour=Group_f), position=position_dodge(0.4), width=0.6)+
    scale_y_continuous(limits=c(0, 20), breaks=c(seq(0, 20, 5)))+
    labs(y=" Mean absolute \n change (%)", x="Parameter")+
    scale_x_discrete(labels=c(expression("K"[h]),expression("Q"[10]["s"]), expression("K"[wue]), expression(lambda[S1]), expression("f"[S1]), expression("T"[opt]),expression("B"[p]),expression("b"[i]), expression("W"[2]^o), expression("K"[s])))+
    scale_color_manual(values=c(GroupCols), " Parameter \n type")+
    annotate("text", y=19, x=0.8, label=cornerLab[i], size=4.5)
  gb<-gb+borderTheme0.5
  print(gb)
  fiName<-paste(varList[i], "Box",".png", sep="")
  #ggsave(filename = fiName, plot=gb, width = 6, height = 2.9, units= "in", device='png', dpi=320)
  
} else if(i==6){
  gb<-ggplot(VarInta, aes(x=ParamNew_f, y=value_abs))+geom_boxplot(aes(colour=Group_f), position=position_dodge(0.4), width=0.6)+
    scale_y_continuous(limits=c(0, 20), breaks=c(seq(0, 20, 5)))+
    labs(y=" Mean absolute \n change (%)", x="Parameter")+
    scale_x_discrete(labels=c(expression("K"[wue]),expression(delta[S1]), expression("f"[S1]), expression("T"[opt]), expression("K"[f]), expression("f"),expression("W"[2]^o),expression("a"[w]),expression("B"[p]), expression("W"[MAX1])))+
    scale_color_manual(values=c(GroupCols), " Parameter \n type")+
    annotate("text", y=19, x=0.8, label=cornerLab[i], size=4.5)
  gb<-gb+borderTheme0.5
  print(gb)
  fiName<-paste(varList[i], "Box",".png", sep="")
  #ggsave(filename = fiName, plot=gb, width = 6, height = 2.9, units= "in", device='png', dpi=320)
  
} else if (i==7){
  gb<-ggplot(VarInta, aes(x=ParamNew_f, y=value_abs))+geom_boxplot(aes(colour=Group_f), position=position_dodge(0.4), width=0.6)+
    scale_y_continuous(limits=c(0, 20), breaks=c(seq(0, 20, 5)))+
    labs(y=" Mean absolute \n change (%)", x="Parameter")+
    scale_x_discrete(labels=c(expression("f"[S1]),expression("K"[wue]), expression(delta[S2]), expression("T"[opt]), expression("K"[f]), expression("f"),expression("W"[2]^o),expression("L"[max]),expression("B"[p]), expression("W"[MAX1])))+
    scale_color_manual(values=c(GroupCols), " Parameter \n type")+
    annotate("text", y=19, x=0.8, label=cornerLab[i], size=4.5)
  gb<-gb+borderTheme0.5
  print(gb)
  fiName<-paste(varList[i], "Box",".png", sep="")
  #ggsave(filename = fiName, plot=gb, width = 6, height = 2.9, units= "in", device='png', dpi=320)
  
} else if (i==8){
  gb<-ggplot(VarInta, aes(x=ParamNew_f, y=value_abs))+geom_boxplot(aes(colour=Group_f), position=position_dodge(0.4), width=0.6)+
    scale_y_continuous(limits=c(0, 20), breaks=c(seq(0, 20, 5)))+
    labs(y=" Mean absolute \n change (%)", x="Parameter")+
    scale_x_discrete(labels=c(expression("Q"[10]["s"]),expression("K"[wue]), expression(delta[S3]), expression(rho[S2]), expression(rho[S1]), expression("T"[opt]),expression(lambda[S1]),expression("B"[p]),expression("K"[f]), expression("f")))+
    scale_color_manual(values=c(GroupCols), " Parameter \n type")+
    annotate("text", y=19, x=0.8, label=cornerLab[i], size=4.5)
  gb<-gb+borderTheme0.5
  print(gb)
  fiName<-paste(varList[i], "Box",".png", sep="")
  #ggsave(filename = fiName, plot=gb, width = 6, height = 2.9, units= "in", device='png', dpi=320)
  
} else if(i==9){
  gb<-ggplot(VarInta, aes(x=ParamNew_f, y=value_abs))+geom_boxplot(aes(colour=Group_f), position=position_dodge(0.4), width=0.6)+
    scale_y_continuous(limits=c(0, 20), breaks=c(seq(0, 20, 5)))+
    labs(y=" Mean absolute \n change (%)", x="Parameter")+
    scale_x_discrete(labels=c(expression("Q"[10]["s"]),expression(delta[S4]), expression("K"[wue]), expression(rho[S2]), expression(rho[S1]), expression("T"[opt]),expression("f"),expression(lambda[S1]),expression("K"[f]), expression("B"[p])))+
    scale_color_manual(values=c(GroupCols), " Parameter \n type")+
    annotate("text", y=19, x=0.8, label=cornerLab[i], size=4.5)
  gb<-gb+borderTheme0.5
  print(gb)
  fiName<-paste(varList[i], "Box",".png", sep="")
  #ggsave(filename = fiName, plot=gb, width = 6, height = 2.9, units= "in", device='png', dpi=320)
}
dfTop1<-data.frame(matrix(ncol=2, nrow=15))
colnames(dfTop1)<-c("variable", "TopParam")
dfTop1$variable<-varList[i]
dfTop1$TopParam<-VarOrder[1:15]
dfTop<-rbind(dfTop, dfTop1)
}

###############FIGURES IN MANUSCRIPT#########
###plots with percent change for LCT
lctList<-c("SumLdoc")
lctParams<-c("Amax", "Topt","aw", "Kwue")
lctLabel<-c("A", "B", "C", "D")
for(i in 1:length(lctList)){
  lctDF1<-change_long[change_long$variable==lctList[i],]
  
  for(x in 1:length(lctParams)){
    lctDF<-lctDF1[lctDF1$Param==lctParams[x],]
    if(x==1){
      gA<-ggplot(lctDF, aes(x=Precip, y=value, colour=Tmax, shape=trait))+geom_point(size=2)+
        scale_shape_manual(values=c(16, 17), labels=c("-10%", "+10%")," ")+
        scale_y_continuous(limits=c(-24, 24), breaks=c(seq(-20, 20, 10)))+
        scale_colour_viridis_c("Tmax (C)", option="B", limits=c(12, 27), end=0.9)+
        labs(y="Change (%)", x=expression('Precipitation'~(cm~yr^-1)))+
        borderTheme0+annotate("text", x=56, y=23.2, label= paste(lctLabel[x], ".", sep=""), size=4.5)
      fiName<-paste("CSI", lctList[i], lctParams[x], ".png", sep="")
      #ggsave(filename = fiName, plot=gA, width = 3.5, height = 2.5, units= "in", device='png', dpi=320)
      
    } else if(x==2){
      gB<-ggplot(lctDF, aes(x=Precip, y=value, colour=Tmax, shape=trait))+geom_point(size=2)+
        scale_shape_manual(values=c(16, 17), labels=c("-10%", "+10%")," ")+
        scale_y_continuous(limits=c(-24, 24), breaks=c(seq(-20, 20, 10)))+
        scale_colour_viridis_c("Tmax (C)", option="B", limits=c(12, 27), end=0.9)+
        labs(y=NULL, x=expression('Precipitation'~(cm~yr^-1)))+
        borderTheme0+annotate("text", x=56, y=23.2, label= paste(lctLabel[x], ".", sep=""), size=4.5)
      fiName<-paste("CSI", lctList[i], lctParams[x], ".png", sep="")
      #ggsave(filename = fiName, plot=gB, width = 3.5, height = 2.5, units= "in", device='png', dpi=320)
      
    } else if(x==3){
      gC<-ggplot(lctDF, aes(x=Precip, y=value, colour=Tmax, shape=trait))+geom_point(size=2)+
        scale_shape_manual(values=c(16, 17), labels=c("-10%", "+10%")," ")+
        scale_y_continuous(limits=c(-24, 24), breaks=c(seq(-20, 20, 10)))+
        scale_colour_viridis_c("Tmax (C)", option="B", limits=c(12, 27), end=0.9)+
        labs(y="Change (%)", x=expression('Precipitation'~(cm~yr^-1)))+
        borderTheme0+annotate("text", x=56, y=23.2, label= paste(lctLabel[x], ".", sep=""), size=4.5)
      fiName<-paste("CSI", lctList[i], lctParams[x], ".png", sep="")
      #ggsave(filename = fiName, plot=gC, width = 3.5, height = 2.5, units= "in", device='png', dpi=320)
    } else if(x==4){
      gD<-ggplot(lctDF, aes(x=Precip, y=value, colour=Tmax, shape=trait))+geom_point(size=2)+
        scale_shape_manual(values=c(16, 17), labels=c("-10%", "+10%")," Parameter \n change")+
        scale_y_continuous(limits=c(-24, 24), breaks=c(seq(-20, 20, 10)))+
        scale_colour_viridis_c("Tmax (C)", option="B", limits=c(12, 27), end=0.9)+
        labs(y=NULL, x=expression('Precipitation'~(cm~yr^-1)))+
        borderTheme0.5+annotate("text", x=56, y=23.2, label= paste(lctLabel[x], ".", sep=""), size=4.5)
      fiName<-paste("CSI", lctList[i], lctParams[x], ".png", sep="")
      #ggsave(filename = fiName, plot=gD, width = 3.5, height = 2.5, units= "in", device='png', dpi=320)
    }}
}

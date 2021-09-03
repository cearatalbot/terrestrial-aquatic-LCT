####Plotting base parameter runs over 60 Precip/Temp
####combindations########
#load packages
library(ggplot2)
library(viridis)  
library(gridExtra)
library(reshape2)

#Temp vs. Precip experiment
setwd("") #set to location of terrestrial-aquatic-LCT folder

#read in annual summary model simulations
df<-read.csv("SensitivityAnalysis/SASoil_base.csv", stringsAsFactors = F)
df$Value<-as.numeric(df$Value)
data_wide <- dcast(df, runNum + Tmax + Precip + trait ~ Variable, value.var = "Value", fun.aggregate = mean)
data_wide<-data_wide[(data_wide$Cw > 0) & (data_wide$NPP > 0), ] #remove unrealistic negative NPP and negative wood C

#my ggplot2 theme
borderTheme0.5<-theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
                      panel.background = element_rect(colour = "black", size=1, fill=NA), 
                      strip.background = element_blank(), 
                      panel.spacing = unit(0.3, "lines"), axis.ticks.length = unit(0.2, "lines"), 
                      text = element_text(size = 12),legend.background=element_blank(), 
                      legend.key = element_rect(fill = NA), aspect.ratio = 1)
#create a color palette
tempCols<-mycols<-viridis(6, alpha = 0.8, begin = 0, end = 0.9, option = "B", direction=1)

#LCT vs. precip, colored by temp
data_wide$Tmax_f<-factor(data_wide$Tmax)
p1<-ggplot(data_wide, aes(x=Precip, y=SumLdoc, group=Tmax_f, colour=Tmax_f))+#geom_boxplot()+
  geom_point(size=1.5)+geom_line(size=1.5)+scale_color_manual(values=tempCols, "Tmax (C)")+
  labs(x=expression('Precipitation'~(cm~yr^-1)), y=expression('LCT'~(g~C~m^-2~yr^-1)))+
  borderTheme0.5+scale_y_continuous(limits=c(0,100), breaks=seq(0,100, 25), expand = c(0, 0))+
  annotate("text", x=59, y=93, label= "A.", size=4.5)
p1 
#ggsave(filename = "TmaxLine_NewS.png", plot=p1, width = 3.5, height = 2.5, units= "in", device='png', dpi=320)

#water export vs. precip
p2<-ggplot(data_wide, aes(x=Precip, y=SumD, group=Tmax_f, colour=Tmax_f))+#geom_boxplot()+
  geom_point(size=1.5)+geom_line(size=1.5)+scale_color_manual(values=tempCols, "Tmax (C)")+
  labs(y=expression('Water drainage'~(cm~yr^-1)), x=expression('Precipitation'~(cm~yr^-1)))+
  borderTheme0.5+#scale_y_continuous(limits=c(0,155), breaks=seq(0,150, 50), expand = c(0, 0))+
  scale_y_continuous(limits=c(0,205), breaks=seq(0,200,50), expand=c(0,0))+
  annotate("text", x=59, y=191, label= "B.", size=4.5)
p2
#ggsave(filename = "drainagePrecip_NewS.png", plot=p2, width = 3.5, height = 2.5, units= "in", device='png', dpi=320)

#NPP vs. precip
p3<-ggplot(data_wide, aes(x=Precip, y=NPP, group=Tmax_f, colour=Tmax_f))+#geom_boxplot()+
  geom_point(size=1.5)+geom_line(size=1.5)+scale_color_manual(values=tempCols, "Tmax (C)")+
  labs(y=expression('NPP'~(g~C~m^-2~yr^-1)), x=expression('Precipitation'~(cm~yr^-1)))+
  borderTheme0.5+
  scale_y_continuous(limits=c(0,1150), breaks=seq(0,1000,200), expand=c(0,0))+
  annotate("text", x=59, y=1050, label= "C.", size=4.5)
p3
#ggsave(filename = "NPPprecip_NewS.png", plot=p3, width = 3.5, height = 2.5, units= "in", device='png', dpi=320)

#Soil respiration vs. precip
p4<-ggplot(data_wide, aes(x=Precip, y=SumRSoil, group=Tmax_f, colour=Tmax_f))+#geom_boxplot()+
  geom_point(size=1.5)+geom_line(size=1.5)+scale_color_manual(values=tempCols, "Tmax (C)")+
  labs(y=expression('Soil respiration'~(g~C~m^-2~yr^-1)), x=expression('Precipitation'~(cm~yr^-1)))+
  borderTheme0.5+#scale_y_continuous(limits=c(0,155), breaks=seq(0,150, 50), expand = c(0, 0))+
  scale_y_continuous(limits=c(0,1050), breaks=seq(0,1000,200), expand=c(0,0))+
  annotate("text", x=59, y=975, label= "D.", size=4.5)
p4
#ggsave(filename = "soilprecip_NewS.png", plot=p4, width = 3.6, height = 2.6, units= "in", device='png', dpi=320)

#Soil C vs. precip
data_wide$CsSum<-data_wide$Cs1 + data_wide$Cs2 + data_wide$Cs3 + data_wide$Cs4
p5<-ggplot(data_wide, aes(x=Precip, y=CsSum, group=Tmax_f, colour=Tmax_f))+
  geom_point(size=1.5)+geom_line(size=1.5)+scale_color_manual(values=tempCols, "Tmax (C)")+
  labs(y=expression('Total soil C'~(g~C~m^-2)), x=expression('Precipitation'~(cm~yr^-1)))+
  borderTheme0.5+
  scale_y_continuous(limits=c(0,5000), breaks=seq(0,5000,1000), expand=c(0,0))+
  annotate("text", x=59, y=4650, label= "E.", size=4.5)
p5
#ggsave(filename = "CsPrecip_NewS.png", plot=p5, width = 3.5, height = 2.5, units= "in", device='png', dpi=320)

#soil DOC vs. precip
data_wide$CdSum<-data_wide$Cdoc1 + data_wide$Cdoc2
p6<-ggplot(data_wide, aes(x=Precip, y=CdSum, group=Tmax_f, colour=Tmax_f))+#geom_boxplot()+
  geom_point(size=1.5)+geom_line(size=1.5)+scale_color_manual(values=tempCols, "Tmax (C)")+
  labs(y=expression('Total soil DOC'~(g~C~m^-2)), x=expression('Precipitation'~(cm~yr^-1)))+
  borderTheme0.5+#scale_y_continuous(limits=c(0,155), breaks=seq(0,150, 50), expand = c(0, 0))+
  scale_y_continuous(limits=c(0,8.5), breaks=seq(0,8,2), expand=c(0,0))+
  annotate("text", x=59, y=7.8, label= "F.", size=4.5)
p6
#ggsave(filename = "CdPrecip_NewS.png", plot=p6, width = 3.5, height = 2.5, units= "in", device='png', dpi=320)


###Plots for only two precips at all temps
#make color palette
mycols2<-viridis(2, alpha = 0.8, begin = 0, end = 0.8, option = "D", direction=1)
data_wide$Precip_f<-factor(data_wide$Precip)
#subset for two precips
SubPrecip<-rbind(data_wide[data_wide$Precip==75,], data_wide[data_wide$Precip==250,])

#LCT vs temp
p1a<-ggplot(SubPrecip, aes(x=Tmax, y=SumLdoc, group=Precip_f, colour=Precip_f))+#geom_boxplot()+
  geom_point(size=1.5)+geom_line(size=1.5)+scale_color_manual(values=mycols2, "Precipitation")+
  labs(x=expression('Tmax'~(C)), y=expression('LCT'~(g~C~m^-2~yr^-1)))+
  borderTheme0.5+scale_y_continuous(limits=c(0,100), breaks=seq(0,100, 25), expand = c(0, 0))+
  scale_x_continuous(limits=c(11,28), breaks=seq(12,27, 3), expand = c(0, 0))+
  annotate("text", x=12.6, y=93, label= "A.", size=4.5)
p1a 
#ggsave(filename = "TmaxLine2a_NewS.png", plot=p1a, width = 3.5, height = 2.5, units= "in", device='png', dpi=320)

#Drainage vs temp
p2a<-ggplot(SubPrecip, aes(x=Tmax, y=SumD, group=Precip_f, colour=Precip_f))+#geom_boxplot()+
  geom_point(size=1.5)+geom_line(size=1.5)+scale_color_manual(values=mycols2, "Precipitation")+
  labs(y=expression('Water drainage'~(cm~yr^-1)), x=expression('Tmax'~(C)))+
  borderTheme0.5+
  scale_x_continuous(limits=c(11,28), breaks=seq(12,27, 3), expand = c(0, 0))+
  scale_y_continuous(limits=c(0,240), breaks=seq(0,200,50), expand=c(0,0))+
  annotate("text", x=12.6, y=222, label= "B.", size=4.5)
p2a
#ggsave(filename = "drainagePrecip2_NewS.png", plot=p2a, width = 3.5, height = 2.5, units= "in", device='png', dpi=320)

#Soil DOC vs temp
p3a<-ggplot(SubPrecip, aes(x=Tmax, y=CdSum, group=Precip_f, colour=Precip_f))+#geom_boxplot()+
  geom_point(size=1.5)+geom_line(size=1.5)+scale_color_manual(values=mycols2, "Precipitation")+
  labs(y=expression('Total soil DOC'~(g~C~m^-2)), x=expression('Tmax'~(C)))+
  borderTheme0.5+
  scale_x_continuous(limits=c(11,28), breaks=seq(12,27, 3), expand = c(0, 0))+
  scale_y_continuous(limits=c(0,9), breaks=seq(0,8,2), expand=c(0,0))+
  annotate("text", x=12.6, y=8.2, label= "C.", size=4.5)
p3a
#ggsave(filename = "MeanCdPrecip2_NewS.png", plot=p3a, width = 3.5, height = 2.5, units= "in", device='png', dpi=320)

#Soil C vs temp
p4a<-ggplot(SubPrecip, aes(x=Tmax, y=CsSum, group=Precip_f, colour=Precip_f))+#geom_boxplot()+
  geom_point(size=1.5)+geom_line(size=1.5)+scale_color_manual(values=mycols2, "Precipitation")+
  labs(y=expression('Total soil C'~(g~C~m^-2)), x=expression('Tmax'~(C)))+
  borderTheme0.5+
  scale_x_continuous(limits=c(11,28), breaks=seq(12,27, 3), expand = c(0, 0))+
  scale_y_continuous(limits=c(0,5300), breaks=seq(0,5000, 1000), expand=c(0,0))+
  annotate("text", x=12.6, y=4920, label= "D.", size=4.5)
p4a
#ggsave(filename = "MeanCsPrecip2_NewS.png", plot=p4a, width = 3.5, height = 2.5, units= "in", device='png', dpi=320)

#NPP vs Temp
p5a<-ggplot(SubPrecip, aes(x=Tmax, y=NPP, group=Precip_f, colour=Precip_f))+#geom_boxplot()+
  geom_point(size=1.5)+geom_line(size=1.5)+scale_color_manual(values=mycols2, "Precipitation")+
  labs(y=expression('NPP'~(g~C~m^-2~yr^-1)), x=expression('Tmax'~(C)))+
  borderTheme0.5+
  scale_y_continuous(limits=c(0,1300), breaks=seq(0,1200,200), expand=c(0,0))+
  scale_x_continuous(limits=c(11,28), breaks=seq(12,27, 3), expand = c(0, 0))+
  annotate("text", x=12.7, y=1200, label= "E.", size=4.5)
p5a
#ggsave(filename = "NPPprecip2_NewS.png", plot=p5a, width = 3.5, height = 2.5, units= "in", device='png', dpi=320)

#soil respiration vs. temp
p6a<-ggplot(SubPrecip, aes(x=Tmax, y=SumRSoil, group=Precip_f, colour=Precip_f))+#geom_boxplot()+
  geom_point(size=1.5)+geom_line(size=1.5)+scale_color_manual(values=mycols2, "Precipitation")+
  labs(y=expression('Soil respiration'~(g~C~m^-2~yr^-1)), x=expression('Tmax'~(C)))+
  borderTheme0.5+
  scale_x_continuous(limits=c(11,28), breaks=seq(12,27, 3), expand = c(0, 0))+
  scale_y_continuous(limits=c(0,1070), breaks=seq(0,1000,200), expand=c(0,0))+
  annotate("text", x=12.8, y=975, label= "F.", size=4.5)
p6a
#ggsave(filename = "soilprecip2_NewS.png", plot=p6a, width = 3.6, height = 2.6, units= "in", device='png', dpi=320)


#extras for supplement
#GPP vs. precip
p1c<-ggplot(data_wide, aes(x=Precip, y=GPP, group=Tmax_f, colour=Tmax_f))+#geom_boxplot()+
  geom_point(size=1.5)+geom_line(size=1.5)+scale_color_manual(values=tempCols, "Tmax (C)")+
  labs(y=expression('GPP'~(g~C~m^-2~yr^-1)), x=expression('Precipitation'~(cm~yr^-1)))+
  borderTheme0.5+
  scale_y_continuous(limits=c(400,2050), breaks=seq(400,2000,400), expand=c(0,0))+
  annotate("text", x=59, y=1900, label= "A.", size=4.5)
p1c
#ggsave(filename = "GPPprecip_NewS.png", plot=p1c, width = 3.6, height = 2.6, units= "in", device='png', dpi=320)

#TER vs precip
data_wide$TER<-(data_wide$GPP+data_wide$NEE) #calculate ecosystem respiration
p2c<-ggplot(data_wide, aes(x=Precip, y=TER, group=Tmax_f, colour=Tmax_f))+#geom_boxplot()+
  geom_point(size=1.5)+geom_line(size=1.5)+scale_color_manual(values=tempCols, "Tmax (C)")+
  labs(y=expression('ER'~(g~C~m^-2~yr^-1)), x=expression('Precipitation'~(cm~yr^-1)))+
  borderTheme0.5+
  scale_y_continuous(limits=c(300,1900), breaks=seq(300,1800,300), expand=c(0,0))+
  annotate("text", x=59, y=1795, label= "B.", size=4.5)
p2c
#ggsave(filename = "ERprecip_NewS.png", plot=p2c, width = 3.6, height = 2.6, units= "in", device='png', dpi=320)


########linear regressions#####
data_wide$WSum<-data_wide$W1 + data_wide$W2 #calculate soil water
data_wide$DOCconc<-data_wide$CdSum/(data_wide$WSum*0.01) #calculate soil DOC concentration

##LCT-precip regression
myfit<-lm(SumLdoc~Precip, data=data_wide)

##standardized effect size of intermediate drivers
#function to calculate z scores
ScoreVals<-function(Vals){
  my_sd<- sd(Vals)*sqrt((length(Vals)-1)/(length(Vals)))
  my_mean <- mean(Vals)
  my_scored<-(Vals-my_mean)/my_sd
  return(my_scored)}

data_wide$CdSum_scored<-ScoreVals(data_wide$CdSum)
data_wide$DOCconc<-ScoreVals(data_wide$DOCconc)
data_wide$WSum_scored<-ScoreVals(data_wide$WSum)
data_wide$SumD_scored<-ScoreVals(data_wide$SumD)
data_wide$NPP_scored<-ScoreVals(data_wide$NPP)
data_wide$SoilResp_scored<-ScoreVals(data_wide$SumRSoil)
data_wide$Precip_scored<-ScoreVals(data_wide$Precip)
data_wide$Tmax_scored<-ScoreVals(data_wide$Tmax)

#standardized effect size (slope) for competing intermediate drivers
myfit<-lm(SumLdoc~CdSum_scored+WSum_scored+SumD_scored,data=data_wide)

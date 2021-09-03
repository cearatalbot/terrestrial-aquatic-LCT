####Function to summarize daily model outputs as annual sums (flux) or means (pool)####
#######################################################################################
library(reshape2)
#df is the daily model output

summaryTAM<-function(df){
    sumStore1<-data.frame(matrix(ncol=2, nrow=37))
    colnames(sumStore1)<-c("Variable", "Value")
    
    sumStore1[,1]<-c("Cw", "Cr", "Cl", "Ccwd", "Cs1","Cs3", "Cdoc1", "Cdoc2", "W1", "W2",
                     "Ca", "Ra", "Rr", "Rs1", "Rs3", "LCT1", "LCT2", "SumLdoc",
                     "Rhdoc", "Rhdoc2", "SumRSoil", "Q2", "Q1", "GPP", "T", 
                     "SumD", "NPP",  "NEE", "LCE/NEE", "Cwmin", "Cs4", "Ds1", "Ds3", 
                     "Ds4", "Ls1", "Ls3", "Cs2")
    sumStore1[1,2]<-mean(df$Cw)
    sumStore1[2,2]<-mean(df$Cr)
    sumStore1[3,2]<-mean(df$Cl)
    sumStore1[4,2]<-mean(df$Ccwd)
    sumStore1[5,2]<-mean(df$Cs1)
    sumStore1[6,2]<-mean(df$Cs3)
    sumStore1[7,2]<-mean(df$Cdoc1)
    sumStore1[8,2]<-mean(df$Cdoc2)
    sumStore1[9,2]<-mean(df$W1)
    sumStore1[10,2]<-mean(df$W2)
    sumStore1[11,2]<-mean(df$Ca)
    sumStore1[12,2]<-sum(df$Ra)
    sumStore1[13,2]<-sum(df$Rr)
    sumStore1[14,2]<-sum(df$Rs1)
    sumStore1[15,2]<-sum(df$Rs3)
    sumStore1[16,2]<-sum(df$LCT1)
    sumStore1[17,2]<-sum(df$LCT2)
    sumStore1[18,2]<-sum(df$LCT1)+sum(df$LCT2)
    sumStore1[19,2]<-sum(df$Rhdoc) 
    sumStore1[20,2]<-sum(df$Rhdoc2)
    sumStore1[21,2]<-sum(df$Rs1)+sum(df$Rs2)+sum(df$Rs3)+sum(df$Rhdoc)+sum(df$Rhdoc2)+sum(df$Ds4)
    sumStore1[22,2]<-sum(df$Q2)
    sumStore1[23,2]<-sum(df$Q1)
    sumStore1[24,2]<-sum(df$GPP)
    sumStore1[25,2]<-sum(df$T)
    sumStore1[26,2]<-sum(df$Q2)+sum(df$Q1)
    sumStore1[27,2]<-sum(df$NPP)
    sumStore1[28,2]<-ifelse(sum(df$NPP) >=0, (sum(df$Rs1)+sum(df$Rs2)+sum(df$Rhdoc)+sum(df$Rhdoc2)+sum(df$Rs3)+sum(df$Ds4))-sum(df$NPP), NA)
    sumStore1[29,2]<-(sum(df$LCT2)+sum(df$LCT1))/((sum(df$Rs1)+sum(df$Rs2)+sum(df$Rhdoc)+sum(df$Rhdoc2)+sum(df$Rs3)+sum(df$Ds4))-sum(df$NPP))
    sumStore1[30,2]<-min(df$Cw)
    sumStore1[31,2]<-mean(df$Cs4)
    sumStore1[32,2]<-sum(df$Ds1)
    sumStore1[33,2]<-sum(df$Ds3)
    sumStore1[34,2]<-sum(df$Ds4)
    sumStore1[35,2]<-sum(df$Ls1)
    sumStore1[36,2]<-sum(df$Ls3)
    sumStore1[37,2]<-mean(df$Cs2)
    return(sumStore1)} 

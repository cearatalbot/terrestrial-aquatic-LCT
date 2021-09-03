########Function to calculate leaf on and leaf off day########
#Usage: Tmean=vector of mean monthly temperature [C] for one year 
#(days within the same month have the same Temp), baseT_on= base temperature
#for growing degree days, latit=latitude in [decimal degrees] 
#leaf on is calculated with growing degree days & leaf off is calculated
#with daylight
#output is a vector of leaf off on day and leaf off day for the input year
gddSum <- function(Tmean, baseT_on, latit) {
  leafOnDF<-data.frame(matrix(ncol=5, nrow=365)) #df to store data
  colnames(leafOnDF)<-c("DOY", "GDD", "sumGDD", "photoOff") 
  leafOnDF$DOY<-1:365 
  #calculate Day On
  for(i in 1:365){
    a <- Tmean[i]-baseT_on
    leafOnDF[i,2]<-ifelse(a >0, a, 0)
    leafOnDF[i,4]<-(daylength(latit, leafOnDF$DOY[i])*60*60)
  }
  
  leafOnDF$sumGDD<-cumsum(leafOnDF$GDD)
  leafOnDay<-as.numeric(head(leafOnDF[leafOnDF$sumGDD >= 100, 1], 1))
  leafOffSub<-leafOnDF[196:365, ]
  #calculate leaf off day
  leafOffDay<-as.numeric(head(leafOffSub[leafOffSub$photoOff < 39300, 1], 1)) #39300
  phenol<-c(leafOnDay, leafOffDay) #first num is leaf on day, second num is leaf off day
  
  return(phenol)}
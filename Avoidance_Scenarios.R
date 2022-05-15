##############################################

#AVOIDANCE ANALYSIS


##############################################
#Code to generate data for Figures 4 & S4 & Tables 4, S3, & S4
source("Data_Functions.R") #need process.model
set.seed(2021)

avoidance.data.function<-function(arena.size=NULL, infection.length=NULL, C=NULL){
  if(arena.size==6){C<-36*5}else if (arena.size==2){C<-4*5}
  
  
  
  df<-main.processes.model(process.of.interest="beh.resist",
                           arena.size=arena.size,infection.length=infection.length,C=C,num.sample.periods = 0)
  
  df<-do.call(rbind.data.frame, df)
  df$arena.size<-paste(paste(arena.size,"x",sep=""), arena.size, sep="")
  df$infection.length<-paste(infection.length)
  df$prop.established<-df$total.parasites/C
  return(df)
}


#360 minutes(6x6)
Avoidance.six.360<-avoidance.data.function(arena.size=6, infection.length=360)
#180 minutes(6x6 and 2x2)
Avoidance.six.180<-avoidance.data.function(arena.size=6, infection.length=180)
Avoidance.two.180<-avoidance.data.function(arena.size=2, infection.length=180)


#135 minutes
Avoidance.six.135<-avoidance.data.function(arena.size=6, infection.length=135)
Avoidance.two.135<-avoidance.data.function(arena.size=2, infection.length=135)


#90 minutes
Avoidance.six.90<-avoidance.data.function(arena.size=6, infection.length=90)
Avoidance.two.90<-avoidance.data.function(arena.size=2, infection.length=90)


#45 minutes
Avoidance.six.45<-avoidance.data.function(arena.size=6, infection.length=45)
Avoidance.two.45<-avoidance.data.function(arena.size=2, infection.length=45)

#5 minutes
Avoidance.six.5<-avoidance.data.function(arena.size=6, infection.length=5)
Avoidance.two.5<-avoidance.data.function(arena.size=2, infection.length=5)


Avoidance.full.df<-rbind(Avoidance.six.360, Avoidance.six.180,Avoidance.two.180,Avoidance.six.135,Avoidance.two.135,Avoidance.six.90,Avoidance.two.90,
                         Avoidance.six.45,Avoidance.two.45, Avoidance.six.5,Avoidance.two.5)

write.csv(Avoidance.full.df,"data/Avoidance_5to360mins.csv")

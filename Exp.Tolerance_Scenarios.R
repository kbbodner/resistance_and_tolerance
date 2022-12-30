#####################################################

#Exposure Tolerance Scenario Generation

#####################################################

#Generates Exposure Tolerance Differences to produce Figure 5 and Tables 5, S5 & S6
source("Data_Functions.R") 


#calculate percentage alive (for fig. 5 data)
percent.alive.func<-function(df,n,t){
  extra.num<-which(!seq(1,t) %in% unique(df$time.died))
  df.withZeros<-data.frame(cbind(extra.num,rep(0,length(extra.num))))
  colnames(df.withZeros)<-c("time.died","frequency")
  
  full.df<-rbind(df,df.withZeros)
  full.df<-full.df[order(full.df$time.died),]
  
  final.df<-data.frame(cbind(full.df$time.died,cumsum(full.df$frequency)))
  colnames(final.df)<-c("Time","cumsum.died")
  
  
  final.df$percent.alive<-100*(n - final.df$cumsum.died)/n
  
  final.df<-rbind(data.frame(Time = 0, cumsum.died = 0,percent.alive=100), final.df)
  
  return(final.df)
}

set.seed(2021)

exp.tol.list<-main.processes.model(process.of.interest="exp.tol", experiment.length=20, sample.size=10)

#data prep for figure 4
Exp.Tolerance.Fig.5<-c()
for (i in 1:length(exp.tol.list)){
  
  df.total<-exp.tol.list[[i]]
  
  
  #CHANGE THIS var.level1<-"higher infection mortality";var.level2<-"higher exposure mortality"
  df.no.exp.mortality<-subset(df.total,var.level=="Higher Exposure Tolerance" & time.died !=0)
  df.exp.mortality<-subset(df.total,var.level=="Higher Infection Tolerance" & time.died !=0)
  
  df.no.exp.mortality.agg<-aggregate(var.level~time.died,df.no.exp.mortality,FUN=length)
  colnames(df.no.exp.mortality.agg)<-c("time.died","frequency")
  
  #for hosts with exposure mortality
  df.exp.mortality.agg<-aggregate(var.level~time.died,df.exp.mortality,FUN=length)
  colnames(df.exp.mortality.agg)<-c("time.died","frequency")
  
  n<-nrow(df.total)/2
  
  #function to create % survived
  only.infect.mort<-percent.alive.func(df=df.no.exp.mortality.agg,n=n,t=20) #set t to experimental length
  only.infect.mort$var.level<-"Higher Exposure Tolerance"
  only.infect.mort$sim.num<-paste(i,"A")
  
  exp.mort<-percent.alive.func(df=df.exp.mortality.agg,n=n,t=20)  #set t to experimental length
  exp.mort$var.level<-"Higher Infection Tolerance"
  exp.mort$sim.num<-paste(i,"B")
  
  Exp.Tolerance.Fig.5<-rbind(Exp.Tolerance.Fig.5,only.infect.mort,exp.mort)
  
  
  
}

#calculates percent died
Exp.Tolerance.Fig.5$percent.died<-100-exposure.tol.df$percent.alive

write.csv(Exp.Tolerance.Fig.5, "data/Exp.Tolerance_Fig5.csv")

#Data for Tables 5, S5, S6

Exp.Tolerance.20days<-do.call(rbind.data.frame, exp.tol.list)
Exp.Tolerance.20days$unique.id<-Exp.Tolerance.20days$sim.num
Exp.Tolerance.20days$experiment.length<-"20"


set.seed(2021)
exp.tol.list.10<-main.processes.model(process.of.interest="exp.tol",experiment.length=10)

Exp.Tolerance.10days<-do.call(rbind.data.frame, exp.tol.list.10)
Exp.Tolerance.10days$unique.id<-Exp.Tolerance.10days$sim.num
Exp.Tolerance.10days$experiment.length<-"10"


set.seed(2021)
exp.tol.list.2<-main.processes.model(process.of.interest="exp.tol",experiment.length=2,sample.size = 20)
Exp.Tolerance.2days<-do.call(rbind.data.frame, exp.tol.list.2)
Exp.Tolerance.2days$unique.id<-Exp.Tolerance.2days$sim.num
Exp.Tolerance.2days$experiment.length<-"2"

Exp.Tolerance.Tables<-rbind(Exp.Tolerance.20days, Exp.Tolerance.10days,Exp.Tolerance.2days)

write.csv(Exp.Tolerance.Tables, "data/Exp.Tolerance_Tables.csv")

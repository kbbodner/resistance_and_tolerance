###############################

#CLEARANCE ANALYSIS


##############################
source("Stats_Functions.R")
clearance.tablesS7<-read.csv("data/Clearance_Tables7.csv",header=T)
clearance.vs.tolerance.diff<-read.csv("data/Clearance_vs_Tolerance_Diff.csv",header=T)
clearance.sampling.Fig6a<-read.csv("data/Clearance_Fig6a.csv",header=T)
clearance.agg.Fig6b<-read.csv("data/Clearance_Fig6b.csv",header=T)

#Table S7

#Function for Poisson GLM Tables
pois.table.function<-function(df=NULL, with.established = F){
  
  #set-up for stats analysis
  df$var.level<-factor(df$var.level, levels=c("Higher Infection Resistance","Higher Clearance"))
  
  #run stats analysis
  df.stats<-pois.glm.iterator(inds.df = df,
                              process.of.interest = "clearance", 
                              clearance.with.time=F, with.established=with.established)
  
  
  
  #create list with table values
  table.values <- list(mean(df.stats$intercept.coef), mean(df.stats$coef.high),
                       sum(df.stats$p.intercept<0.05), sum(df.stats$p.high<0.05))  
  names(table.values) <- c("Avg.coef:Intercept","Avg.coef:High Clearance Group",
                           "p-value<0.05(x/100):Intercept","p-value<0.05(x/100):High Clearance Group") 
  return(table.values)
}

clearance.day1<-subset(clearance.tablesS7,samplingday=="Day1")
clearance.day3<-subset(clearance.tablesS7,samplingday=="Day3")

#Using the known established parasite burden
pois.table.function(df=clearance.day1, with.established = F)
pois.table.function(df=clearance.day3, with.established = F)

#Using the known established parasite burden
pois.table.function(df=clearance.day1, with.established = T)
pois.table.function(df=clearance.day3, with.established = T)

#Table 6

#Clearance differences only
clearance.df<-subset(clearance.vs.tolerance.diff,type.of.diff=="Clearance Diff") 

#run Cox Model
#without interaction
coxph.clearance<-coxph.iterator(inds.df=clearance.df,interaction=F) 
coxph.clearance.na.removed<-coxph.clearance[complete.cases(coxph.clearance),]

apply(coxph.clearance.na.removed[1:4],2,FUN=mean) #coefficients
apply(coxph.clearance.na.removed[5:6],2,FUN=function(x){sum(x<0.05)}) #p-values

#with interaction
coxph.clearance.interaction<-coxph.iterator(inds.df=clearance.df,interaction=T)
coxph.clearance.interaction.na.removed<-coxph.clearance.interaction[complete.cases(coxph.clearance.interaction),]
apply(coxph.clearance.interaction.na.removed[1:6],2,FUN=mean)
apply(coxph.clearance.interaction.na.removed[7:9],2,FUN=function(x){sum(x<0.05)})


#Infection Tolerance with Equal High Clearance
inf.tol.with.high.clearance.df<-subset(clearance.vs.tolerance.diff,type.of.diff=="Tol with High Clearance") 

inf.tol.with.high.clearance.df$var.level<-factor(inf.tol.with.high.clearance.df$var.level,levels=c("high","baseline"))

coxph.temp.df.with.high.clearance<-coxph.iterator(inds.df=inf.tol.with.high.clearance.df,interaction=F)

#without interaction
apply(coxph.temp.df.with.high.clearance[1:4],2,FUN=mean)
apply(coxph.temp.df.with.high.clearance[5:6],2,FUN=function(x){sum(x<0.05)})

#with interaction
coxph.temp.df.with.high.clearance.interaction<-coxph.iterator(inds.df=inf.tol.with.high.clearance.df,interaction=T) #need exp and reg. coef.

apply(coxph.temp.df.with.high.clearance.interaction[1:6],2,FUN=mean)
apply(coxph.temp.df.with.high.clearance.interaction[7:9],2,FUN=function(x){sum(x<0.05)})




# Figure 6

#subset data into different sampling categories
level.one.time<-subset(clearance.sampling.Fig6a, sampling.style=="1 time")
level.two.times.no.post<-subset(clearance.sampling.Fig6a, sampling.style=="2 times w/o post-infection")
level.two.times.with.post<-subset(clearance.sampling.Fig6a, sampling.style=="2 times w/ post-infection")
level.three.times.no.post<-subset(clearance.sampling.Fig6a, sampling.style=="3 times w/o post-infection")
level.three.times.with.post<-subset(clearance.sampling.Fig6a, sampling.style=="3 times w/ post-infection")


# #Run poisson GLM for different sampling scenarios
pois.one.time<-pois.glm.iterator(inds.df=level.one.time,process.of.interest="clearance",interaction=T)
pois.two.times.no.post<-pois.glm.iterator(inds.df=level.two.times.no.post,process.of.interest="clearance",interaction=T)
pois.two.times.with.post<-pois.glm.iterator(inds.df=level.two.times.with.post,process.of.interest="clearance",interaction=T)
 
# PREPARE DATA FOR FIGURE
one<-apply(pois.one.time,2,FUN=mean)
two.no.post<-apply(pois.two.times.no.post,2,FUN=mean)
two.with.post<-apply(pois.two.times.with.post,2,FUN=mean)

time<-seq(0,20,1)
 
#Sample Once (only at end) - Interaction Model
temp.1.1<-data.frame(cbind(time,exp(one["intercept.coef"]+time*one["coef.time"])))
colnames(temp.1.1)<-c("time", "avg.parasites")
temp.1.1$resistance.level<-"Higher Clearance"
temp.1.1$sim.num<-1
temp.1.2<-data.frame(cbind(time,exp(one["intercept.coef"]+time*one["coef.time"] +
                                      one["coef.high"] + time*one["coef.interact"])))
colnames(temp.1.2)<-c("time", "avg.parasites")
temp.1.2$resistance.level<-"Higher Infection Resistance"
temp.1.2$sim.num<-2

one.interact<-rbind(temp.1.1,temp.1.2)
one.interact$sampling.strategy<-"Sample Once (only at end) - Interaction Model"

#Sample Twice (after exposure & end) - Interaction Model
temp.1.1<-data.frame(cbind(time,exp(two.no.post["intercept.coef"]+time*two.no.post["coef.time"])))
colnames(temp.1.1)<-c("time", "avg.parasites")
temp.1.1$resistance.level<-"Higher Clearance"
temp.1.1$sim.num<-3
temp.1.2<-data.frame(cbind(time,exp(two.no.post["intercept.coef"]+
                                      time*two.no.post["coef.time"] + two.no.post["coef.high"] +
                                      time*two.no.post["coef.interact"])))
colnames(temp.1.2)<-c("time", "avg.parasites")
temp.1.2$resistance.level<-"Higher Infection Resistance"
temp.1.2$sim.num<-4

two.no.post<-rbind(temp.1.1,temp.1.2)
two.no.post$sampling.strategy<-"Sample Twice (mid-point & end) - Interaction Model"


temp.1.1<-data.frame(cbind(time,exp(two.with.post["intercept.coef"]+time*two.with.post["coef.time"])))
colnames(temp.1.1)<-c("time", "avg.parasites")
temp.1.1$resistance.level<-"Higher Clearance"
temp.1.1$sim.num<-5
temp.1.2<-data.frame(cbind(time,exp(two.with.post["intercept.coef"]+
                                      time*two.with.post["coef.time"] + two.with.post["coef.high"] +
                                      time*two.with.post["coef.interact"])))

colnames(temp.1.2)<-c("time", "avg.parasites")
temp.1.2$resistance.level<-"Higher Infection Resistance"
temp.1.2$sim.num<-6

two.with.post<-rbind(temp.1.1,temp.1.2)
two.with.post$sampling.strategy<-"Sample Twice (after exposure & end) - Interaction Model"


pois.one.time<-pois.glm.iterator(inds.df=level.one.time,process.of.interest="clearance",interaction=F)
pois.two.times.no.post<-pois.glm.iterator(inds.df=level.two.times.no.post,process.of.interest="clearance",interaction=F)
pois.two.times.with.post<-pois.glm.iterator(inds.df=level.two.times.with.post,process.of.interest="clearance",interaction=F)

#1 sample at end: 
one<-apply(pois.one.time,2,FUN=mean)
two.n.p<-apply(pois.two.times.no.post,2,FUN=mean) #two without post sampling
two.w.p<-apply(pois.two.times.with.post,2,FUN=mean) #two with post sampling

time<-seq(0,20,1)

#Sample Twice (after exposure & end) - Non-interaction Model
temp.1.1<-data.frame(cbind(time,exp(two.w.p["intercept.coef"]+time*two.w.p["coef.time"])))
colnames(temp.1.1)<-c("time", "avg.parasites")
temp.1.1$resistance.level<-"Higher Clearance"
temp.1.1$sim.num<-10
temp.1.2<-data.frame(cbind(time,exp(two.w.p["intercept.coef"]+time*two.w.p["coef.time"] + 
                                      two.w.p["coef.high"])))
colnames(temp.1.2)<-c("time", "avg.parasites")
temp.1.2$resistance.level<-"Higher Infection Resistance"
temp.1.2$sim.num<-11

two.no.interact<-rbind(temp.1.1,temp.1.2)
two.no.interact$sampling.strategy<-"Sample Twice (after exposure & end) - Non-interaction Model"


df.total<-rbind(two.no.interact,one.interact,two.no.post,two.with.post)

df.total$sampling.strategy.2<-factor(df.total$sampling.strategy,levels=c("Sample Once (only at end) - Interaction Model","Sample Twice (mid-point & end) - Interaction Model",
                                                                         "Sample Twice (after exposure & end) - Interaction Model",
                                                                         "Sample Twice (after exposure & end) - Non-interaction Model" ))


#Figure 6

fig.6<-ggplot(data=clearance.agg.Fig6b,aes(x=time,y=avg.parasites,color=resistance.level)) + 
  geom_point() +
  geom_line(data=subset(df.total, time!=0),aes(x=time,y=avg.parasites,group=sim.num,color=resistance.level, 
                                               linetype=sampling.strategy.2),size=1) + 
  scale_color_manual(values = c("purple", "goldenrod3")) + theme_classic() + 
  scale_x_continuous(breaks = seq(1,21,2), limits=c(1,21)) + 
  scale_y_continuous(breaks = seq(0,100,5),limits = c(0,30)) +
  labs(x = "Experiment Time (# of Timesteps)", y = "Avg. Established Parasite Burden", 
       color = "Host Resistance Strategy", linetype="Sampling Design & Model Configuration") +
  scale_linetype_manual(values=c("dotted", "dotdash", "twodash","solid")) + 
  theme_bw() + theme(axis.text=element_text(size=14),axis.title=element_text(size=15),
                     legend.title=element_text(size=15),legend.text=element_text(size=12)) +
  theme(strip.background =element_rect(fill="grey")) +
  theme(strip.text = element_text(colour = 'white',size=19,face="bold"))


# png(filename="Figure6.png",height=2300,width=5200,res=400)
# fig.6
# dev.off()


#Calculate RMSE for average Model Prediction with Average parasite burden
sample.once.error=0
sample.twice.midpoint.error=0
sample.twice.after.exp.error=0
non.interact.error=0

for (i in 1:20){
  predict<-subset(df.total,time==i)
  data<-subset(clearance.agg.Fig6b,time==i)
  
  data.highClear<-data$avg.parasites[data$resistance.level=="Higher Clearance"]
  data.highInfect<-data$avg.parasites[data$resistance.level=="Higher Infection Resistance"]
  
  #Sample Once error (interaction)
  interact.once<-subset(predict,predict$sampling.strategy== "Sample Once (only at end) - Interaction Model")
  error.once.clear<-(interact.once$avg.parasites[interact.once$resistance.level=="Higher Clearance"] - data.highClear)^2
  error.once.resist<-(interact.once$avg.parasites[interact.once$resistance.level=="Higher Infection Resistance"] - data.highInfect)^2
  
  sample.once.error<-sample.once.error+ sum(error.once.clear) + sum(error.once.resist)
  
  #Sample Twice midpoint error (interaction)
  interact.twice.midpoint<-subset(predict,predict$sampling.strategy== "Sample Twice (mid-point & end) - Interaction Model")
  error.twice.midpoint.clear<-(interact.twice.midpoint$avg.parasites[interact.twice.midpoint$resistance.level=="Higher Clearance"] - data.highClear)^2
  error.twice.midpoint.resist<-(interact.twice.midpoint$avg.parasites[interact.twice.midpoint$resistance.level=="Higher Infection Resistance"] - data.highInfect)^2
  
  sample.twice.midpoint.error<-sample.twice.midpoint.error + sum(error.twice.midpoint.clear) + sum(error.twice.midpoint.resist)
  
  
  
  #Sample Twice after exposure error (interaction)
  interact.twice.after.exp<-subset(predict,predict$sampling.strategy== "Sample Twice (after exposure & end) - Interaction Model")
  
  error.twice.after.exp.clear<-(interact.twice.after.exp$avg.parasites[interact.twice.after.exp$resistance.level=="Higher Clearance"] - data.highClear)^2
  error.twice.after.exp.resist<-(interact.twice.after.exp$avg.parasites[interact.twice.after.exp$resistance.level=="Higher Infection Resistance"] - data.highInfect)^2
  
  sample.twice.after.exp.error<-sample.twice.after.exp.error+ sum(error.twice.after.exp.clear) + sum(error.twice.after.exp.resist)
  
  #sample Twice after exposure error (non-interaction)
  non.interact<-subset(predict,predict$sampling.strategy== "Sample Twice (after exposure & end) - Non-interaction Model")
  non.interact.clear<-(non.interact$avg.parasites[non.interact$resistance.level=="Higher Clearance"] - data.highClear)^2
  non.interact.resist<-(non.interact$avg.parasites[non.interact$resistance.level=="Higher Infection Resistance"] - data.highInfect)^2
  
  non.interact.error<-non.interact.error + sum(non.interact.clear) + sum(non.interact.resist)
  
}
n<-nrow(clearance.agg.Fig6b)

RMSE.once<-sqrt(sample.once.error/n)
RMSE.twice.midpoint<-sqrt(sample.twice.midpoint.error/n)
RMSE.twice.after.exp<-sqrt(sample.twice.after.exp.error/n)
RMSE.non.interact<-sqrt(non.interact.error/n)

####ROOT MEAN SQUARED ERROR
#first entry: sampling once at end (interaction model)
#second entry: twice but not immediately post-exposure (interaction model)
#third entry: twice but sampling immediately post-exposure (interaction model)
#fourth entry: no interaction term in model
print(c(RMSE.once,RMSE.twice.midpoint,RMSE.twice.after.exp,RMSE.non.interact))



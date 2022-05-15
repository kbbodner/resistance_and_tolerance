###############################

#CLEARANCE ANALYSIS


##############################
source("Stats_Functions.R")
clearance.tablesS7<-read.csv("data/Clearance_Tables7",header=T)
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
clearance.df<-subset(clearance.vs.tolerance.diff,tolerance.diff=="no") 

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


#Infection Tolerance with Equal Clearance
inf.tol.with.clearance.df<-subset(clearance.vs.tolerance.diff,tolerance.diff=="yes") 

inf.tol.with.clearance.df$var.level<-factor(inf.tol.with.clearance.df$var.level,levels=c("high","baseline"))

coxph.temp.df.with.clearance<-coxph.iterator(inds.df=inf.tol.with.clearance.df,interaction=F)

#without interaction
apply(coxph.temp.df.with.clearance[1:4],2,FUN=mean)
apply(coxph.temp.df.with.clearance[5:6],2,FUN=function(x){sum(x<0.05)})

#with interaction
coxph.temp.df.without.mu.interaction<-coxph.iterator(inds.df=inf.tol.without.mu,interaction=T) #need exp and reg. coef.

apply(coxph.temp.df.without.mu.interaction[1:6],2,FUN=mean)
apply(coxph.temp.df.without.mu.interaction[7:9],2,FUN=function(x){sum(x<0.05)})


#Stats for Figure 6
#data has been read in
##REMINDER: 
clearance.sampling.Fig6a #read.csv("data/Clearance_Fig6a.csv",header=T)
clearance.agg.Fig6b<- #("data/Clearance_Fig6b.csv",header=T)

#Figure 6

ggplot(data=clearance.agg.Fig6,aes(x=time,y=avg.parasites,color=resistance.level)) + 
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





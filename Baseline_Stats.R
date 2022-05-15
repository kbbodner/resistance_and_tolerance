##############################################################################

#Statistics for Baseline Analyses: Infection Resistance and Tolerance

################################################################################
source("Stats_Functions.R")

#INFECTION RESISTANCE

#Read in Baseline Scenario data
baseline.infect.resist_25<-read.csv("data/Baseline_Resist_25%.csv") #for main analysis
baseline.infect.resist_10<-read.csv("data/Baseline_Resist_10%.csv") #for supplementary material
baseline.infect.resist_75<-read.csv("data/Baseline_Resist_75%.csv") #for supplementary material

#Create host grouping variable
baseline.infect.resist_25$var.level<-factor(baseline.infect.resist_25$var.level,levels=c("high","baseline"))
baseline.infect.resist_10$var.level<-factor(baseline.infect.resist_10$var.level,levels=c("high","baseline"))
baseline.infect.resist_75$var.level<-factor(baseline.infect.resist_75$var.level,levels=c("high","baseline"))

#Run Poisson GLMs for scenarios
baseline.resist.temp.25<-pois.glm.iterator.basic.exp(inds.df=baseline.infect.resist_25)
baseline.resist.temp.10<-pois.glm.iterator.basic.exp(inds.df=baseline.infect.resist_10)
baseline.resist.temp.75<-pois.glm.iterator.basic.exp(inds.df=baseline.infect.resist_75)

#write .csv for Figures and Stats
#write.csv(baseline.resist.temp.25,"data/Baseline_Resist_25%_Stats.csv")
#write.csv(baseline.resist.temp.10,"data/Baseline_Resist_10%_Stats.csv")
#write.csv(baseline.resist.temp.75,"data/Baseline_Resist_75%_Stats.csv")

#INFECTION TOLERANCE

#Read in Baseline Scenario data
baseline.infect.tol_5<-read.csv("data/Baseline_Tol_5x.csv") #for main analysis
baseline.infect.tol_2.5<-read.csv("data/Baseline_Tol_2.5x.csv") #for supplementary material
baseline.infect.tol_10<-read.csv("data/Baseline_Tol_10x.csv") #for supplementary material

#Create host grouping variable
baseline.infect.tol_5$var.level<-factor(baseline.infect.tol_5$var.level,levels=c("high","baseline"))
baseline.infect.tol_2.5$var.level<-factor(baseline.infect.tol_2.5$var.level,levels=c("high","baseline"))
baseline.infect.tol_10$var.level<-factor(baseline.infect.tol_10$var.level,levels=c("high","baseline"))

#Run Cox Model for Scenarios
baseline.tol.temp.5<-coxph.iterator(inds.df=baseline.infect.tol_5,interaction=F) 
baseline.tol.temp.2.5<-coxph.iterator(inds.df=baseline.infect.tol_2.5,interaction=F) 
baseline.tol.temp.10<-coxph.iterator(inds.df=baseline.infect.tol_10,interaction=F) 

#write .csv for Figures and Stats
#write.csv(baseline.tol.temp.5,"data/Baseline_Tol_5x_Stats.csv")
#write.csv(baseline.tol.temp.2.5,"data/Baseline_Tol_2.5x_Stats.csv")
#write.csv(baseline.tol.temp.10,"data/Baseline_Tol_10x_Stats.csv")


########################################

#Statistics for Tables 2 & S2

#######################################
source("Data_Functions.R")

####TABLE 2

#INFECTION RESISTANCE
set.seed(2021)

beta.high=0.5*1.25
beta.list<-main.experimental.setup(arena.size=6, 
                                   infection.length=180,
                                   var.of.interest="beta",beta.2=beta.high,C=30)

baseline_infect_resist.stats<-do.call(rbind.data.frame, beta.list)
baseline_infect_resist.stats$resist.diff<-paste("25%")

baseline_infect_resist.stats$unique.id<-paste(baseline_infect_resist.stats$resist.diff,baseline_infect_resist.stats$sim.number)
baseline_infect_resist.stats$var.level<-factor(baseline_infect_resist.stats$var.level,levels=c("high","baseline"))

#Run Poisson GLM
pois.temp.df<-pois.glm.iterator.basic.exp(inds.df=baseline_infect_resist.stats) 

#Mean, min and max for Table 2
apply(pois.temp.df[1:4],2,FUN=mean)
apply(pois.temp.df[1:4],2,FUN=min)
apply(pois.temp.df[1:4],2,FUN=max)
sum(pois.temp.df$p.intercept<0.05)
sum(pois.temp.df$p.high<0.05)

#INFECTION TOLERANCE
set.seed(2021)
alpha.high=0.0005*5
alpha.list<-main.experimental.setup(arena.size=6, 
                                    infection.length=180,num.sample.periods=0, mu=0,
                                    var.of.interest="alpha",alpha.2=alpha.high,C=30)

baseline_infect_tol.stats<-do.call(rbind.data.frame, alpha.list)
baseline_infect_tol.stats$tol.diff<-paste("5x")


baseline_infect_tol.stats$unique.id<-paste(baseline_infect_tol.stats$tol.diff,baseline_infect_tol.stats$sim.number)
baseline_infect_tol.stats$var.level<-factor(baseline_infect_tol.stats$var.level,levels=c("high","baseline"))


#Without Interaction
coxph.temp.df<-coxph.iterator(inds.df=baseline_infect_tol.stats,interaction=F) 

apply(coxph.temp.df[1:6],2,FUN=mean)
apply(coxph.temp.df[1:6],2,FUN=min)
apply(coxph.temp.df[1:6],2,FUN=max)
sum(coxph.temp.df$p.parasites<0.05)
sum(coxph.temp.df$p.high<0.05)

#With Interaction
coxph.temp.df.interact<-coxph.iterator(inds.df=baseline_infect_tol.stats,interaction=T) 

apply(coxph.temp.df.interact[1:6],2,FUN=mean)
apply(coxph.temp.df.interact[1:6],2,FUN=min)
apply(coxph.temp.df.interact[1:6],2,FUN=max)
sum(coxph.temp.df.interact$p.parasites<0.05)
sum(coxph.temp.df.interact$p.high<0.05)
sum(coxph.temp.df.interact$p.high.int<0.05)


####TABLE S2
set.seed(2021)

#INFECTION TOLERANCE
alpha.high=0.0005*50
alpha.list.high<-main.experimental.setup(arena.size=6, 
                                    infection.length=180,num.sample.periods=0, mu=0,
                                    var.of.interest="alpha",alpha.2=alpha.high,C=30)

baseline_infect_tol.stats.high<-do.call(rbind.data.frame, alpha.list.high)
baseline_infect_tol.stats.high$tol.diff<-paste("50x")


baseline_infect_tol.stats.high$unique.id<-paste(baseline_infect_tol.stats.high$tol.diff,baseline_infect_tol.stats.high$sim.number)
baseline_infect_tol.stats.high$var.level<-factor(baseline_infect_tol.stats.high$var.level,levels=c("high","baseline"))


#Without Interaction
coxph.temp.df.high<-coxph.iterator(inds.df=baseline_infect_tol.stats.high,interaction=F) 

apply(coxph.temp.df.high[1:6],2,FUN=mean)
apply(coxph.temp.df.high[1:6],2,FUN=min)
apply(coxph.temp.df.high[1:6],2,FUN=max)
sum(coxph.temp.df.high$p.parasites<0.05)
sum(coxph.temp.df.high$p.high<0.05)

#With Interaction
coxph.temp.df.interact.high<-coxph.iterator(inds.df=baseline_infect_tol.stats.high,interaction=T) 

apply(coxph.temp.df.interact.high[1:6],2,FUN=mean)
apply(coxph.temp.df.interact.high[1:6],2,FUN=min)
apply(coxph.temp.df.interact.high[1:6],2,FUN=max)
sum(coxph.temp.df.interact.high$p.parasites<0.05)
sum(coxph.temp.df.interact.high$p.high<0.05)
sum(coxph.temp.df.interact.high$p.high.int<0.05)




source("Data_Functions.R")

#Data for Table S7

#First sampling on Day 1 (includes both infection resistance and clearance differences)
set.seed(2021)
clearance.list.day1<-main.processes.model(process.of.interest="clearance",
                                          num.sample.periods=0,sample.size=10, 
                                          experiment.length=1, 
                                          infect.resist.and.clearance.diff = T)

clearance.day1<-do.call(rbind.data.frame, clearance.list.day1)
clearance.day1$unique.id<-clearance.day1$sim.num
clearance.day1$samplingday<-"Day1"

#First sampling on Day 5 (includes both infection resistance and clearance)
set.seed(2021)
clearance.list.day3<-main.processes.model(process.of.interest="clearance",
                                               num.sample.periods=0,sample.size=10, 
                                               experiment.length=3,
                                               infect.resist.and.clearance.diff = T)

clearance.day3<-do.call(rbind.data.frame, clearance.list.day3)
clearance.day3$unique.id<-clearance.day3$sim.num
clearance.day3$samplingday<-"Day3"

clearance.tableS7<-rbind(clearance.day1,clearance.day3)

write.csv(clearance.tableS7,"data/Clearance_TableS7.csv")
#####################################################
#Table 6

set.seed(2021)

#Here there are no infection resistance differences (only clearance differences)
clearance.list<-main.processes.model(process.of.interest="clearance",
                                          num.sample.periods=0)

clearance.df<-do.call(rbind.data.frame, clearance.list)
clearance.df$unique.id<-clearance.df$sim.num
clearance.df$type.of.diff<-"Clearance Diff"


#INFECTION TOLERANCE WITH CLEARANCE

set.seed(2021)

#Infection tolerance Difference with Clearance (same as baseline)
inf.tol.with.high.clearance.list<-main.experimental.setup(arena.size=6, 
                                                 infection.length=180,
                                                 var.of.interest="alpha",
                                                 alpha.2=0.0005*5,mu=10*1/100,num.sample.periods = 0)


inf.tol.with.high.clearance.df<-do.call(rbind.data.frame, inf.tol.with.high.clearance.list)
inf.tol.with.high.clearance.df$unique.id<-inf.tol.with.high.clearance.df$sim.number
#to make it match with clearance dataset remove arena.size and infection.length.minutes
inf.tol.with.high.clearance.df<-subset(inf.tol.with.high.clearance.df,select=-c(arena.size,infection.length.minutes))
names(inf.tol.with.high.clearance.df)[15]<-"sim.num" #change sim.number to sim.num to match clearance

#create identifier for both datasets
inf.tol.with.high.clearance.df$type.of.diff<-"Tol with High Clearance"

#with low clearance
set.seed(2021)

#Infection tolerance Difference with Clearance (same as baseline)
inf.tol.with.low.clearance.list<-main.experimental.setup(arena.size=6, 
                                                     infection.length=180,
                                                     var.of.interest="alpha",
                                                     alpha.2=0.0005*5,mu=1/100,num.sample.periods = 0)


inf.tol.with.low.clearance.df<-do.call(rbind.data.frame, inf.tol.with.low.clearance.list)
inf.tol.with.low.clearance.df$unique.id<-inf.tol.with.low.clearance.df$sim.number
#to make it match with clearance dataset remove arena.size and infection.length.minutes
inf.tol.with.low.clearance.df<-subset(inf.tol.with.low.clearance.df,select=-c(arena.size,infection.length.minutes))
names(inf.tol.with.low.clearance.df)[15]<-"sim.num" #change sim.number to sim.num to match clearance

#create identifier for both datasets
inf.tol.with.low.clearance.df$type.of.diff<-"Tol with Low Clearance"





clearance.vs.tolerance.diff<-rbind(clearance.df,inf.tol.with.high.clearance.df, inf.tol.with.low.clearance.df)

write.csv(clearance.vs.tolerance.diff,"data/Clearance_vs_Tolerance_Diff.csv")


##################################################
#Data for Figure 6
clearance.sampling.list<-est.resist.model()


#takes only the sampling style 
clearance.sampling.df<-do.call(rbind.data.frame,clearance.sampling.list[[1]])

clearance.sampling.df$unique.id<-clearance.sampling.df$sim.num

write.csv(clearance.sampling.df, "data/Clearance_Fig6a.csv")

#Aggregated Data (Part b)
clearance.sampling.df.agg<-subset(do.call(rbind.data.frame,clearance.sampling.list[[2]]),alive.or.dead=="alive"& t<21)
colnames(clearance.sampling.df.agg)<-c("time", "avg.parasites","alive.or.dead","resistance.level", "sim.num")
clearance.sampling.df.agg$resistance.level<-factor(clearance.sampling.df.agg$resistance.level,labels=c("Higher Clearance", "Higher Infection Resistance"))

write.csv(clearance.sampling.df.agg, "data/Clearance_Fig6b.csv")



#Create Stats model for Sample Twice (after exposure & end) - Non-interaction Model
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

#Figure  6


#Tables


level.one.time<-subset(all.data.est.mort, sampling.style=="1 time")
level.two.times.no.post<-subset(all.data.est.mort, sampling.style=="2 times w/o post-infection")
level.two.times.with.post<-subset(all.data.est.mort, sampling.style=="2 times w/ post-infection")
level.three.times.no.post<-subset(all.data.est.mort, sampling.style=="3 times w/o post-infection")
level.three.times.with.post<-subset(all.data.est.mort, sampling.style=="3 times w/ post-infection")


pois.one.time<-pois.glm.iterator(inds.df=level.one.time,process.of.interest="est.resist",interaction=T)
pois.two.times.no.post<-pois.glm.iterator(inds.df=level.two.times.no.post,process.of.interest="est.resist",interaction=T)
pois.two.times.with.post<-pois.glm.iterator(inds.df=level.two.times.with.post,process.of.interest="est.resist",interaction=T)
pois.three.times.no.post<-pois.glm.iterator(inds.df=level.three.times.no.post,process.of.interest="est.resist",interaction=T)
pois.three.times.with.post<-pois.glm.iterator(inds.df=level.three.times.with.post,process.of.interest="est.resist",interaction=T)


##########################################
without.dead.day1<-subset(level.three.times.with.post, t==1)
without.dead.day7<-subset(level.three.times.no.post,  t==7)
test<-pois.glm.iterator(inds.df=without.dead.day1,process.of.interest = "beh.resist")
test2<-pois.glm.iterator(inds.df=without.dead.day7,process.of.interest = "beh.resist")

apply(test[1:4],2,FUN=mean)
apply(test2[1:4],2,FUN=mean)
sum(test$p.intercept<0.05)
sum(test$p.high<0.05)
sum(test2$p.intercept<0.05)
sum(test2$p.high<0.05)

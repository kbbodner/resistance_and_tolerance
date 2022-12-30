#Main Scenario Generating Functions
source("IBM.R")

#BASIC EXPERIMENT DESIGN FUNCTION USED BY BASELINE ANLALYSES AND CLEARANCE ANALYSES
#Includes arena size, infection length, sampling periods, experiment length and sample size for research design characteristics
#Can select whether the var.of.interest is beta (higher and lower infection resistance) or alpha (higher and lower infection tolerance)
#runs IBM.R to create scenarios from IBMs
main.experimental.setup<-function(arena.size=6, infection.length=180,mu=0,num.sample.periods=2,var.of.interest=NULL,
                                  beta.2=1,alpha.2=1,resist.and.tol.proc="none",
                                  experiment.length=20, sample.size=10,C=30){
  df.total.list<-list()
  
  
  for (i in 1:100){
    
    print(i)
    #set baseline parameters for experimental design
    C=C; experiment.length=experiment.length;post.infect.sample=T; mu.1=mu;mu.2=mu
    
    #baseline parameters
    alpha<-0.0005; alpha.high=0.0005
    beta<-0.5; beta.high=0.5
    
    #change to variable of interest
    if(var.of.interest=="alpha"){
      alpha.high<-alpha.2 
    }else if(var.of.interest == "beta"){beta.high<-beta.2}#0.625
    
    if(resist.and.tol.proc=="opposite"){ #opposite means higher tolerance and lower immunity
      mu.2=mu*10 #baseline has higher tolerance (alpha is lower) and lower immunity (mu is lower)
      
    }else if (resist.and.tol.proc=="both"){
      mu.1=mu*10  #baseline has higher tolerance (alpha is lower) and higher immunity (mu is higher)
    }
    
    #added experimental length to functions so if things go wrong, remove that.
    baseline<-experimental.model.run(alpha=alpha,mu=mu.1,beta=beta,H=100,C=C,arena.size=arena.size,infection.length=infection.length,
                                     beh.res=F,pent.mort=F,parasite.mort=F, experiment.length=experiment.length) #returns inds_hist (a list)
    high<-experimental.model.run(alpha=alpha.high,mu=mu.2,beta=beta.high,H=100,C=C,arena.size=arena.size,infection.length=infection.length,
                                 beh.res=F,pent.mort=F,parasite.mort=F, experiment.length=experiment.length) #returns inds_hist (a list)
    
    
    baseline.df<-sampling.time(inds_hist=baseline,experiment.length=experiment.length,num.sample.periods=num.sample.periods,
                               post.infect.sample=post.infect.sample, sample.size=sample.size) #returns single dataframe
    baseline.df$var.level<-"baseline"
    
    #sets the number and timing of sampling
    high.df<-sampling.time(inds_hist=high,experiment.length=experiment.length,num.sample.periods=num.sample.periods,
                           post.infect.sample=post.infect.sample, sample.size=sample.size) #returns single dataframe
    high.df$var.level<-"high"
    
    if(class(high.df)!="data.frame"){
      return("Num.sample.periods or num.sampled is too large for the number of hosts available")
    }
    
    
    final.df<-rbind(baseline.df,high.df)
    
    if(C!=0){
    final.df$arena.size<-paste(arena.size)
    final.df$infection.length.minutes<-paste(infection.length)}
    if(C==0){
      final.df$arena.size<-0
      final.df$infection.length.minutes<-0
    }
    
    
    final.df$sim.number<-paste(i)
    
    df.total.list[[i]]<-final.df
    
  }
  
  return(df.total.list)
}


#AVOIDANCE, EXPOSURE TOL. & CLEARANCE FUNCTION USED BY AVOIDANCE_MAIN.R TO MEASURE DIFFERENCES IN THESE THREE PROCESSES
#Allows for differences in arena size, dosage (C), infection length, number of sampling periods, sample size, experiment length and whether there was sampling following exposure
#infect.resist.and.clearance.diff indicates whether this analyses is for clearance (where clearance is tested to see if it can be confounded as tolerance)
main.processes.model<-function(process.of.interest=NULL,arena.size=6,C=30,infection.length=180, num.sample.periods=2,
                               sample.size=10,experiment.length=20,post.infect.sample=T, 
                               infect.resist.and.clearance.diff = F){
  
  df.total.list<-list()
  
  for (i in 1:100){
    
    print(i)
    
    #set baseline parameters for experimental design
    arena.size=arena.size;infection.length=infection.length;
    experiment.length=experiment.length;post.infect.sample=post.infect.sample
    
    
    #baseline parameters
    alpha<-0.0005; alpha.high=0.0005
    beta<-0.5; beta.high=0.5; H=100
    mu=0; mu.high=0
    alpha.pent.mort.high=0; alpha.pent.mort.low=0
    
    #baseline processes
    beh.res=F;pent.mort=F;parasite.mort=F
    
    
    #PROCESS-OF-INTEREST
    if(process.of.interest=="exp.tol"){
      pent.mort=T;
      alpha.high=alpha.high*5 
      alpha.pent.mort.high=0.005; alpha.pent.mort.low=0.005/5
      var.level1<-"Higher Exposure Tolerance";var.level2<-"Higher Infection Tolerance"
      #experiment.length=50
    }else if(process.of.interest=="beh.resist"){
      beta.high=0.5*1.25
      beh.res=T;H=20;experiment.length=1
      num.sample.periods=0
      var.level1<-"Behavioural Resistance";var.level2<-"Higher Infection Resistance"
    }else if (process.of.interest=="clearance"){
      parasite.mort=T;
      if(infect.resist.and.clearance.diff == T){
        beta.high=0.5*1.25;H=20
      }
      mu<-1/100; mu.high<-10*1/100
      var.level1<-"Higher Clearance";var.level2<-"Higher Infection Resistance"
    }
    
    
    baseline<-experimental.model.run(alpha=alpha.high,mu=mu.high,beta=beta.high,alpha.pent=alpha.pent.mort.low,
                                     H=H,C=C,arena.size=arena.size,infection.length=infection.length,
                                     beh.res=beh.res,pent.mort=pent.mort,parasite.mort=parasite.mort) #returns inds_hist (a list)
    
    high<-experimental.model.run(alpha=alpha,mu=mu,beta=beta,alpha.pent=alpha.pent.mort.high,
                                 H=H,C=C,arena.size=arena.size,infection.length=infection.length,
                                 beh.res=F,pent.mort=pent.mort,parasite.mort=parasite.mort) #returns inds_hist (a list)
    
    
    baseline.df<-sampling.time(inds_hist=baseline,experiment.length=experiment.length,num.sample.periods=num.sample.periods,
                               post.infect.sample=post.infect.sample, sample.size=sample.size) #returns single dataframe
    baseline.df$var.level<-var.level1
    
    
    high.df<-sampling.time(inds_hist=high,experiment.length=experiment.length,num.sample.periods=num.sample.periods,
                           post.infect.sample=post.infect.sample, sample.size=sample.size) #returns single dataframe
    high.df$var.level<-var.level2
    
    
    final.df<-rbind(baseline.df,high.df)
    final.df$sim.num<-paste(i)
    
    df.total.list[[i]]<-final.df
    
  }
  
  
  return(df.total.list)
  
  
  
}

#CLEARANCE ANALYSES FOR MEASURING BOTH INFECTION RESISTANCE AND CLEARANCE DIFFERENCES
#Runs number of samples and timing of samples for hosts with clearance differences
#takes in arena size, dosage level (C), length of exposure (infection.length), number of samples (sample size) and experiment length
est.resist.model<-function(arena.size=6,C=30,infection.length=180,
                           sample.size=10,experiment.length=20){
  
  df.total.list<-list()
  df.total.list.beta<-list()
  
  for (i in 1:100){
    
    print(i)
    
    
    #baseline parameters
    alpha<-0.0005; alpha.high=0.0005
    beta<-0.5; beta.high=0.5*1.25
    H=100
    mu<-1/100; mu.high<-10*1/100
    alpha.pent.mort.high=0; alpha.pent.mort.low=0
    var.level1<-"lower parasite immunity";var.level2<-"higher parasite immunity"
    
    
    
    baseline<-experimental.model.run(alpha=alpha.high,mu=mu,beta=beta,alpha.pent=alpha.pent.mort.low,
                                     H=H,C=C,arena.size=arena.size,infection.length=infection.length,
                                     beh.res=F,pent.mort=F,parasite.mort=T) #returns inds_hist (a list)
    
    high<-experimental.model.run(alpha=alpha,mu=mu.high,beta=beta.high,alpha.pent=alpha.pent.mort.high,
                                 H=H,C=C,arena.size=arena.size,infection.length=infection.length,
                                 beh.res=F,pent.mort=F,parasite.mort=T) #returns inds_hist (a list)
    
    #Four different sampling strategies - baseline
    #once
    baseline.one.time<-sampling.time(inds_hist=baseline,experiment.length=experiment.length,num.sample.periods=0,
                                     post.infect.sample=F, sample.size=sample.size) #returns single dataframe
    baseline.one.time$var.level<-var.level1
    baseline.one.time$sampling.style<-"1 time"
    
    #twice
    baseline.two.times.with.post<-sampling.time(inds_hist=baseline,experiment.length=experiment.length,num.sample.periods=1,
                                                post.infect.sample=T, sample.size=sample.size) #returns single dataframe
    baseline.two.times.with.post$var.level<-var.level1
    baseline.two.times.with.post$sampling.style<-"2 times w/ post-infection"
    
    
    baseline.two.times.without.post<-sampling.time(inds_hist=baseline,experiment.length=experiment.length,num.sample.periods=1,
                                                   post.infect.sample=F, sample.size=sample.size) #returns single dataframe
    baseline.two.times.without.post$var.level<-var.level1
    baseline.two.times.without.post$sampling.style<-"2 times w/o post-infection"
    
    
    #three times
    baseline.three.times.with.post<-sampling.time(inds_hist=baseline,experiment.length=experiment.length,num.sample.periods=2,
                                                  post.infect.sample=T, sample.size=sample.size) #returns single dataframe
    baseline.three.times.with.post$var.level<-var.level1
    baseline.three.times.with.post$sampling.style<-"3 times w/ post-infection"
    
    
    baseline.three.times.without.post<-sampling.time(inds_hist=baseline,experiment.length=experiment.length,num.sample.periods=2,
                                                     post.infect.sample=F, sample.size=sample.size) #returns single dataframe
    baseline.three.times.without.post$var.level<-var.level1
    baseline.three.times.without.post$sampling.style<-"3 times w/o post-infection"
    
    
    baseline.df<-rbind(baseline.one.time, baseline.two.times.with.post, baseline.two.times.without.post,
                       baseline.three.times.with.post, baseline.three.times.without.post)
    
    ###############
    #Four different sampling strategies - high
    #once
    high.one.time<-sampling.time(inds_hist=high,experiment.length=experiment.length,num.sample.periods=0,
                                 post.infect.sample=F, sample.size=sample.size) #returns single dataframe
    high.one.time$var.level<-var.level2
    high.one.time$sampling.style<-"1 time"
    
    #twice
    high.two.times.with.post<-sampling.time(inds_hist=high,experiment.length=experiment.length,num.sample.periods=1,
                                            post.infect.sample=T, sample.size=sample.size) #returns single dataframe
    high.two.times.with.post$var.level<-var.level2
    high.two.times.with.post$sampling.style<-"2 times w/ post-infection"
    
    
    high.two.times.without.post<-sampling.time(inds_hist=high,experiment.length=experiment.length,num.sample.periods=1,
                                               post.infect.sample=F, sample.size=sample.size) #returns single dataframe
    high.two.times.without.post$var.level<-var.level2
    high.two.times.without.post$sampling.style<-"2 times w/o post-infection"
    
    
    #three times
    high.three.times.with.post<-sampling.time(inds_hist=high,experiment.length=experiment.length,num.sample.periods=2,
                                              post.infect.sample=T, sample.size=sample.size) #returns single dataframe
    high.three.times.with.post$var.level<-var.level2
    high.three.times.with.post$sampling.style<-"3 times w/ post-infection"
    
    
    high.three.times.without.post<-sampling.time(inds_hist=high,experiment.length=experiment.length,num.sample.periods=2,
                                                 post.infect.sample=F, sample.size=sample.size) #returns single dataframe
    high.three.times.without.post$var.level<-var.level2
    high.three.times.without.post$sampling.style<-"3 times w/o post-infection"
    
    
    high.df<-rbind(high.one.time, high.two.times.with.post, high.two.times.without.post,
                   high.three.times.with.post, high.three.times.without.post)
    
    
    #combine all together
    final.df<-rbind(baseline.df,high.df)
    final.df$sim.num<-paste(i)
    
    df.total.list[[i]]<-final.df
    
    
    
    #ONLY DEAD
    baseline.df.temp<-do.call(rbind.data.frame,baseline)
    baseline.df.dead<-subset(baseline.df.temp,host.mort==1)
    baseline.df.dead.agg<-aggregate(baseline.df.dead$parasites.dead~baseline.df.dead$t,FUN=mean)
    colnames(baseline.df.dead.agg)<-c("t", "total.parasites")
    baseline.df.dead.agg$alive.or.dead<-"dead"
    baseline.df.dead.agg$var.level<-var.level1
    
    high.df.temp<-do.call(rbind.data.frame,high)
    high.df.dead<-subset(high.df.temp,host.mort==1)
    high.df.dead.agg<-aggregate(high.df.dead$parasites.dead~high.df.dead$t,FUN=mean)
    colnames(high.df.dead.agg)<-c("t", "total.parasites")
    high.df.dead.agg$alive.or.dead<-"dead"
    high.df.dead.agg$var.level<-var.level2
    
    #ONLY ALIVE
    baseline.df.beta<-subset(baseline.df.temp,time.died==0)
    baseline.df.beta.agg<-aggregate(baseline.df.beta$No.parasites~baseline.df.beta$t,FUN=mean)
    colnames(baseline.df.beta.agg)<-c("t", "total.parasites")
    baseline.df.beta.agg$alive.or.dead<-"alive"
    baseline.df.beta.agg$var.level<-var.level1
    
    high.df.beta<-subset(high.df.temp,time.died==0)
    high.df.beta.agg<-aggregate(high.df.beta$No.parasites~high.df.beta$t,FUN=mean)
    colnames(high.df.beta.agg)<-c("t", "total.parasites")
    high.df.beta.agg$alive.or.dead<-"alive"
    high.df.beta.agg$var.level<-var.level2
    
    temp.df<-rbind(baseline.df.beta.agg,high.df.beta.agg,baseline.df.dead.agg, high.df.dead.agg)
    temp.df$sim.num<-paste(i)
    
    df.total.list.beta[[i]]<-temp.df
  }
  
  
  
  return(list(df.total.list,df.total.list.beta)) #sampling, total aggregate
  
}


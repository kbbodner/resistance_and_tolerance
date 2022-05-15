#MAIN MODEL FILE


#HOST DEATH -> PARASITE INDUCED MORTALITY
#inds is a matrix, alpha is a single parameter (probability)
#alpha is the probability of host dying PER parasite
#probability is not density dependent

host.mort<-function(inds,alpha=NULL){
  
  pop.size<-nrow(inds)
  
  #binomial distribution where each parasite 
  #has an alpha probability of successfully killing the host
  #is this the same as bernoulli with n*p where n is number of parasites?
  inds[,"host.mort"]<-rbinom(n=pop.size,size=inds[,"No.parasites"],p=alpha)

  #change to 1 or 0 if host died
  temp.dead<-ifelse(inds[,"host.mort"]>0,1,0)
  inds[,"host.mort"]<-temp.dead
  
  #put in parasites that died
  temp.par.dead<-ifelse(inds[,"host.mort"]==1,inds[,"No.parasites"],0)
  inds[,"parasites.dead"]<-temp.par.dead + inds[,"parasites.dead"]
  
  #update number of parasites
  temp.par<-ifelse(inds[,"host.mort"]==1,0,inds[,"No.parasites"])
  inds[,"No.parasites"]<-temp.par
  
  return(inds)
}

#density dependent function
#May be the same as previous - just translating n*p into one dice roll rather than across n
host.dens.mort<-function(inds,alpha=NULL){
  pop.size<-nrow(inds)
  
  #convert rate to probability where alpha is the probability of dying by the next time step (based on parasite load)
  #alpha<-1-exp(-alpha*inds[,"No.parasites"])
  
  
  #here each parasite has a probability of killing you but its ability to do so
  #gets amped up by alot the more parasites you have
  #e.g. 1 parasite vs 10 with alpha=0.1 --> 1*0.1=0.1 chance of dying vs 10*0.1=1
  #this assumes that alpha is a probability
  alpha.dens.dep<-ifelse(inds[,"No.parasites"]*alpha>1,1,inds[,"No.parasites"]*alpha)

  #random variable of whether or not you died
  inds[,"host.mort"]<-rbinom(n=pop.size,size=1,p=alpha.dens.dep)
  
  #put in parasites that died
  temp.par.dead<-ifelse(inds[,"host.mort"]>0,inds[,"No.parasites"],0)
  inds[,"parasites.dead"]<-ifelse(inds[,"parasites.dead"]>0,inds[,"parasites.dead"],temp.par.dead)
  
  #update number of parasites
  temp.par<-ifelse(inds[,"host.mort"]>0,0,inds[,"No.parasites"])
  inds[,"No.parasites"]<-temp.par
  
  return(inds)
}


host.pent.mort<-function(inds,alpha=NULL){
  pop.size<-nrow(inds)
  
  #here each parasite has a probability of killing you but its ability to do so
  alpha.dens.dep<-ifelse(inds[,"penetrated.parasites"]*alpha>1,1,inds[,"penetrated.parasites"]*alpha)
  
  #random variable of whether or not you died
  inds[,"host.pent.mort"]<-rbinom(n=pop.size,size=1,p=alpha.dens.dep)
  
  #put in parasites that died
  temp.par.dead<-ifelse(inds[,"host.pent.mort"]>0,inds[,"No.parasites"],0)
  inds[,"parasites.dead"]<-ifelse(inds[,"parasites.dead"]>0,inds[,"parasites.dead"],temp.par.dead)
  
  #update number of parasites
  temp.par<-ifelse(inds[,"host.pent.mort"]>0,0,inds[,"No.parasites"])
  inds[,"No.parasites"]<-temp.par
  
  return(inds)
}


#PARASITE LOST PARASITE MORTALITY
parasite.mortality<-function(inds,mu=NULL){
  
#density independent
pop.size<-nrow(inds)
inds[,"parasite.mort"]<-rbinom(n=pop.size,size=inds[,"No.parasites"],p=mu)

inds[,"No.parasites"]<-inds[,"No.parasites"] - inds[,"parasite.mort"]
return(inds)
}

#PARASITE LOST PARASITE MORTALITY - need to fix function to either increase parasite mort or decrease 
#depending on whether it is parasite limited or host limited
parasite.dens.mort<-function(inds,mu=NULL){
  
  #density independent
  pop.size<-nrow(inds)
  
  
  #probability of parasites dying multiplied by parasite load but need to go from 1:number of parasite load weights
  #if a host has 10 parasites, prob of 1 dying is mu but the 2nd one is 2*mu, 3 is 3*mu. 
  #this assumes that mu is dependent on the parasite load
  
  #thus need to do a forloop and calculate for each individual and run a bernoulli trial per parasite
  #remember if it is the parasite that is limiting that it would probably kill more per parasite
  #if it is the host that is limiting than it may have a weaker immune response with increasing parasite load
  mu.update<-ifelse(inds[,"No.parasites"]*mu>1,1,inds[,"No.parasites"]*mu)

  
  
  #how many parasites died
  inds[,"parasite.mort"]<-rbinom(n=pop.size,size=inds[,"No.parasites"],p=mu.update)
 
  
  inds[,"No.parasites"]<-inds[,"No.parasites"] - inds[,"parasite.mort"]
  return(inds)
}

#Parasite Establishment
#each day unestablished parasites either live, die or develop
#theta and mu are each vectors
parasite.establishment<-function(inds=NULL,theta.C=NULL,mu.C=0){
  
  #theta - probability of parasites developing
  #mu - probability of parasites dying
  #fi - probability of parasites staying
  fi<-1-mu - theta

  
  #create vector for them
  vec<-c(fi,theta,mu)
 
  #create matrix of how many parasites stay and go
  #column is for individual host and row is status (row1=stay,row2=develop,row3=die)
  parasites.temp<-sapply(inds[,"No.C.parasites"],FUN=function(x){rmultinom(1,size=x,prob=vec)})
 
  
  #row 1 needs to be put back in inds[,"No.C.parasites"]
  inds[,"No.C.parasites"]<-parasites.temp[1,]

  #row 2 needs to be moved into inds[,"No.parasites] (but keep previous parasites there)
  inds[,"No.parasites"]<-inds[,"No.parasites"] + parasites.temp[2,]
  
  return(inds)
  }

#########################################################################################

#EXPERIMENTAL MODEL

#arena.size = size of the infection arena (x by x) so 2 => 2x2; type = int, default = 4 x 4
#infection.length = number of timesteps in the infection in hours; type = int, default = 180 minutes (3 hours)

#PROCESSES
#beh.res = whether there is behavioural resistance; type = bool, default = FALSE
#pent.mort= whether there is exposure tolerance (prevent penetration mortality); type = bool, default=F,
#parasite.mort = whether there is parasite mortality (established parasite resistance); type = bool, default = F

#######################################################################################

#SOURCE INFECTION SCRIPT
source("Avoidance_Movement.R")

#time steps
#infection happens on day 0 and right away parasites establish
experimental.model.run<-function(alpha=0.001,mu=1/5,beta=0.5,alpha.pent=0.01,H=100,C=30,arena.size=5,infection.length=90,
                                 beh.res=F,pent.mort=F,parasite.mort=F){

#create inds matrix, ind_hist
inds<-array(data=0,dim=c(H,11))
rownames(inds)<-paste("ind",seq(1,H),sep="")
colnames(inds)<-c("initial.parasites","infection.length","penetrated.parasites","established.parasites","No.parasites","host.mort","host.pent.mort","parasite.mort","time.died","parasites.dead", "t")

#set initial parasite dose
inds[,"initial.parasites"]<-C

#to keep track of simulation history
inds_hist<-NULL;  
  
#simulation set-up
ts<-1
time_steps<-51; #experimental length (maybe should set?)

num.inds<-nrow(inds);

if(ts==1){
  #1) INFECTION & ESTABLISHMENT (BEHAVIOURAL RESISTANCE)
  
  inds<-infection.function(inds=inds,arena.size=arena.size,infection.length=infection.length, beh.res=beh.res,
                     H=H,C=C,beta=beta)#function that infects hosts

  #fill in infection length
  inds[,"infection.length"]<-infection.length
  
  #2) SET NUMBER OF PARASITES
  inds[,"No.parasites"]<-inds[,"established.parasites"]
  
  
  #3) PENETRATION MORTALITY
  if(pent.mort==T){
    inds<-host.pent.mort(inds=inds, alpha=alpha.pent)
    inds[inds[,"host.pent.mort"]>=1,"time.died"]=1 
  }

  inds[,"t"]<-ts #set day 1 for infection
 inds_hist[[1]]<-inds
  
  ts<-ts + 1
}
  #Established parasite resistance (T/F)
# if(parasite.mort==T){
#     mu=1/100
# }

  
  #3) REMAINDER OF EXPERIMENT
  while(ts<time_steps){ #move this while loop
    #inds<-host.mort(inds=inds, alpha=alpha)
    inds<-host.dens.mort(inds=inds, alpha=alpha)
    inds<-parasite.mortality(inds=inds, mu=mu)
    temp.time.died<-ifelse(inds[,"host.mort"]>0,ts,inds[,"time.died"])
    inds[,"time.died"]<-temp.time.died
    inds[,"t"]<-ts
    inds_hist[[ts]]<-data.frame(inds)
    ts<-ts + 1
    
  }
return(inds_hist)

}

##################################################################

#Experimental time, sampling interval and sample after infection

#INPUTS
#inds_hist = a history of experimental infection; type = LIST
#experiment.length= length of experiment in days, default is 20 days; type = integer 
#num.sample.periods = number of sample periods prior to experiment end.
#If 0 then hosts were only sampled at the end of the experiment; type = integer
#post.infect.sample = if sample was taken right after infection (on day 0)
#sample.size = how many inds sampled at each sample period


##################################################################

#####NOTE FOR ME
#For those that died, will need to change t to be time.died or figure out way to amalgamate
#test sampler properly

sampling.time<-function(inds_hist=NULL,experiment.length=20,num.sample.periods=2,post.infect.sample=T, sample.size=10){
 
  if(num.sample.periods==0 | experiment.length==1){
    
    #take only the end of experiment
    
    samples.df<-data.frame(inds_hist[[experiment.length]])
    
    #create one column with total number of parasites counted (for both alive and dead hosts)
    samples.df$total.parasites<-samples.df$No.parasites + samples.df$parasites.dead
    
    samples.df$total.parasites.detected<-rbinom(n=length(samples.df$total.parasites),size=samples.df$total.parasites,prob=.9)
    
    return(samples.df) #returning dataframe of who has been sampled
  }
  
  
  #create a sampling dataframe
  samples.df<-data.frame() 
  index.sampled<-c() #tracks indexes sampled; remove index sampled if resampling can occur
  #total number of hosts
  total.n<-nrow(inds_hist[[1]])
  
  #baseline check to say sample size is too large
  
  #End of experiment data
  temp.df<-inds_hist[[experiment.length]] 
  
  #if the total alive at the last sampling date is greater than or equal to the num.sample.periods*sample.size
  #then we can continue, if False, return that num.sample.periods or sample.size is too large
  
  sufficient.df<-inds_hist[[round(experiment.length/(num.sample.periods + 1))*num.sample.periods]]
  #sufficient.n<- sum(temp.df[,"time.died"]==0) >= (num.sample.periods+1)*sample.size
  sufficient.n<- sum(sufficient.df[,"time.died"]==0) >= num.sample.periods*sample.size
  
  if(sufficient.n==F){
    return(paste("Num.sample.periods or num.sampled is too large for the number of hosts available"))
  }
  
  #if num.sample*sample size is not too large, then proceed with sampling
  
  #SAMPLING AT ESTABLISHMENT
  if(post.infect.sample==T){
    
    #host indexes that are dead
    index.dead<-which(inds_hist[[1]][,"time.died"] !=0)
    
    #who to sample
    temp.index<-sampler.helper(sample.size=sample.size, index.dead=index.dead, index.sampled=index.sampled,total.n=total.n)
    
    #sample the hosts
    samples.df<-inds_hist[[1]][temp.index,]
    
    #update index sampled so hosts cannot be resampled 
    index.sampled<-c(index.sampled,temp.index)
    
    #note that one sampling period has been completed
    num.sample.periods=num.sample.periods -1
    
  }
  
  #For sampling after infection
  
  #calculate which nth day you should sample
  nth.sampling.day<-round(experiment.length/(num.sample.periods + 1)) #need +1 to account for end of experiment
  
  sampling.interval=1 #counter for iterator
  
  while(num.sample.periods>0){
    
    current.sample.day<-nth.sampling.day*sampling.interval #keep track of sampling interval - what day are we sampling
    temp.inds<-inds_hist[[current.sample.day]] #as inds_hist[1] needs to be t=0 (establishment day)
    
    #host indexes that are dead
    index.dead<-which(temp.inds[,"time.died"] !=0)
    
    #which indexes to sample
    temp.index<-sampler.helper(sample.size=sample.size, index.dead=index.dead, index.sampled=index.sampled,total.n=total.n)
    
    #perform sampling and add to samples dataframe
    samples.df<-rbind(samples.df,temp.inds[temp.index,])
    
    #update index sampled so hosts cannot be resampled 
    index.sampled<-c(index.sampled,temp.index)
    
    #note that one sampling period has been completed
    num.sample.periods=num.sample.periods -1
    
    #update how many times sampled
    sampling.interval=sampling.interval + 1
  }
  
 
  
  #can add all who died in experiment too so only exclude those already sampled
  final.day<-temp.df[-index.sampled,]
  
  samples.df<-rbind(samples.df,final.day)
  
  #create one column with total number of parasites counted (for both alive and dead hosts)
  samples.df$total.parasites<-samples.df$No.parasites + samples.df$parasites.dead
  
  samples.df$total.parasites.detected<-rbinom(n=length(samples.df$total.parasites),size=samples.df$total.parasites,prob=.9)
  
  return(samples.df) #returning dataframe of who has been sampled
  
}


#sampler that excludes those that died and who have already been sampled
#helper function for sampling.time
sampler.helper<-function(sample.size=20, index.dead=NULL, index.sampled=NULL,total.n=200){
  
  all.possible.hosts<-seq(1,total.n)
  all.avail.hosts<-all.possible.hosts[!(all.possible.hosts %in% index.dead)  & !(all.possible.hosts %in% index.sampled)]
  samples<-sample(all.avail.hosts,size=sample.size,replace=F)  
  
  return(samples)
}
  


###########################################################################


#reduces list of dataframes to a middle one
#df <-Reduce(rbind, ind_hist)

#create model run with establishment

model.run.establishment<-function(inds=NA,alpha=0.001, mu.C=1/5, theta.C=1/2, mu.L=1/10){
  ts<-1
  time_steps<-15;
  ind_hist<-NULL;
  necropsied.day7<-data.frame()
  while(ts<time_steps){
    inds<-host.dens.mort(inds=inds, alpha=alpha)
    inds<-parasite.mort(inds=inds, mu=mu.L)
    inds<-parasite.establishment(inds=inds,theta.C=theta.C,mu.C=mu.C)
    temp.time.died<-ifelse(inds[,"host.mort"]>0,ts,inds[,"time.died"])
    inds[,"time.died"]<-temp.time.died
    inds[,"t"]<-ts
    
    if(ts==7){
      #need to remove that number from dataset
      df.inds<-data.frame(inds)
      indexes.1<-which((df.inds$age_group==1 & df.inds$time.died==0))
      indexes.2<-which((df.inds$age_group==2 & df.inds$time.died==0))
      
      #to necropsy
      necropsy.1<-indexes.1[1:floor((length(indexes.1)/2))]
      necropsy.2<-indexes.2[1:floor((length(indexes.2)/2))]
      
      inds[c(necropsy.1,necropsy.2),"necropsied"]<-1
      inds[c(necropsy.1,necropsy.2),"necropsied.parasites"]<-inds[c(necropsy.1,necropsy.2),"No.parasites"]
      inds[c(necropsy.1,necropsy.2),"No.parasites"]<-0
    }
    
    ts<-ts + 1
    ind_hist[[ts]]<-data.frame(inds);
  }
  return(ind_hist)
}


#simulation = short runs for 7 days and has establishment, long runs for 30 days and assumes immediate establishment
#function that takes in blank dataset, runs discrete models with parameters and spits out outputs for 10x, 25x, and 50x higher (or % if for beta)
df.generation<-function(inds=NULL,alpha=1e-03,mu.L=1/10, beta=0.35,var.of.interest="alpha",simulation="long"){

  #change 
  if (var.of.interest=="alpha"){
    alpha.1=c(alpha,10*alpha)
    alpha.2=c(25*alpha,50*alpha)
  } else if (var.of.interest=="mu"){
    mu.L.1=c(mu.L,mu.L*1/10)
    mu.L.2=c(mu.L*1/25,mu.L*1/50)
  } else if (var.of.interest=="beta"){
    beta.1=c(beta, beta*1.10)
  }
  
  if(simulation=="long"){
    #UPTAKE
    #juv
    #inds1[1:top.index.juv,"No.parasites"]<-rpois(n=top.index.juv,beta.1[1])
    inds1[1:top.index.juv,"No.parasites"]<-rbinom(n=top.index.juv,size=30,prob=beta.1[1])
    #adult
    #inds1[bottom.index.adult:top.index.adult,"No.parasites"]<-rpois(n=(top.index.adult - top.index.juv),beta.1[2])
    inds1[bottom.index.adult:top.index.adult,"No.parasites"]<-rbinom(n=(top.index.adult - top.index.juv),
                                                                     size=30,prob=beta.1[2])
    
    #UPTAKE
    #juv
    #inds2[1:top.index.juv,"No.parasites"]<-rpois(n=top.index.juv,beta.2[1])
    inds2[1:top.index.juv,"No.parasites"]<-rbinom(n=top.index.juv,size=30,prob=beta.2[1])
    #adult
    #inds2[bottom.index.adult:top.index.adult,"No.parasites"]<-rpois(n=(top.index.adult - top.index.juv),beta.2[2])
    inds2[bottom.index.adult:top.index.adult,"No.parasites"]<-rbinom(n=(top.index.adult - top.index.juv),
                                                                     size=30,prob=beta.2[2])
    
    df<-model.run(inds=inds1,alpha=alpha.1, mu=mu.L.1)
    df.2<-model.run(inds=inds2,alpha=alpha.2, mu=mu.L.2)
    
    df.total.temp<-data.frame(df[[length(df)]])
    df.total2.temp<-data.frame(df.2[[length(df.2)]])
    
  }else if(simulation=="short"){
    #UPTAKE
    inds1[1:top.index.juv,"No.C.parasites"]<-rbinom(n=top.index.juv,size=30,prob=beta.1[1])
    inds1[bottom.index.adult:top.index.adult,"No.C.parasites"]<-rbinom(n=(top.index.adult - top.index.juv),
                                                                     size=30,prob=beta.1[2])
    inds2[1:top.index.juv,"No.C.parasites"]<-rbinom(n=top.index.juv,size=30,prob=beta.2[1])
    inds2[bottom.index.adult:top.index.adult,"No.C.parasites"]<-rbinom(n=(top.index.adult - top.index.juv),
                                                                     size=30,prob=beta.2[2])
    #run model and get ind_hist and necropsied datasets
    df<-model.run.establishment(inds=inds1,alpha=alpha.1, mu.L=mu.L.1)
    df.2<-model.run.establishment(inds=inds2,alpha=alpha.2, mu.L=mu.L.2)
    
    #dataframe of last timestep in ind_hist
    df.total.temp<-data.frame(df[[length(df)]])
    df.total2.temp<-data.frame(df.2[[length(df.2)]])
    
    #df.total.temp<-df.ind.temp
    #df.total2.temp<-df.ind2.temp
    
    df.total.temp$No.parasites<-df.total.temp$No.parasites + df.total.temp$necropsied.parasites
    df.total2.temp$No.parasites<-df.total2.temp$No.parasites + df.total2.temp$necropsied.parasites
  }
  

df.total.temp$starting.parasites<-inds1[,"No.parasites"]
df.total2.temp$starting.parasites<-inds2[,"No.parasites"]

#Add four age classes
df.total.temp$category<-NA
df.total.temp$category[df.total.temp$age_group==1]<-"baseline"
df.total.temp$category[df.total.temp$age_group==2]<-"10x higher"

df.total2.temp$category<-NA
df.total2.temp$category[df.total2.temp$age_group==1]<-"25x higher"
df.total2.temp$category[df.total2.temp$age_group==2]<-"50x higher"

#bind together
df.total<-rbind(df.total.temp,df.total2.temp)

df.total$var.level<-factor(df.total$category, levels=c("baseline", "10x higher", "25x higher", "50x higher"))

#didn't die due to parasitism then has a zero
df.total$censored<-1
df.total$censored[df.total$time.died==0]<-0
df.total$time.to.death<-df.total$time.died

if(simulation=="long"){
  df.total$time.to.death[df.total$time.died==0]<-20
}else if(simulation=="short"){
  df.total$time.to.death[df.total$time.died==0 & df.total$necropsied==0]<-14
  df.total$time.to.death[df.total$time.died==0 & df.total$necropsied==1]<-7
}

#create parasites
df.total$total.parasites<-df.total$No.parasites + df.total$parasites.dead

return(df.total)
}


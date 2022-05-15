##############################################################





##############################################################
#Packages
library(survival)
library(survminer)

#pois.glm.iterator: Poisson GLM Analyses used in baseline analysis to measure infection resistance
pois.glm.iterator.basic.exp<-function(inds.df=NULL){
  
  unique.id.vec<-unique(inds.df$unique.id)
  #set size of dataframe and number of iterations
  n.iter<-length(unique.id.vec)
  
  
  dataset.all<-data.frame(matrix(NA,nrow=numeric(n.iter),ncol=4))
  colnames(dataset.all)<-c("intercept.coef","coef.high","p.intercept",
                           "p.high")
  
  
  
  
  #iterate through list of dataframes
  for (i in 1:n.iter){
    #take function from above with values
    
    df.total<-subset(inds.df,unique.id==unique.id.vec[i])
    
    
    fit.glm.all<-glm(established.parasites~var.level,data=df.total,family="poisson") #total.parasites
    dataset.all[i,1:2]<-summary(fit.glm.all)$coefficients[,1]
    dataset.all[i,3:4]<-summary(fit.glm.all)$coefficients[,4]
    
    
    
    
    #add arena.size, infection.length, dosage, and resist.diff
    dataset.all$infection.length[i]<-df.total$infection.length[1]
    dataset.all$arena.size[i]<-df.total$arena.size[1]
    dataset.all$initial.parasites[i]<-df.total$initial.parasites[1] 
    dataset.all$resist.diff[i]<-df.total$resist.diff[1] 
    
    
    
    
  }
  return(dataset.all)
  
}


#coxph.iterator: Cox Models Analyses used in baseline analysis to measure infection tolerance and
#used in _________
#inds.df: the scenario data generated from the Scenario scripts
#interaction: whether an interaction should be included between parasites and host group
coxph.iterator<-function(inds.df=NULL,interaction=F){
  
  #set size of dataframe and number of iterations
  unique.id.vec<-unique(inds.df$unique.id)
  #set size of dataframe and number of iterations
  n.iter<-length(unique.id.vec)
  
  #create empty datasets
  dataset.all<-data.frame(matrix(NA,nrow=n.iter,ncol=9))
  colnames(dataset.all)<-c("coef.parasites","coef.high",
                           "exp.coef.parasites","exp.coef.high","p.parasites",
                           "p.high","high.CI_low", "high.CI_high","converge.issue")
  
  dataset.all.interact<-data.frame(matrix(NA,nrow=n.iter,ncol=9))
  colnames(dataset.all.interact)<-c("coef.parasites","coef.high","coef.high.int",
                                    "exp.coef.parasites","exp.coef.high","exp.coef.high.int",
                                    "p.parasites","p.high","p.high.int")
  
  #iterate through list of dataframes
  for (i in 1:n.iter){
    #take function from above with values
    
    df.total<-subset(inds.df,unique.id==unique.id.vec[i])
    
    #create time.to.death & censored variable
    df.total$time.to.death<-NA
    df.total$time.to.death[df.total$time.died==0]<-df.total$t[df.total$time.died==0]
    df.total$time.to.death[df.total$time.died!=0]<-df.total$time.died[df.total$time.died!=0]
    
    #status indicator, normally 0=alive, 1=dead from surv notation
    df.total$censored<-NA
    df.total$censored[df.total$time.died==0]<-0 #alive  
    df.total$censored[df.total$time.died!=0]<-1 #died
    
    
    #surv object for model
    surv_object.all<-Surv(time=df.total$time.to.death,event=df.total$censored,type="right") 
    
    
    #create model
    
    if(interaction==F){
      
      options(warn=2)
      aa <- try(fit <- summary(coxph(surv_object.all~total.parasites+var.level,df.total,control = coxph.control(iter.max = 50))))
      if (class(aa) == "try-error"){
        dataset.all$converge.issue[i]<-"yes"
        print(paste("i =", i, "had error"))
        next
      }else{dataset.all$converge.issue[i]<-"no"}
      
      options(warn=-1)
      fit.coxph.all<-coxph(surv_object.all~total.parasites + var.level,df.total,control = coxph.control(iter.max = 50))
      
      dataset.all[i,1:2]<-summary(fit.coxph.all)$coefficients[,1] #double check
      dataset.all[i,3:4]<-summary(fit.coxph.all)$coefficients[,2]
      dataset.all[i,5:6]<-summary(fit.coxph.all)$coefficients[,5]
      dataset.all[i,7:8]<-summary(fit.coxph.all)$conf.int[2,3:4]
      
      #add arena.size, infection.length and resist.diff
      dataset.all$infection.length[i]<-df.total$infection.length[1]
      dataset.all$arena.size[i]<-df.total$arena.size[1]
      dataset.all$tol.diff[i]<-df.total$tol.diff[1]
      dataset.all$initial.parasites[i]<-df.total$initial.parasites[1]
      
      
    }else{
      
      options(warn=2)
      aa <- try(fit <- summary(coxph(surv_object.all~total.parasites*var.level,df.total,control = coxph.control(iter.max = 50))))
      if (class(aa) == "try-error"){
        dataset.all.interact$converge.issue[i]<-"yes"
        print(paste("i =", i, "had error"))
        next
      }else{dataset.all.interact$converge.issue[i]<-"no"}
      
      options(warn=-1)
      fit.coxph.all.interact<-coxph(surv_object.all~total.parasites*var.level,df.total,control = coxph.control(iter.max = 50))
      
      dataset.all.interact[i,1:3]<-summary(fit.coxph.all.interact)$coefficients[,1] 
      dataset.all.interact[i,4:6]<-summary(fit.coxph.all.interact)$coefficients[,2]
      dataset.all.interact[i,7:9]<-summary(fit.coxph.all.interact)$coefficients[,5]
      
      
      #add arena.size, infection.length and resist.diff
      dataset.all$infection.length[i]<-df.total$infection.length[1]
      dataset.all$arena.size[i]<-df.total$arena.size[1]
      dataset.all$tol.diff[i]<-df.total$tol.diff[1]
      
    }
    
    
  }
  if(interaction==T){
    
    final.df<-dataset.all.interact
    
  }else{
    final.df<-dataset.all
  }
  return(final.df)
  
}

#Pois.glm.iterator function runs a Poisson GLM to find differences in parasite loads across hosts (ie evidence of resistance)
#inds.df: the scenario data generated from Scenarios scripts
#process.of.interest: the additional resistance process in the scenarios (none, beh.resist [ie.avoidance] clearance)
#interaction: whether an interaction should be included in the analysis
#clearance with time: whether time since exposure should be included in the analyses (for clearance analysis)
#with.established: whether the known established parasites should be used in the clearance analyses rather than the parasite load recorded
pois.glm.iterator<-function(inds.df=NULL,process.of.interest="none",interaction=F,
                            clearance.with.time=T, with.established=F){
  
  unique.id.vec<-unique(inds.df$unique.id)
  #set size of dataframe and number of iterations
  n.iter<-length(unique.id.vec)
  
  if(process.of.interest=="clearance"){
    #create empty datasets
    if(interaction==T){
      
      dataset.all<-data.frame(matrix(NA,nrow=numeric(n.iter),ncol=16))
      colnames(dataset.all)<-c("intercept.coef","coef.time","coef.high","coef.interact","p.intercept","p.time",
                               "p.high","p.interact","intercept.CI_low",
                               "intercept.CI_high","time.CI_low",
                               "time.CI_high","group.CI_low",
                               "group.CI_high","interact.CI_low",
                               "interact.CI_high") #2.5 %,97.5 % 
      
    }else if(interaction==F){
      
      if(clearance.with.time==T){
        dataset.all<-data.frame(matrix(NA,nrow=numeric(n.iter),ncol=12))
        colnames(dataset.all)<-c("intercept.coef","coef.time","coef.high","p.intercept","p.time",
                                 "p.high","intercept.CI_low",
                                 "intercept.CI_high","time.CI_low",
                                 "time.CI_high","group.CI_low",
                                 "group.CI_high") #2.5 %,97.5 % 
        }
      
      else if(clearance.with.time==F){
        dataset.all<-data.frame(matrix(NA,nrow=numeric(n.iter),ncol=5))
        colnames(dataset.all)<-c("intercept.coef","coef.high","p.intercept",
                                 "p.high","exp.intercept.high") 
      }
    
    }
    
    
    
  }else if(process.of.interest=="beh.resist"){
    dataset.all<-data.frame(matrix(NA,nrow=numeric(n.iter),ncol=5))
    colnames(dataset.all)<-c("intercept.coef","coef.high","p.intercept",
                             "p.high","exp.intercept.high") 
  }
  
  
  
  #iterate through list of dataframes
  for (i in 1:n.iter){
    #take function from above with values
    
    df.total<-subset(inds.df,unique.id==unique.id.vec[i])
    
    #time.to.death
    df.total$time.to.death<-df.total$t
    df.total$time.to.death[df.total$time.died>0]<-df.total$time.died[df.total$time.died>0]
    
    df.total$prop.established<-df.total$established.parasites/df.total$penetrated.parasites
    #create model with dead and alive
    
    df.alive.only<-subset(df.total,time.died!=0)
    
    
    
    if(process.of.interest=="clearance"){
      
      if(interaction==T){
        
        fit.glm.all<-glm(total.parasites~time.to.death*var.level,data=df.total,family="poisson")
        
        
        dataset.all[i,1:4]<-summary(fit.glm.all)$coefficients[,1]
        dataset.all[i,5:8]<-summary(fit.glm.all)$coefficients[,4]
        dataset.all[i,9:10]<-confint(fit.glm.all,1); #intercept conf int.
        dataset.all[i,11:12]<-confint(fit.glm.all,2); #t conf int.
        dataset.all[i,13:14]<-confint(fit.glm.all,3); #var.level conf int.
        dataset.all[i,15:16]<-confint(fit.glm.all,4); #var.level conf int.
        
        
      }else if (interaction==F){
        
        if(clearance.with.time==T){
        fit.glm.all<-glm(total.parasites~time.to.death + var.level,data=df.total,family="poisson")
        
        
        dataset.all[i,1:3]<-summary(fit.glm.all)$coefficients[,1]
        dataset.all[i,4:6]<-summary(fit.glm.all)$coefficients[,4]
        #dataset.all[i,7]<-exp(summary(fit.glm.all)$coefficients[1,1] + summary(fit.glm.all)$coefficients[3,1])
        dataset.all[i,7:8]<-confint(fit.glm.all,1); #intercept conf int.
        dataset.all[i,9:10]<-confint(fit.glm.all,2); #t conf int.
        dataset.all[i,11:12]<-confint(fit.glm.all,3); #var.level conf int.
        
        
        dataset.all$resist.diff[i]<-df.total$resist.diff[1]
        }
        
        else if(clearance.with.time==F){
          if(with.established==T){
            fit.glm.all<-glm(established.parasites~var.level,data=df.total,family="poisson")}
          else if(with.established==F){ #Using the known established parasites
          fit.glm.all<-glm(total.parasites~var.level,data=df.total,family="poisson") }
          
          dataset.all[i,1:2]<-summary(fit.glm.all)$coefficients[,1]
          dataset.all[i,3:4]<-summary(fit.glm.all)$coefficients[,4]
          dataset.all[i,5]<-exp(summary(fit.glm.all)$coefficients[1,1] + summary(fit.glm.all)$coefficients[2,1])
        }
        
        
      }
      
      
    }else if(process.of.interest=="beh.resist"){
      fit.glm.all<-glm(established.parasites~var.level,data=df.total,family="poisson") 
      dataset.all[i,1:2]<-summary(fit.glm.all)$coefficients[,1]
      dataset.all[i,3:4]<-summary(fit.glm.all)$coefficients[,4]
      dataset.all[i,5]<-exp(summary(fit.glm.all)$coefficients[1,1] + summary(fit.glm.all)$coefficients[2,1])
    }
    
    
    #add arena.size, infection.length and resist.diff
    dataset.all$infection.length[i]<-df.total$infection.length[1]
    dataset.all$arena.size[i]<-df.total$arena.size[1]
    
    
    
    
  }
  return(dataset.all)
  
}



binomial.glm.iterator<-function(inds.df=NULL,process.of.interest="none",parasite.type="established",interaction=F, parasites.coef="yes"){
  unique.id.vec<-unique(inds.df$unique.id)
  #set size of dataframe and number of iterations
  n.iter<-length(unique.id.vec)
  
  
  #create empty datasets
  if(process.of.interest=="beh.resist"){
    dataset.all<-data.frame(matrix(NA,nrow=numeric(n.iter),ncol=4))
    colnames(dataset.all)<-c("intercept.coef","coef.high","p.intercept",
                             "p.high") 
    
  }else if (process.of.interest=="exp.tol"){
    if(interaction==F){
      if (parasites.coef=="no"){
      dataset.all<-data.frame(matrix(NA,nrow=numeric(n.iter),ncol=4))
      colnames(dataset.all)<-c("intercept.coef","coef.high","p.intercept", "p.high")
      }
      
      if (parasites.coef=="yes"){
        dataset.all<-data.frame(matrix(NA,nrow=numeric(n.iter),ncol=6))
        colnames(dataset.all)<-c("intercept.coef","coef.parasites","coef.high","p.intercept",
                                 "p.parasites", "p.high")
      }
     
    }else if(interaction==T){
      
      dataset.all<-data.frame(matrix(NA,nrow=numeric(n.iter),ncol=8))
      colnames(dataset.all)<-c("intercept.coef","coef.parasites", "coef.high", "coef.interaction", "p.intercept",
                               "p.parasites", "p.high", "p.interaction")
    }
    
  }      
  
  
  #iterate through list of dataframes
  for (i in 1:n.iter){
    #take function from above with values
    
    df.total<-subset(inds.df,unique.id==unique.id.vec[i])
    
    
    if(process.of.interest=="beh.resist"){
      
      if(parasite.type=="contacted"){
        
        df.total$prop.contacted<-df.total$established.parasites/df.total$penetrated.parasites
        
        fit.glm.all<-glm(prop.contacted~var.level,data=df.total,family="binomial",
                         weight=penetrated.parasites)
        
        
        
      }else if(parasite.type=="established"){
        
        df.total$prop.established<-df.total$established.parasites/df.total$initial.parasites
        
        fit.glm.all<-glm(prop.established~var.level,data=df.total,family="binomial",
                         weights=initial.parasites) #could add weight but produces the same result
        
        
      }else if (parasite.type=="remaining"){
        
        #df.total$infection.length.cont<-as.numeric(paste(df.total$infection.length))
        df.total$prop.remaining<-(df.total$initial.parasites - df.total$penetrated.parasites)/df.total$initial.parasites
        
        fit.glm.all<-glm(prop.remaining~var.level,data=df.total,
                         family="binomial",weights=initial.parasites)
        
        
        
      }
      
      dataset.all[i,1:2]<-summary(fit.glm.all)$coefficients[,1]
      dataset.all[i,3:4]<-summary(fit.glm.all)$coefficients[,4]
      
    }else if(process.of.interest=="exp.tol"){
      ###survived
      if(parasite.type=="established" & interaction==F){
        
        if(parasites.coef=="yes"){
          fit.glm.all<-glm(mort.bin~established.parasites + var.level,data=df.total,family="binomial")
          dataset.all[i,1:3]<-summary(fit.glm.all)$coefficients[,1]
          dataset.all[i,4:6]<-summary(fit.glm.all)$coefficients[,4]
          
        }else if (parasites.coef=="no"){
          fit.glm.all<-glm(mort.bin~var.level,data=df.total,family="binomial")

          dataset.all[i,1:2]<-summary(fit.glm.all)$coefficients[,1]
          dataset.all[i,3:4]<-summary(fit.glm.all)$coefficients[,4]
        }
        
      }else if(parasite.type=="contacted"& interaction==F){
        
        
        if(parasites.coef=="yes"){
          fit.glm.all<-glm(mort.bin~penetrated.parasites + var.level,data=df.total,family="binomial")
          dataset.all[i,1:3]<-summary(fit.glm.all)$coefficients[,1]
          dataset.all[i,4:6]<-summary(fit.glm.all)$coefficients[,4]
          
        }else if (parasites.coef=="no"){
          fit.glm.all<-glm(mort.bin~var.level,data=df.total,family="binomial")
          
          dataset.all[i,1:2]<-summary(fit.glm.all)$coefficients[,1]
          dataset.all[i,3:4]<-summary(fit.glm.all)$coefficients[,4]
          
        }
        
        
      }else if(parasite.type=="established"& interaction==T){
        
        fit.glm.all<-glm(mort.bin~established.parasites*var.level,data=df.total,family="binomial")
        
        dataset.all[i,1:4]<-summary(fit.glm.all)$coefficients[,1]
        dataset.all[i,5:8]<-summary(fit.glm.all)$coefficients[,4]
        
      }else if(parasite.type=="contacted"& interaction==T){
        
        fit.glm.all<-glm(mort.bin~penetrated.parasites*var.level,data=df.total,family="binomial")
        
        dataset.all[i,1:4]<-summary(fit.glm.all)$coefficients[,1]
        dataset.all[i,5:8]<-summary(fit.glm.all)$coefficients[,4]
        
      }
      
      
    }
    
    
    
    #add arena.size, infection.length and resist.diff
    dataset.all$infection.length[i]<-df.total$infection.length[1]
    dataset.all$arena.size[i]<-df.total$arena.size[1]
    
    
  }
  return(dataset.all)
  
}



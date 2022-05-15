#########################################################

#EXPOSURE TOLERANCE ANALYSIS (STATS, TABLES AND FIGURES)

######################################################
#Code to create Figures:5
#Code to create Tables: 5, S5, & S6
library(ggplot2)

source("Stats_Functions.R")

Fig.5.data<-read.csv("data/Exp.Tolerance_Fig5.csv")
Exp.Tol.Tables.data<-read.csv("data/Exp.Tolerance_Tables.csv")


fig.5<-ggplot(Fig.5.data, aes(x = Time, y = percent.died,group=sim.num)) + 
  geom_line(aes(color=var.level),size=0.5) + 
  scale_color_manual(values = c("#CC79A7", "#009E73")) + theme_classic() + #alpha=0.6
  scale_x_continuous(breaks = seq(0,20,5)) + scale_y_continuous(breaks = seq(0,100,10),
                                                                limits = c(0,60)) +
  labs(x = "Experiment Length (# of Timesteps)", y = "% Hosts Dead", 
       color = "Host Tolerance Strategy") + theme_bw() + theme(axis.text=element_text(size=14),axis.title=element_text(size=15),
                                                               legend.title=element_text(size=15),legend.text=element_text(size=12)) +
  theme(strip.background =element_rect(fill="grey")) +
  theme(strip.text = element_text(colour = 'white',size=19,face="bold"))

# 
# png(filename="Figure5.png",height=2300,width=5200,res=400)
# fig.5
# dev.off()

#Function for Tables

coxph.table.function<-function(df=NULL){
  
  #run stats analysis
  df.stats<-coxph.iterator(inds.df=df)
  
  
  #create list with table values
  table.values <- list(mean(df.stats$coef.parasites), mean(df.stats$coef.high),
                       mean(df.stats$exp.coef.parasites), mean(df.stats$exp.coef.high),
                       sum(df.stats$p.parasites<0.05), sum(df.stats$p.high<0.05))  
  names(table.values) <- c("Avg.coef:Parasites","Avg.coef:HighInfect.Tol.Group.",
                           "Avg.Exp.coef:Parasites","Avg.Exp.coef:High Infect.Tol.Group.",
                           "p-value<0.05(x/100):Parasites","p-value<0.05(x/100):High Infect.Tol.Group") 
  return(table.values)
}

#Table 5

#twenty days
Exp.Tolerance.20<-subset(Exp.Tol.Tables.data,experiment.length=="20")
coxph.table.function(df=Exp.Tolerance.20)

#ten days 
Exp.Tolerance.10<-subset(Exp.Tol.Tables.data,experiment.length=="10")
coxph.table.function(df=Exp.Tolerance.10)

#Table S5
Exp.Tolerance.20.1daygone<-subset(Exp.Tolerance.20, time.died !=1)
coxph.table.function(df=Exp.Tolerance.20.1daygone)

#Table S6
Exp.Tolerance.2<-subset(Exp.Tol.Tables.data,experiment.length=="2")

#Isolate to just those who died/survived on day 1
Exp.Tolerance.2.1daydied<-subset(Exp.Tolerance.2, time.died ==1 | time.died==2 | t==1)
Exp.Tolerance.2.1daydied$mort.bin<-0
Exp.Tolerance.2.1daydied$mort.bin[Exp.Tolerance.2.1daydied$time.died==1]<-1



binomial.table.function<-function(df=NULL, parasite.type=NULL, interaction = NULL, parasites.coef = NULL){
  
  #run stats analysis
  df.stats<-binomial.glm.iterator(inds.df=df,
    process.of.interest = "exp.tol",parasite.type = parasite.type,
    interaction=interaction, parasites.coef = parasites.coef)
  
  #With or without Parasites Coef
  
  if (parasites.coef == "no"){
    #create list with table values
    table.values <- list(mean(df.stats$intercept.coef), mean(df.stats$coef.high),
                         sum(df.stats$p.intercept<0.05), sum(df.stats$p.high<0.05))  
    names(table.values) <- c("Avg.coef:Intercept","Avg.coef:HighExp.Tol.Group.",
                             "p-value<0.05(x/100):Intercept","p-value<0.05(x/100):HighExp.Tol.Group")
  
  }else if (parasites.coef =="yes"){
    
    #Interaction
    if (interaction==F){
      print("yes")
      table.values <- list(mean(df.stats$intercept.coef), mean(df.stats$coef.parasites), mean(df.stats$coef.high),
                           sum(df.stats$p.intercept<0.05), sum(df.stats$p.parasites<0.05), sum(df.stats$p.high<0.05))  
      names(table.values) <- c("Avg.coef:Intercept","Avg.coef:Parasites","Avg.coef:HighExp.Tol.Group.",
                               "p-value<0.05(x/100):Intercept","p-value<0.05(x/100):Parasites",
                               "p-value<0.05(x/100):HighExp.Tol.Group")
      
    }else if (interaction ==T){
      
      table.values <- list(mean(df.stats$intercept.coef), mean(df.stats$coef.parasites),
                           mean(df.stats$coef.high), mean(df.stats$coef.interaction),
                           sum(df.stats$p.intercept<0.05), sum(df.stats$p.parasites<0.05), 
                           sum(df.stats$p.high<0.05), sum(df.stats$p.interaction<0.05))  
      names(table.values) <- c("Avg.coef:Intercept","Avg.coef:Parasites","Avg.coef:HighExp.Tol.Group.",
                               "Avg.coef:Interaction","p-value<0.05(x/100):Intercept",
                               "p-value<0.05(x/100):Parasites","p-value<0.05(x/100):HighExp.Tol.Group",
                               "p-value<0.05(x/100):Interaction" )
      
    }
  }
  
  return(table.values)
}


binomial.table.function(df=Exp.Tolerance.2.1daydied,parasite.type = "contacted",
                        interaction=F, parasites.coef = "no")

binomial.table.function(df=Exp.Tolerance.2.1daydied,parasite.type = "contacted",
                        interaction=F, parasites.coef = "yes")

binomial.table.function(df=Exp.Tolerance.2.1daydied,parasite.type = "established",
                        interaction=F, parasites.coef = "yes")

#with an interaction
binomial.table.function(df=Exp.Tolerance.2.1daydied,parasite.type = "contacted",
                        interaction=T, parasites.coef = "yes")

binomial.table.function(df=Exp.Tolerance.2.1daydied,parasite.type = "established",
                        interaction=T, parasites.coef = "yes")



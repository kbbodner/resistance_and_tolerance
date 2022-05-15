##############################################

#AVOIDANCE ANALYSIS


##############################################
#Code to create Figures:4 & S4
#Code to create Tables: 4, S3, & S4
library(ggplot2)
source("Stats_Functions.R") #need process.model

#Bring in data produced from Avoidance_Scenarios.R
Avoidance.full.df<-read.csv("data/Avoidance_5to360mins.csv") 

Avoidance.180to45<-subset(Avoidance.full.df, infection.length !=5 & infection.length !=360)
Avoidance.180to45$infection.length<-factor(Avoidance.180to45$infection.length,levels=c("45", "90", "135", "180"))
Avoidance.180to45$var.level<-factor(Avoidance.180to45$var.level, levels=c("Higher Infection Resistance", "Behavioural Resistance"))
Avoidance.180to45$prop.given.contact<-Avoidance.180to45$established.parasites/Avoidance.180to45$penetrated.parasites
Avoidance.180to45$prop.remaining<-(Avoidance.180to45$initial.parasites - Avoidance.180to45$penetrated.parasites)/Avoidance.180to45$initial.parasites


#Figure 4

fig.4<-ggplot(data=Avoidance.180to45,aes(x=infection.length,y=prop.established,fill=var.level)) + geom_boxplot()+
  facet_wrap(.~arena.size,nrow=1) + theme_bw() + theme(axis.text=element_text(size=14),axis.title=element_text(size=15),
                                                       legend.title=element_text(size=15),legend.text=element_text(size=12)) +
  scale_fill_manual("Host Resistance Strategy",values=c("goldenrod1", "skyblue1"),
                    labels=c("Higher Infection Resistance", "Avoidance")) +
  theme(strip.background =element_rect(fill="grey"))+
  theme(strip.text = element_text(colour = 'white',size=19,face="bold")) + ylim(0,1) + 
  xlab("Exposure Duration (# of Timesteps)") + ylab("Proportion of Established Parasites")

fig.S3<-ggplot(data=Avoidance.180to45,aes(x=infection.length,y=prop.given.contact,fill=var.level)) + geom_boxplot()+
  facet_wrap(.~arena.size,nrow=1) + theme_bw() + theme(axis.text=element_text(size=14),axis.title=element_text(size=15),
                                                       legend.title=element_text(size=15),legend.text=element_text(size=12)) +
  scale_fill_manual("Host Resistance Strategy",values=c("goldenrod1", "skyblue1"),
                    labels=c("Higher Infection Resistance", "Avoidance")) +
  theme(strip.background =element_rect(fill="grey"))+ 
  theme(strip.text = element_text(colour = 'white',size=19,face="bold")) + ylim(0,1) + 
  xlab("Exposure Duration (# of Timesteps)") + ylab("Proportion of Established Parasites Given Contact")


fig.S4<-ggplot(data=Avoidance.180to45,aes(x=infection.length,y=prop.remaining,fill=var.level)) + geom_boxplot()+
  facet_wrap(.~arena.size,nrow=1) + theme_bw() + theme(axis.text=element_text(size=14),axis.title=element_text(size=15),
                                                       legend.title=element_text(size=15),legend.text=element_text(size=12)) +
  scale_fill_manual("Host Resistance Strategy",values=c("goldenrod1", "skyblue1"),
                    labels=c("Higher Infection Resistance", "Avoidance")) +
  theme(strip.background =element_rect(fill="grey"))+ 
  theme(strip.text = element_text(colour = 'white',size=19,face="bold")) + ylim(0,1) + 
  xlab("Exposure Duration (# of Timesteps)") + ylab("Proportion of Remaining Parasites")

#Uncomment to save images
# png(filename="Figure4.png",height=2300,width=5200,res=400)
#  fig.4
#  dev.off()
# 
# png(filename="FigureS3.png",height=2300,width=5200,res=400)
# fig.S3
# dev.off()
# 
# png(filename="FigureS4.png",height=2300,width=5200,res=400)
# fig.S4
# dev.off()

####Tables 4, S3 & S4

#Function for Binomial GLM

binomial.table.function<-function(df=NULL, parasite.type=NULL){
  
  #set-up for stats analysis
  df$unique.id<-df$sim.num
  df$var.level<-factor(df$var.level, levels=c("Higher Infection Resistance", "Behavioural Resistance"))
  
  #run stats analysis  
  df.stats<-binomial.glm.iterator(inds.df=df,process.of.interest="beh.resist",parasite.type=parasite.type)
  
  #create list with table values
  table.values <- list(mean(df.stats$intercept.coef), mean(df.stats$coef.high),
                       sum(df.stats$p.intercept<0.05), sum(df.stats$p.high<0.05))  
  names(table.values) <- c("Avg.coef:Intercept","Avg.coef:Avoidance Group",
                           "p-value<0.05 (x/100):Intercept","p-value<0.05 (x/100):Avoidance Group") 
  return(table.values)
}

#Function to calculate predictions of differences from the Binomial Model
prediction.function<-function(intercept.coef=NA, group.coef=NA){
  
  Yhat_rank_infection.group <- mean(intercept.coef)
  Yhat_rank_beh.group <- mean(intercept.coef) + mean(group.coef)
  
  # Converting the predicted value from log odds into percentage
  infection.group<-exp( Yhat_rank_infection.group ) / ( 1 + exp( Yhat_rank_infection.group ) )
  beh.group<-exp( Yhat_rank_beh.group) / ( 1 + exp( Yhat_rank_beh.group ) )
  
  prop.diff=beh.group - infection.group
  return(prop.diff)
  
}



#set-up data for tables
Avoidance.two.5<-subset(Avoidance.full.df,infection.length=="5" & arena.size=="2x2")
Avoidance.two.45<-subset(Avoidance.full.df,infection.length=="45" & arena.size=="2x2")
Avoidance.two.90<-subset(Avoidance.full.df,infection.length=="90" & arena.size=="2x2")
Avoidance.two.135<-subset(Avoidance.full.df,infection.length=="135" & arena.size=="2x2")
Avoidance.two.180<-subset(Avoidance.full.df,infection.length=="180" & arena.size=="2x2")

Avoidance.six.45<-subset(Avoidance.full.df,infection.length=="45" & arena.size=="6x6")
Avoidance.six.90<-subset(Avoidance.full.df,infection.length=="90" & arena.size=="6x6")
Avoidance.six.135<-subset(Avoidance.full.df,infection.length=="135" & arena.size=="6x6")
Avoidance.six.180<-subset(Avoidance.full.df,infection.length=="180" & arena.size=="6x6")
Avoidance.six.360<-subset(Avoidance.full.df,infection.length=="360" & arena.size=="6x6")

#Table 4
#PROP of PARASITES THAT ESTABLISHED
#(2x2) for 45, 90, 135, and 180 exposure duration

binomial.table.function(Avoidance.two.45, parasite.type="established")
binomial.table.function(Avoidance.two.90, parasite.type="established")
binomial.table.function(Avoidance.two.135, parasite.type="established")
binomial.table.function(Avoidance.two.180, parasite.type="established")

#PROP of PARASITES THAT ESTABLISHED
#(6x6) for 45, 90, 135, and 180 exposure duration

binomial.table.function(Avoidance.six.45, parasite.type="established")
binomial.table.function(Avoidance.six.90, parasite.type="established")
binomial.table.function(Avoidance.six.135, parasite.type="established")
binomial.table.function(Avoidance.six.180, parasite.type="established")

#PROP of PARASITES REMAINING
#(6x6) for 45, 90, 135, and 180 exposure duration

Avoidance.six.45.remaining<-binomial.table.function(Avoidance.six.45, parasite.type="remaining")
Avoidance.six.90.remaining<-binomial.table.function(Avoidance.six.90, parasite.type="remaining")
Avoidance.six.135.remaining<-binomial.table.function(Avoidance.six.135, parasite.type="remaining")
Avoidance.six.180.remaining<-binomial.table.function(Avoidance.six.180, parasite.type="remaining")

Avoidance.six.45.remaining
Avoidance.six.90.remaining
Avoidance.six.135.remaining
Avoidance.six.180.remaining

#Generate predictions of differences between hosts for parasites remaining
prediction.function(intercept.coef=Avoidance.six.45.remaining$`Avg.coef:Intercept`, group.coef=Avoidance.six.45.remaining$`Avg.coef:Avoidance Group`)
prediction.function(intercept.coef=Avoidance.six.90.remaining$`Avg.coef:Intercept`, group.coef=Avoidance.six.90.remaining$`Avg.coef:Avoidance Group`)
prediction.function(intercept.coef=Avoidance.six.135.remaining$`Avg.coef:Intercept`, group.coef=Avoidance.six.135.remaining$`Avg.coef:Avoidance Group`)
prediction.function(intercept.coef=Avoidance.six.180.remaining$`Avg.coef:Intercept`, group.coef=Avoidance.six.180.remaining$`Avg.coef:Avoidance Group`)


#Table S3
#PROP of PARASITES THAT ESTABLISHED
#(2x2) for 5 exposure duration and (6x6) for 360 exposure duration

binomial.table.function(Avoidance.two.5, parasite.type="established")
binomial.table.function(Avoidance.six.360, parasite.type="established")

#Table S4
#PROP of PARASITES THAT ESTABLISHED GIVEN CONTACT
#(2x2) for 45, 90, 135, and 180 exposure duration

binomial.table.function(Avoidance.two.45, parasite.type="contacted")
binomial.table.function(Avoidance.two.90, parasite.type="contacted")
binomial.table.function(Avoidance.two.135, parasite.type="contacted")
binomial.table.function(Avoidance.two.180, parasite.type="contacted")

#PROP of PARASITES THAT ESTABLISHED GIVEN CONTACT
#(6x6) for 45, 90, 135, and 180 exposure duration

binomial.table.function(Avoidance.six.45, parasite.type="contacted")
binomial.table.function(Avoidance.six.90, parasite.type="contacted")
binomial.table.function(Avoidance.six.135, parasite.type="contacted")
binomial.table.function(Avoidance.six.180, parasite.type="contacted")








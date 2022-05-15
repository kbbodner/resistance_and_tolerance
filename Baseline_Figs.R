############################################################

#Baseline Main File - Infection Resistance and Tolerance

############################################################
#Code to create Figures: 3 & S2
library(ggplot2)

#Source in Poisson GLM Results and Cox Model Results produced from Baseline_Stats.R
baseline_Infect_Resist_25<-read.csv("data/Baseline_Resist_25%_Stats.csv")
baseline_Infect_Tol_5<-read.csv("data/Baseline_Tol_5x_Stats.csv") 


#Function for Infection Resistance Figure 3

baseline_infect_resist_fig<-function(df=NULL,sig=50,dosage="LOW"){
  
  #set-up data for plot
  df$arena.size <- factor(df$arena.size, levels = c("10", "6", "2"))
  df.sig<-subset(df,p.high>=sig)
  df.sig$arena.size <- factor(df.sig$arena.size, levels = c("10", "6", "2"))
  
 fig<-ggplot(df, aes(x=infection.length, y=coef.high.mean)) +
    geom_ribbon(data=df.sig,aes(ymin=coef.high.min,ymax=coef.high.max,fill=arena.size),alpha=0.21) +
    geom_point(aes(x=infection.length, y=coef.high.mean,col=arena.size),size=2.5) +
    geom_line(aes(y=coef.high.mean,
                  col=arena.size),data=df,size=0.9) +
    theme_bw() + 
    geom_point(aes(x=infection.length, y=coef.high.mean),data=subset(df,p.high>=sig),color="black",size=2.5,shape=21)+
    geom_point(aes(x=infection.length, y=coef.high.mean),data=subset(df,p.high<sig),color="gold",size=2.5)+
    geom_point(aes(x=infection.length, y=coef.high.mean),data=subset(df,p.high<sig),color="dark grey",size=2.5,shape=21)+
    theme(axis.text=element_text(size=14),axis.title=element_text(size=15), 
          legend.title=element_text(size=15),legend.text=element_text(size=12))+
    theme(strip.background =element_rect(fill="grey")) + ylim(-0.5,0.02)+
    theme(strip.text = element_text(colour = 'white',size=19,face="bold"))  + 
    xlab("Exposure Duration") + ylab("Host group coef.") + 
   ggtitle(paste(dosage, "DOSAGE")) + theme(plot.title = element_text(hjust = 0.5))
  
  return(fig)
  
}

baseline_infect_tol_fig<-function(df=NULL,sig=50){
  
  df$arena.size <- factor(cox.agg.df$arena.size, levels = c("10", "6", "2"))
  df.sig<-subset(df, p.high>=sig)
  
  fig<-ggplot(df, aes(x=infection.length, y=coef.high.mean)) +
    geom_ribbon(data=df.sig,aes(ymin=coef.high.min,ymax=coef.high.max,fill=arena.size),alpha=0.21) +
    geom_point(aes(x=infection.length, y=coef.high.mean,col=arena.size),size=2.5) +
    geom_line(aes(y=coef.high.mean,
                  col=arena.size),data=df,size=0.9) +
    theme_bw() + 
    geom_point(aes(x=infection.length, y=coef.high.mean),data=subset(df,p.high>=sig),color="black",size=2.5,shape=21)+
    geom_point(aes(x=infection.length, y=coef.high.mean),data=subset(df,p.high<sig),color="gold",size=2.5)+
    geom_point(aes(x=infection.length, y=coef.high.mean),data=subset(df,p.high<sig),color="dark grey",size=2.5,shape=21)+
    theme(axis.text=element_text(size=14),axis.title=element_text(size=15), 
          legend.title=element_text(size=15),legend.text=element_text(size=12))+
    theme(strip.background =element_rect(fill="grey")) + ylim(-4.2,0.5) +
    theme(strip.text = element_text(colour = 'white',size=19,face="bold"))  + 
    xlab("Exposure Duration") + ylab("Host group coef.")
  
  return(fig)
  
}


#FIGURE 3
sig=50 #set 50/100 significant coefficients to be a threshold

baseline_Infect_Resist_25.sig<-do.call(data.frame,aggregate(p.high~infection.length + arena.size + initial.parasites,data=baseline_Infect_Resist_25, FUN=function(x) p.sig=sum(x<0.05)))
baseline_Infect_Resist_25.agg<-do.call(data.frame,aggregate(coef.high~infection.length+ arena.size + initial.parasites,data=baseline_Infect_Resist_25,FUN=function(x) c(mean = mean(x), min = min(x), max=max(x)) ))
baseline_Infect_Resist_25.agg<-merge(baseline_Infect_Resist_25.sig,baseline_Infect_Resist_25.agg)  

#subset into low, medium and high dosage
pois.agg.subset.lowDosage.sig<-subset(baseline_Infect_Resist_25.agg,initial.parasites==10)
pois.agg.subset.medDosage.sig<-subset(baseline_Infect_Resist_25.agg,initial.parasites==30)
pois.agg.subset.highDosage.sig<-subset(baseline_Infect_Resist_25.agg,initial.parasites==90)

#create figures
infect.resist_lowDosage<-baseline_infect_resist_fig(df=pois.agg.subset.lowDosage.sig,sig=sig, dosage="LOW")
infect.resist_medDosage<-baseline_infect_resist_fig(df=pois.agg.subset.medDosage.sig,sig=sig, dosage="MEDIUM")
infect.resist_highDosage<-baseline_infect_resist_fig(df=pois.agg.subset.highDosage.sig,sig=sig,dosage="HIGH")


#Infection Tolerance
cox.converge.removed_5<-subset(baseline_Infect_Tol_5,converge.issue=="no")
cox.sig.df_5<-do.call(data.frame,aggregate(p.high~infection.length + arena.size + initial.parasites,data=cox.converge.removed_5, FUN=function(x) p.sig=sum(x<0.05)))
cox.agg.df_5<-do.call(data.frame,aggregate(coef.high~infection.length+arena.size + initial.parasites,data=cox.converge.removed_5,FUN=function(x) c(mean = mean(x), min = min(x), max=max(x)) ))
baseline_Infect_Tol_5.agg<-merge(cox.sig.df_5,cox.agg.df_5)  


#subset into low, medium and high
cox.agg.subset.lowDosage.sig<-subset(baseline_Infect_Tol_5.agg,initial.parasites==10)
cox.agg.subset.medDosage.sig<-subset(baseline_Infect_Tol_5.agg,initial.parasites==30)
cox.agg.subset.highDosage.sig<-subset(baseline_Infect_Tol_5.agg,initial.parasites==90)

#create figures
infect.tol_lowDosage<-baseline_infect_tol_fig(df=cox.agg.subset.lowDosage.sig,sig=sig)
infect.tol_medDosage<-baseline_infect_tol_fig(df=cox.agg.subset.medDosage.sig,sig=sig)
infect.tol_highDosage<-baseline_infect_tol_fig(df=cox.agg.subset.highDosage.sig,sig=sig)

#final figures
fig.3<-ggarrange(infect.resist_lowDosage,infect.resist_medDosage,infect.resist_highDosage,
          infect.tol_lowDosage,infect.tol_medDosage,infect.tol_highDosage,common.legend = T, legend="bottom")

png(filename="Figure3.png",height=2300,width=5200,res=400)
fig.3
dev.off()



#SUPPLEMENTARY MATERIAL##########
baseline_Infect_Resist_10<-read.csv("data/Baseline_Resist_10%_Stats.csv")
baseline_Infect_Resist_75<-read.csv("data/Baseline_Resist_75%_Stats.csv")

baseline_Infect_Tol_2.5<-read.csv("data/Baseline_Tol_2.5x_Stats.csv")
baseline_Infect_Tol_10<-read.csv("data/Baseline_Tol_10x_Stats.csv")

baseline_resist_supp.mat.fig<-function(df_high=NULL,df_low=NULL, sig=50,dosage="LOW"){

  
  df_high$arena.size <- factor(df_high$arena.size, levels = c("10", "6", "2"))
  df_low$arena.size <- factor(df_low$arena.size, levels = c("10", "6", "2"))
  df.sig_high<-subset(df_high, p.high>=sig)
  df.sig_low<-subset(df_low, p.high>=sig)
  
  #coefficient with confidence intervals
  fig<-ggplot(df_high, aes(x=infection.length, y=coef.high.mean)) +
    geom_ribbon(data=df.sig_high,aes(ymin=coef.high.min,ymax=coef.high.max,fill=arena.size),alpha=0.21) +
    geom_ribbon(data=df.sig_low,aes(ymin=coef.high.min,ymax=coef.high.max,fill=arena.size),alpha=0.21) +
    geom_point(aes(x=infection.length, y=coef.high.mean,col=arena.size),size=2.5) +
    geom_line(aes(y=coef.high.mean, col=arena.size),data=df_high,size=0.9) +
    geom_line(aes(y=coef.high.mean,
                  col=arena.size),data=df_low,size=0.9) +
    theme_bw() + 
    geom_point(aes(x=infection.length, y=coef.high.mean),data=subset(df_high,p.high>=sig),color="black",size=2.5,shape=21)+
    geom_point(aes(x=infection.length, y=coef.high.mean),data=subset(df_high,p.high<sig),color="gold",size=2.5)+
    geom_point(aes(x=infection.length, y=coef.high.mean),data=subset(df_high,p.high<sig),color="dark grey",size=2.5,shape=21)+
    geom_point(aes(x=infection.length, y=coef.high.mean),data=subset(df_low,p.high>=sig),color="black",size=2.5,shape=21)+
    geom_point(aes(x=infection.length, y=coef.high.mean),data=subset(df_low,p.high<sig),color="gold",size=2.5)+
    geom_point(aes(x=infection.length, y=coef.high.mean),data=subset(df_low,p.high<sig),color="dark grey",size=2.5,shape=21)+
    theme(axis.text=element_text(size=14),axis.title=element_text(size=15), 
          legend.title=element_text(size=15),legend.text=element_text(size=12))+
    theme(strip.background =element_rect(fill="grey")) + ylim(-1.15,0.06) +
    theme(strip.text = element_text(colour = 'white',size=19,face="bold"))  + 
    xlab("Exposure Duration") + ylab("Host group coef.") + 
    ggtitle(paste(dosage, "DOSAGE")) + theme(plot.title = element_text(hjust = 0.5))
  
  return(fig)
}

baseline_Tol_supp.mat.fig<-function(df_high=NULL,df_low=NULL, sig=50){
  
  df_high$arena.size <- factor(df_high$arena.size, levels = c("10", "6", "2"))
  df_low$arena.size <- factor(df_low$arena.size, levels = c("10", "6", "2"))
  df.sig_high<-subset(df_high, p.high>=sig)
  df.sig_low<-subset(df_low, p.high>=sig)

fig<-ggplot(df_high, aes(x=infection.length, y=coef.high.mean)) +
  geom_ribbon(data=df.sig_low,aes(ymin=coef.high.min,ymax=coef.high.max,fill=arena.size),alpha=0.21) +
  geom_ribbon(data=df.sig_high,aes(ymin=coef.high.min,ymax=coef.high.max,fill=arena.size),alpha=0.21) +
  geom_point(aes(x=infection.length, y=coef.high.mean,col=arena.size),size=2.5) +
  geom_line(aes(y=coef.high.mean, col=arena.size),data=df_low,size=0.9) +
  geom_line(aes(y=coef.high.mean,
                col=arena.size),data=df_high,size=0.9) +
  theme_bw() + 
  geom_point(aes(x=infection.length, y=coef.high.mean),data=subset(df_high,p.high>=sig),color="black",size=2.5,shape=21)+
  geom_point(aes(x=infection.length, y=coef.high.mean),data=subset(df_high,p.high<sig),color="gold",size=2.5)+
  geom_point(aes(x=infection.length, y=coef.high.mean),data=subset(df_high,p.high<sig),color="dark grey",size=2.5,shape=21)+
  geom_point(aes(x=infection.length, y=coef.high.mean),data=subset(df_low,p.high>=sig),color="black",size=2.5,shape=21)+
  geom_point(aes(x=infection.length, y=coef.high.mean),data=subset(df_low,p.high<sig),color="gold",size=2.5)+
  geom_point(aes(x=infection.length, y=coef.high.mean),data=subset(df_low,p.high<sig),color="dark grey",size=2.5,shape=21)+
  theme(axis.text=element_text(size=14),axis.title=element_text(size=15), 
        legend.title=element_text(size=15),legend.text=element_text(size=12))+
  theme(strip.background =element_rect(fill="grey")) + ylim(-4.5,0.5) +
  theme(strip.text = element_text(colour = 'white',size=19,face="bold"))  + 
  xlab("Exposure Duration") + ylab("Host group coef.")

return(fig)

}


#Infection Resistance
sig=50

#Separate by Dosage (10% higher)
baseline_Infect_Resist_10.sig<-do.call(data.frame,aggregate(p.high~infection.length + arena.size + initial.parasites,data=baseline_Infect_Resist_10, FUN=function(x) p.sig=sum(x<0.05)))
baseline_Infect_Resist_10.agg<-do.call(data.frame,aggregate(coef.high~infection.length+ arena.size + initial.parasites,data=baseline_Infect_Resist_10,FUN=function(x) c(mean = mean(x), min = min(x), max=max(x)) ))
baseline_Infect_Resist_10.agg<-merge(baseline_Infect_Resist_10.sig,baseline_Infect_Resist_10.agg)  

#Separate by Dosage (75% higher)
baseline_Infect_Resist_75.sig<-do.call(data.frame,aggregate(p.high~infection.length + arena.size + initial.parasites,data=baseline_Infect_Resist_75, FUN=function(x) p.sig=sum(x<0.05)))
baseline_Infect_Resist_75.agg<-do.call(data.frame,aggregate(coef.high~infection.length+ arena.size + initial.parasites,data=baseline_Infect_Resist_75,FUN=function(x) c(mean = mean(x), min = min(x), max=max(x)) ))
baseline_Infect_Resist_75.agg<-merge(baseline_Infect_Resist_75.sig,baseline_Infect_Resist_75.agg)  


#Subset into low, medium and high dosage
pois.agg.subset.lowDosage.sig_10<-subset(baseline_Infect_Resist_10.agg,initial.parasites==10)
pois.agg.subset.medDosage.sig_10<-subset(baseline_Infect_Resist_10.agg,initial.parasites==30)
pois.agg.subset.highDosage.sig_10<-subset(baseline_Infect_Resist_10.agg,initial.parasites==90)

pois.agg.subset.lowDosage.sig_75<-subset(baseline_Infect_Resist_75.agg,initial.parasites==10)
pois.agg.subset.medDosage.sig_75<-subset(baseline_Infect_Resist_75.agg,initial.parasites==30)
pois.agg.subset.highDosage.sig_75<-subset(baseline_Infect_Resist_75.agg,initial.parasites==90)

resist.supp.mat_lowDosage<-baseline_resist_supp.mat.fig(df_high=pois.agg.subset.lowDosage.sig_75,df_low=pois.agg.subset.lowDosage.sig_10, sig=sig,dosage="LOW")
resist.supp.mat_medDosage<-baseline_resist_supp.mat.fig(df_high=pois.agg.subset.medDosage.sig_75,df_low=pois.agg.subset.medDosage.sig_10, sig=sig,dosage="MEDIUM")
resist.supp.mat_highDosage<-baseline_resist_supp.mat.fig(df_high=pois.agg.subset.highDosage.sig_75,df_low=pois.agg.subset.highDosage.sig_10, sig=sig,dosage="HIGH")




#Infection Tolerance

cox.converge.removed_2.5<-subset(baseline_Infect_Tol_2.5,converge.issue=="no")
cox.sig.df_2.5<-do.call(data.frame,aggregate(p.high~infection.length + arena.size + initial.parasites,data=cox.converge.removed_2.5, FUN=function(x) p.sig=sum(x<0.05)))
cox.agg.df_2.5<-do.call(data.frame,aggregate(coef.high~infection.length+arena.size + initial.parasites,data=cox.converge.removed_2.5,FUN=function(x) c(mean = mean(x), min = min(x), max=max(x)) ))
baseline_Infect_Tol_2.5.agg<-merge(cox.sig.df_2.5,cox.agg.df_2.5) 

cox.converge.removed_10<-subset(baseline_Infect_Tol_10,converge.issue=="no")
cox.sig.df_10<-do.call(data.frame,aggregate(p.high~infection.length + arena.size + initial.parasites,data=cox.converge.removed_10, FUN=function(x) p.sig=sum(x<0.05)))
cox.agg.df_10<-do.call(data.frame,aggregate(coef.high~infection.length+arena.size + initial.parasites,data=cox.converge.removed_10,FUN=function(x) c(mean = mean(x), min = min(x), max=max(x)) ))
baseline_Infect_Tol_10.agg<-merge(cox.sig.df_10,cox.agg.df_10) 


#subset into low, medium and high
cox.agg.subset.lowDosage.sig.2.5<-subset(baseline_Infect_Tol_2.5.agg,initial.parasites==10)
cox.agg.subset.medDosage.sig.2.5<-subset(baseline_Infect_Tol_2.5.agg,initial.parasites==30)
cox.agg.subset.highDosage.sig.2.5<-subset(baseline_Infect_Tol_2.5.agg,initial.parasites==90)

cox.agg.subset.lowDosage.sig.10<-subset(baseline_Infect_Tol_10.agg,initial.parasites==10)
cox.agg.subset.medDosage.sig.10<-subset(baseline_Infect_Tol_10.agg,initial.parasites==30)
cox.agg.subset.highDosage.sig.10<-subset(baseline_Infect_Tol_10.agg,initial.parasites==90)

#create figures
infect.tol_lowDosage<-baseline_Tol_supp.mat.fig(df_high=cox.agg.subset.lowDosage.sig.10,df_low=cox.agg.subset.lowDosage.sig.2.5, sig=sig)
infect.tol_medDosage<-baseline_Tol_supp.mat.fig(df_high=cox.agg.subset.medDosage.sig.10,df_low=cox.agg.subset.medDosage.sig.2.5, sig=sig)
infect.tol_highDosage<-baseline_Tol_supp.mat.fig(df_high=cox.agg.subset.highDosage.sig.10,df_low=cox.agg.subset.highDosage.sig.2.5, sig=sig)



fig.s2<-ggarrange(resist.supp.mat_lowDosage,resist.supp.mat_medDosage,resist.supp.mat_highDosage,
          infect.tol_lowDosage,infect.tol_medDosage,infect.tol_highDosage.cox,common.legend = T)


png(filename="FigureS2.png",height=2300,width=5200,res=400)
fig.S2
dev.off()









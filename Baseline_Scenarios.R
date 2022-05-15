#########################################################

#Basic infection resistance and tolerance Scenarios
#Generating IBM Scenario Data

##########################################################

#Call in IBM and Main Scenario Functions
source("IBM.R")

#Call in experimental process function

############### GENERATING INFECTION RESISTANCE SCENARIOS FOR BASELINE ANALYSIS #####################

set.seed(2021)

top.k=3 # No. of free-living parasites hosts are exposed to in all scenarios (10, 30, 90 parasites) 
top.j=3 # No. of arena size scenarios (2x2, 6x6, 10x10 arena sizes)
top.i=11 # No. of exposure duration scenarios (10-210 exposure duration) by 20 

#No. of infection resistance difference scenarios (10%, 25%, 75% higher)
#sets resistance difference between two host groups
resist.percent.higher=25 #25% Main Analyses
#resist.percent.higher=10 #10% Supplementary
#resist.percent.higher=75  #75% Supplementary

#sets difference in infection resistance
beta.high=0.5*(1+resist.percent.higher/100)

beta.df<-data.frame()

#iterates through different dosages (10,30,90)
C=10 

for(k in 1:top.k){
  #sets infection resistance for second host group for other scenarios
  if(k==2){C=30}else if(k==3){C=90}
  
  print(paste("k equals",k))
  #iterates through arena size scenarios
  arena.temp=2 #2x2 arena size
  for(j in 1:top.j){
    
    #sets arena size 
    arena.size=arena.temp
    
    arena.temp=arena.temp+4
    
    print(paste("j equals",j))
    #iterates through exposure duration scenarios
    for(i in 1:top.i){
      row.num<-top.i*top.j*(k-1)+top.i*(j-1)+i 
      print(row.num)
      
      #sets exposure duration
      infection.length=10*(2*i-1)
      beta.list<-main.experimental.setup(arena.size=arena.size, 
                                         infection.length=infection.length,
                                         var.of.interest="beta",beta.2=beta.high,C=C)
      
      beta.unlisted<-do.call(rbind.data.frame, beta.list)
      beta.unlisted$resist.diff<-paste(resist.percent.higher,"%",sep="")
      beta.df<-rbind(beta.df,beta.unlisted)}}}

beta.df$unique.id<-paste(beta.df$initial.parasites,beta.df$arena.size,beta.df$infection.length,beta.df$sim.number)


#write.csv(beta.df,"data/Baseline_Resist_25%.csv")
#write.csv(beta.df,"data/Baseline_Resist_10%.csv")
#write.csv(beta.df,"data/Baseline_Resist_75%.csv")

############### GENERATING INFECTION TOLERANCE SCENARIOS #######################
set.seed(2021)

top.k=3 # No. of free-living parasites hosts are exposed to in all scenarios (10, 30, 90 parasites) 
top.j=3 # No. of arena size scenarios (2x2, 6x6, 10x10 arena sizes)
top.i=11 # No. of exposure duration scenarios (10-210 exposure duration) by 20    

#No. of infection tolerance difference scenarios (2.5x, 5x, 10x higher)
#sets tolerance difference between two host groups
tol.xtimes.higher=5 #5x Main Analyses
#tol.xtimes.higher=2.5 #2.5x Supplementary
#tol.xtimes.higher=10  #10x Supplementary


alpha.df<-data.frame()

#sets difference in infection tolerance
alpha.high=0.0005*tol.xtimes.higher

#iterates through different dosages (10,30,90)
C=10
for(k in 1:top.k){
  if(k==2){C=30}else if(k==3){C=90}
  
  print(paste("k equals",k))
  #iterates through arena size scenarios
  arena.temp=2 #2x2 arena size
  for(j in 1:top.j){
    #set arena size
    arena.size=arena.temp
    
    arena.temp=arena.temp+4
    
    print(paste("j equals",j))
    #iterates through exposure duration scenarios
    for(i in 1:top.i){
      row.num<-top.i*top.j*(k-1)+top.i*(j-1)+i 
      print(row.num)
      
      #sets exposure duration
      infection.length=10*(2*i-1)
      alpha.list<-main.experimental.setup(arena.size=arena.size, 
                                          infection.length=infection.length,num.sample.periods=0,
                                          var.of.interest="alpha",alpha.2=alpha.high,C=C)
      
      alpha.unlisted<-do.call(rbind.data.frame, alpha.list)
      alpha.unlisted$tol.diff<-paste(tol.xtimes.higher,"x",sep="")
      alpha.df<-rbind(alpha.df,alpha.unlisted)}}
  
}

alpha.df$unique.id<-paste(alpha.df$initial.parasites,alpha.df$arena.size,alpha.df$infection.length,alpha.df$sim.number)

#write.csv(alpha.df,"data/Baseline_Tol_5x.csv")
#write.csv(alpha.df,"data/Baseline_Tol_2.5x.csv")
#write.csv(alpha.df,"data/Baseline_Tol_10x.csv")
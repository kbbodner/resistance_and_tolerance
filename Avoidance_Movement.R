#discrete stochastic random walk model

movement_avoidance<-function(inds_infect,xmax=NULL,ymax=NULL,avoidance=F){
  if(avoidance==T){
    
    #possible host movement
    possible.spots.matrix<-possible.spots(x.coord=inds_infect[1,"x_loc"],y.coord=inds_infect[1,"y_loc"])
    
    No.parasites.per.spot<-c()
    for(i in 1:nrow(possible.spots.matrix)){
      No.parasites.per.spot[i]<-parasite.counter(inds_infect=inds_infect,coord=possible.spots.matrix[i,])
    }
    
    possible.spots.matrix<-cbind(possible.spots.matrix,No.parasites.per.spot)
    
    #remove indices out of bounds of container
    index.x.out<-(which(possible.spots.matrix[,1]>xmax | possible.spots.matrix[,1]<1))*-1
    index.y.out<-(which(possible.spots.matrix[,2]>ymax | possible.spots.matrix[,2]<1))*-1
    
    if(length(c(index.x.out,index.y.out))>0){
      possible.spots.matrix<-possible.spots.matrix[c(index.x.out,index.y.out),]
    }
    
    #find max parasites to set weight
    max.parasites<-max(possible.spots.matrix[,"No.parasites.per.spot"])
    
    
    #create vector to sample from weights choices based on max parasites
    #NO SCALING - USE THIS INSTEAD FOR BASIC WEIGHT
    weight<-(max.parasites+1) - possible.spots.matrix[,"No.parasites.per.spot"]
    
    #scaling by 2x for the weight
    #weight<-2*((max.parasites+1) - possible.spots.matrix[,"No.parasites.per.spot"])
    indices<-seq(1,length(possible.spots.matrix[,"No.parasites.per.spot"]))
    
    #weighted sample
    weighted.sample<-rep(indices,weight)
    
    new.host.spot.index<-sample(x=weighted.sample,size=1,replace=T);
    
    #host.spot
    inds_infect[1,"x_loc"]<-possible.spots.matrix[new.host.spot.index,1] #take x coordinate at index
    inds_infect[1,"y_loc"]<-possible.spots.matrix[new.host.spot.index,2] #take y coordinate at index
    
    #parasite movement
    x_move<-sample(x=c(-1,0,1),size=sum(inds_infect[,"host"]==0),replace=T);
    y_move<-sample(x=c(-1,0,1),size=sum(inds_infect[,"host"]==0),replace=T);
    
    inds_infect[inds_infect[,"host"]==0,"x_loc"]<-inds_infect[inds_infect[,"host"]==0,"x_loc"] + x_move;
    inds_infect[inds_infect[,"host"]==0,"y_loc"]<-inds_infect[inds_infect[,"host"]==0,"y_loc"] + y_move;
    
  }else{
    x_move<-sample(x=c(-1,0,1),size=nrow(inds_infect),replace=T);
    y_move<-sample(x=c(-1,0,1),size=nrow(inds_infect),replace=T);
    
    inds_infect[,"x_loc"]<-inds_infect[,"x_loc"] + x_move;
    inds_infect[,"y_loc"]<-inds_infect[,"y_loc"] + y_move;
  }
  
  #reflection boundary
  inds_infect[,"x_loc"]<-ifelse(inds_infect[,"x_loc"]>xmax,(xmax - 1),inds_infect[,"x_loc"])
  inds_infect[,"x_loc"]<-ifelse(inds_infect[,"x_loc"]<1,2,inds_infect[,"x_loc"])
  
  inds_infect[,"y_loc"]<-ifelse(inds_infect[,"y_loc"]>ymax,(ymax - 1),inds_infect[,"y_loc"])
  inds_infect[,"y_loc"]<-ifelse(inds_infect[,"y_loc"]<1,2,inds_infect[,"y_loc"])
  
  return(inds_infect)
}

#count who is around you
parasite.counter<-function(inds_infect=NULL,coord=NULL){
  possible.parasite.spots<-possible.spots(x.coord=coord[1], y.coord=coord[2])
  if(nrow(inds_infect)==2){
    parasite.location<-inds_infect[2,c("x_loc","y_loc")]
    parasite.x<-parasite.location["x_loc"]
    parasite.y<-parasite.location["y_loc"]
  }else{parasite.location<-inds_infect[2:nrow(inds_infect),c("x_loc","y_loc")]; 
  parasite.x<-parasite.location[,"x_loc"]; parasite.y<-parasite.location[,"y_loc"] }
  
  #putting the parasites
  number.of.parasites<-c()
  for(i in 1:nrow(possible.parasite.spots)){
    in.x<-which(parasite.x %in% possible.parasite.spots[i,1])
    in.y<-which(parasite.y %in% possible.parasite.spots[i,2])
    number.of.parasites[i]<-sum(in.x %in% in.y)
  }
  total.possible.parasites<-sum(number.of.parasites)
  return(total.possible.parasites)
}


#possible spots matrix
possible.spots<-function(x.coord=NULL,y.coord=NULL){
  #possible x movements
  x.spots<-c(x.coord,(x.coord+1),(x.coord-1))
  y.spots<-c(y.coord,(y.coord+1),(y.coord-1))
  #create matrix
  possible.spots<-cbind(c(rep(x.spots[1],3),rep(x.spots[2],3),rep(x.spots[3],3)),rep(y.spots,3))
  return(possible.spots)
  
}

#infection
infection<-function(inds_infect,beta){
  
  same.x.index<-which(inds_infect[,"x_loc"] %in% inds_infect[inds_infect[,"host"]==1,"x_loc"])
  same.y.index<-which(inds_infect[,"y_loc"] %in% inds_infect[inds_infect[,"host"]==1,"y_loc"])
  
  #number of parasites that penetrated host
  contact.parasites<-sum(same.x.index %in% same.y.index)-1
  
  inds_infect[1,"penetrated.parasites"]<-inds_infect[1,"penetrated.parasites"] + contact.parasites
  
  #number established
  inds_infect[1,"established.parasites"]<-inds_infect[1,"established.parasites"] + rbinom(n=1,size=contact.parasites,p=beta)
  
  #remove parasites
  parasite.row.to.remove<--1*same.x.index[same.x.index %in% same.y.index]
  
  #update matrix but remove host from index
  if (length(parasite.row.to.remove)>1){
    inds_infect<-inds_infect[c(parasite.row.to.remove[-1]),]
  }
  return(inds_infect)
}

############################################################

#Main infection Function


#############################################################


infection.function<-function(inds=NULL,arena.size=NULL,infection.length=NULL, beh.res=NULL,
                   H=NULL,C=NULL,beta=NULL){
  
  
  #assign arena size
  xmax=arena.size
  ymax=arena.size
  
  #create dataset for penetration and establishment
  penetration.df<-c()
  establishment.df<-c()
  
  for (i in 1:H){
    inds_infect<-array(dim=c((C+1),5));
    colnames(inds_infect)<-c("host", "x_loc", "y_loc","penetrated.parasites","established.parasites")
    rownames(inds_infect)<-paste("ind",seq(1,(C+1)),sep="")
    
    #assign host or parasite
    inds_infect[1,"host"]<-1
    inds_infect[2:nrow(inds_infect),"host"]<-0
    
    #zero parasites for everyone
    inds_infect[,"penetrated.parasites"]<-0
    inds_infect[,"established.parasites"]<-0
    
    #assigned a random x,y coordinate
    inds_infect[,"x_loc"]<-sample(x=1:xmax,size=nrow(inds_infect),replace=T) 
    inds_infect[,"y_loc"]<-sample(x=1:ymax,size=nrow(inds_infect),replace=T) 
    
    
    ts<-1;
    time_steps<-infection.length #sets length of infection
    
    while(ts<time_steps){
      if(class(inds_infect)[1]=="matrix"){
        inds_infect<-movement_avoidance(inds_infect,xmax=xmax,ymax=ymax,avoidance=beh.res);
        #suppress warnings
        inds_infect<-infection(inds_infect,beta=beta);}
      
      ts<-ts+1
      
      
    }
    
  #need this here as if all parasites penetrated the matrix changes class from matrix to numeric
  if(class(inds_infect)[1]=="numeric"){
    penetration.df[i]<-inds_infect["penetrated.parasites"]
    establishment.df[i]<-inds_infect["established.parasites"]
  }else{
    penetration.df[i]<-inds_infect[1,"penetrated.parasites"]
    establishment.df[i]<-inds_infect[1,"established.parasites"]
  }
      
    
  } 
  #here just fill in to inds
  inds[,"penetrated.parasites"]<-penetration.df
  inds[,"established.parasites"]<-establishment.df
  
  
  return(inds)
  
}



##############################################################
#create time-step video function to show infection

# inds_infect<-array(dim=c(30,5));
# colnames(inds_infect)<-c("host", "x_loc", "y_loc","penetrated.parasites","time_step")
# rownames(inds_infect)<-paste("ind",seq(1,30),sep="")
# 
# xmax=3
# ymax=3
# 
# #assign host or parasite
# inds_infect[1,"host"]<-1
# inds_infect[2:nrow(inds_infect),"host"]<-0
# 
# #zero parasites
# inds_infect[,"No.parasites"]<-0
# 
# #assign a random location
# inds_infect[,"x_loc"]<-sample(x=1:xmax,size=nrow(inds_infect),replace=T) #+ runif(nrow(inds),min=0,max=0.5)
# inds_infect[,"y_loc"]<-sample(x=1:ymax,size=nrow(inds_infect),replace=T) #+ runif(nrow(inds),min=0,max=0.5)
# 
# inds.no.avoid<-inds_infect
# inds.avoid<-inds_infect
# 
# ts<-1;
# time_steps<-24
# 
# inds.avoid_hist<-NULL
# inds.no.avoid_hist<-NULL
# 
# while(ts<time_steps){
#   if(class(inds.no.avoid)=="matrix"){ #this matrix code will have an error as infection removes a row and then the class will change from matrix to numeric (fix code)
#     inds.no.avoid<-infection(inds.no.avoid,p=0.5);
#     inds.no.avoid<-movement_avoidance(inds.no.avoid,xmax=xmax,ymax=ymax,avoidance=F);
#     inds.no.avoid[,"time_step"]<-ts
#   }else{inds.no.avoid["time_step"]<-ts}
#   
#   
#   if(class(inds.avoid)=="matrix"){
#     inds.avoid<-infection(inds.avoid,p=0.8);
#     inds.avoid<-movement_avoidance(inds.avoid,xmax=xmax,ymax=ymax,avoidance=T);
#     inds.avoid[,"time_step"]<-ts}else{inds.avoid["time_step"]<-ts}
#   
#   
#   inds.avoid_hist[[ts]]<-inds.avoid
#   inds.no.avoid_hist[[ts]]<-inds.no.avoid
#   ts<-ts+1
#   
#   
# }
# 


# ##plot of starting location (red is host)
# #png("C:/Users/kbodn/Desktop/resist and tolerate/host.container.png",units="cm", height=20,width=25,res=600)
# plot(x=inds_infect[,"x_loc"], y=inds_infect[,"y_loc"], pch=20,cex=4, 
#      xlim=c(1,xmax),ylim=c(1,ymax), xlab="x location", mar=c(5,5,1,1), 
#      ylab="y location", cex.lab=1.5, cex.axis=1.5)
# points(x=inds_infect[1,"x_loc"], y=inds_infect[1,"y_loc"],col="red",pch=15,cex=6)
# 
# #dev.off()
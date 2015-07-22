## Page 10 of "Occupancy Estimation and Modeling" points out that it is unwise to include habitat
##  covariates for detection probabilities since they also can influence abundance and occupancy.
##  This makes intuitive sense but I can't find where anyone has formally demonstrated it or measured
##  the bias involved in such methods.  Thats what I attempt to do here.

library(RMark)
require(grid)
set.seed(15)


##
## General simulation parameters: 1000 sites, nvisits=3. Psi=.5, p=.5.
##
generalocc<-function(nsites=1000,nvisits=3,realpsi=.5,realp=0.5){
  truth.occ<- runif(nsites,0,1) <= realpsi
  #mean(truth.occ)
  visits<-matrix(ncol=nvisits,nrow=nsites)
  for (i in 1:nsites){
    for (j in 1:nvisits){
      visits[i,j]<-(runif(1,0,1) <= realp)*truth.occ[i]
    }
  }
  det.hist<-data.frame(ch=apply(visits,1,paste0,collapse=""),freq=1,stringsAsFactors = FALSE)
  model<-mark(det.hist,model="Occupancy",silent=TRUE,output=FALSE)
  return(model)
}
nsim=100
result<-lapply(rep(1000,length.out=nsim),generalocc)
x<-data.frame(p=0,psi=0)
for(i in 1:nsim)  x[i,1:2] <- result[[i]]$results$real[,1]
apply(x,2,mean)



##
## Scenario 1.  Environmental variable influences occupancy but not detection.
##         psi= beta0 +beta1(cov1); beta0=0, beta1=1 (Gives realpsi=0.5 if covmean=0), cov1=rnorm(mean=0, sd=1))
## Since logit(psi)= Beta0+Beta1*cov1, Psi= logit^-1(beta0+beta1*cov1)
## The argument covstruct takes three inputs: "none" "psi" "p" or "both"
cov.occ.fun <- function(nsites=1000,nvisits=3,beta0=0,beta1=1,realp=.5,covstruct="none"){
  
  cov1<-rnorm(nsites)
  realpsi<-plogis(beta0+beta1*cov1)
  
  truth.occ<- runif(nsites,0,1) < (realpsi)
  visits<-matrix(ncol=nvisits,nrow=nsites)
  for (i in 1:nsites){
    for (j in 1:nvisits){
      visits[i,j]<-(runif(1,0,1) <= realp)*truth.occ[i]
    }
  }
  #mean(apply(visits,2,mean))
  det.hist<-data.frame(ch=apply(visits,1,paste0,collapse=""),freq=1,cov1=cov1,stringsAsFactors = FALSE)
  if (covstruct=="none"){
    model<-mark(det.hist,model="Occupancy",silent=TRUE,output=FALSE)
  } else if(covstruct=="psi"){
    model<-mark(det.hist,model="Occupancy",silent=TRUE,output=FALSE,
                model.parameters = list(p=list(formula=~1),Psi=list(formula=~cov1)))
  } else if (covstruct=="p") {
    model<-mark(det.hist,model="Occupancy",silent=TRUE,output=FALSE,
                model.parameters = list(p=list(formula=~cov1),Psi=list(formula=~1)))
  } else if (covstruct=="both") {
    model<-mark(det.hist,model="Occupancy",silent=TRUE,output=FALSE,
                model.parameters = list(p=list(formula=~cov1),Psi=list(formula=~cov1)))
  } 
  model$realpsi<-realpsi
  model$realpsi<-realp
  return(model)
}
nsim=3
result <- lapply(rep(1000,length.out=nsim),cov.occ.fun,covstruct="both")

nbeta<-length(rownames(result[[i]]$results$beta))
betavals <-matrix(nrow=nsim,ncol=1+nbeta)
colnames(betavals) <- c(rownames(result[[i]]$results$beta),"MeanCov")
for (i in 1:nsim){
  betavals[i,1:nbeta] <-  result[[i]]$results$beta[,1]
  betavals[i,(nbeta+1)] <- mean(result[[i]]$data$data$cov1)
}
apply(betavals,2,mean)
mean(plogis(betavals[,2]+betavals[,3]*betavals[,4]))

realvals<-data.frame(phat=0,p.se=0,psihat=0,psi.se=0,realpsi=0)
for (i in 1:nsim){
  realvals[i,1:2] <- result[[i]]$results$real[1,1:2]
  realvals[i,3:4] <- result[[i]]$results$real[2,1:2] 
  realvals$realpsi[i] <- mean(result[[i]]$realpsi)
}
apply(realvals,2,mean)


##Loops Simulation
nsim=2
for (i in 1:11) result[[i]] <- lapply(rep(1000,length.out=nsim),cov.occ.fun,covstruct="both",beta1=(i-1)/5)
realvals  <- data.frame(phat=rep(0,length.out=nsim),p.se=0,psihat=0,psi.se=0)
truevals  <- data.frame(truepsi=rep(0,length.out=nsim),truep=0,truebeta1=0)
nbeta     <- length(rownames(result[[1]][[1]]$results$beta))
betavals  <- matrix(nrow=nsim,ncol=1+nbeta)
colnames(betavals) <- c(rownames(result[[1]][[1]]$results$beta),"MeanCov")
graphdata <- cbind(realvals,truevals,betavals)
for (j in 1:11){
  for (i in 1:nsim){
    realvals[i,1:2]       <- result[[j]][[i]]$results$real[1,1:2]
    realvals[i,3:4]       <- result[[j]][[i]]$results$real[2,1:2] 
    truevals$truepsi[i]   <- result[[j]][[i]]$realpsi
    truevals$truep[i]     <- mean(result[[j]][[i]]$realp)
    truevals$truebeta1[i] <- ((j-1)/5)
    betavals[i,1:nbeta]   <- result[[j]][[i]]$results$beta[,1]
    betavals[i,(nbeta+1)] <- mean(result[[j]][[i]]$data$data$cov1)
    duh<-cbind(realvals,truevals,betavals)
    graphdata[j,]<-apply(duh,2,mean)
  }
}

ggplot(graphdata,aes(x=truebeta1))+
  geom_pointrange(aes(y=phat,ymin=phat-p.se,ymax=phat+p.se,color="Red"))+
  geom_pointrange(aes(y=psihat,ymin=psihat-psi.se,ymax=psihat+psi.se,color="Blue"))+theme_classic()+ylim(0,1)

##
## Scenario 2.  Environmental variable influences detection but not occupancy.
##         psi= beta0 +beta1(cov1); beta0=0, beta1=1 (Gives realpsi=0.5 if covmean=0), cov1=rnorm(mean=0, sd=1))
## Since logit(psi)= Beta0+Beta1*cov1, Psi= logit^-1(beta0+beta1*cov1)
## The argument covstruct takes three inputs: "none" "psi" "p" or "both"
cov.det.fun <- function(nsites=1000,nvisits=3,beta0.p=0,beta1.p=1,realpsi=.5,covstruct="none"){
  
  cov1<-rnorm(nsites)
  realp<-plogis(beta0.p+beta1.p*cov1)
  
  truth.occ<- runif(nsites,0,1) < (realpsi)
  visits<-matrix(ncol=nvisits,nrow=nsites)
  for (i in 1:nsites){
    for (j in 1:nvisits){
      visits[i,j]<-(runif(1,0,1) <= realp[i])*truth.occ[i]
    }
  }
  #mean(apply(visits,2,mean))
  det.hist<-data.frame(ch=apply(visits,1,paste0,collapse=""),freq=1,cov1=cov1,stringsAsFactors = FALSE)
  if (covstruct=="none"){
    model<-mark(det.hist,model="Occupancy",silent=TRUE,output=FALSE)
  } else if(covstruct=="psi"){
    model<-mark(det.hist,model="Occupancy",silent=TRUE,output=FALSE,
                model.parameters = list(p=list(formula=~1),Psi=list(formula=~cov1)))
  } else if (covstruct=="p") {
    model<-mark(det.hist,model="Occupancy",silent=TRUE,output=FALSE,
                model.parameters = list(p=list(formula=~cov1),Psi=list(formula=~1)))
  } else if (covstruct=="both") {
    model<-mark(det.hist,model="Occupancy",silent=TRUE,output=FALSE,
                model.parameters = list(p=list(formula=~cov1),Psi=list(formula=~cov1)))
  } 
  model$realpsi<-realpsi
  model$realpsi<-realp
  return(model)
}


nsim=3
result <- lapply(rep(1000,length.out=nsim),cov.det.fun,covstruct="both")

nbeta<-length(rownames(result[[i]]$results$beta))
betavals <-matrix(nrow=nsim,ncol=1+nbeta)
colnames(betavals) <- c(rownames(result[[i]]$results$beta),"MeanCov")
for (i in 1:nsim){
  betavals[i,1:nbeta] <-  result[[i]]$results$beta[,1]
  betavals[i,(nbeta+1)] <- mean(result[[i]]$data$data$cov1)
}
apply(betavals,2,mean)
mean(plogis(betavals[,2]+betavals[,3]*betavals[,4]))

realvals<-data.frame(phat=0,p.se=0,psihat=0,psi.se=0,realpsi=0)
for (i in 1:nsim){
  realvals[i,1:2] <- result[[i]]$results$real[1,1:2]
  realvals[i,3:4] <- result[[i]]$results$real[2,1:2] 
  realvals$realpsi[i] <- mean(result[[i]]$realpsi)
}
apply(realvals,2,mean)


##Loops Simulation
nsim=2
for (i in 1:11) result[[i]] <- lapply(rep(1000,length.out=nsim),cov.occ.fun,covstruct="both",beta1=(i-1)/5)
realvals  <- data.frame(phat=rep(0,length.out=nsim),p.se=0,psihat=0,psi.se=0)
truevals  <- data.frame(truepsi=rep(0,length.out=nsim),truep=0,truebeta1=0)
nbeta     <- length(rownames(result[[1]][[1]]$results$beta))
betavals  <- matrix(nrow=nsim,ncol=1+nbeta)
colnames(betavals) <- c(rownames(result[[1]][[1]]$results$beta),"MeanCov")
graphdata <- cbind(realvals,truevals,betavals)
for (j in 1:11){
  for (i in 1:nsim){
    realvals[i,1:2]       <- result[[j]][[i]]$results$real[1,1:2]
    realvals[i,3:4]       <- result[[j]][[i]]$results$real[2,1:2] 
    truevals$truepsi[i]   <- result[[j]][[i]]$realpsi
    truevals$truep[i]     <- mean(result[[j]][[i]]$realp)
    truevals$truebeta1[i] <- ((j-1)/5)
    betavals[i,1:nbeta]   <- result[[j]][[i]]$results$beta[,1]
    betavals[i,(nbeta+1)] <- mean(result[[j]][[i]]$data$data$cov1)
    duh<-cbind(realvals,truevals,betavals)
    graphdata[j,]<-apply(duh,2,mean)
  }
}

ggplot(graphdata,aes(x=truebeta1))+
  geom_pointrange(aes(y=phat,ymin=phat-p.se,ymax=phat+p.se,color="Red"))+
  geom_pointrange(aes(y=psihat,ymin=psihat-psi.se,ymax=psihat+psi.se,color="Blue"))+theme_classic()+ylim(0,1)


##
## Scenario 3.  Environmental variable influences occupancy and detection equally.  Used for detection.
##         psi= beta0 +beta1(cov1); beta0=0, beta1=1 (Gives realpsi=0.5 if covmean=0), cov1=rnorm(mean=0, sd=1))
## Since logit(psi)= Beta0+Beta1*cov1, Psi= logit^-1(beta0+beta1*cov1)
## The argument covstruct takes three inputs: "none" "psi" or "p"
cov.occ.det.fun <- function(nsites=1000,nvisits=3,beta0.psi=0,beta0.p=0, beta1=1,covstruct="none"){
  
  cov1<-rnorm(nsites)
  realpsi<-plogis(beta0.psi+beta1*cov1)
  realp<-plogis(beta0.p+beta1*cov1)
  truth.occ<- runif(nsites,0,1) < (realpsi)
  visits<-matrix(ncol=nvisits,nrow=nsites)
  for (i in 1:nsites){
    for (j in 1:nvisits){
      visits[i,j]<-(runif(1,0,1) <= realp[i])*truth.occ[i]
    }
  }
  #mean(apply(visits,2,mean))
  det.hist<-data.frame(ch=apply(visits,1,paste0,collapse=""),freq=1,cov1=cov1,stringsAsFactors = FALSE)
  if (covstruct=="none"){
    model<-mark(det.hist,model="Occupancy",silent=TRUE,output=FALSE)
  } else if(covstruct=="psi"){
    model<-mark(det.hist,model="Occupancy",silent=TRUE,output=FALSE,
                model.parameters = list(p=list(formula=~1),Psi=list(formula=~cov1)))
  } else if (covstruct=="p") {
    model<-mark(det.hist,model="Occupancy",silent=TRUE,output=FALSE,
                model.parameters = list(p=list(formula=~cov1),Psi=list(formula=~1)))
  } else if (covstruct=="both") {
    model<-mark(det.hist,model="Occupancy",silent=TRUE,output=FALSE,
                model.parameters = list(p=list(formula=~cov1),Psi=list(formula=~cov1)))
  }
  model$realpsi<-realpsi
  model$realp <- realp
  return(model)
}

nsim=3
confused<-list()
confused <- lapply(rep(1000,length.out=nsim),cov.occ.det.fun,covstruct="none",beta1=(i-1)/2)

nbeta<-length(rownames(result[[i]]$results$beta))
betavals <- matrix(nrow=nsim,ncol=1+nbeta)
colnames(betavals) <- c(rownames(result[[i]]$results$beta),"MeanCov")
for (i in 1:nsim){
  betavals[i,1:nbeta] <-  result[[i]]$results$beta[,1]
  betavals[i,(nbeta+1)] <- mean(result[[i]]$data$data$cov1)
}
apply(betavals,2,mean)
mean(plogis(betavals[,2]+betavals[,3]*betavals[,4]))

realvals<-data.frame(phat=0,p.se=0,psihat=0,psi.se=0,realpsi=0)
for (i in 1:nsim){
  realvals[i,1:2] <- result[[i]]$results$real[1,1:2]
  realvals[i,3:4] <- result[[i]]$results$real[2,1:2] 
  realvals$realpsi[i] <- mean(result[[i]]$realpsi)
}
apply(realvals,2,mean)


##Loops Simulation
nsim=2
for (i in 1:11) result[[i]] <- lapply(rep(1000,length.out=nsim),cov.occ.det.fun,covstruct="both",beta1=(i-1)/5)
realvals  <- data.frame(phat=rep(0,length.out=nsim),p.se=0,psihat=0,psi.se=0)
truevals  <- data.frame(truepsi=rep(0,length.out=nsim),truep=0,truebeta1=0)
nbeta     <- length(rownames(result[[1]][[1]]$results$beta))
betavals  <- matrix(nrow=nsim,ncol=1+nbeta)
colnames(betavals) <- c(rownames(result[[1]][[1]]$results$beta),"MeanCov")
graphdata <- cbind(realvals,truevals,betavals)
for (j in 1:11){
  for (i in 1:nsim){
    realvals[i,1:2]       <- result[[j]][[i]]$results$real[1,1:2]
    realvals[i,3:4]       <- result[[j]][[i]]$results$real[2,1:2] 
    truevals$truepsi[i]   <- mean(result[[j]][[i]]$realpsi)
    truevals$truep[i]     <- mean(result[[j]][[i]]$realp)
    truevals$truebeta1[i] <- ((j-1)/5)
    betavals[i,1:nbeta]   <- result[[j]][[i]]$results$beta[,1]
    betavals[i,(nbeta+1)] <- mean(result[[j]][[i]]$data$data$cov1)
    duh<-cbind(realvals,truevals,betavals)
    graphdata[j,]<-apply(duh,2,mean)
  }
}

ggplot(graphdata,aes(x=truebeta1))+
  geom_pointrange(aes(y=phat,ymin=phat-p.se,ymax=phat+p.se,color="phat"))+
  geom_pointrange(aes(y=psihat,ymin=psihat-psi.se,ymax=psihat+psi.se,color="Psihat"))+theme_classic()+ylim(0,1)+
  geom_line(aes(y=truepsi,color="Real value"))

###CODE FOR OVERNIGHT RUN OF SIMULATION 1, ALL FOUR COV Structures.  Output from each cov structre saved in graphdatanone-both
###
###
###

OvernightSim1Run <- function(){
  nsim=200
  print("Sim 1 Starting none at:")
  print(Sys.time())
  for (i in 1:21) result[[i]] <- lapply(rep(1000,length.out=nsim),cov.occ.fun,covstruct="none",beta1=(i-1)/10)
  realvals  <- data.frame(phat=rep(0,length.out=nsim),p.se=0,psihat=0,psi.se=0)
  truevals  <- data.frame(truepsi=rep(0,length.out=nsim),truep=0,truebeta1=0)
  nbeta     <- length(rownames(result[[1]][[1]]$results$beta))
  betavals  <- matrix(nrow=nsim,ncol=1+nbeta)
  colnames(betavals) <- c(rownames(result[[1]][[1]]$results$beta),"MeanCov")
  graphdata.none <- cbind(realvals,truevals,betavals)
  alldata.none<- list()
  for (j in 1:21){
    for (i in 1:nsim){
      realvals[i,1:2]       <- result[[j]][[i]]$results$real[1,1:2]
      realvals[i,3:4]       <- result[[j]][[i]]$results$real[2,1:2] 
      truevals$truepsi[i]   <- mean(result[[j]][[i]]$realpsi)
      truevals$truep[i]     <- mean(result[[j]][[i]]$realp)
      truevals$truebeta1[i] <- ((j-1)/5)
      betavals[i,1:nbeta]   <- result[[j]][[i]]$results$beta[,1]
      betavals[i,(nbeta+1)] <- mean(result[[j]][[i]]$data$data$cov1)
      duh1<-cbind(realvals,truevals,betavals)
    }
    graphdata.none[j,]<-apply(duh1,2,mean)
    alldata.none[[j]]<-duh1
  }
  
  
  print("Starting p at:")
  print(Sys.time())
  for (i in 1:21) result[[i]] <- lapply(rep(1000,length.out=nsim),cov.occ.fun,covstruct="p",beta1=(i-1)/10)
  realvals  <- data.frame(phat=rep(0,length.out=nsim),p.se=0,psihat=0,psi.se=0)
  truevals  <- data.frame(truepsi=rep(0,length.out=nsim),truep=0,truebeta1=0)
  nbeta     <- length(rownames(result[[1]][[1]]$results$beta))
  betavals  <- matrix(nrow=nsim,ncol=1+nbeta)
  colnames(betavals) <- c(rownames(result[[1]][[1]]$results$beta),"MeanCov")
  graphdata.p <- cbind(realvals,truevals,betavals)
  alldata.p<- list()
  for (j in 1:21){
    for (i in 1:nsim){
      realvals[i,1:2]       <- result[[j]][[i]]$results$real[1,1:2]
      realvals[i,3:4]       <- result[[j]][[i]]$results$real[2,1:2] 
      truevals$truepsi[i]   <- mean(result[[j]][[i]]$realpsi)
      truevals$truep[i]     <- mean(result[[j]][[i]]$realp)
      truevals$truebeta1[i] <- ((j-1)/5)
      betavals[i,1:nbeta]   <- result[[j]][[i]]$results$beta[,1]
      betavals[i,(nbeta+1)] <- mean(result[[j]][[i]]$data$data$cov1)
      duh2<-cbind(realvals,truevals,betavals)
    }
    graphdata.p[j,]<-apply(duh2,2,mean)
    alldata.p[[j]]<-duh2
  }
  
  print("Starting psi at:")
  print(Sys.time())
  for (i in 1:21) result[[i]] <- lapply(rep(1000,length.out=nsim),cov.occ.fun,covstruct="psi",beta1=(i-1)/10)
  realvals  <- data.frame(phat=rep(0,length.out=nsim),p.se=0,psihat=0,psi.se=0)
  truevals  <- data.frame(truepsi=rep(0,length.out=nsim),truep=0,truebeta1=0)
  nbeta     <- length(rownames(result[[1]][[1]]$results$beta))
  betavals  <- matrix(nrow=nsim,ncol=1+nbeta)
  colnames(betavals) <- c(rownames(result[[1]][[1]]$results$beta),"MeanCov")
  graphdata.psi <- cbind(realvals,truevals,betavals)
  alldata.psi<- list()
  for (j in 1:21){
    for (i in 1:nsim){
      realvals[i,1:2]       <- result[[j]][[i]]$results$real[1,1:2]
      realvals[i,3:4]       <- result[[j]][[i]]$results$real[2,1:2] 
      truevals$truepsi[i]   <- mean(result[[j]][[i]]$realpsi)
      truevals$truep[i]     <- mean(result[[j]][[i]]$realp)
      truevals$truebeta1[i] <- ((j-1)/5)
      betavals[i,1:nbeta]   <- result[[j]][[i]]$results$beta[,1]
      betavals[i,(nbeta+1)] <- mean(result[[j]][[i]]$data$data$cov1)
      duh3<-cbind(realvals,truevals,betavals)
    }
    graphdata.psi[j,]<-apply(duh3,2,mean)
    alldata.psi[[j]]<-duh3
  }
  
  print("Starting both at:")
  print(Sys.time())
  for (i in 1:21) result[[i]] <- lapply(rep(1000,length.out=nsim),cov.occ.fun,covstruct="both",beta1=(i-1)/10)
  realvals  <- data.frame(phat=rep(0,length.out=nsim),p.se=0,psihat=0,psi.se=0)
  truevals  <- data.frame(truepsi=rep(0,length.out=nsim),truep=0,truebeta1=0)
  nbeta     <- length(rownames(result[[1]][[1]]$results$beta))
  betavals  <- matrix(nrow=nsim,ncol=1+nbeta)
  colnames(betavals) <- c(rownames(result[[1]][[1]]$results$beta),"MeanCov")
  graphdata.both <- cbind(realvals,truevals,betavals)
  alldata.both<- list()
  for (j in 1:21){
    for (i in 1:nsim){
      realvals[i,1:2]       <- result[[j]][[i]]$results$real[1,1:2]
      realvals[i,3:4]       <- result[[j]][[i]]$results$real[2,1:2] 
      truevals$truepsi[i]   <- mean(result[[j]][[i]]$realpsi)
      truevals$truep[i]     <- mean(result[[j]][[i]]$realp)
      truevals$truebeta1[i] <- ((j-1)/5)
      betavals[i,1:nbeta]   <- result[[j]][[i]]$results$beta[,1]
      betavals[i,(nbeta+1)] <- mean(result[[j]][[i]]$data$data$cov1)
      duh4<-cbind(realvals,truevals,betavals)
      
    }
    graphdata.both[j,]<-apply(duh4,2,mean)
    alldata.both[[j]]<-duh4
  }
  
  print("Finished at:")
  print(Sys.time())
  
  #return(duh1,duh2,duh3,duh4)
  return(list(graphdata.none,graphdata.p,graphdata.psi,graphdata.both,alldata.none,alldata.p,alldata.psi,alldata.both))
}

TotalResult1<-OvernightSim1Run()
save(TotalResult1, file = "TotalResult1.rda")
load("TotalResult1.rda")

se(TotalResult1[[1]]$phat)
mean(TotalResult1[[1]]$p.se)


p1 <-TotalPlot(TotalResult1[[1]][1:21,], "p(.) psi(.)")
p2 <-TotalPlot(TotalResult1[[2]][1:21,], "p(cov) psi(.)")
p3 <-TotalPlot(TotalResult1[[3]][1:21,], "p(.) psi(cov)")
p4 <-TotalPlot(TotalResult1[[4]][1:21,], "p(cov) psi(cov)")
gridExtra::grid.arrange(p1,p2,p3,p4)


TotalPlot<- function(x,plottitle="TITLE") {
  ggplot(x,aes(x=truebeta1))+
    geom_pointrange(aes(y=phat,ymin=phat-se(phat),ymax=phat+se(phat),color="phat"),size=0.60)+
    geom_pointrange(aes(y=psihat,ymin=psihat-psi.se,ymax=psihat+psi.se,color="Psihat"))+theme_classic()+ylim(0,1)+
    geom_line(aes(y=truepsi,color="Real value"))+labs(title=plottitle, x="Beta1")+
    theme(legend.position=c(0,0), legend.justification=c(0,0),text=element_text(size=15),
          plot.margin=unit(c(0,.5,0.2,.5),"cm"), legend.title=element_blank(), legend.background=element_blank(),legend.key=element_blank())
}



###CODE FOR OVERNIGHT RUN OF SIMULATION 1, ALL FOUR COV Structures.  Output from each cov structre saved in graphdatanone-both
###
###
###

OvernightSim2Run <- function(){
  nsim=200
  print("Sim 2 Starting none at:")
  print(Sys.time())
  for (i in 1:21) result[[i]] <- lapply(rep(1000,length.out=nsim),cov.det.fun,covstruct="none",beta1.p=(i-1)/10)
  realvals  <- data.frame(phat=rep(0,length.out=nsim),p.se=0,psihat=0,psi.se=0)
  truevals  <- data.frame(truepsi=rep(0,length.out=nsim),truep=0,truebeta1=0)
  nbeta     <- length(rownames(result[[1]][[1]]$results$beta))
  betavals  <- matrix(nrow=nsim,ncol=1+nbeta)
  colnames(betavals) <- c(rownames(result[[1]][[1]]$results$beta),"MeanCov")
  graphdata.none <- cbind(realvals,truevals,betavals)
  alldata.none<- list()
  for (j in 1:21){
    for (i in 1:nsim){
      realvals[i,1:2]       <- result[[j]][[i]]$results$real[1,1:2]
      realvals[i,3:4]       <- result[[j]][[i]]$results$real[2,1:2] 
      truevals$truepsi[i]   <- mean(result[[j]][[i]]$realpsi)
      truevals$truep[i]     <- mean(result[[j]][[i]]$realp)
      truevals$truebeta1[i] <- ((j-1)/5)
      betavals[i,1:nbeta]   <- result[[j]][[i]]$results$beta[,1]
      betavals[i,(nbeta+1)] <- mean(result[[j]][[i]]$data$data$cov1)
      duh1<-cbind(realvals,truevals,betavals)
    }
    graphdata.none[j,]<-apply(duh1,2,mean)
    alldata.none[[j]]<-duh1
  }
  
  
  print("Starting p at:")
  print(Sys.time())
  for (i in 1:21) result[[i]] <- lapply(rep(1000,length.out=nsim),cov.det.fun,covstruct="p",beta1.p=(i-1)/10)
  realvals  <- data.frame(phat=rep(0,length.out=nsim),p.se=0,psihat=0,psi.se=0)
  truevals  <- data.frame(truepsi=rep(0,length.out=nsim),truep=0,truebeta1=0)
  nbeta     <- length(rownames(result[[1]][[1]]$results$beta))
  betavals  <- matrix(nrow=nsim,ncol=1+nbeta)
  colnames(betavals) <- c(rownames(result[[1]][[1]]$results$beta),"MeanCov")
  graphdata.p <- cbind(realvals,truevals,betavals)
  alldata.p<- list()
  for (j in 1:21){
    for (i in 1:nsim){
      realvals[i,1:2]       <- result[[j]][[i]]$results$real[1,1:2]
      realvals[i,3:4]       <- result[[j]][[i]]$results$real[2,1:2] 
      truevals$truepsi[i]   <- mean(result[[j]][[i]]$realpsi)
      truevals$truep[i]     <- mean(result[[j]][[i]]$realp)
      truevals$truebeta1[i] <- ((j-1)/5)
      betavals[i,1:nbeta]   <- result[[j]][[i]]$results$beta[,1]
      betavals[i,(nbeta+1)] <- mean(result[[j]][[i]]$data$data$cov1)
      duh2<-cbind(realvals,truevals,betavals)
    }
    graphdata.p[j,]<-apply(duh2,2,mean)
    alldata.p[[j]]<-duh2
  }
  
  print("Starting psi at:")
  print(Sys.time())
  for (i in 1:21) result[[i]] <- lapply(rep(1000,length.out=nsim),cov.det.fun,covstruct="psi",beta1.p=(i-1)/10)
  realvals  <- data.frame(phat=rep(0,length.out=nsim),p.se=0,psihat=0,psi.se=0)
  truevals  <- data.frame(truepsi=rep(0,length.out=nsim),truep=0,truebeta1=0)
  nbeta     <- length(rownames(result[[1]][[1]]$results$beta))
  betavals  <- matrix(nrow=nsim,ncol=1+nbeta)
  colnames(betavals) <- c(rownames(result[[1]][[1]]$results$beta),"MeanCov")
  graphdata.psi <- cbind(realvals,truevals,betavals)
  alldata.psi<- list()
  for (j in 1:21){
    for (i in 1:nsim){
      realvals[i,1:2]       <- result[[j]][[i]]$results$real[1,1:2]
      realvals[i,3:4]       <- result[[j]][[i]]$results$real[2,1:2] 
      truevals$truepsi[i]   <- mean(result[[j]][[i]]$realpsi)
      truevals$truep[i]     <- mean(result[[j]][[i]]$realp)
      truevals$truebeta1[i] <- ((j-1)/5)
      betavals[i,1:nbeta]   <- result[[j]][[i]]$results$beta[,1]
      betavals[i,(nbeta+1)] <- mean(result[[j]][[i]]$data$data$cov1)
      duh3<-cbind(realvals,truevals,betavals)
    }
    graphdata.psi[j,]<-apply(duh3,2,mean)
    alldata.psi[[j]]<-duh3
  }
  
  print("Starting both at:")
  print(Sys.time())
  for (i in 1:21) result[[i]] <- lapply(rep(1000,length.out=nsim),cov.det.fun,covstruct="both",beta1.p=(i-1)/10)
  realvals  <- data.frame(phat=rep(0,length.out=nsim),p.se=0,psihat=0,psi.se=0)
  truevals  <- data.frame(truepsi=rep(0,length.out=nsim),truep=0,truebeta1=0)
  nbeta     <- length(rownames(result[[1]][[1]]$results$beta))
  betavals  <- matrix(nrow=nsim,ncol=1+nbeta)
  colnames(betavals) <- c(rownames(result[[1]][[1]]$results$beta),"MeanCov")
  graphdata.both <- cbind(realvals,truevals,betavals)
  alldata.both<- list()
  for (j in 1:21){
    for (i in 1:nsim){
      realvals[i,1:2]       <- result[[j]][[i]]$results$real[1,1:2]
      realvals[i,3:4]       <- result[[j]][[i]]$results$real[2,1:2] 
      truevals$truepsi[i]   <- mean(result[[j]][[i]]$realpsi)
      truevals$truep[i]     <- mean(result[[j]][[i]]$realp)
      truevals$truebeta1[i] <- ((j-1)/5)
      betavals[i,1:nbeta]   <- result[[j]][[i]]$results$beta[,1]
      betavals[i,(nbeta+1)] <- mean(result[[j]][[i]]$data$data$cov1)
      duh4<-cbind(realvals,truevals,betavals)
      
    }
    graphdata.both[j,]<-apply(duh4,2,mean)
    alldata.both[[j]]<-duh4
  }
  
  print("Finished at:")
  print(Sys.time())
  
  #return(duh1,duh2,duh3,duh4)
  return(list(graphdata.none,graphdata.p,graphdata.psi,graphdata.both,alldata.none,alldata.p,alldata.psi,alldata.both))
}

TotalResult2<-OvernightSim2Run()
save(TotalResult2, file = "TotalResult2.rda")

se(TotalResult1[[1]]$phat)
mean(TotalResult1[[1]]$p.se)


p1 <-TotalPlot(TotalResult2[[1]][1:21,], "p(.) psi(.)")
p2 <-TotalPlot(TotalResult2[[2]][1:21,], "p(cov) psi(.)")
p3 <-TotalPlot(TotalResult2[[3]][1:21,], "p(.) psi(cov)")
p4 <-TotalPlot(TotalResult2[[4]][1:21,], "p(cov) psi(cov)")
gridExtra::grid.arrange(p1,p2,p3,p4)


TotalPlot<- function(x,plottitle="TITLE") {
  ggplot(x,aes(x=truebeta1))+
    geom_pointrange(aes(y=phat,ymin=phat-se(phat),ymax=phat+se(phat),color="phat"),size=0.60)+
    geom_pointrange(aes(y=psihat,ymin=psihat-psi.se,ymax=psihat+psi.se,color="Psihat"))+theme_classic()+ylim(0,1)+
    geom_line(aes(y=truepsi,color="Real value"))+labs(title=plottitle, x="Beta1")+
    theme(legend.position=c(0,0), legend.justification=c(0,0),text=element_text(size=15),
          plot.margin=unit(c(0,.5,0.2,.5),"cm"), legend.title=element_blank(), legend.background=element_blank(),legend.key=element_blank())
}



###CODE FOR OVERNIGHT RUN OF SIMULATION 3, ALL FOUR COV Structures.  Output from each cov structre saved in graphdatanone-both
###
###
###

OvernightSim3Run <- function(){
  nsim=200
  print("Sim 3 starting none at:")
  print(Sys.time())
  for (i in 1:21) result[[i]] <- lapply(rep(1000,length.out=nsim),cov.occ.det.fun,covstruct="none",beta1=(i-1)/10)
  realvals  <- data.frame(phat=rep(0,length.out=nsim),p.se=0,psihat=0,psi.se=0)
  truevals  <- data.frame(truepsi=rep(0,length.out=nsim),truep=0,truebeta1=0)
  nbeta     <- length(rownames(result[[1]][[1]]$results$beta))
  betavals  <- matrix(nrow=nsim,ncol=1+nbeta)
  colnames(betavals) <- c(rownames(result[[1]][[1]]$results$beta),"MeanCov")
  graphdata.none <- cbind(realvals,truevals,betavals)
  for (j in 1:21){
    for (i in 1:nsim){
      realvals[i,1:2]       <- result[[j]][[i]]$results$real[1,1:2]
      realvals[i,3:4]       <- result[[j]][[i]]$results$real[2,1:2] 
      truevals$truepsi[i]   <- mean(result[[j]][[i]]$realpsi)
      truevals$truep[i]     <- mean(result[[j]][[i]]$realp)
      truevals$truebeta1[i] <- ((j-1)/5)
      betavals[i,1:nbeta]   <- result[[j]][[i]]$results$beta[,1]
      betavals[i,(nbeta+1)] <- mean(result[[j]][[i]]$data$data$cov1)
      duh<-cbind(realvals,truevals,betavals)
      graphdata.none[j,]<-apply(duh,2,mean)
    }
  }
  
  
  print("Starting p at:")
  print(Sys.time())
  for (i in 1:21) result[[i]] <- lapply(rep(1000,length.out=nsim),cov.occ.det.fun,covstruct="p",beta1=(i-1)/10)
  realvals  <- data.frame(phat=rep(0,length.out=nsim),p.se=0,psihat=0,psi.se=0)
  truevals  <- data.frame(truepsi=rep(0,length.out=nsim),truep=0,truebeta1=0)
  nbeta     <- length(rownames(result[[1]][[1]]$results$beta))
  betavals  <- matrix(nrow=nsim,ncol=1+nbeta)
  colnames(betavals) <- c(rownames(result[[1]][[1]]$results$beta),"MeanCov")
  graphdata.p <- cbind(realvals,truevals,betavals)
  for (j in 1:21){
    for (i in 1:nsim){
      realvals[i,1:2]       <- result[[j]][[i]]$results$real[1,1:2]
      realvals[i,3:4]       <- result[[j]][[i]]$results$real[2,1:2] 
      truevals$truepsi[i]   <- mean(result[[j]][[i]]$realpsi)
      truevals$truep[i]     <- mean(result[[j]][[i]]$realp)
      truevals$truebeta1[i] <- ((j-1)/5)
      betavals[i,1:nbeta]   <- result[[j]][[i]]$results$beta[,1]
      betavals[i,(nbeta+1)] <- mean(result[[j]][[i]]$data$data$cov1)
      duh<-cbind(realvals,truevals,betavals)
      graphdata.p[j,]<-apply(duh,2,mean)
    }
  }
  
  print("Starting psi at:")
  print(Sys.time())
  for (i in 1:21) result[[i]] <- lapply(rep(1000,length.out=nsim),cov.occ.det.fun,covstruct="psi",beta1=(i-1)/10)
  realvals  <- data.frame(phat=rep(0,length.out=nsim),p.se=0,psihat=0,psi.se=0)
  truevals  <- data.frame(truepsi=rep(0,length.out=nsim),truep=0,truebeta1=0)
  nbeta     <- length(rownames(result[[1]][[1]]$results$beta))
  betavals  <- matrix(nrow=nsim,ncol=1+nbeta)
  colnames(betavals) <- c(rownames(result[[1]][[1]]$results$beta),"MeanCov")
  graphdata.psi <- cbind(realvals,truevals,betavals)
  for (j in 1:21){
    for (i in 1:nsim){
      realvals[i,1:2]       <- result[[j]][[i]]$results$real[1,1:2]
      realvals[i,3:4]       <- result[[j]][[i]]$results$real[2,1:2] 
      truevals$truepsi[i]   <- mean(result[[j]][[i]]$realpsi)
      truevals$truep[i]     <- mean(result[[j]][[i]]$realp)
      truevals$truebeta1[i] <- ((j-1)/5)
      betavals[i,1:nbeta]   <- result[[j]][[i]]$results$beta[,1]
      betavals[i,(nbeta+1)] <- mean(result[[j]][[i]]$data$data$cov1)
      duh<-cbind(realvals,truevals,betavals)
      graphdata.psi[j,]<-apply(duh,2,mean)
    }
  }
  
  print("Starting both at:")
  print(Sys.time())
  for (i in 1:21) result[[i]] <- lapply(rep(1000,length.out=nsim),cov.occ.det.fun,covstruct="both",beta1=(i-1)/10)
  realvals  <- data.frame(phat=rep(0,length.out=nsim),p.se=0,psihat=0,psi.se=0)
  truevals  <- data.frame(truepsi=rep(0,length.out=nsim),truep=0,truebeta1=0)
  nbeta     <- length(rownames(result[[1]][[1]]$results$beta))
  betavals  <- matrix(nrow=nsim,ncol=1+nbeta)
  colnames(betavals) <- c(rownames(result[[1]][[1]]$results$beta),"MeanCov")
  graphdata.both <- cbind(realvals,truevals,betavals)
  for (j in 1:21){
    for (i in 1:nsim){
      realvals[i,1:2]       <- result[[j]][[i]]$results$real[1,1:2]
      realvals[i,3:4]       <- result[[j]][[i]]$results$real[2,1:2] 
      truevals$truepsi[i]   <- mean(result[[j]][[i]]$realpsi)
      truevals$truep[i]     <- mean(result[[j]][[i]]$realp)
      truevals$truebeta1[i] <- ((j-1)/5)
      betavals[i,1:nbeta]   <- result[[j]][[i]]$results$beta[,1]
      betavals[i,(nbeta+1)] <- mean(result[[j]][[i]]$data$data$cov1)
      duh<-cbind(realvals,truevals,betavals)
      graphdata.both[j,]<-apply(duh,2,mean)
    }
  }
  
  print("Finished at:")
  print(Sys.time())
  
  return(list(graphdata.none,graphdata.p,graphdata.psi,graphdata.both,alldata.none,alldata.p,alldata.psi,alldata.both))
}

TotalResult3<-OvernightSim3Run()
save(TotalResult3, file = "TotalResult3.rda")

graphdata<-TotalResult[[3]][1:21,]

p1 <-TotalPlot(TotalResult3[[1]][1:21,], "p(.) psi(.)")
p2 <-TotalPlot(TotalResult3[[2]][1:21,], "p(cov) psi(.)")
p3 <-TotalPlot(TotalResult3[[3]][1:21,], "p(.) psi(cov)")
p4 <-TotalPlot(TotalResult3[[4]][1:21,], "p(cov) psi(cov)")
gridExtra::grid.arrange(p1,p2,p3,p4,)


TotalPlot<- function(x,plottitle="TITLE") {
  ggplot(x,aes(x=truebeta1))+
    geom_pointrange(aes(y=phat,ymin=phat-se(phat),ymax=phat+se(phat),color="phat"),size=0.60)+
    geom_pointrange(aes(y=psihat,ymin=psihat-psi.se,ymax=psihat+psi.se,color="Psihat"))+theme_classic()+ylim(0,1)+
    geom_line(aes(y=truepsi,color="Real value"))+labs(title=plottitle, x="Beta1")+
    theme(legend.position=c(0,0), legend.justification=c(0,0),text=element_text(size=15),
          plot.margin=unit(c(0,.5,0.2,.5),"cm"), legend.title=element_blank(), legend.background=element_blank(),legend.key=element_blank())
}


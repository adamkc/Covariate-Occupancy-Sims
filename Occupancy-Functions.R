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
  model<-mark(det.hist,model="Occupancy",silent=TRUE,output=FALSE,delete=TRUE)
  return(model)
}
nsim=10
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
    model<-mark(det.hist,model="Occupancy",silent=TRUE,output=FALSE,delete=TRUE)
  } else if(covstruct=="psi"){
    model<-mark(det.hist,model="Occupancy",silent=TRUE,output=FALSE,delete=TRUE,
                model.parameters = list(p=list(formula=~1),Psi=list(formula=~cov1)))
  } else if (covstruct=="p") {
    model<-mark(det.hist,model="Occupancy",silent=TRUE,output=FALSE,delete=TRUE,
                model.parameters = list(p=list(formula=~cov1),Psi=list(formula=~1)))
  } else if (covstruct=="both") {
    model<-mark(det.hist,model="Occupancy",silent=TRUE,output=FALSE,delete=TRUE,
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
    model<-mark(det.hist,model="Occupancy",silent=TRUE,output=FALSE,delete=TRUE)
  } else if(covstruct=="psi"){
    model<-mark(det.hist,model="Occupancy",silent=TRUE,output=FALSE,delete=TRUE,
                model.parameters = list(p=list(formula=~1),Psi=list(formula=~cov1)))
  } else if (covstruct=="p") {
    model<-mark(det.hist,model="Occupancy",silent=TRUE,output=FALSE,delete=TRUE,
                model.parameters = list(p=list(formula=~cov1),Psi=list(formula=~1)))
  } else if (covstruct=="both") {
    model<-mark(det.hist,model="Occupancy",silent=TRUE,output=FALSE,delete=TRUE,
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
    model<-mark(det.hist,model="Occupancy",silent=TRUE,output=FALSE,delete=TRUE)
  } else if(covstruct=="psi"){
    model<-mark(det.hist,model="Occupancy",silent=TRUE,output=FALSE,delete=TRUE,
                model.parameters = list(p=list(formula=~1),Psi=list(formula=~cov1)))
  } else if (covstruct=="p") {
    model<-mark(det.hist,model="Occupancy",silent=TRUE,output=FALSE,delete=TRUE,
                model.parameters = list(p=list(formula=~cov1),Psi=list(formula=~1)))
  } else if (covstruct=="both") {
    model<-mark(det.hist,model="Occupancy",silent=TRUE,output=FALSE,delete=TRUE,
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
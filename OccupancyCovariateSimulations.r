## Page 10 of "Occupancy Estimation and Modeling" points out that it is unwise to include habitat
##  covariates for detection probabilities since they also can influence abundance and occupancy.
##  This makes intuitive sense but I can't find where anyone has formally demonstrated it or measured
##  the bias involved in such methods.  Thats what I attempt to do here.

library(RMark)
library(ggplot2)
library(grid)
nsim <- 50

###CODE FOR OVERNIGHT RUN OF SIMULATION 0, ALL FOUR COV Structures. 
###
###
###

OvernightSim0Run <- function(){

  output <- list()
  print("Sim 0 Starting none at:")
  print(Sys.time())
  set.seed(17)
  output[[1]] <- LoopedSimulation(covstruct = "none",occ.fun = "general.occ",nsim = nsim)
  
  print("Starting p at:")
  print(Sys.time())
  set.seed(17)
  output[[2]] <- LoopedSimulation(covstruct = "p",occ.fun = "general.occ",nsim = nsim)
  
  print("Starting psi at:")
  print(Sys.time())
  set.seed(17)
  output[[3]] <- LoopedSimulation(covstruct = "psi",occ.fun = "general.occ",nsim = nsim)
  
  print("Starting both at:")
  print(Sys.time())
  set.seed(17)
  output[[4]] <- LoopedSimulation(covstruct = "both",occ.fun = "general.occ",nsim = nsim)
  
  print("Finished at:")
  print(Sys.time())
  
  return(output)
}

TotalResult0<-OvernightSim0Run()
save(TotalResult0, file = "TotalResult0.rda")
#load("TotalResult0.rda")

AICcPlot(TotalResult0,"Simulation 0: No Covariance")

p1 <-TotalPlot(TotalResult0[[1]][[1]][1:21,], "p(.) psi(.)")
p2 <-TotalPlot(TotalResult0[[2]][[1]][1:21,], "p(cov) psi(.)")
p3 <-TotalPlot(TotalResult0[[3]][[1]][1:21,], "p(.) psi(cov)")
p4 <-TotalPlot(TotalResult0[[4]][[1]][1:21,], "p(cov) psi(cov)")
gridExtra::grid.arrange(p1,p2,p3,p4,main=textGrob("Simulation 0: No Covariance",gp=gpar(font=4)))


###CODE FOR OVERNIGHT RUN OF SIMULATION 1, ALL FOUR COV Structures.  Output from each cov structre saved in graphdatanone-both
###
###
###

OvernightSim1Run <- function(){
  output <- list()
  print("Sim 1 Starting none at:")
  print(Sys.time())
  set.seed(17)
  output[[1]] <- LoopedSimulation(covstruct = "none",occ.fun = "cov.occ.fun",nsim = nsim)
  
  print("Starting p at:")
  print(Sys.time())
  set.seed(17)
  output[[2]] <- LoopedSimulation(covstruct = "p",occ.fun = "cov.occ.fun",nsim = nsim)
  
  print("Starting psi at:")
  print(Sys.time())
  set.seed(17)
  output[[3]] <- LoopedSimulation(covstruct = "psi",occ.fun = "cov.occ.fun",nsim = nsim)

  print("Starting both at:")
  print(Sys.time())
  set.seed(17)
  output[[4]] <- LoopedSimulation(covstruct = "both",occ.fun = "cov.occ.fun",nsim = nsim)

  print("Finished at:")
  print(Sys.time())
  
  return(output)
}

TotalResult1<-OvernightSim1Run()
save(TotalResult1, file = "TotalResult1.rda")
#load("TotalResult1.rda")

AICcPlot(TotalResult1,"Simulation 1: Psi~cov")


p1 <-TotalPlot(TotalResult1[[1]][[1]][1:21,], "p(.) psi(.)")
p2 <-TotalPlot(TotalResult1[[2]][[1]][1:21,], "p(cov) psi(.)")
p3 <-TotalPlot(TotalResult1[[3]][[1]][1:21,], "p(.) psi(cov)")
p4 <-TotalPlot(TotalResult1[[4]][[1]][1:21,], "p(cov) psi(cov)")
gridExtra::grid.arrange(p1,p2,p3,p4,main=textGrob("Simulation 1: Psi~cov",gp=gpar(font=4)))




###CODE FOR OVERNIGHT RUN OF SIMULATION 2, ALL FOUR COV Structures.  Output from each cov structre saved in graphdatanone-both
###
###
###

OvernightSim2Run <- function(){

  output <- list()
  print("Sim 2 Starting none at:")
  print(Sys.time())
  set.seed(17)
  output[[1]] <- LoopedSimulation(covstruct = "none",occ.fun = "cov.det.fun",nsim = nsim)
  
  print("Starting p at:")
  print(Sys.time())
  set.seed(17)
  output[[2]] <- LoopedSimulation(covstruct = "p",occ.fun = "cov.det.fun",nsim = nsim)
  
  print("Starting psi at:")
  print(Sys.time())
  set.seed(17)
  output[[3]] <- LoopedSimulation(covstruct = "psi",occ.fun = "cov.det.fun",nsim = nsim)
  
  print("Starting both at:")
  print(Sys.time())
  set.seed(17)
  output[[4]] <- LoopedSimulation(covstruct = "both",occ.fun = "cov.det.fun",nsim = nsim)
  
  print("Finished at:")
  print(Sys.time())
  
  return(output)
}

TotalResult2<-OvernightSim2Run()
save(TotalResult2, file = "TotalResult2.rda")
#load("TotalResult2.rda")


AICcPlot(TotalResult2,"Simulation 2: Psi~cov")



p1 <-TotalPlot(TotalResult2[[1]][[1]][1:21,], "p(.) psi(.)")
p2 <-TotalPlot(TotalResult2[[2]][[1]][1:21,], "p(cov) psi(.)")
p3 <-TotalPlot(TotalResult2[[3]][[1]][1:21,], "p(.) psi(cov)")
p4 <-TotalPlot(TotalResult2[[4]][[1]][1:21,], "p(cov) psi(cov)")
gridExtra::grid.arrange(p1,p2,p3,p4,main=textGrob("Simulation 2: p~cov",gp=gpar(font=4)))


###CODE FOR OVERNIGHT RUN OF SIMULATION 3, ALL FOUR COV Structures.  Output from each cov structre saved in graphdatanone-both
###
###
###

OvernightSim3Run <- function(){

  output <- list()
  print("Sim 3 Starting none at:")
  print(Sys.time())
  set.seed(17)
  output[[1]] <- LoopedSimulation(covstruct = "none",occ.fun = "cov.occ.det.fun",nsim = nsim)
  
  print("Starting p at:")
  print(Sys.time())
  set.seed(17)
  output[[2]] <- LoopedSimulation(covstruct = "p",occ.fun = "cov.occ.det.fun",nsim = nsim)
  
  print("Starting psi at:")
  print(Sys.time())
  set.seed(17)
  output[[3]] <- LoopedSimulation(covstruct = "psi",occ.fun = "cov.occ.det.fun",nsim = nsim)
  
  print("Starting both at:")
  print(Sys.time())
  set.seed(17)
  output[[4]] <- LoopedSimulation(covstruct = "both",occ.fun = "cov.occ.det.fun",nsim = nsim)
  
  print("Finished at:")
  print(Sys.time())
  
  return(output)
}


TotalResult3<-OvernightSim3Run()
save(TotalResult3, file = "TotalResult3.rda")
#load("TotalResult3.rda")

AICcPlot(TotalResult3,"Simulation 3: Psi~cov")

p1 <-TotalPlot(TotalResult3[[1]][[1]][1:21,], "p(.) psi(.)")
p2 <-TotalPlot(TotalResult3[[2]][[1]][1:21,], "p(cov) psi(.)")
p3 <-TotalPlot(TotalResult3[[3]][[1]][1:21,], "p(.) psi(cov)")
p4 <-TotalPlot(TotalResult3[[4]][[1]][1:21,], "p(cov) psi(cov)")
gridExtra::grid.arrange(p1,p2,p3,p4,main=textGrob("Simulation 3: Psi~cov,p~cov",gp=gpar(font=4)))



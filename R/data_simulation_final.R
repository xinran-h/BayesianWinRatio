data.simulation<- function(N.sim,N.max,mu.trt,Sigma.trt,
                           mu.ctrl,Sigma.ctrl,cens_upper
){
  myData<-array(dim=c(N.max,7,N.sim));   
  
  for (j in 1:N.sim){
    
    # Simulate randomization groups
    rand.list = sample(c(0,1), N.max, replace = T)
    
    # Simulatetime to recurrence and death
    trt.events <- 	as.data.frame(MASS::mvrnorm(sum(rand.list),
                                         mu= mu.trt,
                                         Sigma= Sigma.trt))     
    trt.events$id =  which(rand.list==1); trt.events$grp = 1; 
    
    ctrl.events <- as.data.frame(MASS::mvrnorm(sum(rand.list==0),
                                         mu= mu.ctrl,
                                         Sigma= Sigma.ctrl))     
    
    ctrl.events$id = which(rand.list==0); ctrl.events$grp = 0
    
    # simulate time to censor
    time.censor = stats::runif(N.max, 0, cens_upper); 
    
    eventsall = rbind(trt.events, ctrl.events); eventsall = eventsall[order(eventsall$id),]
    
    myData[,1,j]<- exp(eventsall[,1])                 # time to recurrence 
    myData[,2,j]<- exp(eventsall[,2])                      # time to death
    myData[,3,j]<- time.censor
    myData[,4,j]<- eventsall[,4]                      # treatment group assnments
    myData[,5,j]<-  eventsall[,3]                     # id
  }
  
  # Create dimension names
  dim_names <- list(
    NULL,
    Variables = c("recurrence_t", "death_t", "censor_t",   "group", "id", "delta1", "delta2"),
    Simulations = paste("Simulation", 1:N.sim)
  )
  
  # Assign the dimension names to the array
  dimnames(myData) <- dim_names
  return(myData)
}
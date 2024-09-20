#' data.simulation
#' 
#' This function simulates data. 
#' 
#' @param N.sim Number of simulations.
#' @param N.max   Maximum number of patients to enroll.
#' @param mu.trt A vector containing the mean time to each event (logarithm) for the treatment arm.
#' @param Sigma.trt A vector representing the variance-covariance matrix of the time to each event (logarithm) for the treatment arm. 
#' @param mu.ctrl A vector containing the mean time to each event (logarithm) for the control arm.
#' @param Sigma.ctrl A vector representing the variance-covariance matrix of the time to each event (logarithm) for the control arm. 
#' @param cens_upper Upper limit for the censoring time, assuming that the censoring time is generated from Uniform(0, cens_upper).
#' @returns A list with the following components:\tabular{ll}{
#'    \code{trial.stop} \tab A value of 1 or 0, 1 = trial stopped and 0 = not stopped.  \cr
#'    \tab \cr
#'    \code{trialER.stop} \tab A value of 1 or 0, 1 = trial stopped early and 0 = not stopped early. \cr
#'    \tab \cr
#'    \code{pts.stop} \tab A numeric value represents the actual sample size used. \cr
#'    \tab \cr
#'    \code{probs} \tab A numeric value, which is the posterior probability of \eqn{(\widehat{WR} > eta)}.  \cr
#'    \tab \cr
#'    \code{WR} \tab A numeric value, which is the estimated posterior win ratio. \cr    
#' }
#' @examples
#' \dontrun{
#' data.simulation(N.sim = 1000, N.max = 20,
#'                mu.trt = c(0.2,0.3), Sigma.trt = matrix(c(1,0.5,0.5,1), nrow = 2, byrow = T),
#'                mu.ctrl = c(0.2,0.3),Sigma.ctrl = matrix(c(1,0.5,0.5,1), nrow = 2, byrow = T),
#'                cens_upper = 5)
#' }

#' @export
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
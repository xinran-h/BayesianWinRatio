#' OCC.Table
#' 
#' This is the main function of BayesWinRatio. 
#' This function perform Bayesian futility monitoring based on one simulation data or one real data,  and returns the operating characteristics.
#' 
#' @useDynLib BayesianWinRatio, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @exportPattern "^[[:alpha:]]+"
#' 
#' @param myData A matrix with N.max rows and 7 columns. The first five columns represent: time to recurrence, time to death, time to censor, treatment arm, id. Leave the last two columns blank. Note that time to recurrence and time to death cannot be NA. If these are NA due to censoring, please use a large number to represent the time. For the treatment arm, use 1 for the treatment arm and 0 for the control arm.
#' @param N.max   Maximum number of patients to enroll.
#' @param design A numeric value indicating the type of design. 1 = proposed design, 2 = time to recurrence design, 3= time to death design, 4 = time to first event design.
#' @param cohort  Interim cohort, which is a numeric vector of the number of patients enrolled at each interim look.
#' @param recruit.int Recruitment interval.
#' @param m0  Prior mean for mu
#' @param L0  Prior variance for mu 
#' @param v0  For the proposed design, this is the prior degrees of freedom for Sigma. For the traditional designs, v0/2 is the prior shape for Sigma.
#' @param S0  For the proposed design, this is the prior scale matrix for Sigma. For the traditional designs, v0*S0/2 is the prior scale for Sigma.
#' @param time_max The upper limit for the recurrence and death time sampled from truncated normal. This will set the upper limit to to time_max rather than Inf.
#' @param eta A pre-specified lower bound of acceptable performance based on historical information.
#' @param lambda  Cutoff parameter.
#' @param thin_MCMC  Thinning degree.
#' @param Niter Number of iterations for gibbs sampler.
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
#' OCC.Table(myData = data,
#' N.max = 100,
#' design = 1, 
#' cohort = c(40,60,80),
#' recruit.int  = 0.25,
#' m0 = c(0,0),L0 = diag(10^6, 2), 
#' v0 = 4,S0 = diag(10^(-6), 2),
#' time_max = 10,eta = 1.5,lambda = 0.25,
#'  thin_MCMC = 5,Niter = 100000)
#' }




#' @export
OCC.Table<- function(myData,N.max,design, cohort, recruit.int,
                     m0,L0, v0, S0,time_max,eta, lambda, thin_MCMC,
                     Niter){
  WR = NA
    
    tryCatch({
  Time.entry <-Time.current <-n.current <-j.cohort<-0;
  trial.stop<- trialER.stop<-pts.stop<-  0;
  
  
  while (n.current < N.max &&  trial.stop !=1 && j.cohort < length(cohort))
  {
    j.cohort<-j.cohort+1;
    n.current<- cohort[j.cohort];
    Time.entry <- recruit.int*(c(1:n.current)-1)        ; ## take into account recruit interval;
    Time.current <- recruit.int*(n.current-1)  ; 
    
    temp = myData[1:n.current-1,] 
    if (design == 1) {
      
      currentData =  temp
      
      result = winratio(currentdd = currentData,
                        n_current = n.current,
                        Time_current = Time.current,
                        Time_entry = Time.entry,
                        m0 = m0,
                        L0 = L0,
                        v0 = v0,
                        S0 = S0,
                        time_max = time_max,
                        eta = eta,
                        lambda = lambda,
                        N.max = N.max,
                        thin_MCMC = thin_MCMC,
                        Niter = Niter)
      probs = result$probs;
      cutoff = result$cutoff;
      WR = result$WR;
      if (!is.na(probs) && probs <= cutoff) { trial.stop<-1; trialER.stop <-1 }
    }else if (design %in% c(2,3,4)){
      # update actual censoring time at this interim look
      temp[,3] = pmin(Time.current-Time.entry[-n.current], temp[,3], temp[,2])
      if(design %in% c(2,3)){
        currentData = temp[,c(design-1, 3:6)]
      }else{
        currentData = data.frame(evnt.time = pmin(temp[,1], temp[,2]), 
                                 temp[,c(3:6)])
        currentData = as.matrix(currentData )
      }
      
      currentData[,5] <- as.integer(currentData[,1] <= currentData[,2])
      
      # update event time to reflect the observed data at this interim look 
      currentData[,1] = ifelse(currentData[,5]==1, currentData[,1], NA)
      
      # number enrolled to each group
      n.current.trt = sum(currentData[,3])
      n.current.ctrl =  n.current - 1 - n.current.trt 
      
      # separate data to trt and control
      currentData.trt = currentData[currentData[,3]==1,]
      currentData.ctrl = currentData[currentData[,3]==0,]

      # update theta
      trt.post = update_theta_univariate(N_iter = Niter, dd = currentData.trt, 
                                                     n = n.current.trt, L0 = L0,m0 = m0, v0 = v0, S0 = S0,
                                                     time_max = time_max)
      ctrl.post = update_theta_univariate(N_iter = Niter, dd = currentData.ctrl, 
                                                      n = n.current.ctrl, L0 = L0, m0 = m0, v0 = v0, S0 = S0,
                                                      time_max = time_max)
      burn_MCMC = as.integer(0.3*Niter)
      idxs <- seq(burn_MCMC, Niter, by = thin_MCMC)
      
      probs = mean(trt.post$theta[ idxs,1] - ctrl.post$theta[ idxs,1] > eta)
      cutoff =  (n.current/N.max)**lambda
      if (!is.na(probs) && probs <= cutoff) { trial.stop<-1; trialER.stop <-1 }
   }
  }
  
  pts.stop<-n.current;
  if (pts.stop==N.max) {trialER.stop<-0}   
  ### early stop for the last cohort; trial.stop =1, while not early stoped.
  if (trialER.stop==0) {pts.stop <- N.max}
  
  # Return the results
  list(trial.stop = trial.stop, trialER.stop = trialER.stop, pts.stop = pts.stop,
       probs = probs, cutoff = cutoff, WR = WR)
}, error = function(e) {
  # Handle the error (e.g., print a message or log it)
  cat("Error:", conditionMessage(e), "\n")
   })
}

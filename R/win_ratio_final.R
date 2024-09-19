

#' est_cens_time
#' 
#' This function is an intermediate function used inside the winratio function to generate censoring time from estimated censoring distribution.
#' 
#' @param fit   Estimated censoring distribution.
#' @param n   A numeric value indicating number of censoring time to generate. Either n.current.trt or n.current.ctrl.
#' @returns Censoring time for each subject.
#' @examples
#' \dontrun{est_cens_time(trt_cens.est,n.current.trt)}


est_cens_time = function(fit,n){
  Cens.generated = as.numeric(stats::quantile(fit, probs = stats::runif(n) )$quantile)
  Cens.generated[is.na(Cens.generated)] = Inf # replace NAs with Inf
  return (Cens.generated)
}


#' winratio
#' 
#' This function is an intermediate function used inside the OCC.Table function to obtain the posterior probability of (WR > eta).


#' @param currentdd   The data at a interim analysis point.
#' @param n_current   Sample size at a interim analysis point.
#' @param Time_current Time at a interim analysis point.
#' @param Time_entry  Time at entry for each patient at a interim analysis point.
#' @param m0  Prior mean for mu.
#' @param L0  Prior covariance for mu.
#' @param v0  Prior degrees of freedom for Sigma.
#' @param S0  Prior scale matrix for Sigma.
#' @param time_max The upper limit for the recurrence and death time sampled from truncated normal. This will set the upper limit to to time_max rather than Inf.
#' @param eta A pre-specified lower bound of acceptable performance based on historical information.
#' @param lambda  Cutoff parameter.
#' @param N.max   Maximum number of patients to enroll.
#' @param thin_MCMC  Thinning degree.
#' @param Niter Number of iterations for gibbs sampler.
#' @returns A list with the following components:\tabular{ll}{
#'    \code{probs} \tab A numeric value, which is the posterior probability of \eqn{(\widehat{WR} > eta)}.  \cr
#'    \tab \cr
#'    \code{cutoff} \tab A numeric value. cutoff = (n_current/N.max)**lambda.  \cr
#'    \tab \cr
#'    \code{WR} \tab A numeric value. cutoff = the average posterior win ratio.  \cr
#' }
#' @examples
#' \dontrun{
#' winratio(currentdd = currentData,n_current = n.current, Time_current = Time.current,
#' Time_entry = Time.entry,m0 = c(0,0),L0 = diag(10^6, 2),v0 = 4,S0 = diag(10^(-6), 2),
#' time_max = 20,eta = 1,lambda = 0.25,N.max = 100, thin_MCMC = 5,
#' Niter = 100000)
#' } 

winratio = function(currentdd,n_current,Time_current,
                    Time_entry, m0,L0,v0,S0,
                    time_max,eta,lambda,N.max, thin_MCMC, Niter){
  
  # update actual censoring time at this interim look
  currentdd[,3] = pmin(Time_current-Time_entry[-n_current], currentdd[,3])
  currentdd[,6] <- as.integer(currentdd[,1] <= pmin(currentdd[,3],currentdd[,2]))
  currentdd[,7] <-  as.integer(currentdd[,2] <= currentdd[,3])
 
 
  # update recurrence and death  and the censoring time to reflect the observed data at this interim look    
  currentdd[,1] = ifelse(currentdd[,6]==1, currentdd[,1], NA)
  currentdd[,3]  = pmin(currentdd[,2] , currentdd[,3] )
  currentdd[,2] = ifelse(currentdd[,7]==1, currentdd[,2], NA)
  
  # number enrolled to each group
  n.current.trt = sum(currentdd[,4])
  n.current.ctrl =  n_current - 1 - n.current.trt 
  
  # separate data to trt and control
  currentData.trt = currentdd[currentdd[,4]==1,]
  currentData.ctrl = currentdd[currentdd[,4]==0,]
  
  # Step 1:Update posterior theta 
  trt.post = update_theta(N_iter = Niter, dd = currentData.trt, 
                          n = n.current.trt, m0 = m0, L0 = L0, S0 = S0, v0 = v0, 
                          time_max = time_max)

  ctrl.post = update_theta(N_iter = Niter, dd = currentData.ctrl, 
                           n = n.current.ctrl,m0 = m0, L0 = L0, S0 = S0, v0 = v0, 
                           time_max = time_max)
  burn_MCMC = as.integer(0.3*Niter)
  idxs <- seq(burn_MCMC, Niter, by = thin_MCMC)
  M_iter = length(idxs)
  
  # Step 2: estimating censoring distribution
  trt_cens.est = survival::survfit(Surv(currentData.trt[,"censor_t"], 1 - currentData.trt[,"delta2"]) ~ 1 )
  ctrl_cens.est =survival::survfit(Surv(currentData.ctrl[,"censor_t"], 1 - currentData.ctrl[,"delta2"]) ~ 1 )
  
  # Step 3: sample R and D for each group, and generate censoring time for each patient
  ### initialize posterior data
  postData <-array(dim=c(n_current-1,6, M_iter));  
  # Create dimension names
  dim_names <- list(
    NULL,
    Variables = c("R", "D", "C", "group",  "delta1", "delta2"),
    Iterations = paste("Iteration", 1:M_iter)
  )
  # Assign the dimension names to the array
  dimnames(postData) <- dim_names
  
  postData[1:n.current.trt,1:2,] = trt.post$post_data[,,idxs]
  postData[(n.current.trt+1):(n_current-1),1:2,] = ctrl.post$post_data[,,idxs]
  
  res1 = replicate(M_iter, est_cens_time(trt_cens.est,n.current.trt), simplify = F)
  res2 = replicate(M_iter, est_cens_time(ctrl_cens.est, n.current.ctrl), simplify = F)
  postData[1:n.current.trt,3,]  = matrix(unlist(res1), nrow = n.current.trt, ncol = M_iter, byrow = F)
  postData[(n.current.trt+1):(n_current-1),3,] =  matrix(unlist(res2), nrow = n.current.ctrl, ncol = M_iter, byrow = F)
  
  postData[1:n.current.trt,4,] = 1
  postData[(n.current.trt+1):(n_current-1),4,] = 0
  postData[, 5,] <- as.integer(postData[, 1,] <= pmin(postData[, 2,], postData[, 3,]))
  postData[, 6,] <- as.integer(postData[, 2,] <= postData[, 3,])
  postData[, 1,] <- ifelse( postData[, 5,]==1,  postData[, 1,], NA)
  postData[, 3,] <- pmin(postData[, 2,], postData[, 3,])
  postData[, 2,] <- ifelse( postData[, 6,]==1,  postData[, 2,], NA)
  

  # Step 4: calculate the IPCW-adjusted WR
  time.trt <- time.ctrl <- surv.trt <- surv.ctrl <- vector("list", M_iter)
  for (j.iter in 1:M_iter){
    C = postData[,3,j.iter] 
    delta2 = postData[,6,j.iter]
    R = postData[,1,j.iter] 
    D = postData[,2,j.iter] 
    
    fit.trt = survival::survfit(Surv(C[1:n.current.trt], 1-delta2[1:n.current.trt]) ~ 1)
    fit.ctrl = survival::survfit(Surv(C[(n.current.trt+1):(n_current-1)], 
                            1-delta2[(n.current.trt+1):(n_current-1)]) ~ 1)
    
    time.trt[[j.iter]] = summary(fit.trt)$time
    time.ctrl[[j.iter]] = summary(fit.ctrl)$time
    surv.trt[[j.iter]] = summary(fit.trt)$surv
    surv.ctrl[[j.iter]] = summary(fit.ctrl)$surv
  }
  
  ## getting the survival probability for event times, from the censoring distribution from both groups
  ## calculate WR
  result = compare(M_iter = M_iter, n_current_ctrl = n.current.ctrl, 
                                n_current_trt = n.current.trt, postData = postData,
                                time_trt = time.trt, time_ctrl = time.ctrl, surv_trt = surv.trt, surv_ctrl = surv.ctrl)
  WR = mean(result$WR, na.rm = T)
  probs = mean(result$WR > eta, na.rm  = T)  
  cutoff = (n_current/N.max)**lambda
  
  rm(result,postData, res1, res2)
  gc()
  
  return(list(probs = probs, cutoff = cutoff, WR = WR))
}



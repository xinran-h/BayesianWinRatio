% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/OCandCutoff_final.R
\name{OCC.Table}
\alias{OCC.Table}
\title{OCC.Table}
\usage{
OCC.Table(
  myData,
  N.max,
  design,
  cohort,
  recruit.int,
  m0,
  L0,
  v0,
  S0,
  time_max,
  eta,
  lambda,
  thin_MCMC,
  Niter
)
}
\arguments{
\item{myData}{A matrix with N.max rows and 7 columns. The first five columns represent: time to recurrence, time to death, time to censor, treatment arm, id. Leave the last two columns blank. Note that time to recurrence and time to death cannot be NA. If these are NA due to censoring, please use a large number to represent the time. For the treatment arm, use 1 for the treatment arm and 0 for the control arm.}

\item{N.max}{Maximum number of patients to enroll.}

\item{design}{A numeric value indicating the type of design. 1 = proposed design, 2 = time to recurrence design, 3= time to death design, 4 = time to first event design.}

\item{cohort}{Interim cohort, which is a numeric vector of the number of patients enrolled at each interim look.}

\item{recruit.int}{Recruitment interval.}

\item{m0}{Prior mean for mu}

\item{L0}{Prior variance for mu}

\item{v0}{For the proposed design, this is the prior degrees of freedom for Sigma. For the traditional designs, v0/2 is the prior shape for Sigma.}

\item{S0}{For the proposed design, this is the prior scale matrix for Sigma. For the traditional designs, v0*S0/2 is the prior scale for Sigma.}

\item{time_max}{The upper limit for the recurrence and death time sampled from truncated normal. This will set the upper limit to to time_max rather than Inf.}

\item{eta}{A pre-specified lower bound of acceptable performance based on historical information.}

\item{lambda}{Cutoff parameter.}

\item{thin_MCMC}{Thinning degree.}

\item{Niter}{Number of iterations for gibbs sampler.}
}
\value{
A list with the following components:\tabular{ll}{
\code{trial.stop} \tab A value of 1 or 0, 1 = trial stopped and 0 = not stopped.  \cr
\tab \cr
\code{trialER.stop} \tab A value of 1 or 0, 1 = trial stopped early and 0 = not stopped early. \cr
\tab \cr
\code{pts.stop} \tab A numeric value represents the actual sample size used. \cr
\tab \cr
\code{probs} \tab A numeric value, which is the posterior probability of \eqn{(\widehat{WR} > eta)}.  \cr
\tab \cr
\code{WR} \tab A numeric value, which is the estimated posterior win ratio. \cr
}
}
\description{
This is the main function of BayesWinRatio.
This function perform Bayesian futility monitoring based on one simulation data or one real data,  and returns the operating characteristics.
}
\examples{
\dontrun{
OCC.Table(myData = data,
N.max = 100,
design = 1, 
cohort = c(40,60,80),
recruit.int  = 0.25,
m0 = c(0,0),L0 = diag(10^6, 2), 
v0 = 4,S0 = diag(10^(-6), 2),
time_max = 10,eta = 1.5,lambda = 0.25,
 thin_MCMC = 5,Niter = 100000)
}
}

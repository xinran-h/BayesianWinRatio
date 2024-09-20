
<!-- README.md is generated from README.Rmd. Please edit that file -->

# BayesianWinRatio

<!-- badges: start -->
<!-- badges: end -->

This package performs Bayesian monitoring using the win ratio approach
based on the manuscript “The Win Ratio Approach in Bayesian Monitoring
for Two-Arm Phase II Clinical Trial Designs with Multiple Time-to-Event
Endpoints”.

## Installation and Load libraries

You can install the development version of BayesianWinRatio from
[GitHub](https://github.com/) with:

``` r
install.packages("devtools")
devtools::install_github("xinran-h/BayesianWinRatio")
```

## BayesianWinRatio in a nutshell

The function OCC.Table is the core function of BayesianWinRatio. This
function perform Bayesian futility monitoring based on one simulation
data or one real data, and returns the operating characteristics.

- **myData**: a matrix with N.max rows and 7 columns. The first five
  columns represent: time to recurrence, time to death, time to censor,
  treatment arm, id. Leave the last two columns blank. Note that time to
  recurrence and time to death cannot be NA. If these are NA due to
  censoring, please use a large number to represent the time. For the
  treatment arm, use 1 for the treatment arm and 0 for the control arm.
- **N.max**: maximum number of patients to enroll.
- **design**: a numeric value indicating the type of design. 1 =
  proposed design, 2 = time to recurrence design, 3= time to death
  design, 4 = time to first event design.
- **cohort**: interim cohort, which is a numeric vector of the number of
  patients enrolled at each interim look.
- **recruit.int**: recruitment interval.
- **m0**: prior mean for mu.
- **L0**: prior variance for mu.
- **v0**: for the proposed design, this is the prior degrees of freedom
  for Sigma. For the traditional designs, v0/2 is the prior shape for
  Sigma.
- **S0**:for the proposed design, this is the prior scale matrix for
  Sigma. For the traditional designs, v0\*S0/2 is the prior scale for
  Sigma.
- **time_max**: the upper limit for the recurrence and death time
  sampled from truncated normal. This will set the upper limit to to
  time_max rather than Inf.
- **eta**: a pre-specified lower bound of acceptable performance based
  on historical information.
- **lambda**: cutoff parameter.
- **thin_MCMC**: thinning degree.
- **Niter**: number of iterations for gibbs sampler.

The output is a list with the following components:

- **trial.stop**: a value of 1 or 0, 1 = trial stopped and 0 = not
  stopped.
- **trialER.stop**: a value of 1 or 0, 1 = trial stopped early and 0 =
  not stopped early.
- **pts.stop**: a numeric value represents the actual sample size used.
- **probs**: a numeric value, which is the posterior probability of
  $(\widehat{WR} > \eta)$.
- **WR**: a numeric value, which is the estimated posterior win ratio.

## An example

Suppose we conduct a trial with a maximum sample size of 20 patients,
who are randomized to two arms (control = 1, treatment = 2). Patients
are followed for two events, recurrence and death. We want to monitor
the trial to decide if we should stop the trial early when the the
treatment arm therapeutic effect is not satisfactory. The interim look
occurs when 5, 10, and 15 patients are enrolled, with a total of 3
interim looks. We use a composite endpoint, the win ratio, to monitor
the trial. The futility monitoring rule at the interim look time k, k =
1, 2, 3, is

$P(\widehat{WR} > 1.5 \vert Y_{k1}, Y_{k0}) > c(n_k)$,

where $c(n_k) = (\frac{n_k}{N.max})^{\lambda}$ is the cutoff parameter,
$\lambda$ is a tuning parameter tuned to achieve desired operating
characteristics, and $Y_{k1}$ and $Y_{k0}$ are observed data from both
arms. This monitoring rule specifies that, if there is a low probability
that the win ratio is greater than 1.5, we will stop the trial early.
The interim look time is inferred by the recruitment interval, which is
the time between enrolling two consecutive patients. For example, if we
enroll 4 patients over 1 month, the recruitment interval is 0.25.

### demo data

The demo data can be loaded by the following code.

``` r
library(BayesianWinRatio)
test
```

### Tuning lambda

In the design phase, we first determine the tuning parameter $\lambda$
by simulation. We tune $\lambda$ to achieve a percentage of early
termination (PET) at a desirable level (say, 0.1) under the null
scenario, where the two arms have identical distribution.Assuming that
the two event time follow a bivariate lognormal distribution, we
estimate the mean time to each event (log) and the variance-covariance
matrix between the two events using historical data. Suppose the
recruitment interval is 0.25, and the censoring time follows a uniform
distribution of $U(0,2)$ we generate 1000 simulation data, by calling
the function data.simulation, as follows:

``` r
library(MASS) 
library(survival) 
library(parallel) 
data = data.simulation(N.sim = 1000, N.max = 20,
                       mu.trt = c(0.2,0.3), Sigma.trt = matrix(c(1,0.5,0.5,1), nrow = 2, byrow = T),
                       mu.ctrl = c(0.2,0.3),Sigma.ctrl = matrix(c(1,0.5,0.5,1), nrow = 2, byrow = T),
                       cens_upper = 2)
```

The parameter $\lambda$ is calibrated by calling the function OCC.Table.
The OCC.Table uses data from one simulation, and returns the operating
characteristics of this simulation data. We run the simulation 10,000
times and calculate the percentage of early termination (PET). The
following code uses a bisection method to find the lambda:

``` r
Cutoff2Prob <- 0.1 # PET = 0.1 under the null scenario

# Define the bisection parameters
a <- 0  # The lowest value of mylambda, update based on the simulation results.
b <- 30 # The highest value of mylambda, update based on the simulation results.
tolerance <- 0.01   # Tolerance for stopping the bisection

N.sim = 10000
# Initialize the bisection loop
while (b - a > tolerance) {
  c <- (a + b) / 2

  
  ### run simulation in parallel, can be changed to other parallel packages such as doParallel, future, etc.
  sim_results <- parallel::mclapply(1:N.sim, function(x){OCC.Table(
    myData = data[,,x],
    N.max = 20,
    design = 1,
    cohort = c(5, 10, 15),
    recruit.int  = 0.25,
    m0 = c(0,0),
    L0 = diag(10^6, 2),
    v0 = 4,
    S0 = diag(10^(-6), 2),
    time_max = 100,
    eta = 1.5,
    lambda = c,
    thin_MCMC = 5,Niter = 100000)}, mc.cores = 100)
  
  ### summarize results from N.sim simulations
  stop.all<-stop.early <-pts.all <-rep(0, N.sim);
  stop.all <- sapply(sim_results, function(res) res$trial.stop)
  stop.early <- sapply(sim_results, function(res) res$trialER.stop)
  pts.all <- sapply(sim_results, function(res) res$pts.stop)
  MyRaw <- data.frame(PRN=1-stop.all, PEN=stop.early, EN=pts.all);
  mysim1<-(apply(MyRaw,2,mean)); 
  
  rm(MyRaw, sim_results, stop.all, stop.early, pts.all)
  gc()
  
  if (mysim1[2] == Cutoff2Prob) {
    break
  } 
  else if (mysim1[2] < Cutoff2Prob) {
    b <- c
  } else {
    a <- c
  }
}
# The final result
  c
```

### Trial Implementation

After determining the tuning parameter $\lambda$, we can implement the
trial. We first preprocess the real data to make it compatible with the
function OCC.Table. The following code shows how to do the data
manipulation:

``` r

# data preprocessing
dd = test
## recode arm: control - 0, trt - 1
dd$arm = as.numeric(dd$arm == 2)
## recode unobserved event time with a very large value, say 99999
dd$recurrence[is.na(dd$recurrence)] = 99999
dd$death[is.na(dd$death)] = 99999
myData = matrix(NA, nrow = nrow(dd), ncol = 7)
myData[,1:5] = as.matrix(dd[, c(3:5,2,1)])
colnames(myData) = c("recurrence_t", "death_t", "censor_t",   "group", "id", "delta1", "delta2")
```

After the data preprocessing, we can run the function OCC.Table to get
the operating characteristics of the trial. At the first interim look
when 5 patients are enrolled, we run the function OCC.Table as follows:

``` r
# run OCC.Table using myData
OCC.Table(
   myData = myData,
    N.max = 20,
    design = 1,
    cohort = c(5),
    recruit.int  = 0.25,
    m0 = c(0,0),
    L0 = diag(10^6, 2),
    v0 = 4,
    S0 = diag(10^(-6), 2),
    time_max = 100,
    eta = 1.5,
    lambda = 10, # using the lambda calibrated above
    thin_MCMC = 5,Niter = 100000)
```

If trialER.stop = 1, we stop the trial early. Otherwise, we continue to
enroll patients. At the second interim look when 10 patients are
enrolled, we run the function OCC.Table as follows:

``` r
# run OCC.Table using myData
OCC.Table(
   myData = myData,
    N.max = 20,
    design = 1,
    cohort = c(10),
    recruit.int  = 0.25,
    m0 = c(0,0),
    L0 = diag(10^6, 2),
    v0 = 4,
    S0 = diag(10^(-6), 2),
    time_max = 100,
    eta = 1.5,
    lambda = 10, # using the lambda calibrated above
    thin_MCMC = 5,Niter = 100000)
```

If trialER.stop = 1, we stop the trial early. Otherwise, we continue to
enroll patients. At the third interim look when 15 patients are
enrolled, we run the function OCC.Table:

``` r
# run OCC.Table using myData
OCC.Table(
   myData = myData,
    N.max = 20,
    design = 1,
    cohort = c(15),
    recruit.int  = 0.25,
    m0 = c(0,0),
    L0 = diag(10^6, 2),
    v0 = 4,
    S0 = diag(10^(-6), 2),
    time_max = 100,
    eta = 1.5,
    lambda = 10, # using the lambda calibrated above
    thin_MCMC = 5,Niter = 100000)
```

If trialER.stop = 1, we stop the trial early. Otherwise, we continue to
enroll patients until maximum sample size is reached.

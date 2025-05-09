% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/flgi_cut_off_unknown_var.r
\name{flgi_cut_off_unknown_var}
\alias{flgi_cut_off_unknown_var}
\title{Cut-off Value of the Forward-looking Gittins Index rule in Continuous Endpoint with Unknown Variances}
\usage{
flgi_cut_off_unknown_var(
  Gittinstype,
  df,
  gittins = NULL,
  Pats,
  nMax,
  TimeToOutcome,
  enrollrate,
  K,
  noRuns2,
  Tsize,
  block,
  rule,
  prior_n,
  prior_mean1,
  prior_sd1,
  mean,
  sd,
  side
)
}
\arguments{
\item{Gittinstype}{type of Gittins indices, should be set to 'UNKV' in this function}

\item{df}{discount factor which is the multiplier for loss at each additional patient in the future.
Available values are 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.99 and 0.995. The maximal sample size can be up to 10000.}

\item{gittins}{user specified Gittins indices for calculation in this function. If \code{gittins} is provided,
\code{Gittinstype} and \code{df} should be NULL.}

\item{Pats}{the number of patients accrued within a certain time frame indicates the
count of individuals who have been affected by the disease during that specific period,
for example, a month or a day. If this number is 10, it represents that
10 people have got the disease within the specified time frame.}

\item{nMax}{the assumed maximum accrued number of patients with the disease in the population, this number
should be chosen carefully to ensure a sufficient number of patients are simulated,
especially when considering the delay mechanism.}

\item{TimeToOutcome}{the distribution of delayed response times or a fixed delay time for responses.
The delayed time could be a month, a week or any other time frame. When the unit changes,
the number of TimeToOutcome should also change. It can be in the format
of expression(rnorm( length( vStartTime ),30, 3)), representing delayed responses
with a normal distribution, where the mean is 30 days and the standard deviation is 3 days.}

\item{enrollrate}{probability that patients in the population can enroll in the trial.
This parameter is related to the number of people who have been affected by the disease in the population,
following an exponential distribution.}

\item{K}{number of total arms in the trial.}

\item{noRuns2}{number of simulations for simulated allocation probabilities within each block. Default value is
set to 100 times, which is recommended in Villar et al., 2015.}

\item{Tsize}{maximal sample size for the trial.}

\item{block}{block size.}

\item{rule}{rules can be used in this function, with values 'FLGI PM', 'FLGI PD' or 'CFLGI'.
'FLGI PM' stands for making decision based on posterior mean;
'FLGI PD' stands for making decision based on posterior distribution;
'CFLGI' stands for controlled forward-looking Gittins Index.}

\item{prior_n}{a vector representing the number of observations assumed in prior distribution, eg: c(1,1) for a two-armed trial.}

\item{prior_mean1}{a vector representing mean of observations assumed in prior distributions, eg: c(0,0,0) for a three-armed trial,
rep(0,K) can be used to simplify the process. If a negative effect is expected, adjust the mean to a negative value.}

\item{prior_sd1}{a vector representing the standard deviation of observations assumed in prior distribution, eg: rep(1,3) for a three-armed trial.}

\item{mean}{a vector of mean hypotheses, for example, as c(0.1,0.1) where 0.1 stands for the mean
for both groups. Another example is c(0.1,0.3) where 0.1 and 0.3 stand for the mean for the control and
a treatment group, respectively.}

\item{sd}{a vector of standard deviation hypotheses, for example, as c(0.64,0.64) where 0.64 stands for the standard deviation
for both groups. Another example is c(0.64,0.4) where 0.64 and 0.4 stand for the standard deviation for the control and
a treatment group, respectively.}

\item{side}{direction of one-sided test with the values of 'upper' or 'lower'.}
}
\value{
Value of T test statistics for one trial.
}
\description{
Function for simulating cut-off values at the final stage using the forward-looking Gittins Index rule
and the controlled forward-looking Gittins Index rule for continuous outcomes with known variance in trials with
2-5 arms. The prior distributions follow Normal-Inverse-Gamma (NIG) (\eqn{(\mu,\sigma^2) \sim NIG({\sf mean}=m,{\sf variance}=V \times \sigma^2,{\sf shape}=a,{\sf rate}=b)})
distributions and should be the same for each arm.
}
\details{
This function simulates trials using the forward-looking Gittins Index rule and the
controlled forward-looking Gittins Index rule under both no delay and delayed scenarios to obtain
cut-off values at the final stage, with control of type I error. The user is expected to run this function
multiple times to determine a reasonable cut-off value for statistical inference.
Considering the delay mechanism, \code{Pats} (the number of patients accrued within a certain time frame),
\code{nMax} (the assumed maximum accrued number of patients with the disease in the population) and 
\code{TimeToOutcome} (the distribution of delayed response times or a fixed delay time for responses) 
are parameters in the functions adapted from \url{https://github.com/kwathen/IntroBayesianSimulation}.
Refer to the website for more details.
}
\examples{
#The forward-looking Gittins Index rule with delayed responses follow a normal 
#distribution with a mean of 60 days and a standard deviation of 3 days
#One can run the following command 20000 times to obtain the selected cut-off 
#value around -1.9298 with an overall lower one-sided type I error 0.025
\donttest{
stopbound1<-lapply(1:20000,function(x){ 
flgi_cut_off_unknown_var(Gittinstype='UNKV',df=0.5,Pats=10,nMax=50000,
TimeToOutcome=expression(rnorm( length( vStartTime ),60, 3)),enrollrate=0.9,
K=3,noRuns2=100,Tsize=852,block=20,rule='FLGI PM',prior_n=rep(2,3),
prior_mean1=rep(9/100,3),prior_sd1=rep(0.006324555,3),
mean=c(9.1/100,9.1/100,9.1/100),sd=c(0.009,0.009,0.009),side='lower')})
stopbound1a<-do.call(rbind,stopbound1)
sum(stopbound1a<(-1.9298) )/20000
#The selected cut-off value is around -1.9298 with an overall lower one-sided 
#type I error of 0.025, based on 20000 simulations.
}
}
\references{
\insertRef{Williamson2019}{RARtrials}
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/brar_select_au_known_var.r
\name{brar_select_au_known_var}
\alias{brar_select_au_known_var}
\title{Select au in Bayesian Response-Adaptive Randomization with a Control Group for Continuous Endpoint with Known Variances}
\usage{
brar_select_au_known_var(
  Pats,
  nMax,
  TimeToOutcome,
  enrollrate,
  N1,
  armn,
  N2,
  tp,
  armlabel,
  blocksize,
  mean,
  sd,
  minstart,
  deltaa,
  tpp,
  deltaa1,
  mean10 = 0,
  mean20 = mean10,
  mean30 = mean10,
  mean40 = mean10,
  mean50 = mean10,
  sd10 = 1,
  sd20 = sd10,
  sd30 = sd10,
  sd40 = sd10,
  sd50 = sd10,
  n10 = 1,
  n20 = n10,
  n30 = n10,
  n40 = n10,
  n50 = n10,
  side,
  output = NULL,
  ...
)
}
\arguments{
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

\item{N1}{number of participants with equal randomization in the 'initialization' period.
Recommend using 10 percent of the total sample size.}

\item{armn}{number of total arms in the trial.}

\item{N2}{maximal sample size for the trial.}

\item{tp}{tuning parameter. Some commonly used numbers are 0.5, 1 and n/2N.}

\item{armlabel}{a vector of treatment labels with an example of c(1, 2), where 1 and 2 describe
how each arm is labeled in a two-armed trial.}

\item{blocksize}{size of block used for equal randomization regarding participants in the 'initialization' period.
Recommend to be an even multiple of the number of total arms.}

\item{mean}{a vector of means in hypotheses, for example, as c(10,10) where 10 stands for the mean
in both groups. Another example is c(10,12) where 10 and 12 stand for the mean
for the control and the treatment group, respectively.}

\item{sd}{a vector of standard deviations in hypotheses, for example, as c(2,2) where 2 stands for the standard deviation
in both groups. Another example is c(1,2) where 1 and 2 stand for the standard deviation
for the control and the treatment group, respectively.}

\item{minstart}{a specified number of participants when one starts to check decision rules.}

\item{deltaa}{a vector of minimal effect expected to be observed for early futility stopping in
each arm is approximately \eqn{1\%}. The length of this parameter is \code{armn}-1.}

\item{tpp}{indicator of \code{tp} equals to n/2N. When \code{tp} is n/2N, \code{tpp} should be assigned 1. Default value is set to 0.}

\item{deltaa1}{a vector of pre-specified minimal effect size expected to be observed at the final stage
for each arm. The length of this parameter is \code{armn}-1.}

\item{mean10, sd10}{prior mean and sd in \eqn{N(mean,sd)} of arm 1 in the trial, which stands for the control. Default value is set to 1.}

\item{mean20, sd20}{prior mean and sd in \eqn{N(mean,sd)} of arm 2 in the trial. Default value is set to \code{mean10} and \code{sd10}.}

\item{mean30, sd30}{prior mean and sd in \eqn{N(mean,sd)} of arm 3 in the trial. Default value is set to \code{mean10} and \code{sd10}.}

\item{mean40, sd40}{prior mean and sd in \eqn{N(mean,sd)} of arm 4 in the trial. Default value is set to \code{mean10} and \code{sd10}.}

\item{mean50, sd50}{prior mean and sd in \eqn{N(mean,sd)} of arm 5 in the trial. Default value is set to \code{mean10} and \code{sd10}.}

\item{n10}{explicit prior n of arm 1 in the trial, which stands for the control. Default value is set to 1.}

\item{n20}{explicit prior n of arm 2 in the trial. Default value is set to \code{n10}.}

\item{n30}{explicit prior n of arm 3 in the trial. Default value is set to \code{n10}.}

\item{n40}{explicit prior n of arm 4 in the trial. Default value is set to \code{n10}.}

\item{n50}{explicit prior n of arm 5 in the trial. Default value is set to \code{n10}.}

\item{side}{direction of a one-sided test, with values 'upper' or 'lower'.}

\item{output}{control the output of brar_select_au_binary. If the user does not specify anything, the function returns
the entire dataset used to select the stopping boundary for each iteration. If the user specifies 'B', the function
only returns the selected stopping boundary for each iteration.}

\item{...}{additional arguments to be passed to \code{\link[stats]{integrate}} (such as rel.tol) from this function.}
}
\value{
A list of results generated from formula \eqn{Pr(\mu_k>\mu_{control}+\delta|data_{t-1})} at each step.
Note that before final stage of the trial, test statistics is calculated from \code{deltaa}, and test statistics is
calculated from \code{deltaa1} at the final stage.
}
\description{
\code{brar_select_au_known_var} involves selecting au in Bayesian Response-Adaptive Randomization with a control group
for continuous endpoints with known variance in trials with two to five arms. The conjugate prior distributions follow
Normal (\eqn{N(mean,sd)}) distributions and can be specified individually for each arm.
}
\details{
This function generates a data set or a value in one iteration for selecting the appropriate au using Bayesian
response-adaptive randomization with a control group under null hypotheses with no delay and delayed scenarios.
The function can handle trials with up to 5 arms for continuous outcomes with known variances. This function uses the formula
\eqn{\frac{Pr(\mu_k={\sf max}\{\mu_1,...,\mu_K\})^{tp}} {\sum_{k=1}^{K}{Pr(\mu_k={\sf max}\{\mu_1,...,\mu_K\})^{tp}}}} with \code{side} equals to 'upper',
and \eqn{\frac{Pr(\mu_k={\sf min}\{\mu_1,...,\mu_K\})^{tp}} {\sum_{k=1}^{K}{Pr(\mu_k={\sf min}\{\mu_1,...,\mu_K\}){tp}}}} 
with \code{side} equals to 'lower', utilizing available data at each step.
Considering the delay mechanism, \code{Pats} (the number of patients accrued within a certain time frame),
\code{nMax} (the assumed maximum accrued number of patients with the disease in the population) and 
\code{TimeToOutcome} (the distribution of delayed response times or a fixed delay time for responses) 
are parameters in the functions adapted from \url{https://github.com/kwathen/IntroBayesianSimulation}.
Refer to the website for more details.
}
\examples{
#brar_select_au_known_var with delayed responses follow a normal distribution with
#a mean of 30 days and a standard deviation of 3 days, where mean=c(8.9/100,8.74/100,8.74/100),
#sd=c(0.009,0.009,0.009), tp=0.5 and the minimal effect size is 0.
set.seed(789)
stopbound1<-lapply(1:10,function(x){
brar_select_au_known_var(Pats=10,nMax=50000,TimeToOutcome=expression(
rnorm(length( vStartTime ),30, 3)),enrollrate=0.1, N1=21,armn=3,
N2=189,tp=0.5,armlabel=c(1,2,3),blocksize=6,mean=c((8.9/100+8.74/100+8.74/100)/3,
(8.9/100+8.74/100+8.74/100)/3,(8.9/100+8.74/100+8.74/100)/3),
sd=c(0.009,0.009,0.009),minstart=21,deltaa=c(0,0.001),tpp=0,deltaa1=c(0,0),
mean10=0.09,mean20=0.09,mean30=0.09, sd10=0.01,sd20=0.01,sd30=0.01,n10=1,n20=1,
n30=1,side='lower')})
simf<-list()
simf1<-list()
for (xx in 1:10){
 if (any(stopbound1[[xx]][21:188,2]<0.01)){
      simf[[xx]]<-NA
   }  else{
      simf[[xx]]<-stopbound1[[xx]][189,2]
 }
 if (any(stopbound1[[xx]][21:188,3]<0.01)){
      simf1[[xx]]<-NA
   }  else{
      simf1[[xx]]<-stopbound1[[xx]][189,3]
 }
}
simf2<-do.call(rbind,simf)
sum(is.na(simf2)) #1, achieve around 10\% futility
simf3<-do.call(rbind,simf1)
sum(is.na(simf3)) #1, achieve around 10\% futility
stopbound1a<-cbind(simf2,simf3)
stopbound1a[is.na(stopbound1a)] <- 0
sum(stopbound1a[,1]>0.973 | stopbound1a[,2]>0.973)/10 #0.1
#the selected stopping boundary is 0.973 with an overall lower one-sided type I
#error of 0.1, based on 10 simulations. Because it is under the permutation null hypothesis,
#the selected deltaa should be an average of 0 and 0.001 which is 0.0005, although
#deltaa could be close to each other with larger simulation numbers.
#It is recommended to conduct more simulations (i.e., 1000) to obtain an accurate deltaa and au.
#As the simulation number increases, the choice of deltaa could be consistent for comparisons
#of each arm to the control.

}
\references{
\insertRef{Wathen2017}{RARtrials}
}

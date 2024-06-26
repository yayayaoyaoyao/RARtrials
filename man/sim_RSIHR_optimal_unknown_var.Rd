% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sim_RSIHR_optimal_unknown_var.r
\name{sim_RSIHR_optimal_unknown_var}
\alias{sim_RSIHR_optimal_unknown_var}
\title{Simulate a Trial Using Generalized RSIHR Allocation for Continuous Endpoint with Unknown Variances}
\usage{
sim_RSIHR_optimal_unknown_var(
  Pats,
  nMax,
  TimeToOutcome,
  enrollrate,
  N1,
  N2,
  armn,
  mean,
  sd,
  alphaa = 0.025,
  armlabel,
  cc,
  side
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

\item{N2}{maximal sample size for the trial.}

\item{armn}{number of total arms in the trial.}

\item{mean}{a vector of hypotheses of mean, with the first one serving as the control group.}

\item{sd}{a vector of hypotheses of standard deviation with the first one serving as the control group.}

\item{alphaa}{the overall type I error to be controlled for the one-sided test. Default value is set to 0.025.}

\item{armlabel}{a vector of arm labels with an example of c(1, 2), where 1 and 2 describe
how each arm is labeled in a two-armed trial.}

\item{cc}{value in the formula of measure of treatment effectiveness, usually take the average
of mean responses in the hypotheses. \code{cc} is the same as C in the details.}

\item{side}{direction of a one-sided test, with values 'upper' or 'lower'.}
}
\value{
\code{sim_RSIHR_optimal_unknown_var} returns an object of class "RSIHRoptimal". An object of class "RSIHRoptimal" is a list containing 
final decision based on the T test statistics with 1 stands for selected and 0 stands for not selected,
T test statistics, the simulated data set and participants accrued for each arm at the time of termination of that group in one trial.
The simulated data set includes 5 columns: participant ID number, enrollment time, observed time of results,
allocated arm, and participants' result.
}
\description{
\code{sim_RSIHR_optimal_unknown_var} simulates a trial for continuous endpoints with unknown variances,
and the allocation probabilities change based on results of accumulated participants in the trial.
}
\details{
This function aims to minimize the criteria \eqn{\sum_{i=1}^{K}n_i\Psi_i}
with constraints \eqn{\frac{\sigma_1^2}{n_1}+\frac{\sigma_k^2}{n_k}\leq C}, where \eqn{k=2,...,K}
for some fixed C. It is equivalent to generalized RSIHR allocation for continuous endpoints with unknown variances.
With more than two arms the one-sided nominal level of each test is \code{alphaa} divided 
by \code{arm*(arm-1)/2}; a Bonferroni correction.
Considering the delay mechanism, \code{Pats} (the number of patients accrued within a certain time frame),
\code{nMax} (the assumed maximum accrued number of patients with the disease in the population) and 
\code{TimeToOutcome} (the distribution of delayed response times or a fixed delay time for responses) 
are parameters in the functions adapted from \url{https://github.com/kwathen/IntroBayesianSimulation}.
Refer to the website for more details.
}
\examples{
#Run the function with delayed responses follow a normal distribution with
#a mean of 30 days and a standard deviation of 3 days under null hypothesis
#in a three-armed trial
sim_RSIHR_optimal_unknown_var(Pats=10,nMax=50000,TimeToOutcome=expression(
rnorm( length( vStartTime ),30, 3)),enrollrate=0.9,N1=8,N2=88,armn=2,
mean=c(9.1/100,9.1/100),sd=c(0.009,0.009),alphaa=0.025,armlabel=c(1,2),
cc=mean(c(9.1/100,9.1/100)),side='upper')

#Run the function with delayed responses follow a normal distribution with
#a mean of 30 days and a standard deviation of 3 days under alternative hypothesis
#in a three-armed trial
sim_RSIHR_optimal_unknown_var(Pats=10,nMax=50000,TimeToOutcome=expression(
rnorm( length( vStartTime ),30, 3)),enrollrate=0.9,N1=8,N2=88,armn=2,
mean=c(9.1/100,8.47/100),sd=c(0.009,0.009),alphaa=0.025,armlabel=c(1,2),
cc=mean(c(9.1/100,8.47/100)),side='upper')
}
\references{
\insertRef{Biswas2011}{RARtrials}
}

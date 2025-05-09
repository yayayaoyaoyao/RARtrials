% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sim_dabcd_max_power.r
\name{sim_dabcd_max_power}
\alias{sim_dabcd_max_power}
\title{Simulate a Trial Using Doubly Adaptive Biased Coin Design with Maximal Power Strategy for Binary Endpoint}
\usage{
sim_dabcd_max_power(
  Pats,
  nMax,
  TimeToOutcome,
  enrollrate,
  N1,
  N2,
  armn,
  armlabel,
  h,
  BB,
  type,
  gamma = 2,
  alphaa = 0.025,
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

\item{N1}{number of participants with equal randomization in the burn-in period.
Recommend using 10 percent of the total sample size.}

\item{N2}{maximal sample size for the trial.}

\item{armn}{number of total arms in the trial.}

\item{armlabel}{a vector of arm labels with an example of c(1, 2), where 1 and 2 describes
how each arm is labeled in a two-armed trial.}

\item{h}{a vector of success probabilities in hypotheses, for example, as c(0.1,0.1) where 0.1 stands for the success probability
for both groups. Another example is c(0.1,0.3) where 0.1 and 0.3 stand for the success probabilities
for the control and the treatment group, respectively.}

\item{BB}{the minimal allocation probabilities for each arm, which is within the
range of \eqn{[0,1/armn]}.}

\item{type}{allocation type, with choices from 'RSIHR' and 'Neyman'.}

\item{gamma}{tuning parameter in Hu & Zhang's formula. When dabcd=0, this parameter does not need
to be specified. Default value is set to 2.}

\item{alphaa}{the overall type I error to be controlled for the one-sided test. Default value is set to 0.025.}

\item{side}{direction of a one-sided test, with values 'upper' or 'lower'.}
}
\value{
\code{sim_dabcd_max_power} returns an object of class "dabcd". An object of class "dabcd" is a list containing 
final decision based on the Z test statistics with 1 stands for selected and 0 stands for not selected,
Z test statistics, the simulated data set and participants accrued for each arm at the time of termination of that group in one trial.
The simulated data set includes 5 columns: participant ID number, enrollment time, observed time of results,
allocated arm, and participants' result.
}
\description{
\code{sim_dabcd_max_power} can be used for doubly adaptive biased coin design with maximal power
strategy for binary outcomes, targeting generalized Neyman allocation and generalized RSIHR allocation.
}
\details{
The function simulates a trial for doubly adaptive biased coin design with maximal power strategy targeting
generalized Neyman allocation with 2-5 arms which is provided in \insertCite{Tymofyeyev2007}{RARtrials} and
generalized RSIHR allocation with 2-3 arms which is provided in \insertCite{Jeon2010}{RARtrials}, with modifications for typos
in \insertCite{Sabo2016}{RARtrials}. All of those methods are not smoothed. The output of this function is based on Hu \code{\&} Zhang's formula \insertCite{Hu2004}{RARtrials}.
With more than two armd the one-sided nominal level of each test is \code{alphaa} divided by \code{arm*(arm-1)/2}; a Bonferroni correction.
Considering the delay mechanism, \code{Pats} (the number of patients accrued within a certain time frame),
\code{nMax} (the assumed maximum accrued number of patients with the disease in the population) and 
\code{TimeToOutcome} (the distribution of delayed response times or a fixed delay time for responses) 
are parameters in the functions adapted from \url{https://github.com/kwathen/IntroBayesianSimulation}.
Refer to the website for more details.
}
\examples{
sim_dabcd_max_power(Pats=10,nMax=50000,TimeToOutcome=expression(rnorm(
length( vStartTime ),30, 3)),enrollrate=0.9,N1=30,N2=300,armn=3,
armlabel=c(1,2,3),h=c(0.2,0.3,0.2),BB=0.1,type='Neyman',
side='upper')
sim_dabcd_max_power(Pats=10,nMax=50000,TimeToOutcome=expression(rnorm(
length( vStartTime ),60, 3)),enrollrate=0.1,N1=50,N2=500,armn=3,
armlabel=c(1,2,3),h=c(0.2,0.3,0.3),BB=0.15,type='RSIHR',
side='upper')

}
\references{
\insertRef{Hu2004}{RARtrials}

\insertRef{Tymofyeyev2007}{RARtrials}

\insertRef{Jeon2010}{RARtrials}

\insertRef{Sabo2016}{RARtrials}
}

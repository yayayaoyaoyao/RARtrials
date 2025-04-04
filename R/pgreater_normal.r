#' @title Calculate the Futility Stopping Probability for Continuous Endpoint with Known Variances Using Normal Distribution
#' @description Calculate the futility stopping probability in Bayesian response-adaptive randomization with
#' a control group using the Thall \eqn{\&} Wathen method for continuous outcomes with known variances. The conjugate prior distributions
#' follow Normal (\eqn{N(mean,sd)}) distributions and can be specified individually for each treatment group.
#' @details This function calculates the results of \eqn{Pr(\mu_k>\mu_{{\sf control}}+\delta|{\sf data})} for \code{side} equals to
#' 'upper' and the results of \eqn{Pr(\mu_{{\sf control}}>\mu_k+\delta|{\sf data})} for \code{side} equals to 'lower'.
#' The result indicates the posterior probability of stopping a treatment group due to futility around \eqn{1\%} in Bayesian
#' response-adaptive randomization with a control arm using Thall \eqn{\&} Wathen method, with accumulated results
#' during the conduct of trials. 
#' @aliases pgreater_normal
#' @export pgreater_normal
#' @param mean1,sd1  mean and sd in \eqn{N({\sf mean},{\sf sd})}, current estimated mean and sd for the control group.
#' @param mean2,sd2  mean and sd in \eqn{N({\sf mean},{\sf sd})}, current estimated mean and sd for the treatment group which is compared to the control group.
#' @param delta pre-specified minimal effect size expected to be observed between the control group and the compared treatment group.
#' @param side direction of a one-sided test, with values 'upper' or 'lower'.
#' @param ... additional arguments to be passed to stats::integrate() (such as rel.tol) from this function.
#' @return a posterior probability of \eqn{Pr(\mu_k>\mu_{{\sf control}}+\delta|{\sf data})} with \code{side} equals to 'upper';
#' a posterior probability of \eqn{Pr(\mu_{{\sf control}}>\mu_k+\delta|{\sf data})} with \code{side} equals to 'lower'.
#' @importFrom stats pt
#' @importFrom stats dt
#' @examples
#' pgreater_normal(mean1=0.091,sd1=0.09,mean2=0.097,sd2=0.08,delta=0,side='upper')
#' pgreater_normal(mean1=0.091,sd1=0.09,mean2=0.087,sd2=0.1,delta=0,side='lower')
#' @references 
#' \insertRef{Wathen2017}{RARtrials}
#' \insertRef{Kevin2007}{RARtrials}

pgreater_normal<-function(mean1=NULL,sd1=NULL,mean2=NULL,sd2=NULL,delta=0,side,...){
  if (side=='lower'){
    result<-pnorm(delta, mean1-mean2,
                       sqrt(sd1^2+sd2^2),lower.tail=FALSE)
    
  }else if (side=='upper'){
    result<-pnorm(delta, mean2-mean1,
                  sqrt(sd1^2+sd2^2),lower.tail=FALSE)
    
  }
  
  return(result)
}


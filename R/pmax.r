#' @title pmax
#' @description Calculate probability that a particular arm is the best in a trial using Bayesian response-adaptive randomization with
#' a control group using Thall \& Wathen method. The prior distributions follow Beta (\eqn{beta(\alpha,\beta)}) distributions
#' for binary outcomes, Normal (\eqn{N(mean,sd)}) distributions for continuous outcomes with known variance, and 
#' Normal-Inverse-Gamma (NIG) (\eqn{NIG(V,m,a,b)}) distributions which are equivalent to Normal-Inverse-Chi-Squared distributions 
#' for continuous outcomes with unknown variances for each arm and can be specified individually.
#' @details This function calculates the results of formula \eqn{Pr(p_k=max\{p_1,...,p_K\})} for
#' \code{side} equals to 'upper' and the results of formula \eqn{Pr(p_k=min\{p_1,...,p_K\})} for
#' \code{side} equals to 'lower'. This function returns the probability that the posterior probability of arm
#' \eqn{k} is maximal or minimal in trials with up to five arms.
#' @aliases pmax
#' @author Chuyao Xu, Thomas Lumley, Alain Vandal
#' @export pmax
#' @param outcome type of outcomes, with values of choices from 'binary', 'UNKV' and 'KV', representing binary outcomes,
#' continuous outcomes with known variance and continuous outcomes with unknown variance respectively.
#' When 'binary' is specified, \code{mean1} to \code{mean5}, \code{sd1} to \code{sd5}, \code{par1} to \code{par5} should be NULL.
#' When 'UNKV' is specified, \code{a1} to \code{a5}, \code{b1} to \code{b5}, \code{mean1} to \code{mean5}, \code{sd1} to \code{sd5} should be NULL.
#' When 'KV' is specified, \code{a1} to \code{a5}, \code{b1} to \code{b5}, \code{par1} to \code{par5} should be NULL.
#' @param armn number of arms in the trial with values up to 5. When \code{armn}=2,
#' only \code{a1} to \code{a2} and \code{b1} to \code{b2} OR \code{mean1} to \code{mean2} and \code{sd1} to \code{sd2} 
#' OR \code{par1} to \code{par2} need to be specified.
#' When \code{armn}=3, only \code{a1} to \code{a3} and  \code{b1} to \code{b3} OR \code{mean1} to \code{mean3} and \code{sd1} to \code{sd3} 
#' OR \code{par1} to \code{par3} need to be specified.
#' When \code{armn}=4, only \code{a1} to \code{a4} and  \code{b1} to \code{b4} OR \code{mean1} to \code{mean4} and \code{sd1} to \code{sd4} 
#' OR \code{par1} to \code{par4} need to be specified.
#' When \code{armn}=5, \code{a1} to \code{a5} and  \code{b1} to \code{b5} OR \code{mean1} to \code{mean5} and \code{sd1} to \code{sd5} 
#' OR \code{par1} to \code{par5} need to be specified.
#' @param a1 \eqn{\alpha} in prior \eqn{beta(\alpha,\beta)} for the arm to calculate the allocation probability of.
#' @param b1 \eqn{\beta} in prior \eqn{beta(\alpha,\beta)} for the arm to calculate the allocation probability of.
#' @param a2 \eqn{\alpha} in prior \eqn{beta(\alpha,\beta)} for one of the remaining arms.
#' @param a3 \eqn{\alpha} in prior \eqn{beta(\alpha,\beta)} for one of the remaining arms.
#' @param a4 \eqn{\alpha} in prior \eqn{beta(\alpha,\beta)} for one of the remaining arms.
#' @param a5 \eqn{\alpha} in prior \eqn{beta(\alpha,\beta)} for one of the remaining arms.
#' @param b2 \eqn{\beta} in prior \eqn{beta(\alpha,\beta)} for one of the remaining arms.
#' @param b3 \eqn{\beta} in prior \eqn{beta(\alpha,\beta)} for one of the remaining arms.
#' @param b4 \eqn{\beta} in prior \eqn{beta(\alpha,\beta)} for one of the remaining arms.
#' @param b5 \eqn{\beta} in prior \eqn{beta(\alpha,\beta)} for one of the remaining arms.
#' @param mean1 mean in prior Normal(mean,sd) for the arm to calculate the allocation probability of.
#' @param sd1 sd in prior Normal(mean,sd) for the arm to calculate the allocation probability of.
#' @param mean2 mean in prior Normal(mean,sd) for one of the remaining arms.
#' @param sd2 sd in prior Normal(mean,sd) for one of the remaining arms.
#' @param mean3 mean in prior Normal(mean,sd) for one of the remaining arms.
#' @param sd3 sd in prior Normal(mean,sd) for one of the remaining arms.
#' @param mean4 mean in prior Normal(mean,sd) for one of the remaining arms.
#' @param sd4 sd in prior Normal(mean,sd) for one of the remaining arms.
#' @param mean5 mean in prior Normal(mean,sd) for one of the remaining arms.
#' @param sd5 sd in prior Normal(mean,sd) for one of the remaining arms.
#' @param par1 a vector of parameters for the arm with a Normal-Inverse-Chi-Squared prior
#' to calculate the allocation probability of.
#' @param par2 a vector of parameters for one of the remaining arms with a Normal-Inverse-Chi-Squared prior.
#' @param par3 a vector of parameters for one of the remaining arms with a Normal-Inverse-Chi-Squared prior.
#' @param par4 a vector of parameters for one of the remaining arms with a Normal-Inverse-Chi-Squared prior.
#' @param par5 a vector of parameters for one of the remaining arms with a Normal-Inverse-Chi-Squared prior.
#' @param side direction of a one-sided test, with values 'upper' or 'lower'.
#' @param ... additional arguments to be passed to \code{\link[stats]{integrate}} (such as rel.tol) from this function.
#' @return a probability that a particular arm is the best in trials up to five arms.
#' @importFrom stats integrate
#' @importFrom stats dbeta
#' @importFrom stats pbeta
#' @importFrom stats dnorm
#' @examples
#' pmax(outcome='binary',armn=5,a1=8,b1=10,a2=5,b2=19,a3=8,b3=21,
#' a4=6, b4=35, a5=15, b5=4, side='upper')
#' pmax(outcome='KV',armn=4,mean1=8, sd1=2,mean2=5,sd2=1,mean3=8,
#' sd3=1.5,mean4=6,sd4=2,side='upper')
#' para<-list(V=1/2,a=0.5,m=9.1/100,b=0.00002)
#' par<-convert_gamma_to_chisq(para)
#' set.seed(123451)
#' y1<-rnorm(100,0.091,0.009)
#' par11<-update_par_nichisq(y1, par)
#' set.seed(123452)
#' y2<-rnorm(90,0.09,0.009)
#' par22<-update_par_nichisq(y2, par)
#' set.seed(123453)
#' y3<-rnorm(110,0.0892,0.009)
#' par33<-update_par_nichisq(y3, par)
#' pmax(outcome='UNKV',armn=3,par1=par11,par2=par22,par3=par33,side='upper')
#' pmax(outcome='UNKV',armn=3,par1=par11,par2=par22,par3=par33,side='lower')


pmax<-function(outcome,armn,a1=NULL,b1=NULL,a2=NULL,b2=NULL,a3=NULL,b3=NULL,a4=NULL,b4=NULL,a5=NULL,b5=NULL,
                mean1=NULL,sd1=NULL,mean2=NULL,sd2=NULL,mean3=NULL,sd3=NULL,mean4=NULL,sd4=NULL,mean5=NULL,sd5=NULL,
                par1=NULL,par2=NULL,par3=NULL,par4=NULL,par5=NULL,side,...){
if (outcome=='binary'){
  if (side=='upper'){

    if (armn==2){
      f1b<-function(y){
          pbeta(y,a2,b2,lower.tail=TRUE)*dbeta(y,a1,b1)
      }
      integral<-integrate(f1b,0,1, rel.tol = 1e-6, stop.on.error = FALSE,...)$value
    } else if (armn==3){
      f2b<-function(y){
          pbeta(y,a2,b2,lower.tail=TRUE)*pbeta(y,a3,b3,lower.tail=TRUE)*
          dbeta(y,a1,b1)
      }
      integral<-integrate(f2b,0,1, rel.tol = 1e-6, stop.on.error = FALSE,...)$value
    }else if (armn==4){
      f3b<-function(y){
          pbeta(y,a2,b2,lower.tail=TRUE)*pbeta(y,a3,b3,lower.tail=TRUE)*
          pbeta(y,a4,b4,lower.tail=TRUE)*
          dbeta(y,a1,b1)
      }
      integral<-integrate(f3b,0,1, rel.tol = 1e-6, stop.on.error = FALSE,...)$value
    }else if (armn==5){
      f4b<-function(y){
          pbeta(y,a2,b2,lower.tail=TRUE)*pbeta(y,a3,b3,lower.tail=TRUE)*
          pbeta(y,a4,b4,lower.tail=TRUE)*pbeta(y,a5,b5,lower.tail=TRUE)*
          dbeta(y,a1,b1)
      }
      integral<-integrate(f4b,0,1, rel.tol = 1e-6, stop.on.error = FALSE,...)$value
    }


  }else if (side=='lower'){

    if (armn==2){
      f5b<-function(y){
        pbeta(y,a2,b2,lower.tail=FALSE)*dbeta(y,a1,b1)
      }
      integral<-integrate(f5b,0,1, rel.tol = 1e-6, stop.on.error = FALSE,...)$value
    }else if (armn==3){
      f6b<-function(y){
        pbeta(y,a2,b2,lower.tail=FALSE)*pbeta(y,a3,b3,lower.tail=FALSE)*
          dbeta(y,a1,b1)
      }
      integral<-integrate(f6b,0,1, rel.tol = 1e-6, stop.on.error = FALSE,...)$value
    }else if (armn==4){
      f7b<-function(y){
        pbeta(y,a2,b2,lower.tail=FALSE)*pbeta(y,a3,b3,lower.tail=FALSE)*
          pbeta(y,a4,b4,lower.tail=FALSE)*
          dbeta(y,a1,b1)
      }
      integral<-integrate(f7b,0,1, rel.tol = 1e-6, stop.on.error = FALSE,...)$value
    }else if (armn==5){
       f8b<-function(y){
          pbeta(y,a2,b2,lower.tail=FALSE)*pbeta(y,a3,b3,lower.tail=FALSE)*
          pbeta(y,a4,b4,lower.tail=FALSE)*pbeta(y,a5,b5,lower.tail=FALSE)*
          dbeta(y,a1,b1)
       }
       integral<-integrate(f8b,0,1, rel.tol = 1e-6, stop.on.error = FALSE,...)$value
    }

  }
  

} else if  (outcome=='KV'){

  if (side=='upper'){

    if (armn==2){
      f1c<-function(y){
        pnorm(y,mean=mean2,sd=sd2,lower.tail=T)*
          dnorm(y,mean=mean1,sd=sd1)
      }
      integral<-integrate(f1c,mean1-5*sd1,mean1+5*sd1, rel.tol = 1e-6, stop.on.error = FALSE,...)$value
    } else if (armn==3){
      f2c<-function(y){
        pnorm(y,mean=mean2,sd=sd2,lower.tail=T)*
          pnorm(y,mean=mean3,sd=sd3,lower.tail=T)*
          dnorm(y,mean=mean1,sd=sd1)
      }
      integral<-integrate(f2c,mean1-5*sd1,mean1+5*sd1, rel.tol = 1e-6, stop.on.error = FALSE,...)$value
    } else if (armn==4){
      f3c<-function(y){
        pnorm(y,mean=mean2,sd=sd2,lower.tail=T)*
          pnorm(y,mean=mean3,sd=sd3,lower.tail=T)*
          pnorm(y,mean=mean4,sd=sd4,lower.tail=T)*
          dnorm(y,mean=mean1,sd=sd1)
      }
      integral<-integrate(f3c,mean1-5*sd1,mean1+5*sd1, rel.tol = 1e-6, stop.on.error = FALSE,...)$value
    } else if (armn==5){
      f4c<-function(y){
        pnorm(y,mean=mean2,sd=sd2,lower.tail=T)*
          pnorm(y,mean=mean3,sd=sd3,lower.tail=T)*
          pnorm(y,mean=mean4,sd=sd4,lower.tail=T)*
          pnorm(y,mean=mean5,sd=sd5,lower.tail=T)*
          dnorm(y,mean=mean1,sd=sd1)
      }
      integral<-integrate(f4c,mean1-5*sd1,mean1+5*sd1, rel.tol = 1e-6, stop.on.error = FALSE,...)$value
    }

  }else if (side=='lower'){

    if (armn==2){
      f5c<-function(y){
          pnorm(y,mean=mean2,sd=sd2,lower.tail=F)*
          dnorm(y,mean=mean1,sd=sd1)
      }
      integral<-integrate(f5c,mean1-5*sd1,mean1+5*sd1, rel.tol = 1e-6, stop.on.error = FALSE,...)$value
    } else if (armn==3){
      f6c<-function(y){
         pnorm(y,mean=mean2,sd=sd2,lower.tail=F)*
         pnorm(y,mean=mean3,sd=sd3,lower.tail=F)*
         dnorm(y,mean=mean1,sd=sd1)
    }
      integral<-integrate(f6c,mean1-5*sd1,mean1+5*sd1, rel.tol = 1e-6, stop.on.error = FALSE,...)$value
   } else if (armn==4){
     f7c<-function(y){
        pnorm(y,mean=mean2,sd=sd2,lower.tail=F)*
        pnorm(y,mean=mean3,sd=sd3,lower.tail=F)*
        pnorm(y,mean=mean4,sd=sd4,lower.tail=F)*
        dnorm(y,mean=mean1,sd=sd1)
     }
     integral<-integrate(f7c,mean1-5*sd1,mean1+5*sd1, rel.tol = 1e-6, stop.on.error = FALSE,...)$value
  } else if (armn==5){
    f8c<-function(y){
          pnorm(y,mean=mean2,sd=sd2,lower.tail=F)*
          pnorm(y,mean=mean3,sd=sd3,lower.tail=F)*
          pnorm(y,mean=mean4,sd=sd4,lower.tail=F)*
          pnorm(y,mean=mean5,sd=sd5,lower.tail=F)*
          dnorm(y,mean=mean1,sd=sd1)
    }
    integral<-integrate(f8c,mean1-5*sd1,mean1+5*sd1, rel.tol = 1e-6, stop.on.error = FALSE,...)$value
    }

  }
  

} else if  (outcome=='UNKV'){
  if (armn==2){
    f1a1<-function(x) pnichisq_mu1(x,par2,side)*dnichisq_mu(x,par1)
    m<-par1$mu
    max_mu<-max(par1$mu,par2$mu)
    min_mu<-min(par1$mu,par2$mu)
    s<-sqrt(par1$sigsq/par1$kappa)+abs(max_mu-min_mu)
    integral<-integrate(f1a1,m-4*s, m+4*s,rel.tol = 1e-6, stop.on.error = FALSE,...)$value
    
  }else if (armn==3){
    f1a2<-function(x) pnichisq_mu1(x,par2,side)*pnichisq_mu1(x,par3,side)*dnichisq_mu(x,par1)
    m<-par1$mu
    max_mu<-max(par1$mu,par2$mu,par3$mu)
    min_mu<-min(par1$mu,par2$mu,par3$mu)
    s<-sqrt(par1$sigsq/par1$kappa)+abs(max_mu-min_mu)
    integral<-integrate(f1a2,m-4*s, m+4*s,rel.tol = 1e-6, stop.on.error = FALSE,...)$value
    
  }else if (armn==4){
    f1a3<-function(x) pnichisq_mu(x,par2,side)*pnichisq_mu(x,par3,side)*pnichisq_mu(x,par4,side)*dnichisq_mu(x,par1)
    m<-par1$mu
    max_mu<-max(par1$mu,par2$mu,par3$mu,par4$mu)
    min_mu<-min(par1$mu,par2$mu,par3$mu,par4$mu)
    s<-sqrt(par1$sigsq/par1$kappa)+abs(max_mu-min_mu)
    integral<-integrate(f1a3,m-4*s, m+4*s,rel.tol = 1e-6, stop.on.error = FALSE,...)$value
    
  }else if (armn==5){
    f1a4<-function(x) pnichisq_mu(x,par2,side)*pnichisq_mu(x,par3,side)*pnichisq_mu(x,par4,side)*pnichisq_mu(x,par5,side)*dnichisq_mu(x,par1)
    m<-par1$mu
    max_mu<-max(par1$mu,par2$mu,par3$mu,par4$mu,par5$mu)
    min_mu<-min(par1$mu,par2$mu,par3$mu,par4$mu,par5$mu)
    s<-sqrt(par1$sigsq/par1$kappa)+abs(max_mu-min_mu)
    integral<-integrate(f1a4,m-4*s, m+4*s,rel.tol = 1e-6, stop.on.error = FALSE,...)$value
    
  }

}
 return(integral)
}

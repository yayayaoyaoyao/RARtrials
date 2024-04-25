#' @title pgreater
#' @description Calculate the futility stopping probability in Bayesian response-adaptive randomization with
#' a control group using Thall \& Wathen method for binary outcomes. The prior distributions follow
#' Beta (\eqn{beta(\alpha,\beta)}) distributions and can be specified individually for each treatment group.
#' @details This function calculates the results of \eqn{Pr(p_k>p_{control}+\delta|data)} for \code{side} equals to
#' 'upper' and the results of \eqn{Pr(p_{control}>p_k+\delta|data)} for \code{side} equals to 'lower'.
#' The result indicates the posterior probability of stopping a treatment group due to futility around \eqn{1\%} in Bayesian
#' response-adaptive randomization with a control arm using Thall \& Wathen method, with accumulated results
#' during the conduct of trials. 
#' @aliases pgreater
#' @author Chuyao Xu, Thomas Lumley, Alain Vandal
#' @export pgreater
#' @param A  \eqn{\alpha} in \eqn{beta(\alpha,\beta)}, current accumulated \eqn{\alpha} for the control group.
#' @param B  \eqn{\beta} in \eqn{beta(\alpha,\beta)}, current accumulated \eqn{\beta} for the control group.
#' @param a  \eqn{\alpha} in \eqn{beta(\alpha,\beta)}, current accumulated \eqn{\alpha} for the treatment group which is compared to the control group.
#' @param b  \eqn{\beta} in \eqn{beta(\alpha,\beta)}, current accumulated \eqn{\beta} for the treatment group which is compared to the control group.
#' @param delta  expected difference in success probabilities between the control group and the treatment group.
#' @param side direction of a one-sided test, with values 'upper' or 'lower'.
#' @param ... additional arguments to be passed to stats::integrate() (such as rel.tol) from this function.
#' @return a posterior probability of \eqn{Pr(p_k>p_{control}+\delta|data)} with \code{side} equals to 'upper';
#' a posterior probability of \eqn{Pr(p_{control}>p_k+\delta|data)} with \code{side} equals to 'lower'.
#' @importFrom stats pbeta
#' @importFrom stats dbeta
#' @examples
#' pgreater(a=5, b=19, A=8, B=10, delta=0.1, side='upper')
#' pgreater(a=58, b=68, A=65, B=79, delta=0, side='lower')
#' @references 
#' \insertRef{Wathen2017}{RARtrials}


pgreater<-function(a,b,A,B,delta,side,...){
  if (side=='upper'){
    f<-function(y){
      pbeta(y+delta,a,b,lower.tail=FALSE)*dbeta(y,A,B)
    }
  }else if (side=='lower'){
    f<-function(y){
      pbeta(y+delta,a,b,lower.tail=TRUE)*dbeta(y,A,B)
    }
  }
  integrate(f,0,1, rel.tol = 1e-6, stop.on.error = FALSE,...)$value
}

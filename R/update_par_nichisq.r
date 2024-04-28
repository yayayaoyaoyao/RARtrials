#' @title update_par_nichisq
#' @description Update hyperparameters of a Normal-Inverse-Chi-Squared distribution
#  with available data.
#' @details This function updates hyperparameters of a Normal-Inverse-Chi-Squared distribution
#' with available data to hyperparameters of a posterior Normal-Inverse-Chi-Squared distribution.
#' Those updated hyperparameters can be converted to hyperparameters in a Normal-Inverse-Gamma distribution
#' for continuous outcomes with unknown variances. 
#' @aliases update_par_nichisq
#' @author Chuyao Xu, Thomas Lumley, Alain Vandal
#' @export update_par_nichisq
#' @param y observed data.
#' @param par a vector of current parameters from a Normal-Inverse-Chi-Squared distribution.
#' @return a list of parameters for a posterior Normal-Inverse-Chi-Squared distribution incorporating available data.
#' @examples
#' para<-list(V=1/2,a=0.5,m=9.1/100,b=0.00002)
#' par<-convert_gamma_to_chisq(para)
#' set.seed(123451)
#' y1<-rnorm(100,0.091,0.009)
#' update_par_nichisq(y1, par)
#' set.seed(123452)
#' y2<-rnorm(90,0.09,0.009)
#' update_par_nichisq(y2, par)
#' @references 
#' \insertRef{Kevin2007}{RARtrials}

update_par_nichisq<-function(y, par){
  oldpar<-par
  n<-length(y)
  ybar<-mean(y)
  yss<-sum((y-ybar)^2)
  par$kappa<-oldpar$kappa+n
  par$nu<-oldpar$nu+n
  par$mu<-(oldpar$mu*oldpar$kappa+n*ybar)/par$kappa
  par$sigsq<- (1/par$nu)*(oldpar$nu*oldpar$sigsq+yss+((n*oldpar$kappa)/(n+oldpar$kappa))*(ybar-oldpar$mu)^2)
  par
}


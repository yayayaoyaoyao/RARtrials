#' @title Convert parameters from a Normal-Inverse-Chi-Squared Distribution to a
#' Normal-Inverse-Gamma Distribution 
#' @description Convert parameters from a Normal-Inverse-Chi-Squared distribution to a
#' Normal-Inverse-Gamma distribution.
#' @details This function convert parameters from a Normal-Inverse-Chi-Squared
#' (\eqn{(\mu,\sigma^2) \sim NIX({\sf mean}=\mu,{\sf effective sample size}=\kappa,{\sf degrees of freedom}=\nu,{\sf variance}=\sigma^2/\kappa)}) 
#' distribution to a Normal-Inverse-Gamma 
#' (\eqn{(\mu,\sigma^2) \sim NIG({\sf mean}=m,{\sf variance}=V \times \sigma^2,{\sf shape}=a,{\sf rate}=b)}) 
#' distribution.
#' @aliases convert_chisq_to_gamma
#' @export convert_chisq_to_gamma
#' @param cpar a list of parameters including mu, kappa, nu, sigsq from a Normal-Chi-Squared distribution. 
#' @return A list of parameters including m, V, a, b from a Normal-Inverse-Gamma distribution.
#' @examples
#' convert_chisq_to_gamma(list(mu=0.091,kappa=2,nu=1,sigsq=4e-05))
#' @references 
#' \insertRef{Kevin2007}{RARtrials}

convert_chisq_to_gamma<-function(cpar){
  list(
    m=cpar$mu,
    V=1/cpar$kappa,
    a=cpar$nu/2,
    b=cpar$nu*cpar$sigsq/2
  )
}




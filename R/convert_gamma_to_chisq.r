#' @title Convert parameters from a Normal-Inverse-Gamma Distribution to a
#' Normal-Inverse-Chi-Squared Distribution
#' @description Convert parameters from a Normal-Inverse-Gamma distribution to a
#' Normal-Inverse-Chi-Squared distribution.
#' @details This function convert parameters from a Normal-Inverse-Gamma 
#' (\eqn{(\mu,\sigma^2) \sim NIG({\sf mean}=m,{\sf variance}=V \times \sigma^2,{\sf shape}=a,{\sf rate}=b)}) 
#' distribution to a Normal-Inverse-Chi-Squared 
#' (\eqn{(\mu,\sigma^2) \sim NIX({\sf mean}=\mu,{\sf effective sample size}=\kappa,{\sf degrees of freedom}=\nu,{\sf variance}=\sigma^2/\kappa)}) 
#' distribution.
#' @aliases convert_gamma_to_chisq
#' @export convert_gamma_to_chisq
#' @param gpar a list of parameters including m, V, a, b from a Normal-Inverse-Gamma distribution.
#' @return A list of parameters including mu, kappa, nu, sigsq from a Normal-Inverse-Chi-Squared distribution.
#' @examples
#' convert_gamma_to_chisq(list(V=1/2,a=0.5,m=9.1/100,b=0.00002))
#' @references 
#' \insertRef{Kevin2007}{RARtrials}

convert_gamma_to_chisq<-function(gpar){
  list(
    mu=gpar$m,
    kappa=1/gpar$V,
    nu=2*gpar$a,
    sigsq=gpar$b/gpar$a
  )
}




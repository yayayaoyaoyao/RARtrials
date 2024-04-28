#' @title convert_gamma_to_chisq
#' @description Convert parameters from a normal-inverse-gamma distribution to a
#' normal-inverse-chi-squared distribution.
#' @details This function convert parameters from a normal-inverse-gamma distribution to a
#' normal-inverse-chi-squared distribution.
#' @aliases convert_gamma_to_chisq
#' @author Chuyao Xu, Thomas Lumley, Alain Vandal
#' @export convert_gamma_to_chisq
#' @param gpar a list of parameters including m, V, a, b from a normal-inverse-gamma distribution.
#' @return a list of parameters of a normal-inverse-gamma distribution.
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




% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/convert_gamma_to_chisq.r
\name{convert_gamma_to_chisq}
\alias{convert_gamma_to_chisq}
\title{Convert parameters from a Normal-Inverse-Gamma Distribution to a
Normal-Inverse-Chi-Squared Distribution}
\usage{
convert_gamma_to_chisq(gpar)
}
\arguments{
\item{gpar}{a list of parameters including m, V, a, b from a Normal-Inverse-Gamma distribution.}
}
\value{
A list of parameters including mu, kappa, nu, sigsq from a Normal-Inverse-Chi-Squared distribution.
}
\description{
Convert parameters from a Normal-Inverse-Gamma distribution to a
Normal-Inverse-Chi-Squared distribution.
}
\details{
This function convert parameters from a Normal-Inverse-Gamma 
(\eqn{(\mu,\sigma^2) \sim NIG({\sf mean}=m,{\sf variance}=V \times \sigma^2,{\sf shape}=a,{\sf rate}=b)}) 
distribution to a Normal-Inverse-Chi-Squared 
(\eqn{(\mu,\sigma^2) \sim NIX({\sf mean}=\mu,{\sf effective sample size}=\kappa,{\sf degrees of freedom}=\nu,{\sf variance}=\sigma^2/\kappa)}) 
distribution.
}
\examples{
convert_gamma_to_chisq(list(V=1/2,a=0.5,m=9.1/100,b=0.00002))
}
\references{
\insertRef{Kevin2007}{RARtrials}
}

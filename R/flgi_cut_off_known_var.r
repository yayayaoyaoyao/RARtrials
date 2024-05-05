#' @title Cut-off Value of Forward-looking Gittins index rule in Continuous Endpoint with Known Variances
#' @description Function for simulating cut-off values at the final stage using the forward-looking Gittins index
#' and the controlled forward-looking Gittins index algorithm for continuous outcomes with known variance in trials with
#' 2-5 arms. The conjugate prior distributions follow Normal (\eqn{N(mean,sd)}) distributions and should be the same for each arm.
#' @details This function simulates trials using the forward-looking Gittins index and the
#' controlled forward-looking Gittins index algorithm under both no delay and delay scenarios to obtain
#' cut-off values at the final stage, with control of type I error. The user is expected to run this function
#' multiple times to determine a reasonable cut-off value for statistical inference.
#' @aliases flgi_cut_off_known_var
#' @author Chuyao Xu, Thomas Lumley, Alain Vandal
#' @export flgi_cut_off_known_var
#' @param Gittinstype type of Gittins indices, should be set to 'KV' in this function.
#' @param df discount factor which is the multiplier for loss at each additional patient in the future.
#' Available values are 0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 0.99 and 0.995. The maximal sample size can be up to 10000.
#' @param gittins user specified Gittins indices for calculation in this function. If \code{gittins} is provided,
#' \code{Gittinstype} and \code{df} should be NULL.
#' @param Pats the number of patients accrued within a certain time frame indicates the
#' count of individuals who have been affected by the disease during that specific period,
#' for example, a month or a day. If this number is 10, it represents that
#' 10 people have got the disease within the specified time frame.
#' @param nMax the maximal accrued number of patients with the disease, this number
#' should be chosen carefully to ensure a sufficient number of patients are simulated,
#' especially when considering the delay mechanism.
#' @param TimeToOutcome Representation of the time distribution of delayed responses. The accrual times
#' could be a month, a week or any other time frame. When the unit changes,
#' the number of TimeToOutcome should also change. It can be in the format
#' of expression(rnorm( length( vStartTime ),30, 3)), representing delayed responses
#' with a normal distribution, where the mean is 30 days and the standard deviation is 3 days.
#' These related functions are adapted from \url{http://github.com/kwathen/IntroBayesianSimulation}.
#' Refer to the website for more details.
#' @param enrollrate probability that patients in the population can enroll in the trial.
#' This parameter is related to the number of people who have been affected by the disease in the population,
#' following an exponential distribution.
#' @param K number of total arms in the trial.
#' @param noRuns2 number of simulations for simulated allocation probabilities within each block. Default value is
#' set to 100, which is recommended in \insertCite{Villar2015}{RARtrials}.
#' @param Tsize maximal sample size for the trial.
#' @param block block size.
#' @param rule rules can be used in this function, with values 'FLGI PM', 'FLGI PD' or 'CFLGI'.
#' 'FLGI PM' stands for making decision based on posterior mean;
#' 'FLGI PD' stands for making decision based on posterior distribution;
#' 'CFLGI' stands for controlled forward-looking Gittins index.
#' @param prior_n a vector representing the number of observations assumed in prior distributions, eg: c(1,1) for a two-armed trial.
#' @param prior_mean a vector representing mean of observations assumed in prior distributions, eg: c(0,0,0) for a three-armed trial,
#' rep(0,K) can be used to simplify the process. If a negative effect is expected, adjust the mean to a negative value.
#' @param mean a vector of mean hypotheses, for example, c(0.1,0.1) where 0.1 stands for the mean
#' for both groups. Another example is c(0.1,0.3) where 0.1 and 0.3 stand for the mean for the control and
#' a treatment group, respectively.
#' @param sd a vector of standard deviation in hypotheses, for example, as c(0.64,0.64) where 0.64 stands for the standard deviation
#' for both groups. Another example is c(0.64,0.4) where 0.64 and 0.4 stand for the standard deviation for the control and
#' a treatment group, respectively.
#' @param side direction of one-sided test with the values of 'upper' or 'lower'.
#' @return Value of Z test statistics for one trial.
#' @importFrom stats runif
#' @importFrom stats rnorm
#' @examples
#' #forward-looking Gittins index rule with delayed responses follow a normal distribution
#' #with a mean of 30 days and a standard deviation of 3 days
#' set.seed(12345)
#' stopbound1<-lapply(1:10,function(x){
#' flgi_cut_off_known_var(Gittinstype='KV',df=0.995,Pats=10,nMax=50000,
#' TimeToOutcome=expression(rnorm( length( vStartTime ),30, 3)),enrollrate=0.5,
#' K=3,noRuns2=100,Tsize=120,block=8,rule='FLGI PM',prior_n=rep(1,3),
#' prior_mean=rep(0,3),mean=c(-0.05,-0.05,-0.05),sd=c(0.346,0.346,0.346),side='upper')})
#' stopbound1a<-do.call(rbind,stopbound1)
#' sum(stopbound1a[,1]>1.84 | stopbound1a[,2]>1.84 )/10 #0.1
#' #the selected cut-off value is 1.84 with an overall upper one-sided type I
#' #error of 0.1, based on 10 simulations.
#' #It is recommended to conduct more simulations to obtain an accurate cut-off value.
#'
#' #forward-looking Gittins index rule with no delay responses
#' set.seed(12345)
#' stopbound2<-replicate(10,flgi_cut_off_known_var(Gittinstype='KV',
#' df=0.995,Pats=10,nMax=50000,TimeToOutcome=0,enrollrate=1,K=2,noRuns2=100,
#' Tsize=72,block=9,rule='FLGI PM',prior_n=rep(2,3),prior_mean=rep(0,3),
#' mean=c(0.155,0.155),sd=c(0.64,0.64),side='lower'))
#' stopbound2a<-matrix(rbind(stopbound2),ncol=1)
#' sum(stopbound2a[,1]<=-1.656 )/10 #0.1
#' #the selected cut-off value is -1.656 with an overall lower one-sided type I
#' #error of 0.1, based on 10 simulations.
#' #It is recommended to conduct more simulations to obtain an accurate cut-off value.
#'
#' @references
#' \insertRef{Williamson2019}{RARtrials}


flgi_cut_off_known_var<-function(Gittinstype,df,gittins=NULL,Pats,nMax,TimeToOutcome,enrollrate,K,noRuns2,Tsize,block,rule,
                             prior_n,prior_mean,mean,sd,side ){

  if (is.null(gittins)){
    GI_Normal_known <- Gittins(Gittinstype,df)
  }else{
    GI_Normal_known <- gittins
  }

  index<-matrix(0,nrow=K,1)
  meanhat<-matrix(0,nrow=1,K)
  sigmahat<-matrix(0,nrow=1,K)
  GI_Std<-rep(0,K)
  zs1<-matrix(0,nrow=1,K-1)
  ap1<-matrix(0,nrow=1,K-1)

  popdat<-pop(Pats,nMax,enrollrate)
  vStartTime<-sort(popdat[[3]][1:Tsize], decreasing = FALSE)
  vOutcomeTime<-SimulateOutcomeObservedTime(vStartTime,TimeToOutcome)

  data1<-matrix(NA_real_,nrow=Tsize,ncol=5)
  data1[,1]<-1:Tsize
  data1[,2]<-vStartTime
  data1[,3]<-vOutcomeTime
  n=matrix(NA,nrow=K,1)
  nn<-rep(0,K)
  sample_mean<-rep(NA,K)

  for (t in 0:((Tsize/block)-1)){

    alp=allocation_probabilities_kn_var(GI_Normal_known=GI_Normal_known,tt=t,data1=data1,arms=K,b=block,runs=noRuns2,
                                        prior_mean=prior_mean,prior_n=prior_n,sd1=sd,side=side)

    if (rule=='Controlled FLGI'  ){
      alp[1]=1/(K-1)
      elp_e=allocation_probabilities_kn_var1(GI_Normal_known=GI_Normal_known,tt=t,data1=data1,arms=K,b=block,runs=noRuns2,
                                             prior_mean=prior_mean,prior_n=prior_n,sd1=sd,side=side)
      c=alp[1]+sum(elp_e)
      alp=(1/c)*c(alp[1],elp_e)
    }

    alp=cumsum(c(0,alp))
    Pob<-rep(0,block)
    Pos<-rep(0,block)
    for (p in 1:block){
      Pob[p]<-runif(1)
      for (k in 1:K){
        if (Pob[p]>alp[k] & Pob[p]<=alp[k+1]){
          Pos[p]=rnorm(1, mean[k], sd[k])
          data1[t*block+p,4]=k
          data1[t*block+p,5]=Pos[p]
        }
      }

    }
  }

  #if (floor(Tsize/block)*block!=Tsize){
  if ((Tsize %% block)!=0){
    Pob<-rep(0,Tsize %% block)
    Posi<-rep(0,Tsize %% block)
    for (p in 1:(Tsize %% block)){
      Pob[p]<-runif(1)
      for (k in 1:K){
        if (Pob[p]>alp[k] & Pob[p]<=alp[k+1]){
          Posi[p]=rnorm(1, mean[k], sd[k])
          data1[floor(Tsize/block)*block+p,4]=k
          data1[floor(Tsize/block)*block+p,5]=Posi[p]
        }
      }
    }
  }

  sdd<-rep(NA,K)
  sddd<-rep(NA,K)
  for (k in 1:K){
    n[k,1]=nrow(data1[data1[,4]==k ,,drop=FALSE])
    meanhat[1,k]=(sum(data1[data1[,4]==k,5])) / n[k,1]
    indexs<-n[k,1]+prior_n[k]
    GI_Std[k] <- GI_Normal_known[ indexs ]
  }

  pc<-matrix(0,1,K-1)


  for (k in 1:(K-1)){
    zs1[1,k]=(meanhat[1,k+1]-meanhat[1,1])/sqrt((sd[1])^2/n[1,1]+(sd[k+1])^2/n[k+1,1])
  }

  return(zs1)

}



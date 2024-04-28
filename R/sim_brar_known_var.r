#' @title sim_brar_known_var
#' @description \code{sim_brar_known_var} simulate a trial with two to five arms using Bayesian Response-Adaptive 
#' Randomization with a control group for continuous outcomes with known variances. The prior distributions follow
#' Normal (\eqn{N(mean,sd)}) distributions and can be specified individually for each arm.
#' @details This function generates a designed trial using Bayesian response-adaptive randomization with
#' a control group under no delay and delay scenarios for continuous outcomes with known variances. The function can handle trials with up to
#' 5 arms. This function uses the formula
#' \eqn{\frac{Pr(p_k=max\{p_1,...,p_K\})^tp} {\sum_{k=1}^{K}{Pr(p_k=max\{p_1,...,p_K\})^tp}}} with \code{side} equals to 'upper'
#' and \eqn{\frac{Pr(p_k=min\{p_1,...,p_K\})^tp} {\sum_{k=1}^{K}{Pr(p_k=min\{p_1,...,p_K\})^tp}}}
#' with \code{side} equals to 'lower', utilizing available data at each step.
#' @aliases sim_brar_known_var
#' @author Chuyao Xu, Thomas Lumley, Alain Vandal
#' @export sim_brar_known_var
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
#' @param N1 number of participants with equal randomization in the 'initialization' period.
#' Recommend using 10 percent of the total sample size.
#' @param armn number of total arms in the trial.
#' @param au a vector of cut-off values in the final selection at the end of the trial,
#' with a length equal to the number of arms minus 1.
#' @param N2 maximal sample size for the trial.
#' @param tp tuning parameter. Some commonly used numbers are 0.5, 1 and n/2N.
#' @param armlabel a vector of treatment labels with an example of c(1, 2), where 1 and 2 describe
#' how each arm is labeled in a two-armed trial.
#' @param blocksize size of block used for equal randomization regarding participants in the 'initialization' period.
#' Recommend to be an even multiple of the number of total arms.
#' @param mean a vector of means in hypotheses, for example, as c(10,10) where 10 stands for the mean
#' in both groups. Another example is c(10,12) where 10 and 12 stand for the mean
#' for the control and the treatment group, respectively.
#' @param sd a vector of standard deviations in hypotheses, for example, as c(2,2) where 2 stands for the standard deviation
#' in both groups. Another example is c(1,2) where 1 and 2 stand for the standard deviation
#' for the control and the treatment group, respectively.
#' @param minstart a specified number of participants when one starts to check decision rules.
#' @param deltaa a vector of minimal effect expected to be observed for early futility stopping in
#' each arm is approximately \eqn{1\%}. The length of this parameter is \code{armn}-1.
#' @param tpp indicator of \code{tp} equals to n/2N. When \code{tp} is n/2N, \code{tpp} should be assigned 1. Default value is set to 0.
#' @param deltaa1 a vector of pre-specified minimal effect size expected to be observed at the final stage
#' for each arm. The length of this parameter is \code{armn}-1.
#' @param mean10 prior mean of arm 1 in the trial, which stands for the control. Default value is set to 1.
#' @param mean20 prior mean of arm 2 in the trial. Default value is set to \code{mean10}.
#' @param mean30 prior mean of arm 3 in the trial. Default value is set to \code{mean10}.
#' @param mean40 prior mean of arm 4 in the trial. Default value is set to \code{mean10}.
#' @param mean50 prior mean of arm 5 in the trial. Default value is set to \code{mean10}.
#' @param sd10 prior standard deviation of arm 1 in the trial, which stands for the control. Default value is set to 0.
#' @param sd20 prior standard deviation of arm 2 in the trial. Default value is set to \code{sd10}.
#' @param sd30 prior standard deviation of arm 3 in the trial. Default value is set to \code{sd10}.
#' @param sd40 prior standard deviation of arm 4 in the trial. Default value is set to \code{sd10}.
#' @param sd50 prior standard deviation of arm 5 in the trial. Default value is set to \code{sd10}.
#' @param n10 explicit prior n of arm 1 in the trial, which stands for the control. Default value is set to 1.
#' @param n20 explicit prior n of arm 2 in the trial. Default value is set to \code{n10}.
#' @param n30 explicit prior n of arm 3 in the trial. Default value is set to \code{n10}.
#' @param n40 explicit prior n of arm 4 in the trial. Default value is set to \code{n10}.
#' @param n50 explicit prior n of arm 5 in the trial. Default value is set to \code{n10}.
#' @param side direction of a one-sided test, with values 'upper' or 'lower'.
#' @param ... additional arguments to be passed to \code{\link[stats]{integrate}} (such as rel.tol) from this function.
#' @return A list of results, including final decision, test statistics, the simulated data set
#' and participants accrued for each arm at the time of termination of that group in one trial.
#' The simulated data set includes 5 columns: participant ID number, enrollment time, observed time of results,
#' allocated arm, and participants' results.
#' In the final decision, 'Superiorityfinal' refers to the selected arm, while 'Not Selected' indicates the arm stopped due to
#' futility, and 'Control Selected' denotes the control arm chosen because other arms did not meet futility criteria before the 
#' final stage or were not deemed effective at the final stage. 
#' Note that before final stage of the trial, test statistics is calculated from \code{deltaa}, and test statistics is
#' calculated from \code{deltaa1} at the final stage.
#' @importFrom stats pnorm
#' @importFrom stats rnorm
#' @importFrom Rdpack reprompt
#' @examples
#' #sim_brar_known_var with delayed responses follow a normal distribution with
#' #a mean of 30 days and a standard deviation of 3 days, where mean=c(8.74/100,8.74/100,8.74/100),
#' #sd=c(0.009,0.009,0.009), tp=0.5 and the minimal effect size is 0.
#' sim_brar_known_var(Pats=10,nMax=50000,TimeToOutcome=expression(rnorm(
#' length(vStartTime),30, 3)),enrollrate=0.1, N1=21,armn=3,au=c(0.901,0.901),
#' N2=189,tp=0.5,armlabel=c(1,2,3),blocksize=6,mean=c(8.74/100,8.74/100,8.74/100),
#' sd=c(0.009,0.009,0.009),minstart=21,deltaa=c(0.00048,0.00048),tpp=0,deltaa1=c(0,0),
#' mean10=0.09,mean20=0.09,mean30=0.09,sd10=0.01,sd20=0.01,sd30=0.01,n10=1,n20=1,n30=1,side='lower')
#' @references 
#' \insertRef{Wathen2017}{RARtrials}

sim_brar_known_var<-function(Pats,nMax,TimeToOutcome,enrollrate,N1,armn,au,N2,tp,armlabel,blocksize,
                             mean,sd,minstart,deltaa,tpp,deltaa1,mean10=0,mean20=mean10,mean30=mean10,mean40=mean10,mean50=mean10,
                             sd10=1,sd20=sd10,sd30=sd10,sd40=sd10,sd50=sd10,n10=1,n20=n10,n30=n10,n40=n10,n50=n10,side,...){

  popdat<-pop(Pats,nMax,enrollrate)
  vStartTime<-sort(popdat[[3]][1:N2], decreasing = FALSE)
  vOutcomeTime<-SimulateOutcomeObservedTime(vStartTime,TimeToOutcome)
  assign1<-blockrand(blocksize,N1,armn,armlabel)

  data1<-matrix(NA_real_,nrow=N2,ncol=5)
  data1[,1]<-1:N2
  data1[,2]<-vStartTime
  data1[,3]<-vOutcomeTime
  data1[1:N1,4]<-assign1$arm[1:N1]

  for (i in 1:(N1)){
    for (j in 1:armn) {
      if (data1[i, 4]==j ){
        data1[i,5]<-rnorm(1,mean=mean[j],
                          sd=sd[j])
      }
    }
  }

  armleft<-c(1:armn)
  decision<-rep(NA,armn )
  phi<-rep(NA,armn )
  stopp<-rep(NA,armn )

  mean0<-list(mean10,mean20,mean30,mean40,mean50)
  sd0<-list(sd10,sd20,sd30,sd40,sd50)
  n0<-list(n10,n20,n30,n40,n50)
  for (jjj in minstart:N2){


    if (jjj>minstart){
      treat<-sample(armleft,size =1, prob = as.vector(pii))
      data1[jjj,4]<-treat
      data1[jjj,5]<-rnorm(1,mean=mean[treat],
                          sd=sd[treat])

    }

    if (jjj<N2){
      total<-sum (as.numeric(data1[1:jjj,3])<=as.numeric(data1[jjj,2]))
    }else if (jjj==N2){
      total<-N2
    }

    result<-vector("list",length(armleft))
    mat<-vector("list",armn)

    for (j in 1:length(armleft)) {
      
      
      if (jjj!=N2){
        data2<-matrix(data1[which(as.numeric(data1[1:jjj,3])<=as.numeric(data1[jjj,2])),],ncol=5)
      }else if (jjj==N2){
        data2<-data1
      }
      tot<-as.numeric(data2[which(data2[,4]==as.numeric(armleft[j])),5])
      
      if (identical(tot, numeric(0))){
        mat[[armleft[j]]]<-matrix(c( mean0[[1]],   
                                     0,
                                     sd0[[armleft[j]]]^2),
                                  nrow=1)
      }
      
      mat[[armleft[j]]]<-matrix(c( (1/((length(tot)/(sd[armleft[j]]^2))+
                                         n0[[armleft[j]]]/(sd0[[armleft[j]]]^2)))*
                                     (mean0[[armleft[j]]]/(sd0[[armleft[j]]]^2)+sum(tot)/(sd[armleft[j]]^2)),   
                                   length(tot),
                                   sqrt(1/((length(tot)/(sd[armleft[j]]^2))+
                                             n0[[armleft[j]]]/(sd0[[armleft[j]]]^2))))
                                ,nrow=1)
      
      
      
      
    }

    if (length(armleft)>1){
      for (j in 1:length(armleft)){
        #if (total>0){
        if (j>1){
          if (side=='lower'){
            result[[j]]<-pnorm(deltaa[armleft[j]-1], mat[[1]][1,1]-mat[[armleft[j]]][1,1],
                               sqrt(mat[[1]][1,3]^2+mat[[armleft[j]]][1,3]^2),lower.tail=FALSE)
            
          }else if (side=='upper'){
            result[[j]]<-pnorm(deltaa[armleft[j]-1], mat[[j]][1,1]-mat[[armleft[1]]][1,1],
                               sqrt(mat[[1]][1,3]^2+mat[[armleft[j]]][1,3]^2),lower.tail=FALSE)
            
          }
        }else if (j==1){
          result[[1]]<-0
        }
        
        
        
      }
      
      aloo<-vector("list",length(armleft))
      aloo<-alofun_kn_var(mat=mat,total=total,armleft=armleft,side=side)
      
    }



    if (jjj==N2){
      resultt<-vector("list",length(armleft))
      if (length(armleft)>1){
        for (j in 1:length(armleft)){
          # if (total>0){
          if (j>1){
            if (side=='lower'){
              resultt[[j]]<-pnorm(deltaa1[armleft[j]-1], mat[[1]][1,1]-mat[[armleft[j]]][1,1],
                                  sqrt(mat[[1]][1,3]^2+mat[[armleft[j]]][1,3]^2),lower.tail=FALSE)
            }else if (side=='upper'){
              resultt[[j]]<-pnorm(deltaa1[armleft[j]-1], mat[[armleft[j]]][1,1]-mat[[1]][1,1],
                                  sqrt(mat[[1]][1,3]^2+mat[[armleft[j]]][1,3]^2),lower.tail=FALSE)
            }
          }else if (j==1){
            resultt[[1]]<-0
          }
          
          
          
        }
      }
      
    }


      posteriorp<-vector("list",length(armleft))
      for (j in 1:length(armleft)){
        if (total>0){
          if (j>1){
            if (side=='lower'){
              posteriorp[[j]]<-pnorm(deltaa1[armleft[j]-1], mat[[1]][1,1]-mat[[armleft[j]]][1,1],
                                     sqrt(mat[[1]][1,3]^2+mat[[armleft[j]]][1,3]^2),lower.tail=FALSE)

            }else if (side=='upper'){
              posteriorp[[j]]<-pnorm(deltaa1[armleft[j]-1], mat[[armleft[j]]][1,1]-mat[[1]][1,1],
                                     sqrt(mat[[1]][1,3]^2+mat[[armleft[j]]][1,3]^2),lower.tail=FALSE)

            }

          }else if (j==1){
            posteriorp[[j]]<-0
          }

        }else if (total==0 ){
          if (j>1){
            if (side=='lower'){
              posteriorp[[j]]<-pnorm(deltaa1[armleft[j]-1], mean0[[1]]-mean0[[armleft[j]]],
                                     sqrt(sd[1]^2+sd[armleft[j]]^2),lower.tail=FALSE)

            }else  if (side=='upper'){
              posteriorp[[j]]<-pnorm(deltaa1[armleft[j]-1], mean0[[armleft[j]]]-mean0[[1]],
                                     sqrt(sd[1]^2+sd[armleft[j]]^2),lower.tail=FALSE)

            }

          }else if (j==1){
            posteriorp[[j]]<-0
          }
        }
      }





    posteriorp1<-do.call(cbind,posteriorp)

    pii<-as.data.frame(do.call(cbind,aloo))
    colnames(pii)<-armleft
    if (jjj<N2){
      sim111<-do.call(cbind,result)
      colnames(sim111)<-armleft
    } else if (jjj==N2) {
      sim111t<-do.call(cbind,resultt)
      colnames(sim111t)<-armleft
    }

    if ( jjj<N2){
      for (k in 2:length(armleft)) {

        if (sim111[1,k]<0.01 & is.na(decision[armleft[k]])){
          decision[armleft[k]]<-'Futility'
          stopp[armleft[k]]<-jjj
          phi[armleft[k]]<-posteriorp1[1,k]

        }
      }

      if ( 'Futility' %in% decision ){
        armleft<-armleft[! armleft %in%  which (decision %in%  'Futility')]
        pii<- pii[,sprintf("%s",armleft)]

      }

      if((length(armleft)==1 & length( which (decision %in%  'Futility'))==(armn-1))){
        stopp[ which (is.na(decision))]<-jjj#decision %in%  NA
        data11<-data1[1:jjj,]
        if (is.na(decision[1])==TRUE) {#decision[1] %in% NA
          decision[1]<-'Control Selected'
          phi[1]<-posteriorp1[1,1]

        }

        return(list(decision,phi,data11,stopp))

        break
      }

      for (yy in 1:length(armleft)){
        if ( pii[colnames(pii) %in%  armleft[yy]]<0.1){
          pii[colnames(pii) %in%  armleft[yy]]=0.1
        }else if ( pii[colnames(pii) %in%  armleft[yy]]>0.9){
          pii[colnames(pii) %in%  armleft[yy]]=0.9
        }
      }


      if (tpp==1){
        pii<-(pii^(jjj/(2*N2)))^tp
      }else if (tpp==0){
        pii<-pii^tp
      }
      pii<-pii/sum(pii)
      pii<-pii[ , order(names(pii))]
      armleft<-sort(armleft,decreasing = FALSE)

    }else if (jjj==N2){

      for (k in 2:length(armleft)) {
        if (sim111t[1,k-1]>au[armleft[k]-1] & is.na(decision[armleft[k]])){
          decision[armleft[k]]<-'Superiorityfinal'
          stopp[armleft[k]]<-jjj
          phi[armleft[k]]<-posteriorp1[1,k]

        }
      }
      for (k in 1:length(armleft)) {
        if (is.na(decision[armleft[k]] )) {
          decision[armleft[k]]<-'Not Selected'
          stopp[armleft[k]]<-jjj
          phi[armleft[k]]<-posteriorp1[1,k]

        }
      }
      if (length( which (decision %in%  'Futility'))==(armn-1) &
          is.na(decision[1])){
        decision[1]<-'Control Selected'
        phi[1]<-posteriorp1[1,1]

      }
      data11<-data1[1:N2,]
      return(list(decision,phi,data11,stopp))

    }
  }

}

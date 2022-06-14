#' @name rmst2
#' @aliases rmst2
#' @title Comparing restricted mean survival time
#' @description Performs two-sample comparisons using the restricted mean survival time (RMST) as a summary measure of the survival time distribution.
#' Three kinds of between-group contrast metrics (i.e., the difference in RMST, the ratio of RMST and the ratio of the restricted mean time lost (RMTL)) are computed.
#' The Greenwood plug-in estimator is used for the asymptotic variance. It performs ANCOVA-type adjusted analyses when covariates are passed to it as an argument.
#' @usage  rmst2(time, status, arm, tau = NULL, covariates = NULL, alpha = 0.05)
#' @param time The follow-up time for right censored data.
#' @param status The status indicator, 1=event, and 0=right censored.
#' @param arm The group indicator for comparison. The elements of this vector take either 1 or 0. Normally, 0=control group, 1=active treatment group.
#' @param tau A scaler value to specify the truncation time point for the RMST calculation.
#' When \code{tau = NULL}, the default value is used. See Details for the definition of the default tau.
#' @param covariates This specifies covariates to be used for the adjusted analyses. When NULL, unadjusted analyses are performed.
#' When non NULL, the ANCOVA-type adjusted analyses are performed using those variables passed as \code{covariates}.
#' This can be one variable (vector) or more than one variables (matrix).
#' @param alpha The default is 0.05. (1-\code{alpha}) confidence intervals are reported.
#' @details The definition of the default tau.
#' Let x1 and x0 be the maximum observed time in Group 1 and Group 0, respectively.
#' Case 1: if the last observations in Group 1 and Group 0 are "event," then
#' tau = max(x1, x0).
#' Case 2-1: if the last observation in Group 1 is "event," the last observation in Group 0 is "censor," and x1 <= x0,
#' tau = max(x1, x0) = x0.
#' Case 2-2: if the last observation in Group 0 is "event," the last observation in Group 1 is "censor," and x1 > x0,
#' tau = max(x1, x0) = x1.
#' Case 3-1: if the last observation in Group 1 is "event," the last observation in Group 0 is "censor," and x1 > x0,
#' tau = min(x1, x0) = x0.
#' Case 3-2: if the last observation in Group 0 is "event," the last observation in Group 1 is "censor," and x1 <= x0,
#' tau = min(x1, x0) = x1.
#' Case 4: the last observations in Group 1 and Group 0 are "censor," then
#' tau = min(x1, x0).
#'
#' @return an object of class rmst2.
#' @return \item{tau}{the truncation time used in the analyses}
#' @return \item{note}{a note regarding the truncation time}
#' @return \item{RMST.arm1}{RMST results in arm 1. This is generated only when \code{covariates} is not specified.}
#' @return \item{RMST.arm0}{RMST results in arm 0. This is generated only when \code{covariates} is not specified.}
#' @return \item{unadjusted.result}{Results of the unadjusted analyses. This is generated only when \code{covariates} is not specified.}
#' The values below are generated when some covariates are passed to the function.
#' @return \item{adjusted.result}{Results of the adjusted analyses.}
#' @return \item{RMST.difference.adjusted}{Results of the parameter estimates with the model to derive an adjusted difference in RMST.}
#' @return \item{RMST.ratio.adjusted}{Results of the parameter estimates with the model to derive an adjusted ratio of RMST.}
#' @return \item{RMTL.ratio.adjusted}{Results of the parameter estimates with the model to derive an adjusted ratio of RMTL.}
#' @references
#' Uno H, Claggett B, Tian L, Inoue E, Gallo P, Miyata T, Schrag D,
#' Takeuchi M, Uyama Y, Zhao L, Skali H, Solomon S, Jacobus S, HughesM,
#' Packer M, Wei LJ. Moving beyond the hazard ratio in quantifying the between-group difference in survival analysis.
#' Journal of clinical Oncology 2014, 32, 2380-2385. doi:10.1200/JCO.2014.55.2208.
#'
#' Tian L, Zhao L,  Wei LJ. Predicting the restricted mean event time with the subject's baseline covariates in survival analysis. Biostatistics 2014, 15, 222-233. doi:10.1093/biostatistics/kxt050.
#' @author Hajime Uno, Lu Tian, Miki Horiguchi, Angel Cronin, Chakib Battioui, James Bell
#' @examples
#' #--- sample data ---#
#' D=rmst2.sample.data()
#' time=D$time
#' status=D$status
#' arm=D$arm
#' tau=NULL
#' x=D[,c(4,6,7)]
#' #--- without covariates ----
#' a=rmst2(time, status, arm, tau=10)
#' print(a)
#' plot(a, xlab="Years", ylab="Probability", density=60)
#' #--- with covariates ----
#' a=rmst2(time, status, arm, tau=10, covariates=x)
#' print(a)


#'@export
#########################################
# rmst2 (2-arm) contrast (main function)
#########################################
#rmst2=function(time, status, arm, tau=NULL, covariates=NULL, adjust.method="reg", alpha=0.05){
rmst2=function(time, status, arm, tau=NULL, covariates=NULL,                      alpha=0.05){
  #-- time
  #-- statuts
  #-- arm (1 or 0)
  #-- covariates (matrix)
  #-- adjust = "reg"-- regression ("aug" -- augumentation)
  #-- alpha=0.05

  #==================================
  #  initial check
  #==================================

  #===== tau =====
  idx=arm==0; tt=time[idx]; tt0max=max(tt); ss=status[idx]; ss0max=min(ss[tt==tt0max]);
  idx=arm==1; tt=time[idx]; tt1max=max(tt); ss=status[idx]; ss1max=min(ss[tt==tt1max]);

  ttmax = max(tt0max, tt1max)
  ttmin = min(tt0max, tt1max)

  #--case 1: the last obs (smaller one)=event, the last obs (longer one)=event
  if(ss0max==1 & ss1max==1){
    if(!is.null(tau)){
      if(tau>ttmax){stop(paste("The truncation time, tau, needs to be shorter than or equal to ", round(ttmax, digits=2)))}
      if(tau<=ttmax){tau=tau; NOTE=paste("The truncation time: tau =", tau, " was specified.")}
    }else{
      tau = ttmax
      NOTE=(paste("The truncation time, tau, was not specified. Thus, the default tau ", round(ttmax, digits=2)," is used."))
    }
  }

  #--case 2: the last obs (smaller one)=event, the last obs (longer one)=censor
  if((ss0max==0 & ss1max==1 & tt0max>=tt1max) | (ss0max==1 & ss1max==0 & tt1max>tt0max)){
    if(!is.null(tau)){
      if(tau>ttmax){stop(paste("The truncation time, tau, needs to be shorter than or equal to ", round(ttmax, digits=2)))}
      if(tau<=ttmax){tau=tau; NOTE=paste("The truncation time: tau =", tau, " was specified.")}
    }else{
      tau = ttmax
      NOTE=paste("The truncation time, tau, was not specified. Thus, the default tau ", round(ttmax, digits=2)," is used.")
    }
  }

  #--case 3: the last obs (smaller one)=censor, the last obs (longer one)=event
  if((ss0max==1 & ss1max==0 & tt0max>=tt1max) | (ss0max==0 & ss1max==1 & tt1max>tt0max)){
    if(!is.null(tau)){
      if(tau>ttmin){stop(paste("The truncation time, tau, needs to be shorter than or equal to ", round(ttmin, digits=2)))}
      if(tau<=ttmin){tau=tau; NOTE=paste("The truncation time: tau =", tau, " was specified.")}
    }else{
      tau = ttmin
      NOTE=(paste("The truncation time, tau, was not specified. Thus, the default tau ", round(ttmin, digits=2)," is used."))
    }
  }

  #--case 4: the last obs (smaller one)=censor, the last obs (longer one)=censor
  if(ss0max==0 & ss1max==0){
    if(!is.null(tau)){
      if(tau<=ttmin){
        NOTE=paste("The truncation time: tau =", tau, " was specified.")
      }
      if(tau>ttmin){
        stop(paste("The truncation time, tau, needs to be shorter than or equal to the minimum of the largest observed time on each of the two groups: ", round(ttmin, digits=2)))
      }
    }else{
      tau = ttmin
      NOTE=(paste("The truncation time, tau, was not specified. Thus, the default tau ", round(ttmin, digits=2)," is used."))
    }
  }


  Z=list()
  Z$tau=tau
  Z$note=NOTE

  #==================================
  #  unadjusted analysis
  #==================================
  if(is.null(covariates)){

    wk1=rmst1(time[arm==1], status[arm==1], tau, alpha)
    wk0=rmst1(time[arm==0], status[arm==0], tau, alpha)

    Z$RMST.arm1=wk1
    Z$RMST.arm0=wk0


    #--- contrast (RMST difference) ---
    rmst.diff.10     = wk1$rmst[1]-wk0$rmst[1]
    rmst.diff.10.se  = sqrt(wk1$rmst.var + wk0$rmst.var)
    rmst.diff.10.low = rmst.diff.10 - qnorm(1-alpha/2)*rmst.diff.10.se
    rmst.diff.10.upp = rmst.diff.10 + qnorm(1-alpha/2)*rmst.diff.10.se
    rmst.diff.pval   = pnorm(-abs(rmst.diff.10)/rmst.diff.10.se)*2
    rmst.diff.result = c(rmst.diff.10, rmst.diff.10.low, rmst.diff.10.upp, rmst.diff.pval)

    #--- contrast (RMST ratio) ---
    rmst.log.ratio.10     = log(wk1$rmst[1]) - log(wk0$rmst[1])
    rmst.log.ratio.10.se  = sqrt(wk1$rmst.var/wk1$rmst[1]/wk1$rmst[1] + wk0$rmst.var/wk0$rmst[1]/wk0$rmst[1])
    rmst.log.ratio.10.low = rmst.log.ratio.10 - qnorm(1-alpha/2)*rmst.log.ratio.10.se
    rmst.log.ratio.10.upp = rmst.log.ratio.10 + qnorm(1-alpha/2)*rmst.log.ratio.10.se
    rmst.log.ratio.pval   = pnorm(-abs(rmst.log.ratio.10)/rmst.log.ratio.10.se)*2
    rmst.ratio.result     = c(exp(rmst.log.ratio.10), exp(rmst.log.ratio.10.low), exp(rmst.log.ratio.10.upp),rmst.log.ratio.pval)

    #--- contrast (RMTL ratio  0/1) ---
    # rmtl.log.ratio.01     = log(wk0$rmtl[1]) - log(wk1$rmtl[1])
    # rmtl.log.ratio.01.se  = sqrt(wk1$rmst.var/wk1$rmtl[1]/wk1$rmtl[1] + wk0$rmst.var/wk0$rmtl[1]/wk0$rmtl[1])
    # rmtl.log.ratio.01.low = rmtl.log.ratio.01 - qnorm(1-alpha/2)*rmtl.log.ratio.01.se
    # rmtl.log.ratio.01.upp = rmtl.log.ratio.01 + qnorm(1-alpha/2)*rmtl.log.ratio.01.se
    # rmtl.log.ratio.pval   = pnorm(-abs(rmtl.log.ratio.01)/rmtl.log.ratio.01.se)*2
    # rmtl.ratio.result     = c(exp(rmtl.log.ratio.01), exp(rmtl.log.ratio.01.low), exp(rmtl.log.ratio.01.upp),rmtl.log.ratio.pval)

    #--- contrast (RMTL ratio  1/0) ---
    rmtl.log.ratio.10     = log(wk1$rmtl[1]) - log(wk0$rmtl[1])
    rmtl.log.ratio.10.se  = sqrt(wk1$rmst.var/wk1$rmtl[1]/wk1$rmtl[1] + wk0$rmst.var/wk0$rmtl[1]/wk0$rmtl[1])
    rmtl.log.ratio.10.low = rmtl.log.ratio.10 - qnorm(1-alpha/2)*rmtl.log.ratio.10.se
    rmtl.log.ratio.10.upp = rmtl.log.ratio.10 + qnorm(1-alpha/2)*rmtl.log.ratio.10.se
    rmtl.log.ratio.pval   = pnorm(-abs(rmtl.log.ratio.10)/rmtl.log.ratio.10.se)*2
    rmtl.ratio.result     = c(exp(rmtl.log.ratio.10), exp(rmtl.log.ratio.10.low), exp(rmtl.log.ratio.10.upp),rmtl.log.ratio.pval)



    #--- results ---
    out=rbind(rmst.diff.result, rmst.ratio.result , rmtl.ratio.result )
    rownames(out)=c("RMST (arm=1)-(arm=0)","RMST (arm=1)/(arm=0)","RMTL (arm=1)/(arm=0)")
    colnames(out)=c("Est.", paste("lower .",round((1-alpha)*100, digits=0), sep=""), paste("upper .",round((1-alpha)*100, digits=0), sep=""), "p")


    #--- output ---
    Z$unadjusted.result = out
    # Z$RMST.difference=out[1,]
    # Z$RMST.ratio=out[2,]
    # Z$RMTL.ratio=out[3,]
    # Z
  }

  #==================================
  #  Adjusted analysis
  #==================================
  if (!is.null(covariates)){

    ## if (adjust.method=="reg"){

    aa=rmst2reg(time, status, as.matrix(cbind(arm, covariates)), arm, tau, type="difference", alpha=alpha)
    bb=rmst2reg(time, status, as.matrix(cbind(arm, covariates)), arm, tau, type="ratio", alpha=alpha)
    cc=rmst2reg(time, status, as.matrix(cbind(arm, covariates)), arm, tau, type="lossratio", alpha=alpha)

    #--- output ---
    out.adj=rbind(aa[2,c(1,5,6,4)], bb[2, c(5,6,7,4)], cc[2, c(5,6,7,4)])
    rownames(out.adj)=c("RMST (arm=1)-(arm=0)","RMST (arm=1)/(arm=0)","RMTL (arm=1)/(arm=0)")
    colnames(out.adj)=c("Est.", paste("lower .",round((1-alpha)*100, digits=0), sep=""), paste("upper .",round((1-alpha)*100, digits=0), sep=""), "p")


    #--- output ---
    Z$adjusted.result = out.adj

    Z$RMST.difference.adjusted = aa
    Z$RMST.ratio.adjusted      = bb
    Z$RMTL.ratio.adjusted      = cc

    ## }else{
    ##  stop "Please sepcify adjust.method"
    ## }

  }


  class(Z)="rmst2"

  Z

}
NULL

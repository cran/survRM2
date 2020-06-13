#####################
# rmst2reg -- hidden
#####################
rmst2reg=function(y, delta, x, arm, tau, type="difference", alpha=0.05){

  if(type!="difference" && type!="ratio" && type!="lossratio")
    print("Type must be difference, ratio or lossratio.")

  if(type=="difference" || type=="ratio" || type=="lossratio"){

    n=length(y)
    x=cbind(1, x)
    p=length(x[1,])

    y0=pmin(y, tau)
    d0=delta
    d0[y0==tau]=1

    d10=d0[arm==1]
    d00=d0[arm==0]
    y10=y0[arm==1]
    y00=y0[arm==0]
    x1=x[arm==1,]
    x0=x[arm==0,]
    n1=length(d10)
    n0=length(d00)

    id1=order(y10)
    y10=y10[id1]
    d10=d10[id1]
    x1=x1[id1,]

    id0=order(y00)
    y00=y00[id0]
    d00=d00[id0]
    x0=x0[id0,]

    fitc1=func_surv(y10, 1-d10)
    fitc0=func_surv(y00, 1-d00)

    weights1=d10/rep(fitc1$surv, table(y10))
    weights0=d00/rep(fitc0$surv, table(y00))

    weights=c(weights1, weights0)

    if(type=="difference")
    {fitt=lm(c(y10,y00)~rbind(x1, x0)-1, weights=weights)
    beta0=fitt$coef

    error1=y10-as.vector(x1%*%beta0)
    score1=x1*weights1*error1

    error0=y00-as.vector(x0%*%beta0)
    score0=x0*weights0*error0
    }

    if(type=="ratio")
    {fitt=glm(c(y10,y00)~rbind(x1, x0)-1, family="quasipoisson", weights=weights)
    beta0=fitt$coef

    error1=y10-exp(as.vector(x1%*%beta0))
    score1=x1*weights1*error1

    error0=y00-exp(as.vector(x0%*%beta0))
    score0=x0*weights0*error0
    }


    if(type=="lossratio")
    {fitt=glm(c(tau-y10,tau-y00)~rbind(x1, x0)-1, family="quasipoisson", weights=weights)
    beta0=fitt$coef

    error1=tau-y10-exp(as.vector(x1%*%beta0))
    score1=x1*weights1*error1

    error0=tau-y00-exp(as.vector(x0%*%beta0))
    score0=x0*weights0*error0
    }

    #---kappa.arm1---
    kappa1   <- score1
    y10table <- table(y10)
    y10revcumsum <- (n1+1)-as.vector(rev(cumsum(rev(y10table))))
    y10pos <- rep(y10revcumsum,y10table)
    y10loc <- n1+1-y10pos
    y10forward <- rep(as.vector(cumsum(y10table)),y10table)

    tab <- matrix(rep(NA,n1*p),ncol=p)
    for(i in 1:p){
      temp   <- rev(cumsum(rev(score1[,i])))
      tab[,i]<- temp[y10pos]
    }

    kappa2  <- tab*(1-d10)/y10loc
    kappa3a <- kappa2/y10loc
    kappa3b <- matrix(rep(NA,p*n1),ncol=p)
    for(i in 1:p){
      kappa3b[,i] <- cumsum(kappa3a[,i])
    }
    kappa3c    <- kappa3b[y10forward,]
    kappa.arm1 <- kappa1+kappa2-kappa3c


    #---kappa.arm0---
    kappa1   <- score0
    y00table <- table(y00)
    y00revcumsum <- (n0+1)-as.vector(rev(cumsum(rev(y00table))))
    y00pos <- rep(y00revcumsum,y00table)
    y00loc <- n0+1-y00pos
    y00forward <- rep(as.vector(cumsum(y00table)),y00table)

    tab <- matrix(rep(NA,n0*p),ncol=p)
    for(i in 1:p){
      temp   <- rev(cumsum(rev(score0[,i])))
      tab[,i]<- temp[y00pos]
    }

    kappa2  <- tab*(1-d00)/y00loc
    kappa3a <- kappa2/y00loc
    kappa3b <- matrix(rep(NA,p*n0),ncol=p)
    for(i in 1:p){
      kappa3b[,i] <- cumsum(kappa3a[,i])
    }
    kappa3c    <- kappa3b[y00forward,]
    kappa.arm0 <- kappa1+kappa2-kappa3c


    #----
    if(type=="difference")
    {gamma=t(kappa.arm1)%*%kappa.arm1+t(kappa.arm0)%*%kappa.arm0
    A=t(x)%*%x
    varbeta=solve(A)%*%gamma%*%solve(A)
    }

    if(type=="ratio" || type=="lossratio")
    {gamma=t(kappa.arm1)%*%kappa.arm1+t(kappa.arm0)%*%kappa.arm0
    A=t(x*exp(as.vector(x%*%beta0)))%*%x
    varbeta=solve(A)%*%gamma%*%solve(A)
    }


    if(type=="difference")
    {beta0=beta0
    se0=sqrt(diag(varbeta))
    z0=beta0/se0
    p0=1-pchisq(z0^2, 1)
    cilow=beta0-se0*qnorm(1-alpha/2)
    cihigh=beta0+se0*qnorm(1-alpha/2)
    result=cbind(beta0, se0, z0, p0, cilow, cihigh)
    colnames(result)=c("coef", "se(coef)", "z","p",paste("lower .",round((1-alpha)*100, digits=0), sep=""), paste("upper .",round((1-alpha)*100, digits=0), sep=""))
    }


    if(type=="ratio" || type=="lossratio")
    {beta0=beta0
    se0=sqrt(diag(varbeta))
    z0=beta0/se0
    p0=1-pchisq(z0^2, 1)
    r0=exp(beta0)
    cilow=exp(beta0-se0*qnorm(1-alpha/2))
    cihigh=exp(beta0+se0*qnorm(1-alpha/2))
    result=cbind(beta0, se0, z0, p0, exp(beta0), cilow, cihigh)
    colnames(result)=c("coef", "se(coef)", "z","p","exp(coef)",paste("lower .",round((1-alpha)*100, digits=0), sep=""), paste("upper .",round((1-alpha)*100, digits=0), sep=""))
    }

    if(p==2)
      rownames(result)=c("intercept", "x")

    if(p>2)
      rownames(result)=c("intercept", colnames(x[,-1]))

    return(result)
  }
}

#' @name plot.rmst2
#' @aliases plot.rmst2
#' @title plot.rmst2
#' @description S3 method for class 'rmst2'
#' @param x Results of the unadjusted analyses.
#' @param xlab x label.
#' @param ylab y label.
#' @param col Color for line. Default is red.
#' @param col.RMST Color for areas of RMST. Default is pink.
#' @param col.RMTL Color for areas of RMTL. Default is orange.
#' @param density Density of shading lines, in lines per inch. Default is 80.
#' @param angle Slope of shading lines, given as an angle in degrees (counter-clockwise). Default is 85.
#' @param ... Further arguments ignored in this function.
#' @return returns a plot
#' @export
######################################
# plot.rmst2 (hidden)
######################################
plot.rmst2=function(x, xlab="", ylab="", col="red", col.RMST="pink", col.RMTL="orange", density=80, angle=85, ...){

  if(is.null(x$unadjusted.result)) stop("Please run rmst2 without covariates")

  if(!is.null(x$unadjusted.result)){

    ft1=x$RMST.arm1$fit
    ft0=x$RMST.arm0$fit
    tau=x$tau

    par(mfrow=c(1,2))

    #=== arm 1 ===
    fit=ft1

    tmp.xx=c(0, fit$time); tmp.yy=c(1, fit$surv) ;
    idx=tmp.xx<=tau
    y.tau = min(tmp.yy[idx])
    xx=c(tmp.xx[idx],   tau)
    yy=c(tmp.yy[idx], y.tau)
    x.step=sort(c(0, xx, xx))
    y.step=rev(sort(c(1,1,yy, yy[-length(yy)])))

    #--------
    plot(fit, mark.time=F, conf.int=F, lwd=2, main="arm=1", xlab=xlab, ylab=ylab, col=col, sub=paste("RMST:",round(x$RMST.arm1$rmst[1], digits=2)), ...)

    for (i in seq(1, length(x.step), by=2)){
      polygon(c(x.step[i], x.step[i+1], x.step[i+1], x.step[i]), c(0, 0, y.step[i+1], y.step[i]), col= col.RMST, density=density, angle=angle, lwd=2)
    }
    for (i in seq(1, length(x.step), by=2)){
      polygon(c(x.step[i], x.step[i+1], x.step[i+1], x.step[i]), c(y.step[i], y.step[i+1], 1,1), col= col.RMTL, density=density, angle=angle, lwd=2)
    }

    x.step=sort(c(0, tmp.xx, tmp.xx))
    y.step=rev(sort(c(1,1,tmp.yy, tmp.yy[-length(tmp.yy)])))
    lines(x.step, y.step, col=col, lwd=3)
    # text(5,0.4, paste(round(rmst$rmst[1], digits=2),"years"), cex=0.9)


    #=== arm 0 ===
    fit=ft0

    tmp.xx=c(0, fit$time); tmp.yy=c(1, fit$surv);
    idx=tmp.xx<=tau
    y.tau = min(tmp.yy[idx])
    xx=c(tmp.xx[idx],   tau)
    yy=c(tmp.yy[idx], y.tau)
    x.step=sort(c(0, xx, xx))
    y.step=rev(sort(c(1,1,yy, yy[-length(yy)])))

    #--------
    plot(fit, mark.time=F, conf.int=F, lwd=2, main="arm=0", xlab=xlab, ylab=ylab, col=col, sub=paste("RMST:",round(x$RMST.arm0$rmst[1], digits=2)), ...)
    for (i in seq(1, length(x.step), by=2)){
      polygon(c(x.step[i], x.step[i+1], x.step[i+1], x.step[i]), c(0, 0, y.step[i+1], y.step[i]), col= col.RMST, density=density, angle=angle, lwd=2)
    }
    for (i in seq(1, length(x.step), by=2)){
      polygon(c(x.step[i], x.step[i+1], x.step[i+1], x.step[i]), c(y.step[i], y.step[i+1], 1,1), col= col.RMTL, density=density, angle=angle, lwd=2)
    }

    x.step=sort(c(0, tmp.xx, tmp.xx))
    y.step=rev(sort(c(1,1,tmp.yy, tmp.yy[-length(tmp.yy)])))
    lines(x.step, y.step, col=col, lwd=3)
    # text(5,0.4, paste(round(rmst$rmst[1], digits=2),"years"), cex=0.9)

  }


}
NULL

#' @name print.rmst2
#' @aliases print.rmst2
#' @title print.rmst2
#' @description S3 method for class 'rmst2'
#' @param x Object to be printed.
#' @param digits Integer indicating the number of decimal places.
#' @param ... Further arguments ignored in this function.
#' @return returns summary output for class 'rmst2'
#' @export
######################################
# print.rmst2 (hidden)
######################################
print.rmst2=function(x, digits=3, ...){

  cat("\n")

  cat(x$note,"\n\n")



  #--- unadjusted analysis --
  if(!is.null(x$unadjusted.result)){

    RMST=rbind(x$RMST.arm1$result[1,], x$RMST.arm0$result[1,])
    RMTL=rbind(x$RMST.arm1$result[2,], x$RMST.arm0$result[2,])
    rownames(RMST)=c("RMST (arm=1)","RMST (arm=0)")
    rownames(RMTL)=c("RMTL (arm=1)","RMTL (arm=0)")

    cat ("Restricted Mean Survival Time (RMST) by arm \n")

    prmatrix(round(RMST , digits=digits))

    cat("\n\n")

    cat ("Restricted Mean Time Lost (RMTL) by arm \n")

    prmatrix(round(RMTL, digits=digits))

    cat("\n\n")

    cat ("Between-group contrast \n")

    prmatrix(round(x$unadjusted.result, digits=digits))

  }


  #--- Adjusted analysis --
  if(!is.null(x$adjusted.result)){

    cat ("Summary of between-group contrast (adjusted for the covariates) \n")
    prmatrix(round(x$adjusted.result, digits=digits))

    cat("\n\n")

    cat ("Model summary (difference of RMST) \n")
    prmatrix(round(x$RMST.difference.adjusted, digits=digits))

    cat("\n\n")

    cat ("Model summary (ratio of RMST) \n")
    prmatrix(round(x$RMST.ratio.adjusted, digits=digits))

    cat("\n\n")

    cat ("Model summary (ratio of time-lost) \n")
    prmatrix(round(x$RMTL.ratio.adjusted, digits=digits))


  }

  invisible(x)
}
NULL

#' @name survRM2-package
#' @aliases  survRM2-package
#' @docType  package
#' @title Comparing Restricted Mean Survival Time
#' @description
#' Performs two-sample comparisons using the restricted mean survival time (RMST) as a summary measure of the survival time distribution.
#' Three kinds of between-group contrast metrics (i.e., the difference in RMST, the ratio of RMST and the ratio of the restricted mean time lost (RMTL)) are computed.
#' The package has a function to perform an ANCOVA-type covariate adjustment as well as unadjusted analyses for those measures.
#' @author Hajime Uno, Lu Tian, Miki Horiguchi, Angel Cronin, Chakib Battioui, James Bell
#'
#' Maintainer: Hajime Uno <huno@jimmy.harvard.edu>
#' @references
#' Uno H, Claggett B, Tian L, Inoue E, Gallo P, Miyata T, Schrag D,
#' Takeuchi M, Uyama Y, Zhao L, Skali H, Solomon S, Jacobus S, HughesM,
#' Packer M, Wei LJ. Moving beyond the hazard ratio in quantifying the between-group difference in survival analysis.
#' Journal of clinical Oncology 2014, 32, 2380-2385. doi:10.1200/JCO.2014.55.2208.
#'
#' Tian L, Zhao L,  Wei LJ. Predicting the restricted mean event time with the subject's baseline covariates in survival analysis. Biostatistics 2014, 15, 222-233. doi:10.1093/biostatistics/kxt050.
#' @keywords
#' survival
#' @seealso
#' survival
#' @import survival
#' @importFrom graphics lines par plot polygon
#' @importFrom stats glm lm pchisq pnorm qnorm
#' @importFrom utils data
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
NULL

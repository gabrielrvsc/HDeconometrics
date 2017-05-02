#' Estimates bagging coefficients for a given pre-testing procedure.
#'
#' Estimates the Inuoe and Kilian (2008) bagging for a given pre-testing procedure.
#'
#' @details This function returns the pre-testing coefficients for all bootstrap samples. This coefficients may then be used to calculate forecasts.
#'
#' There are three types of pre-testing:
#' \itemize{
#' \item{joint}{Tests all variables in one shot. This is the pre-testing used by Inuoe and Kilian (2008). It is not indicated when the number of variables is too close to the number of observations and it is unfeasible if the number of variables is bigger than the number of observations.}
#' \item{group-joint}{The pre-testing is performed on random groups of variables. It is feasible even when the number of variables is bigger than the number of observations. The number of groups may be chosen by the user (default=10).}
#' \item{individual}{Tests each variable individually. It is feasible when the number of observations is bigger than the number of variables but it consumes more time than group-joint.}
#' }
#'
#' @param x Matrix of independent variables. Each row is an observation and each column is a variable.
#' @param y Response variable equivalent to the function.
#' @param fn If pre-testing="personal" the user must define the pre-testing function through this argument. This function must return a vector of coefficients where the first coefficient is the intercept. The first argument must be a matrix where the first column is the y and the remaining columns are x.
#' @param R Number of bootstrap replucations.
#' @param l lenght of the blocks for the block-boostrap.
#' @param sim tsboot argument.
#' @param pre.testing The type of pre-testing (see details).
#' @param fixed.controls numeric or character vector indicating variables that must be used as fixed controls in the pre-testing. These variables are always selected.
#' @param ... Other arguments passed to tsboot and to personal pre-testing
#' @return An object with S3 class bagging.
#' \item{coefficients}{Boostrap coefficients on each sample.}
#' \item{orig.coef}{Coefficients on the original sample.}
#' \item{fitted.values}{In-sample fitted values.}
#' \item{residuals}{Model residuals.}
#' \item{pre.testing}{The pre-testing used.}
#' \item{call}{The matched call.}
#' @keywords Bagging
#' @export
#' @importFrom stats lm
#' @importFrom utils combn
#' @examples
#' ## == This example uses the Brazilian inflation data from
#' #Garcia, Medeiros and Vasconcelos (2017) == ##
#' data("BRinf")
#'
#' ## == Data preparation == ##
#' ## == The model is yt = a + Xt-1'b + ut == ##
#' ## == The autorregressive is a fixed control == ##
#' aux = embed(BRinf,2)
#' y=aux[,1]
#' x=aux[,-c(1:ncol(BRinf))]
#'
#' model=bagging(x,y,pre.testing = "group-joint")
#' model$orig.coef
#'
#' ## == check selection frequency == ##
#' coef=coef(model)
#' coef[coef!=0]=1
#' frequency=(colSums(coef))[-1] # remove intercept
#' barplot(frequency)
#'
#' ## == see fitted values == ##
#' plot(y,type="l")
#' lines(fitted(model),col=2)
#'
#' @references
#' Inoue, Atsushi, and Lutz Kilian. "How useful is bagging in forecasting economic time series? A case study of US consumer price inflation." Journal of the American Statistical Association 103.482 (2008): 511-522.
#'
#' Garcia, Medeiros and Vasconcelos (2017).
#'
#' @seealso \code{\link{predict.bagging}}


bagging=function(x,y,fn=NULL,R=100,l=3,sim="fixed",pre.testing=c("group-joint","joint","individual","personal"),fixed.controls=NULL,...){

  pre.testing=match.arg(pre.testing)
  tsd=cbind(y,x)
  if(pre.testing=="personal"){
    bag=boot::tsboot(tseries=tsd,statistic=fn,R=R,l=l
               ,sim=sim, ...)
  }else{
    bag=boot::tsboot(tseries=tsd,statistic=baggit,R=R,l=l
               ,sim=sim, pre.testing=pre.testing, fixed.controls=fixed.controls,... )
  }

  original=bag$t0
  bootres=bag$t
  colnames(bootres)=names(original)

  yhat=rowMeans(cbind(1,x)%*%t(bootres),na.rm=TRUE)
  residuals=y-yhat

  result=list("coefficients"=bootres,"orig.coef"=original,"fitted.values"=yhat,"residuals"=residuals,"pre.testing"=pre.testing,"call"=match.call())
  class(result)="bagging"
  return(result)
}

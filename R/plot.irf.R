#' VAR ipulse response functions
#'
#' Estimates impulse response coefficients of a HDvar or lbvar model h steps ahead. Confidence bands may be computed by bootstrap.
#'
#'
#' @param x A irf object.
#' @param impulse Impulse variable name or index.
#' @param response Response Variable name or index.
#' @param alpha Significance level used only in case the irf was estimated using bootstrap. May be more than one value.
#' @param lty.cb Confidence bands graphic control.
#' @param lwd.cb Confidence bands graphic control.
#' @param ... Other graphical parameters.
#' @export
#' @keywords VAR, irf, High-dimension, Bayesian models
#' @examples
#' ## == This example uses the Brazilian inflation data from
#' #Garcia, Medeiros and Vasconcelos (2017) == ##
#'
#' # = This is an ilustrative example = #
#' # = The identification ignores which variables are more exogenous = #
#' data("BRinf")
#' Y=BRinf[,1:59]# remove expectation variables
#' modelB=lbvar(Y,p=3)
#' identB=identification(modelB)
#' irfB=irf(modelB,identB,h=12,boot = TRUE,M=100)
#' plot(irfB,1,2,alpha=0.1)
#'
#'
#' @references
#' Garcia, Medeiros and Vasconcelos (2017).
#' @seealso \code{\link{predict}}, \code{\link{lbvar}}, \code{\link{identification}}, \code{\link{irf}}


plot.irf=function(x,impulse,response,alpha=0.05,lty.cb=2,lwd.cb=1,...){
  bootirf=x$density[[impulse]][[response]]
  ir=x$point.irf[[impulse]][,response]
  graphics::plot(ir,type="l",...)
  graphics::abline(h=0,col="yellow",lty=2)
  if(!is.null(bootirf)){
    for(i in 1:length(alpha)){
      aux=round(stats::quantile(1:ncol(bootirf),probs=c(alpha[i]/2,1-alpha[i]/2)))
      graphics::lines(bootirf[,aux[1]],col=i+1,lty=lty.cb,lwd=lwd.cb)
      graphics::lines(bootirf[,aux[2]],col=i+1,lty=lty.cb,lwd=lwd.cb)
    }
  }
}

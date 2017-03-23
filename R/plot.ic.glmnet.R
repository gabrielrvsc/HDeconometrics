#' Plot the information criterion curve through the regularization path
#'
#' Plots the information criterion curve as a function of log(lambda) through the regularization path.
#'
#' @param x A ic.glmnet object estimated using the ic.glmnet function.
#' @param ... Arguments to be passed to other methods.
#' @export
#' @examples
#' ## == This example uses the Brazilian inflation data from
#' #Garcia, Medeiros and Vasconcelos (2017) == ##
#' data("BRinf")
#'
#' ## == Data preparation == ##
#' ## == The model is yt = a + Xt-1'b + ut == ##
#' aux = embed(BRinf,2)
#' y=aux[,1]
#' x=aux[,-c(1:ncol(BRinf))]
#'
#' ## == LASSO == ##
#' lasso=ic.glmnet(x,y,crit = "bic")
#' plot(lasso)
#'

plot.ic.glmnet=function(x,...){
  n=x$glmnet$df
  llambda=log(x$glmnet$lambda)
  ic=x$ic.range
  ylab=names(x$ic[which(x$ic==ic[which.min(ic)])])

  graphics::plot(llambda, ic, xlab = "log(Lambda)", ylab = substr(ylab, 1, 3),...)
  graphics::abline(v = llambda[which.min(ic)], lty = 2)
  graphics::axis(3, at = llambda, labels = n)
}

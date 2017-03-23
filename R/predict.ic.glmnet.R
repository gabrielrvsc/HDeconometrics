#' Predict method for ic.glmnet objects
#'
#' Predicted values based on ic.glmnet object.
#'
#' @param object A ic.glmnet object estimated using the ic.glmnet function.
#' @param newdata An optional data to look for the explanatory variables used to predict. If omitted, the fitted values are used.
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
#' ## == Break the data into in-sample and out-of-sample
#' y.in=y[1:100]; y.out=y[-c(1:100)]
#' x.in=x[1:100,]; x.out=x[-c(1:100),]
#'
#' ## == LASSO == ##
#' lasso=ic.glmnet(x.in,y.in,crit = "bic")
#'
#' pred=predict(lasso,newdata=x.out)
#' plot(y.out, type="l")
#' lines(pred, col=2)
#'

predict.ic.glmnet=function(object,newdata,...){
  if(missing(newdata)){
    return(fitted(object))
  }
  coefficients = coef(object)
  if (is.vector(newdata)) {
    pred = c(1, newdata) %*% coefficients
  }
  else {
    pred = cbind(1, as.matrix(newdata)) %*% coefficients
  }
  return(pred)
}

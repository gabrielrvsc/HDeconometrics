#' Predict method for bagging objects
#'
#' Predicted values based on bagging object.
#'
#' @param object A bagging object estimated using the ic.glmnet function.
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
#' ## == break data (in-sample and out-of-sample)
#' yin=y[1:120]
#' yout=y[-c(1:120)]
#' xin=x[1:120,]
#' xout=x[-c(1:120),]
#'
#' model=bagging(xin,yin,pre.testing = "group-joint")
#' pred=predict(model,xout)
#'
#' plot(yout,type="l")
#' lines(pred,col=2)
#'

predict.bagging=function(object,newdata=NULL,...){
  if(is.null(newdata)){
    return(fitted(object))
  }

  parameters = coef(object)
  parameters[is.na(parameters)]=0
  if (is.vector(newdata)) {
    individual.prediction = c(1, newdata) %*% t(parameters)
    final.prediction = mean(individual.prediction)
  }
  else {
    individual.prediction = cbind(1, newdata) %*% t(parameters)
    final.prediction = rowMeans(individual.prediction)
  }
  return(final.prediction)
}

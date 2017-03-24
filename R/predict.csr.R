#' Predict method for csr objects
#'
#' Predicted values based on csr object.
#'
#' @param object A csr object estimated using the ic.glmnet function.
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
#' ## == The autorregressive is a fixed control == ##
#' aux = embed(BRinf,2)
#' y=aux[,1]
#' x=aux[,-c(1:ncol(BRinf))]
#'
#' model=csr(x,y,K=20,k=4,fixed.controls = 1)
#' ## == Break the data into in-sample and out-of-sample
#' y.in=y[1:100]; y.out=y[-c(1:100)]
#' x.in=x[1:100,]; x.out=x[-c(1:100),]
#'
#' model=csr(x.in,y.in,K=20,k=4,fixed.controls = 1)
#' plot(y.out,type="l")
#' lines(predict(model,newdata = x.out),col=2)
#'

predict.csr=function(object,newdata,...){
  if(missing(newdata)){
    return(fitted(object))
  }
  coefficients = coef(object)
  if (is.vector(newdata)) {
    individual.prediction = c(1, newdata) %*% t(coefficients)
    pred = mean(individual.prediction)
  }
  else {
    individual.prediction = cbind(1, newdata) %*% t(coefficients)
    pred = rowMeans(individual.prediction)
  }
  return(pred)
}


#' @export

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

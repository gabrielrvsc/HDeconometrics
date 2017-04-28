#' @export

predict.boosting=function(object,newdata=NULL,...){
  if(is.null(newdata)){
    return(fitted(object))
  }

  parameters = coef(object)
  ybar = mean(object$y)
  final.prediction = ybar + newdata %*% parameters
}


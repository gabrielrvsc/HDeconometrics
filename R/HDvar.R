#' Estimate High-dimensional VARs
#'
#' Estimate High-dimensional Vector Autorregressive models using the model chosen by the user. The model must be parametric.
#'
#'
#' @param Y Time-series matrix or data.frame with the VAR endogenous variables.
#' @param p Lag order (default = 1).
#' @param fn A function that receives x (independent variables) and y (dependent variable) and returns the model coefficients (default is OLS).
#' @param xreg Exogenous controls.
#' @return An object with S3 class "HDeconometricsVAR", "HDvar".
#' \item{coef.by.equation}{Coefficients listed by each VAR equation.}
#' \item{coef.by.block}{Coefficients separated by blocks (intercepts, lags, exogenous).}
#' \item{fitted.values}{In-sample fitted values.}
#' \item{residuals}{The residuals.}
#' \item{Y}{Supplied endogenous data.}
#' \item{p}{VAR lag order chosen by the user.}
#' \item{N}{Number of endogenous variables.}
#' \item{covmat}{Residuals covariance matrix.}
#' \item{xreg}{Exogenous controls supplied by the user.}
#' \item{T}{Number of observations.}
#' \item{fn}{The estimation function.}
#' \item{call}{The matched call.}
#' @keywords VAR, LASSO, High-dimension, Bayesian models
#' @export
#' @examples
#' ## == This example uses the Brazilian inflation data from
#' #Garcia, Medeiros and Vasconcelos (2017) == ##
#' data("BRinf")
#' Y=BRinf[,1:59]# remove expectation variables
#'
## Estimate by LASSO ##
#' fn=function(x,y){
#'   model=ic.glmnet(x,y)
#'   coef(model)
#' }
#'
#' modelH=HDvar(Y,p=4,fn=fn)
#'
#' # take a look at the coefficients
#' eq=coef(modelH,type="equation")
#' block=coef(modelH,type="block")
#' block$Lag1
#'
#'
#' @references
#' Garcia, Medeiros and Vasconcelos (2017).
#' @seealso \code{\link{predict}}, \code{\link{lbvar}}, \code{\link{irf}}


HDvar=function(Y,p,fn=NULL,xreg=NULL){

  if(is.null(fn)){
    fn=function(x,y){
      model=lm(y~x)
      return(coef(model))
    }
  }


  if (!is.matrix(Y)) {
    Y = as.matrix(Y)
  }

  if (length(xreg) != 0) {
    if (!is.matrix(xreg)) {
      xreg = as.matrix(xreg)
    }
    if(is.null(colnames(xreg))){
      colnames(xreg)=paste("xreg",1:ncol(xreg),sep="")
    }
  }
  N = ncol(Y)

  if(is.null(colnames(Y))){
    colnames(Y)=paste("V",1:N,sep="")
  }

  aux = stats::embed(Y, p + 1)
  Yreg = aux[, 1:N]
  Xreg = aux[, -c(1:N)]

  if(length(xreg)!=0){
    xreg=as.matrix(tail(xreg,nrow(Yreg)))
    Xreg=cbind(Xreg,xreg)
  }
  coef.by.equation=t(apply(Yreg,2,fn,x=Xreg))
  fitted=cbind(1,Xreg)%*%t(coef.by.equation)
  residuals=tail(Y,nrow(fitted))-fitted

  ## == variable names == ##
  rownames(coef.by.equation)=colnames(Y)
  namvar=rep(colnames(Y),p)
  namlag=sort(rep(1:p,ncol(Y)))
  cnams=c("cons.", paste(namvar," lag", namlag,sep=""))
  if(length(xreg)!=0){
    cnams=c(cnams,colnames(xreg))
  }
  colnames(coef.by.equation) = cnams

  ## == coef in VAR blocks == ##
  auxcoef=coef.by.equation
  coef.by.block = list(intersect = auxcoef[, 1])
  auxcoef = auxcoef[, -1]
  for (i in 1:p) {
    coef.by.block[[i + 1]] = auxcoef[, (ncol(Y) *
                                          i - ncol(Y) + 1):(ncol(Y) * i)]
    colnames(coef.by.block[[i+1]])=colnames(coef.by.block[[i+1]])=colnames(Y)
  }
  if (length(xreg) != 0) {
    aux = ncol(auxcoef)
    coef.by.block[[length(coef.by.block) + 1]] =as.matrix( auxcoef[,
                                                         (aux - ncol(xreg) + 1):aux])
    colnames(coef.by.block[[length(coef.by.block)]])=colnames(xreg)
    rownames(coef.by.block[[length(coef.by.block)]])=colnames(Y)
  }

  # = names blocks = #
  namblock=c("cons",paste("Lag",1:p,sep=""))
  if(length(xreg)!=0){
    namblock=c(namblock,"xreg")
  }
  names(coef.by.block)=namblock

  colnames(fitted)=colnames(residuals)
  rownames(fitted)=rownames(residuals)
  sigmae=(t(residuals)%*%residuals)/(nrow(Yreg)-ncol(Yreg))
  colnames(sigmae)=rownames(sigmae)=colnames(Y)

  result=list(coef.by.equation = coef.by.equation, coef.by.block = coef.by.block,
              fitted.values = fitted , residuals=residuals, Y = Y, p = p, N=N, covmat = sigmae,
              xreg = xreg, T = nrow(Xreg),fn=fn ,call=match.call())
  class(result)=c("HDeconometricsVAR","HDvar")
  return(result)
}

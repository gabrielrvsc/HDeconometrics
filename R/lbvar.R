#' Estimate Large Bayesian VARs
#'
#' Estimate large Bayesian Vector Autorregressive models from Banbura et al. (2010)
#'
#'
#' @param Y Time-series matrix or data.frame with the VAR endogenous variables.
#' @param p Lag order (default = 1).
#' @param delta Numeric vector indicating the prior for the autorregressive coefficients (default = 0 for all variables). If the prior is the same for all variables the user may supply a single number. Otherwise the vector must have one element for each variable.
#' @param lambda Constant that regulates the importance given to the priors (default = 0.05). If lambda = 0 the model ignores the data and the posterior equal the prior. For bigger lambda the model converges to the OLS  estimates.
#' @param xreg Exogenous controls.
#' @param ps If TRUE the priors on the sum of the coefficients will be included.
#' @param tau Controls the shrinkage in the priors on the sum of the coefficients.
#' @return An object with S3 class "HDeconometricsVAR", "lbvar".
#' \item{coef.by.equation}{Coefficients listed by each VAR equation.}
#' \item{coef.by.block}{Coefficients separated by blocks (intercepts, lags, exogenous).}
#' \item{fitted.values}{In-sample fitted values.}
#' \item{residuals}{The residuals.}
#' \item{Y}{Supplied endogenous data.}
#' \item{p}{VAR lag order chosen by the user.}
#' \item{N}{Number of endogenous variables.}
#' \item{covmat}{Residuals covariance matrix.}
#' \item{xreg}{Exogenous controls supplied by the user.}
#' \item{Ts}{Number of real observations, number of dummy observations and the sum of both.}
#' \item{delta}{The delta chosen.}
#' \item{lambda}{The lambda chosen.}
#' \item{call}{The matched call.}
#' @keywords VAR, BVAR, High-dimension, Bayesian models
#' @export
#' @examples
#' ## == This example uses the Brazilian inflation data from
#' #Garcia, Medeiros and Vasconcelos (2017) == ##
#' data("BRinf")
#' Y=BRinf[,1:59]# remove expectation variables
#' modelB=lbvar(Y,p=4)
#'
#' # take a look at the coefficients
#' eq=coef(modelB,type="equation")
#' block=coef(modelB,type="block")
#' block$Lag1
#'
#' @references
#' Banbura, M., Giannone, D., & Reichlin, L. (2010). Large Bayesian vector autoregressions. Journal of Applied Econometrics, 25, 71–92.
#'
#' Garcia, Medeiros and Vasconcelos (2017).
#' @seealso \code{\link{predict}}, \code{\link{HDvar}}, \code{\link{irf}}, \code{\link{fitLambda}}


lbvar=function (Y, p = 1, delta = 0, lambda = 0.05, xreg = NULL,ps=FALSE,tau=10*lambda)
{
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
  Sig = stats::cov(Yreg)
  sig = sqrt(diag(Sig))

  ## == Create Dummy Variables == ##
  aux1 = delta * diag(sig, N)/lambda
  aux2 = matrix(0, N * (p - 1), N)
  aux3 = diag(sig, N)
  aux4 = rep(0, N)
  Yd = rbind(aux1, aux2, aux3, aux4)
  aux1 = diag(1:p, p)
  aux2 = diag(sig, N)/lambda
  aux3 = kronecker(aux1, aux2)
  aux4 = matrix(0, N, N * p)
  aux5 = rep(0, N * p)
  aux6 = rbind(aux3, aux4, aux5)
  if (length(xreg) > 0) {
    aux7 = matrix(0, nrow(aux6), ncol(xreg))
    aux6 = cbind(aux6, aux7)
    Xreg = cbind(Xreg, tail(xreg, nrow(Xreg)))
  }
  Xd = cbind(c(rep(0, nrow(aux6) - 1), 0.1), aux6)

  if(ps==TRUE){
    YDps=diag(colMeans(Yreg))/tau
    XDps=cbind(0,kronecker(t(1:p),YDps))
    if(length(xreg)>0){
      aux=matrix(0,nrow(XDps),ncol(Xd)-ncol(XDps))
      XDps=cbind(XDps,aux)
    }

    Xd=rbind(Xd,XDps)
    Yd=rbind(Yd,YDps)
  }

  # == Variable Star == #
  Ystar = rbind(Yd, Yreg)
  Xstar = rbind(Xd, cbind(1, Xreg))

  # == betas and fitted == #
  betas=stats::coef(stats::lm(Ystar~-1+Xstar))
  fitted = cbind(1, Xreg) %*% betas
  sigmae = (1/nrow(Ystar)) * t(Ystar - cbind(Xstar) %*% betas) %*%
    (Ystar - cbind(Xstar) %*% betas)
  coef.by.equation = t(betas)

  ## == variable names == ##
  rownames(coef.by.equation)=colnames(Y)
  namvar=rep(colnames(Y),p)
  namlag=sort(rep(1:p,ncol(Y)))
  cnams=c("cons.", paste(namvar," lag", namlag,sep=""))
  if(length(xreg)!=0){
    cnams=c(cnams,colnames(xreg))
  }
  colnames(coef.by.equation) = cnams
  colnames(sigmae)=rownames(sigmae)=colnames(Y)

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

  residuals=tail(Y,nrow(fitted))-fitted
  colnames(fitted)=colnames(residuals)
  rownames(fitted)=rownames(residuals)

  result=list(coef.by.equation = coef.by.equation, coef.by.block = coef.by.block,
              fitted.values = fitted , residuals=residuals, Y = Y, p = p, N=N, covmat = sigmae,
              xreg = xreg, Ts = c(T = nrow(Xreg), Td = nrow(Xd), Tstar = nrow(Xstar)),delta=delta,lambda=lambda,ps=ps,tau=tau,call=match.call())
  class(result)=c("HDeconometricsVAR","lbvar")

  return(result)
}

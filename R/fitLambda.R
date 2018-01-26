#' Fit the parameter lambda of the Bayesian VAR
#'
#' Fit the parameter lambda of the Bayesian VAR. This parameter controls the importance given to the priors. If lambda=0 the model is the same as the OLS case. For bigger values of lambda more importance is given to the priors and less importance to the data.
#'
#' @details The choice of lambda is arbitrary. However, Banbura et al. (2010) uses the fit of a smaller model with just a few variables as a target to the Bayesian VAR. In other words, this function chooses the lambda that matches the fit of a smaller model chosen by the user on the chosen Bayesian VAR. If lambda = 0 the model ignores the data and the posterior equal the prior. For bigger lambda the model converges to the OLS  estimates.
#'
#' @param Y Time-series matrix or data.frame with the VAR endogenous variables.
#' @param variables Either a numeric vector indication the position of the variables to be included in the small model or a characted vector with the variable names.
#' @param lambdaseq Sequence of lambdas to be tested.
#' @param p Lag order (default = 1).
#' @param p.reduced Lag order of the small model.
#' @param delta Numeric vector indicating the prior for the autorregressive coefficients (default = 0 for all variables). If the prior is the same for all variables the user may supply a single number. Otherwise the vector must have one element for each variable.
#' @param xreg Exogenous controls.
#' @param scale If TRUE the variables are centered with variance equal 1 (default is TRUE).
#' @param ps If TRUE the priors on the sum of the coefficients will be included.
#' @param tau Controls the shrinkage in the priors on the sum of the coefficients.
#' @keywords VAR, BVAR, High-dimension, Bayesian models
#' @export
#' @examples
#' ## == This example uses the Brazilian inflation data from
#' #Garcia, Medeiros and Vasconcelos (2017) == ##
#' data("BRinf")
#' Y=BRinf[,1:59]# remove expectation variables
#' lambda=fitLambda(Y,variables=c(1,4,10),
#'                 lambdaseq = seq(0,0.1,0.005),
#'                 p=2,p.reduced = 2)
#'model=lbvar(Y,p=2,lambda=lambda)
#'
#' @references
#' Banbura, M., Giannone, D., & Reichlin, L. (2010). Large Bayesian vector autoregressions. Journal of Applied Econometrics, 25, 71–92.
#'
#' Garcia, Medeiros and Vasconcelos (2017).
#' @seealso \code{\link{lbvar}}


fitLambda=function(Y,variables=1,lambdaseq=seq(0,1,0.1),p=1,p.reduced=p,delta=0,xreg=NULL,scale=TRUE,ps=FALSE,tau=10*lambdaseq){
  if(scale==TRUE){
    Y=scale(Y)
  }
  Yreduced=Y[,variables]
  model=HDvar(Yreduced,p.reduced)
  resid=stats::residuals(model)
  rmse=sum(sqrt(colMeans(resid^2)))

  save.rmseV=rep(0,length(lambdaseq))
  for(i in 1:length(lambdaseq)){
    lbv=lbvar(Y,p,lambda = lambdaseq[i],delta=delta,xreg=xreg,ps=ps,tau=tau[i])
    residV=stats::residuals(lbv)[,variables]
    rmseV=sum(sqrt(colMeans(residV^2)))
    save.rmseV[i]=rmseV
  }
  aux=abs(save.rmseV-rmse)
  graphics::plot(lambdaseq,aux,type="l",xlab="lambda",ylab="RMSE abs difference")
  return(lambdaseq[which.min(aux)])
}

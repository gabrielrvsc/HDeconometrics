#' Estimate the Complete Subset Regressions (CSR) with pre-testing
#'
#' Estimate the Complete Subset Regressions Elliott, Gargano and Timmermann (2013) with the possibility of pre-testing the variables as in Garcia, Medeiros and Vasconcelos (2017) and Medeiros and Vasconcelos(2016).
#'
#' @details The Complete Subset Regressions estimates, for a given set of K variables, all possible models with k variables. The number of regressions to be fitted grow very fast with K, therefore, in some cases a pre-testing to select only the most relevant variables may be the only way to make the model computationally feasible.
#'
#' The pre-testing is activated automatically if the total number of variables in x is bigger than K. If the user chooses K=ncol(x) the procedure is exactly the one presented in Elliott, Gargano and Timmermann (2013).
#'
#' If the user chooses to include fixed controls in the model they will included in all models and each model will have k+length(fixed.controls) variables.
#'
#' @param x Matrix of independent variables. Each row is an observation and each column is a variable.
#' @param y Response variable equivalent to the function.
#' @param K Number of variables to be selected after the pre-testing. If K=ncol(x) the pre-testing is redundant.
#' @param k Number of variables in each subset. Must be smaller than K.
#' @param fixed.controls A vector indicatin which variables in x are fixed controls. May be a character vector with the variable names or a numeric vector with the variables position in x.
#' @return An object with S3 class csr.
#' \item{coefficients}{Coefficients on each model of the subset.}
#' \item{fitted.values}{In-sample fitted values.}
#' \item{residuals}{Model residuals.}
#' \item{call}{The matched call.}
#' @keywords Complete subset regression
#' @export
#' @importFrom stats lm
#' @importFrom utils combn
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
#' plot(y,type="l")
#' lines(fitted(model),col=2)
#'
#' @references
#' Elliott, Graham, Antonio Gargano, and Allan Timmermann. "Complete subset regressions." Journal of Econometrics 177.2 (2013): 357-373.
#'
#' Garcia, Medeiros and Vasconcelos (2017).
#'
#' Medeiros, Marcelo C., and Gabriel FR Vasconcelos. "Forecasting macroeconomic variables in data-rich environments." Economics Letters 138 (2016): 50-52.
#'
#' @seealso \code{\link{predict.csr}}


csr=function (x, y, K = min(20,ncol(x)), k = 4, fixed.controls = NULL)
{

  if(!is.matrix(x)){
    x=as.matrix(x)
  }
  if(ncol(x)<2){
    stop("Only one variable in x, csr is senseless in this case.")
  }

  if (length(fixed.controls) != 0) {
    if(is.character(fixed.controls)) fixed.controls=match(fixed.controls,colnames(x))
    w = x[, fixed.controls]
    nonw = setdiff(1:ncol(x), fixed.controls)
  }else {
    w = rep(0, nrow(x))
    nonw = 1:ncol(x)
  }
  save.stat = matrix(NA, ncol(x), 3)
  save.stat[fixed.controls, 3] = fixed.controls
  save.stat[fixed.controls, -3] = 0
  for (i in nonw) {
    if(length(fixed.controls)!=0){
      save.stat[i, ] = c(abs(summary(lm(y ~ x[, i] + w))$coefficients[2,c(3, 4)]), i)
    }else{
      save.stat[i, ] = c(abs(summary(lm(y ~ x[, i]))$coefficients[2,c(3, 4)]), i)
    }
  }
  t.ord = order(save.stat[, 1], decreasing = TRUE)
  save.stat = save.stat[t.ord, ]
  aux = setdiff(1:ncol(x), fixed.controls)
  if (K > length(aux)) {
    stop("K bigger than the number of possible variables. Choose a value for K that is equal \n         smaller than ncol(X)-length(fixed.controls)")
  }
  selected = save.stat[1:K, 3]
  aux = combn(selected, k)
  m = ncol(aux)
  final.coef = matrix(0, m, ncol(x))
  final.const = rep(0, m)
  if (length(fixed.controls) != 0) {
    for (i in 1:m) {
      model = coef(lm(y ~ w + x[, aux[, i]]))
      final.const[i] = model[1]
      model = model[-1]
      final.coef[i, fixed.controls] = model[1:length(fixed.controls)]
      final.coef[i, aux[, i]] = model[-c(1:length(fixed.controls))]
    }
  }else {
    for (i in 1:m) {
      model = coef(lm(y ~ x[, aux[, i]]))
      final.const[i] = model[1]
      final.coef[i, aux[, i]] = model[-1]
    }
  }
  colnames(final.coef) = colnames(x)
  final.coef[is.na(final.coef)]=0
  final.coef = cbind(intersect = final.const, final.coef)

  fitted.values.i = cbind(1, x) %*% t(final.coef)
  fitted.values=rowMeans(fitted.values.i)
  residuals=y-fitted.values

  result=list(coefficients = final.coef,fitted.values=fitted.values,residuals=residuals,call=match.call())
  class(result)="csr"

  return(result)
}

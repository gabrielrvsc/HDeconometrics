#' Recursive identification of a VAR model
#'
#' Estimates the contemporaneous coefficients and the structural shocks covariance matrix using recursive identification.
#'
#' @param object A HDvar or lbvar object.
#' @return A list with the following objects:
#' \item{A}{Matrix of contemporaneous coefficients.}
#' \item{sigma2u}{Structural shocks covariance matrix.}
#' @export
#' @examples
#' ## == This example uses the Brazilian inflation data from
#' #Garcia, Medeiros and Vasconcelos (2017) == ##
#'
#' # = This is an ilustrative example = #
#' # = The identification ignores which variables are more exogenous = #
#'
#' data("BRinf")
#' Y=BRinf[,1:59]# remove expectation variables
#' modelB=lbvar(Y,p=4)
#' identB=identification(modelB)
#'
#' @references
#' Garcia, Medeiros and Vasconcelos (2017).
identification=function(object){
  covmat_e=object$covmat
  choles=chol(covmat_e)
  covmat_u_sqrt=diag(diag(choles))

  A=solve(t(solve(covmat_u_sqrt)%*%choles))
  sigma2u=covmat_u_sqrt^2
  colnames(A)=colnames(sigma2u)=rownames(sigma2u)=rownames(A)=rownames(covmat_e)
  A[upper.tri(A)]=0
  return(list("A"=A,"covmatu"=sigma2u))
}

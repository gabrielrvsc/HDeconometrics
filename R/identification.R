#' Recursive identification of a VAR model
#'
#' Estimates the contemporaneous coefficients and the structural shocks covariance matrix using recursive identification or using a identity matrix.
#'
#' @param object A HDvar or lbvar object.
#' @param type "recursive" or "identity". The last is just using a identity matrix as contemporaneous shocks.
#' @return A list with the following objects:
#' \item{A}{Matrix of contemporaneous coefficients.}
#' \item{sigma2u}{Structural shocks covariance matrix.}
#' \item{type}{Matched type.}
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
identification=function(object,type=c("recursive","identity")){
  type=match.arg(type)
  if(type=="recursive"){
    covmat_e=object$covmat
    choles=chol(covmat_e)
    covmat_u_sqrt=diag(diag(choles))

    A=solve(t(solve(covmat_u_sqrt)%*%choles))
    sigma2u=covmat_u_sqrt^2
    colnames(A)=colnames(sigma2u)=rownames(sigma2u)=rownames(A)=rownames(covmat_e)
    A[upper.tri(A)]=0
  }
  if(type=="identity"){
    covmat_e=object$covmat
    sigma2u=diag(diag(covmat_e))
    A=diag(nrow(sigma2u))
    colnames(A)=colnames(sigma2u)=rownames(sigma2u)=rownames(A)=rownames(covmat_e)
  }

  return(list("A"=A,"covmatu"=sigma2u,type=type))
}

#' Estimate component-wise boosting for dynamic models.
#'
#' Estimates the Bai & Ng (2009) component-wise boosting for dynamic models.
#'
#' @details This is an implementation of time-series component-wise boosting using the results from Bai and Ng (2009). The algorithm has its own way of determining when to stop. Keep ic.break=TRUE if you want to use the standard stop criterion based on information criterion.
#'
#' Note that the information criterion automaticaly adjusts the degrees of freedom of the model considering that the boosting may select the same variable more than once.
#'
#' @param x Matrix of independent variables. Each row is an observation and each column is a variable.
#' @param y Response variable equivalent to the function.
#' @param v Algorithm step size.
#' @param minIt Minimum number of iterations in case ic.break=TRUE.
#' @param maxIt Maximum number of iterations.
#' @param ic.break If TRUE, algorithm breaks when the minimum information criteria is likely to be found. If FALSE the algorithm stops only when maxIt is reached (default=TRUE).
#' @return An object with S3 class boosting.
#' \item{coefficients}{Boosting coefficients for the model with the smallest information criteria.}
#' \item{fitted.values}{In-sample fitted values.}
#' \item{residuals}{Model residuals.}
#' \item{best.crit}{The smalles information criterion found.}
#' \item{crit}{The information criterion in each iteration.}
#' \item{df}{Degrees of freedom.}
#' \item{coef.selection.count}{How many times each variable was selected.}
#' \item{y}{The supplied y.}
#' \item{call}{The matched call.}
#' @keywords Boosting, Dynamic models
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
#' aux = embed(BRinf,2)
#' y=aux[,1]
#' x=aux[,-c(1:ncol(BRinf))]
#'
#' ## == Use factors == ##
#' factors=prcomp(x,scale. = TRUE)
#' xfact=factors$x[,1:10]
#'
#' model=boosting(xfact,y)
#' coef(model)
#'
#' plot(y,type="l")
#' lines(fitted(model),col=2)
#'
#' @references
#' Bai, Jushan, and Serena Ng. "Boosting diffusion indices." Journal of Applied Econometrics 24.4 (2009): 607-629.
#'
#' Garcia, Medeiros and Vasconcelos (2017).
#'
#' @seealso \code{\link{predict.boosting}}


boosting=function(x,y,v=0.2,minIt=ncol(x)/2,maxIt=10*ncol(x),ic.break=TRUE){

  phi=rep(mean(y),length(y))
  B=rep(0,ncol(x))
  BB=(1/nrow(x))*matrix(1,nrow(x),nrow(x))
  df.final=Inf
  M=maxIt
  save.B=matrix(NA,ncol(x),maxIt)
  save.crit=rep(Inf,M)
  coef.selection.count=rep(0,ncol(x))

  for(m in 1:M){

    u=y-phi
    b=rep(NA,ncol(x))
    e=matrix(NA,nrow(x),ncol(x))
    ssr=rep(NA,ncol(x))


    iter=apply(x,2,auxboost,u=u)
    b=iter[1,]
    ssr=iter[2,]


    best.ssr=which(ssr==min(ssr))
    best.coef=b[best.ssr]
    coef.selection.count[best.ssr]=coef.selection.count[best.ssr]+1
    phi=phi+v*x[,best.ssr]*best.coef
    aux=rep(0,ncol(x))
    aux[best.ssr]=best.coef
    B=B+v*aux

    BB=BB+v*x[,best.ssr]%*%solve(t(x[,best.ssr])
                                 %*%x[,best.ssr])%*%t(x[,best.ssr])%*%(diag(1,nrow(x))-BB)
    df=sum(diag(BB))
    crit=log(sum((y-phi)^2)) + (log(nrow(x))*df)/nrow(x)

    if(crit<min(save.crit,na.rm = TRUE)){
      df.final=df
    }

    save.crit[m]=crit

    save.B[,m]=B
    if(ic.break==TRUE){
      if(m>=minIt){
        if(min(save.crit[1:(m/2)])<crit){
          break
        }
      }
    }

  }

  coef.final=save.B[,which(save.crit==min(save.crit,na.rm=TRUE))]
  names(coef.final)=colnames(x)
  fitted=mean(y)+x%*%coef.final
  save.crit=save.crit[1:m]

  residuals=y-fitted

  result=list("coefficients"=coef.final,"fitted.values"=fitted,"residuals"=residuals,"best.crit"=min(save.crit,na.rm=TRUE),"crit"=save.crit,"df"=df.final,"y"=y,"coef.selection.count"=coef.selection.count,"call"=match.call())
  class(result)="boosting"
  return(result)
}

#' @export

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

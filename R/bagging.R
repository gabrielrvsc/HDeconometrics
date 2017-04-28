#' @export

bagging=function(x,y,fn=NULL,R=100,l=3,sim="fixed",pre.testing=c("group-joint","joint","individual","personal"),fixed.controls=NULL,...){

  pre.testing=match.arg(pre.testing)
  tsd=cbind(y,x)
  if(pre.testing=="personal"){
    bag=boot::tsboot(tseries=tsd,statistic=fn,R=R,l=l
               ,sim=sim, ...)
  }else{
    bag=boot::tsboot(tseries=tsd,statistic=baggit,R=R,l=l
               ,sim=sim, pre.testing=pre.testing, fixed.controls=fixed.controls,... )
  }

  original=bag$t0
  bootres=bag$t
  colnames(bootres)=names(original)

  yhat=rowMeans(cbind(1,x)%*%t(bootres),na.rm=TRUE)
  residuals=y-yhat

  result=list("coefficients"=bootres,"orig.coef"=original,"fitted.values"=yhat,"residuals"=residuals,"pre.testing"=pre.testing,"call"=match.call())
  class(result)="bagging"
  return(result)
}

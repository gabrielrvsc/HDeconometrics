
ic.glmnet = function (x, y, crit=c("bic","aic","aicc","hqc"),...)
{
  if (is.matrix(x) == FALSE) {
    x = as.matrix(x)
  }
  if (is.vector(y) == FALSE) {
    y = as.vector(y)
  }
  crit=match.arg(crit)
  n=length(y)
  model = glmnet(x = x, y = y, ...)
  coef = coef(model)
  lambda = model$lambda
  df = model$df

  yhat=cbind(1,x)%*%coef
  residuals = (y- yhat)
  mse = colMeans(residuals^2)
  sse = colSums(residuals^2)

  nvar = df + 1
  bic = n*log(mse)+nvar*log(n)
  aic = n*log(mse)+2*nvar
  aicc = aic+(2*nvar*(nvar+1))/(n-nvar-1)
  hqc = n*log(mse)+2*nvar*log(log(n))

  sst = (n-1)*var(y)
  r2 = 1 - (sse/sst)
  adjr2 = (1 - (1 - r2) * (n - 1)/(nrow(x) - nvar - 1))

  crit=switch(crit,bic=bic,aic=aic,aicc=aicc,hqc=hqc)

  selected=best.model = which(crit == min(crit))

  ic=c(bic=bic[selected],aic=aic[selected],aicc=aicc[selected],hqc=hqc[selected])

  result=list(coef=coef[,selected],ic=ic,lambda = lambda[selected], nvar=nvar[selected],
              glmnet=model,residuals=residuals[,selected],fitted=yhat[,selected],call = match.call())

  class(result)="ic.glmnet"
  return(result)
}



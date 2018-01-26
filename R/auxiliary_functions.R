auxboost=function(x,u){
  b=solve(t(x)%*%x)%*%(t(x)%*%u)
  mod=lm(u~-1+x)
  b=coef(mod)
  e=u-x%*%matrix(b)
  ssr=t(e)%*%e
  return(c("b"=b,"ssr"=ssr))
}


irfaux=function(object,ident,h,unity.shock=TRUE){
  p=object$p
  k=ncol(object$fitted)
  J=matrix(0,k,k*p)
  diag(J[1:k,1:k])=1

  aux=diag(k)
  if(unity.shock==FALSE){
    aux=aux%*%sqrt(ident$covmatu)
  }
  phi0=solve(ident$A)%*%aux

  A=unlist(object$coef.by.block[2:(p+1)])

  companion <- matrix(0, nrow = k * p, ncol = k * p)
  companion[1:k, 1:(k * p)] = A
  if (p > 1) {
    j <- 0
    for (i in (k + 1):(k * p)) {
      j <- j + 1
      companion[i, j] <- 1
    }
  }
  store.phi=list()
  store.phi[[1]]=diag(k)
  aux=companion
  for(i in 2:(h+1)){
    store.phi[[i]]=J%*%aux%*%t(J)
    aux=companion%*%aux
  }
  irmat=lapply(store.phi,function(x) x %*% phi0)

  aux=matrix(NA,h+1,k)
  colnames(aux)=colnames(object$Y)
  ir=list()
  for(i in 1:k){
    ir[[i]]=aux
  }

  for(i in 1:k){
    for(j in 1:length(irmat)){
      ir[[i]][j,]=irmat[[j]][,i]
    }
  }
  names(ir)=colnames(object$Y)
  return(ir)
}



##############

HDvarboot=function(Y,p,fn,xreg=NULL,databoot){
  N = ncol(Y)

  Yreg=databoot[,1:N]
  Xreg=databoot[,-c(1:N)]

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
              fitted.values = fitted , residuals=residuals, Y = Y, p = p, N=N, covmat = sigmae, type = "var",
              xreg = xreg, T = nrow(Xreg),fn=fn ,call=match.call())
  class(result)=c("HDeconometricsVAR","HDvar")
  return(result)
}



lbvarboot=function (Y, p = 1, delta = 0, lambda = 0.05, xreg = NULL,databoot, ps=TRUE, tau=10*lambda)
{

  N = ncol(Y)
  Yreg=databoot[,1:N]
  Xreg=databoot[,-c(1:N)]
  if(length(xreg)!=0){
    xregnam=colnames(xreg)
    xreg=as.matrix(xreg)
    aux=(ncol(Xreg)-ncol(xreg)+1):ncol(Xreg)
    xreg=as.matrix(Xreg[,aux])
    Xreg=Xreg[,-aux]
    colnames(xreg)=xregnam
  }

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
  betas=coef(lm(Ystar~-1+Xstar))
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
              fitted.values = fitted , residuals=residuals, Y = Y, p = p, N=N, covmat = sigmae, type = "var",
              xreg = xreg, Ts = c(T = nrow(Xreg), Td = nrow(Xd), Tstar = nrow(Xstar)),delta=delta,lambda=lambda,call=match.call())
  class(result)=c("HDeconometricsVAR","lbvar")

  return(result)
}

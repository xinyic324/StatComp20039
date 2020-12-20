#' @title Principal component analysis using R
#' @description Use R package \code{pca} to achieve data dimensionality reduction of sample matrix
#' @import qpdf
#' @import DT
#' @import bootstrap
#' @import graphics
#' @import DAAG
#' @import boot
#' @import energy
#' @import kedd
#' @import MASS
#' @import RANN
#' @import Ball
#' @import FNN
#' @import stats
#' @param X the Sample matrix(numeric)
#' @return a random sample list
#' @export
pca=function(X){
  x0=apply(X,1,mean)
  sigma=(X-x0)%*%t(X-x0)
  svdw=eigen(sigma)
  percentage2=function(X,per){
    arraysum=sum(X)
    m=length(X)
    temarray=0
    num=0
    for (i in 1:m) {
      temarray=temarray+X[i]
      num=num+1
      if(temarray>=arraysum*per)break
    }
    return(num)
  }
  d=percentage2(svdw$values,per=0.90)
  w1=svdw$vectors
  W=w1[,(1:d)]
  Y=t(W)%*%(X-x0)
  reY=W%*%Y+x0
  return(list(Y=Y,reY=reY,x0=x0,num=d))
}

#' @title Linear Discriminant Analysis using R
#' @description Use R package \code{lda} to achieve data dimensionality reduction
#' @param X the Sample matrix(numeric)
#' @param Y the Sample matrix(numeric)
#' @param d the number of dimension of sample matrix
#' @return a random sample list
#' @export
lda=function(X,Y,d){
  m1=nrow(X)
  m2=nrow(Y)
  m3=ncol(X)
  m4=ncol(Y)
  XY=rbind(X,Y)
  u=apply(XY,2,mean)
  u1=apply(X,2,mean)
  u2=apply(Y,2,mean)
  sb=m1*((u1-u)%*%t(u1-u))+m2*((u2-u)%*%t(u2-u))
  sw1=matrix(0,nrow = m3,ncol = m3)
  sw2=matrix(0,nrow = m4,ncol = m4)
  for (i in 1:m1) {
    sw1=as.matrix((u1-X[i,]))%*%t(as.matrix((u1-X[i,])))+sw1
  }
  for (j in 1:m2) {
    sw2=as.matrix((u2-Y[j,]))%*%t(as.matrix((u2-Y[j,])))+sw2
  }
  sw=sw1+sw2
  w=solve(sw)%*%sb
  wl=eigen(w)
  M=wl$vectors
  w0=M[,c(1:d)]
  X1=X%*%w0
  Y1=Y%*%w0
  return(list(x1=X1,y1=Y1,vectors=M))
} 

#' @title Locally linear embedding using R
#' @description Use R package \code{lle} to achieve data dimensionality reduction
#' @import FNN
#' @importFrom FNN knn.index
#' @param X the Sample matrix(numeric)
#' @param d the number of dimension of sample matrix
#' @param k the number of nearest neighbors
#' @return a random sample matrix
#' @export
lle=function(X,d,k){
  X=t(X)
  l=knn.index(X,k)
  n=nrow(X)
  w=matrix(nrow=n,ncol=k)
  wx=matrix(0,nrow=n,ncol = n)
  for (i in 1:n) {
    z=X[c(l[i,]),]
    lk=apply(-z,1,FUN = "+",X[i,])
    q=t(lk)%*%lk
    lambda=5
    em=solve(q+lambda*diag(k))
    w[i,]=apply(em,1,sum)/sum(em) 
    sh=l[i,]
    for (j in 1:k) {
      wx[i,sh[j]]=w[i,j]
    }
  }
  M=t(diag(n)-wx)%*%(diag(n)-wx)
  NM=eigen(M)$vectors
  D=d+1
  Y=t(NM[,(2:D)])*sqrt(n)
  return(Y)
}

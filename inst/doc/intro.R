## ----eval=FALSE---------------------------------------------------------------
#  function(X){
#    x0=apply(X,1,mean)
#    sigma=(X-x0)%*%t(X-x0)
#    svdw=eigen(sigma)
#    percentage2=function(X,per){
#      arraysum=sum(X)
#      m=length(X)
#      temarray=0
#      num=0
#      for (i in 1:m) {
#        temarray=temarray+X[i]
#        num=num+1
#        if(temarray>=arraysum*per)break
#      }
#      return(num)
#    }
#    d=percentage2(svdw$values,per=0.90)
#    w1=svdw$vectors
#    W=w1[,(1:d)]
#    Y=t(W)%*%(X-x0)
#    reY=W%*%Y+x0
#    return(list(Y=Y,reY=reY,x0=x0,num=d))
#  }

## ----eval=TRUE----------------------------------------------------------------
set.seed(1234)
options(warn =-1)
library(DT)
library(StatComp20039)
X=matrix(runif(7,min=165,max=180),nrow = 1,ncol = 7)
Y=matrix(runif(7,min=55,max=75),nrow = 1,ncol = 7)
Z=matrix(runif(7,min = 17,max = 20),nrow = 1,ncol = 7)
sigma=c(1:ncol(Z))
Acc = rbind(X,Y,Z)
Acc = round(Acc)
datatable(Acc,rownames = c("Heights","Weights","age"),
colnames = paste("S",sigma))

## ----eval=TRUE----------------------------------------------------------------
library(StatComp20039)
X=matrix(runif(7,min=165,max=180),nrow = 1,ncol = 7)
Y=matrix(runif(7,min=55,max=75),nrow = 1,ncol = 7)
Z=matrix(runif(7,min = 17,max = 20),nrow = 1,ncol = 7)
sigma=c(1:ncol(Z))
sigma=c(1:ncol(Z))
Acc = rbind(X,Y,Z)
Acc = round(Acc)
rY=StatComp20039::pca(Acc)$Y
rY=round(rY,3)
datatable(rY,rownames = paste("PC",(1:nrow(rY))),
colnames = paste("S",sigma))

## ----eval=FALSE---------------------------------------------------------------
#  function(X,Y,d){
#    m1=nrow(X)
#    m2=nrow(Y)
#    m3=ncol(X)
#    m4=ncol(Y)
#    XY=rbind(X,Y)
#    u=apply(XY,2,mean)
#    u1=apply(X,2,mean)
#    u2=apply(Y,2,mean)
#    sb=m1*((u1-u)%*%t(u1-u))+m2*((u2-u)%*%t(u2-u))
#    sw1=matrix(0,nrow = m3,ncol = m3)
#    sw2=matrix(0,nrow = m4,ncol = m4)
#    for (i in 1:m1) {
#      sw1=as.matrix((u1-X[i,]))%*%t(as.matrix((u1-X[i,])))+sw1
#    }
#    for (j in 1:m2) {
#      sw2=as.matrix((u2-Y[j,]))%*%t(as.matrix((u2-Y[j,])))+sw2
#    }
#    sw=sw1+sw2
#    w=solve(sw)%*%sb
#    wl=eigen(w)
#    M=wl$vectors
#    w0=M[,c(1:d)]
#    X1=X%*%w0
#    Y1=Y%*%w0
#    return(list(x1=X1,y1=Y1,vectors=M))
#  }

## ----eval=TRUE----------------------------------------------------------------
library(StatComp20039)
X=matrix(c(2.95,6.63,2.53,7.79,3.57,5.65,3.16,5.47),nrow = 4,ncol = 2,byrow = TRUE)
Y=matrix(c(2.58,4.46,2.16,6.22,3.27,3.52),nrow = 3,ncol = 2,byrow = TRUE)
ex=apply(X,2,mean)
ey=apply(Y,2,mean)
E=apply(rbind(X,Y),2,mean)
ee=rbind(ex,ey,E)
par(mfrow=c(1,1))
plot(X[,1],X[,2], pch = 1, xlab = "x", ylab = "y", main = "x vs y",xlim = c(2,3.6),ylim = c(3.5,8))
points(Y[,1],Y[,2],pch = 2,col="red")  
points(ee[,1],ee[,2],pch=8,col="blue")
par(mfrow=c(1,1))
m=StatComp20039::lda(X=X,Y=Y,d=1)$vectors
m1=m[,1]
m2=m[,2]
k1=m1[2]/m1[1]
b1=0
k2=m2[2]/m2[1]
b2=E[2]-k2*E[1]
plot(X[,1],X[,2], pch = 1, xlab = "x", ylab = "y", main = "x vs y",xlim = c(2,6),ylim = c(0,8))
points(Y[,1],Y[,2],pch = 2,col="red")  
points(ee[,1],ee[,2],pch=8,col="blue")
abline(0,k1,col="red")
abline(b2,k2,col="blue")
newx1=length(nrow(X))
newy1=length(nrow(X))
newx2=length(nrow(Y))
newy2=length(nrow(Y))
for (i in 1:4) {
  tem=X[i,]
  newx1[i]=(tem[1]+tem[2]*k1)/(k1*k1+1)
  newy1[i]=k1*newx1[i]
}
for (j in 1:3) {
  temp=Y[j,]
  newx2[j]=(temp[1]+k1*temp[2])/(k1*k1+1)
  newy2[j]=k1*newx2[j]
}
points(newx1,newy1,type = "p",pch=19,col="blueviolet")
points(newx2,newy2,type = "p",col="blue",pch=8)

## ----eval=FALSE---------------------------------------------------------------
#  function(X,d,k){
#    X=t(X)
#    l=knn.index(X,k)
#    n=nrow(X)
#    w=matrix(nrow=n,ncol=k)
#    wx=matrix(0,nrow=n,ncol = n)
#    for (i in 1:n) {
#      z=X[c(l[i,]),]
#      lk=apply(-z,1,FUN = "+",X[i,])
#      q=t(lk)%*%lk
#      lambda=5
#      em=solve(q+lambda*diag(k))
#      w[i,]=apply(em,1,sum)/sum(em)
#      sh=l[i,]
#      for (j in 1:k) {
#        wx[i,sh[j]]=w[i,j]
#      }
#    }
#    M=t(diag(n)-wx)%*%(diag(n)-wx)
#    NM=eigen(M)$vectors
#    D=d+1
#    Y=t(NM[,(2:D)])*sqrt(n)
#    return(Y)
#  }

## ----eval=TRUE----------------------------------------------------------------
library(DT)
library(StatComp20039)
X=matrix(runif(7,min=165,max=180),nrow = 1,ncol = 7)
Y=matrix(runif(7,min=55,max=75),nrow = 1,ncol = 7)
Z=matrix(runif(7,min = 17,max = 20),nrow = 1,ncol = 7)
sigma=c(1:ncol(Z))
Acc = rbind(X,Y,Z)
Acc = round(Acc)
datatable(Acc,rownames = c("Heights","Weights","age"),
colnames = paste("S",sigma))

## ----eval=TRUE----------------------------------------------------------------
options(warn =-1)
library(DT)
library(StatComp20039)
X=matrix(runif(7,min=165,max=180),nrow = 1,ncol = 7)
Y=matrix(runif(7,min=55,max=75),nrow = 1,ncol = 7)
Z=matrix(runif(7,min = 17,max = 20),nrow = 1,ncol = 7)
sigma=c(1:ncol(Z))
Acc = rbind(X,Y,Z)
Acc = round(Acc)
rY=StatComp20039::lle(Acc,d=1,k=2)
rY=round(rY,3)
datatable(rY,rownames = paste("PC",(1:nrow(rY))),
colnames = paste("S",sigma))


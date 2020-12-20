## ----iris---------------------------------------------------------------------
plot(iris)

## ----cars---------------------------------------------------------------------
knitr::kable(head(cars))

## ----eval=FALSE---------------------------------------------------------------
#  n <- 1000
#  u <- runif(n)
#  x <- 2*((1-u)^(-1/2))
#  hist(x, prob = TRUE, main = expression(Pareto(2, 2)))
#  y <- seq(0, 100, .01)
#  lines(y,8/(y^3) )

## -----------------------------------------------------------------------------
n <- 1000;
u <- runif(n);
y <- 2*((1-u)^(-1/4))-2;
m <- 1000;
r <- 4; 
beta <- 2;
lambda <- rgamma(m, r, beta)
x <- rexp(m, lambda) 
par(mfrow=c(1,2))
hist(y, prob = TRUE, main = expression("Empirical Pareto distribution"))
y <- seq(0, 10, .01)
lines(y, r*(beta^r)*((beta+y)^(-r-1)),col="blue")
hist(x, prob = TRUE, main = expression("Theoretical Pareto distribution"))
x <- seq(0, 10, .01)
lines(x, r*(beta^r)*((beta+x)^(-r-1)),col="red")

## -----------------------------------------------------------------------------
festim<-function(a){
x<-runif(a,min=0,max=pi/3)
theta.hat<-mean(sin(x))*(pi/3)##theta.hat表示monte carlo方法生成的估计值
print(c(theta.hat,1-cos(pi/3)))##theta.hat表示由monte carlo方法生成的估计值与理论值
}
festim(1e4)

## -----------------------------------------------------------------------------
mcesit<-function(x,antithetic=TRUE){
u<-runif(x/2)
if(!antithetic)v<-runif(x/2)
else
  v<-1-u
u<-c(u,v)
xdf<-mean(exp(u))

}
set.seed(123)
mc1<-mcesit(x=1e4)##mc1表示有对偶变量法求出的估计值
print(mc1)
set.seed(123)
mc2<-mcesit(x=1e4,antithetic = FALSE)##
print(mc2)##mc2表示有简单蒙特卡罗方法求出的估计值

## -----------------------------------------------------------------------------
n<-1e3
mc1<-numeric(n)
mc2<-numeric(n)
for(i in 1:n){
  mc1[i]<-mcesit(x=1e3)
  mc2[i]<-mcesit(x=1e3,antithetic = FALSE)
}
print(sd(mc1))
print(sd(mc2))
cdl<-(var(mc1)-var(mc2))/var(mc1)
print(cdl)##cdl表示由使用对偶变量法减少的方差
cdl2<-(-exp(2)-3*exp(1)-1)/(-1.5*exp(2)-5*exp(1)-2.5)
print(cdl2)##cdl2表示理论上应该减少的方差

## ----eval=FALSE---------------------------------------------------------------
#  x<-seq(1,5,.005)
#  z<-2
#  f1<-exp(-(x-1)/2)*0.5
#  f2<- x*exp(-(x^2)/2+1/2)
#  L<-(x^2/sqrt(2*pi))*exp(-0.5*x^2)##L表示给定的g(x)
#  plot(x,L,col="black",lwd=z)
#  lines(x,f1,lty=2,col="red",lwd=z)
#  lines(x,f2,lty=3,col="blue",lwd=z)
#  legend("topright", legend = c("L","f1","f2"),
#  col = c("black","red","blue"),lty = 1:3, lwd = z, inset = 0.05)

## ----eval=FALSE---------------------------------------------------------------
#  n<-1e4
#  lambda.hat<-numeric(1)
#  se<-numeric(1)
#  mhh<-numeric(n)
#  HZ<-function(x){
#    (x^2/sqrt(2*pi))*exp(-0.5*x^2)*(x>1)
#  }
#  for(i in 1:n){
#  x<-rexp(1,1/2)+1##这里使用importance function f1
#  L<-function(x){
#  (exp(-(x-1)/2)*0.5)
#  }
#  mhh[i]<- HZ(x)/L(x)
#  }
#  lambda.hat[1]<-mean(mhh)
#  se[1]<-sd(mhh)
#  print(lambda.hat)##估计结果
#  print(se)##估计方差

## -----------------------------------------------------------------------------
n<-1e4
lambda.hat<-numeric(1)
se<-numeric(1)
AS<-numeric(n)
HZ<-function(x){
  (x^2/sqrt(2*pi))*exp(-0.5*x^2)*(x>1)
}
for(i in 1:n){
u<-runif(1)
x<- sqrt(1-2*log(1-u)) ##这里使用importance function f2
S<-function(x){
  x*exp(-(x^2)/2+1/2)
}
AS[i]<-HZ(x)/S(x)
}
lambda.hat<-mean(AS)
se<-sd(AS)
print(lambda.hat)##估计结果
print(se)##估计方差

## -----------------------------------------------------------------------------
a<-1e4
b<-5 ##分层数
k<-a/b ##每层重复次数
n<-30
T2<-numeric(b)
xesitx<-matrix(0,n,2)
g <- function(x) {
exp(-x - log(1+x^2))*(x > 0)*(x < 1)
}

c<-numeric(b)
for(i in 1:b)
{
  c[i]=(1-exp(-1))/(exp((-i+1)/b)-exp(-i/b))
}
print(c)

for(i in 1:n){
u <- runif(a) #f3
x <- - log(1 - u * (1 - exp(-1)))
fg <- g(x) / (exp(-x) / (1 - exp(-1)))
xesitx[i,1]=mean(fg)
for (j in 1:b) {
  u<-runif(a/b,(j-1)/b,j/b)
  x<- -log(exp((-j+1)/b) - (u * (1 - exp(-1))/c[j]))
  T2[j]=mean( g(x) / (c[j]*exp(-x) / (1 - exp(-1))))
}
xesitx[i,2]=sum(T2)
}
xesitx##所得的矩阵

apply(xesitx,2,mean)##估计的均值

apply(xesitx,2,var)##估计的方差



## -----------------------------------------------------------------------------
mu<-2
beta<-1
m<-1e4
n<-20
beta.hat1<-numeric(m)
beta.hat2<-numeric(m)
for(i in 1:m){
x<-numeric(n)
x<-rlnorm(n,2,1)  
y<-log(x)  
z<-sum(log(x))/n
hl<-var(log(x))
beta.hat1[i]<-z-1.96*sqrt(hl/n)##置信上限
beta.hat2[i]<-z+1.96*sqrt(hl/n)##置信下限
}
estib1<-mean(beta.hat1)
estib2<-mean(beta.hat2)
print(c(estib1,estib2))##这是我们求出来的置信区间
indictor<-function(j){
for(i in 1:m){
  if(beta.hat1[i]<=2&&beta.hat2[i]>=2)
  j=j+1
  else{j=j}
}
  j
}
print(indictor(0)/m)##这是我们估计的置信水平

## -----------------------------------------------------------------------------
n <- 20
alpha <- .05
UCL <- replicate(1000, expr = {
x <- rnorm(n, mean = 0, sd = 2)
(n-1) * var(x) / qchisq(alpha, df = n-1)
} )
mean(UCL>4)

## -----------------------------------------------------------------------------
m<-20
n<-1e4
alpha<-0.05
theta.hat1<-numeric(n)
theta.hat2<-numeric(n)
for(i in 1:n){
x<-rchisq(m,2)
theta.hat1[i]<-mean(x)-qt(0.975,19)*sd(x)/sqrt(20)##95%的置信区间下限
theta.hat2[i]<-mean(x)+qt(0.975,19)*sd(x)/sqrt(20)##95%的置信区间上限
}
indictor<-function(j){
for(i in 1:n){
  if(theta.hat1[i]<=2&&theta.hat2[i]>=2)
  j=j+1
  else{j=j}
}
  j
}
print(indictor(0)/n)

## -----------------------------------------------------------------------------
set.seed(1234)
m<-c(10,20,30,50,100)#sample sizes
k<-qnorm(.975,0, sqrt(6*(m-2) / ((m+1)*(m+3))))##k denotes the vector of critical values of the normal distribution
skew<-function(x){
  bar1<-mean(x)
  me1<-mean((x-bar1)^3)
  me2<-mean((x-bar1)^2)
  return(me1/(me2)^(3/2))
}
p1<-length(m)
p2<-length(m)
n<-1e4
xesits1<-length(n)
xesits2<-length(n)
for(i in 1:length(m)){
  for (j in 1:n) {
  x<-rbeta(m[i],1,1)
  y<-rt(m[i],3)##Here we assume v=3
  xesits1[j]<-as.integer(abs(skew(x))>=k[i])
  xesits2[j]<-as.integer(abs(skew(y))>=k[i])
  }
  p1[i]<-mean(xesits1)
  p2[i]<-mean(xesits2)
}
p1##power of the skewness test of normality against symmetric Beta(α, α) distributions 
p2##power of the skewness test of normality against symmetric t(v) distributions



## ----eval=FALSE---------------------------------------------------------------
#  set.seed(123)
#  count5test <- function(x, y) {
#  X <- x - mean(x)
#  Y <- y - mean(y)
#  outx <- sum(X > max(Y)) + sum(X < min(Y))
#  outy <- sum(Y > max(X)) + sum(Y < min(X))
#  return(as.integer(max(c(outx, outy)) > 5))
#  }
#  sigma1 <- 1
#  sigma2 <- 1.5
#  x<-length(m)
#  y<-length(m)
#  m<-1e4
#  n<-c(20,100,1e4)## the small, medium, and large sample sizes
#  power<-length(n)
#  for (i in 1:length(n)) {
#   power[i]<- mean(replicate(m, expr={
#  x <- rnorm(n[i], 0, sigma1)
#  y <- rnorm(n[i], 0, sigma2)
#  count5test(x, y)
#  }))
#  }
#  print(power)

## ----eval=FALSE---------------------------------------------------------------
#  set.seed(123)
#  alpha<-0.055
#  M<-c(20,100,1e4)## the small, medium, and large sample sizes
#  sigma1<-1
#  sigma2<-1.5
#  m<-1000
#  resit<-numeric(m)
#  pw3<-length(M)
#  for(j in 1:length(M)){
#  for(i in 1:m){
#  Y1<-rnorm(M[j],0,sigma1)
#  Y2<-rnorm(M[j],0,sigma2)
#  res.est<-var.test(Y1, Y2, alternative = "two.sided")
#  resit[i]<-as.integer(res.est$p.value<=alpha)
#  }
#    pw3[j]<-mean(resit)
#  }
#  print(pw3)

## ----eval=FALSE---------------------------------------------------------------
#  set.seed(1234)
#  library(MASS)
#  m<-100
#  d<-2
#  M<-c(20,30,50,100)
#  k<-qchisq(0.95,d*(d+1)*(d+2)/6)
#  xskes<-function(x,y,n){
#    z<-cbind(x,y)
#    zbar<-apply(z,2,mean)
#    iz<-z-zbar
#    hatsig<-solve(cov(z))
#    po<-matrix(0,n,n)
#    for(i in 1:n){
#      g<-as.matrix(iz[i,])
#      for(j in 1:n){
#        f<-as.matrix(iz[j,])
#        po[i,j]<-(t(g)%*%hatsig%*%f)^3
#      }
#    }
#    return(sum(po)/(n^2))
#  }## 求多元偏度的值
#  p2<-numeric(length(M))
#  for(i in 1:length(M)){
#    xesite<-numeric(m)
#    for(j in 1:m){
#      mean1<-c(0,0)
#      sigma4<-matrix(c(1,0,0,1),2,2)
#      hats<-mvrnorm(M[i],mean1,sigma4)##生成二维独立正态随机变量
#      xhat<-as.matrix(hats)
#      x<-xhat[,1]
#      y<-xhat[,2]
#      xesite[j]<-as.integer(abs(M[i]*xskes(x,y,M[i])/6)>=k)
#    }
#    p2[i]<-mean(xesite)
#  }
#  p2

## ----eval=FALSE---------------------------------------------------------------
#  set.seed(1234)
#  xskes<-function(z,n){
#    zbar<-apply(z,2,mean)
#    iz<-z-zbar
#    hatsig<-solve(cov(z))
#    sum<-0
#    for(i in 1:n){
#      for(j in 1:n){
#        sum=sum+(t(iz[i,])%*%hatsig%*%iz[j,])^3
#      }
#    }
#    return(sum/n^2)
#  }##计算偏度的函数
#  library(MASS)
#  n<-30
#  m<-100
#  e<-c(seq(0,.15,.01),seq(.15,1,.05))
#  d<-2
#  h<-qchisq(0.95,d*(d+1)*(d+2)/6)
#  pwr1<-length(e)
#  xhat<-matrix(nrow=n,ncol=2)
#  sigma9<-numeric(m)
#  for (i in 1:length(e)) {
#    for (j in 1:m) {
#      sigma<-sample(c(1,100), replace = TRUE, size = n, prob = c(1-e[i], e[i]))
#      for (k in 1:n)
#      { xhat[k,]<-mvrnorm(1,rep(0,2),sigma[k]*diag(2))}
#      sigma9[j]<-as.integer(abs(xskes(xhat,n)*5) >= h)
#      }
#    pwr1[i]<-mean(sigma9)
#  }
#  print(pwr1)
#  ##图像
#  plot(e, pwr1, type = "b",
#  xlab = bquote(epsilon), ylim = c(0,1))
#  abline(h = .1, lty = 3)
#  se <- sqrt(pwr1 * (1-pwr1) / m) #add standard errors
#  lines(e, pwr1+se, lty = 3)
#  lines(e, pwr1-se, lty = 3)

## -----------------------------------------------------------------------------
set.seed(1234)
data(law, package = "bootstrap")
x<-law$LSAT
y<-law$GPA
theta.hat<-cor(x,y)
n<-nrow(law)
theta.jack<-nrow(n)
for (i in 1:n) {
  theta.jack[i]<-cor(x[-i],y[-i])
}
bias<-(n-1)*(mean(theta.jack)-theta.hat)
se<-sqrt((n-1)*
  mean((theta.jack-mean(theta.jack))^2))
print(list(bias=bias,se=se))

## ----eval=FALSE---------------------------------------------------------------
#  set.seed(1234)
#  data(aircondit,package = "boot")
#  y<-as.matrix(aircondit)
#  x<-mean(y)
#  lambda.hat<-1/x
#  n<-nrow(aircondit)
#  B<-200
#  x1<-numeric(B)
#  lambda1<-numeric(B)
#  
#  for (i in 1:B) {
#    i<-sample(1:n,size = n,replace = TRUE)
#    x1[i]<-mean(y[i])
#    lambda1[i]<-1/x1[i]
#  }
#  print(list(lambda.hat=lambda.hat,sd=sd(lambda1),bias=mean(lambda1-lambda.hat)))

## ----eval=FALSE---------------------------------------------------------------
#  library(boot)
#  data(aircondit, package = "bootstrap")
#  lambda.boot<-function(x,i){
#    mean(x[i,])
#  }
#  x<-as.matrix(aircondit)
#  boot.lambda<-boot(x,statistic = lambda.boot,R=2000)
#  print(boot.lambda)
#  print(boot.ci(boot.lambda,type=c("norm","basic","perc","bca")))

## -----------------------------------------------------------------------------
data(scor, package = "bootstrap")
sigma1<-cov(scor)
ev<-eigen(sigma1)
lambda1<-ev$values[1]
theta.hat<-lambda1/sum(ev$values)
n<-nrow(scor)
theta.jack <- numeric(n)
for (i in 1:n) {
sigma1<-cov(scor[-i,])
ev<-eigen(sigma1)
lambda1<-ev$values[1]
theta.jack[i] <-lambda1/sum(ev$values)
}
bias <- (n - 1) * (mean(theta.jack) - theta.hat)
se <- sqrt((n-1) *
mean((theta.jack - mean(theta.jack))^2))
print(list(bias=bias,se=se))

## -----------------------------------------------------------------------------
data(ironslag,package = "DAAG")
n<-length(ironslag$magnetic)
m<-(n-1)*n/2
mt<-as.matrix(ironslag$magnetic)
ct<-as.matrix(ironslag$chemical)
theta1.hat<-theta2.hat<-theta3.hat<-theta4.hat<-matrix(nrow=m,ncol=1)
l<-1
for (i in 1:n) {
  for(j in 1:n){
    if(i!=j&&i<j){
    yt<-mt[-c(i,j),]
    xt<-ct[-c(i,j),]
    
    L1<-lm(yt ~ xt)
    y1.hat1<-L1$coef[1]+L1$coef[2]*ironslag$chemical[i]
    y2.hat1<-L1$coef[1]+L1$coef[2]*ironslag$chemical[j]
    theta1.hat[l,]<-(ironslag$magnetic[i]-y1.hat1)^2+(ironslag$magnetic[j]-y2.hat1)^2
        
    L2 <- lm(yt ~ xt + I(xt^2))
    y1.hat2 <- L2$coef[1] + L2$coef[2] * ironslag$chemical[i] +L2$coef[3] * ironslag$chemical[i]^2
    y2.hat2 <- L2$coef[1] + L2$coef[2] * ironslag$chemical[j] +L2$coef[3] * ironslag$chemical[j]^2
    theta2.hat[l,] <- (ironslag$magnetic[i] - y1.hat2)^2+(ironslag$magnetic[j] - y2.hat2)^2
    
    L3 <- lm(log(yt) ~ xt)
    logy1.hat3 <- L3$coef[1] + L3$coef[2] * ironslag$chemical[i]
    logy2.hat3 <- L3$coef[1] + L3$coef[2] * ironslag$chemical[j]
     y1.hat3 <- exp(logy1.hat3)
     y2.hat3<-exp(logy2.hat3)
    theta3.hat[l,]<- (ironslag$magnetic[i] - y1.hat3)^2+(ironslag$magnetic[j] - y2.hat3)^2
    
    L4 <- lm(log(yt) ~ log(xt))
    logy1.hat4 <- L4$coef[1] + L4$coef[2] * log(ironslag$chemical[i])
    logy2.hat4 <- L4$coef[1] + L4$coef[2] * log(ironslag$chemical[j])
    y1.hat4 <- exp(logy1.hat4)
    y2.hat4 <- exp(logy2.hat4)
    theta4.hat[l,]<- (ironslag$magnetic[i] - y1.hat4)^2+(ironslag$magnetic[j] - y2.hat4)^2
    
    if(l<=m-1) l=l+1 else(break)
    }
    
  }
}
print(list(se1=sum(theta1.hat)/(2*m),se2=sum(theta2.hat)/(2*m),se3=sum(theta3.hat)/(2*m),se4=sum(theta4.hat)/(2*m)))

## -----------------------------------------------------------------------------
data("ironslag",package = "DAAG")
model2<-lm(ironslag$magnetic~ironslag$chemical+I(ironslag$chemical^2))
model2

## ----eval=FALSE---------------------------------------------------------------
#  set.seed(1234)
#  R<-999
#  n1<-20
#  n2<-30
#  mu1<-mu2<-0
#  sigma1<-sigma2<-1
#  outx1<-numeric(R)
#  outy1<-numeric(R)
#  m<-200
#  P<-numeric(m)
#  for (j in 1:m) {
#  x<-rnorm(n1,mu1,sigma1)
#  y<-rnorm(n2,mu2,sigma2)
#  z<-c(x,y)
#  D<-numeric(R)
#  X <- x - mean(x)
#  Y <- y - mean(y)
#  outx <- sum(X > max(Y)) + sum(X < min(Y))
#  outy <- sum(Y > max(X)) + sum(Y < min(X))
#  D0<-max(c(outx, outy))
#   for (i in 1:R) {
#    k<-sample(1:50,size = 20,replace = FALSE)
#    x1<-z[k]-mean(z[k])
#    y1<-z[-k]-mean(z[-k])
#    outx1[i]<- sum(x1 > max(y1)) + sum(x1 < min(y1))
#    outy1[i]<- sum(y1 > max(x1)) + sum(y1 < min(x1))
#    D[i]=max(c(outx1[i],outy1[i]))
#  }
#  P[j]<-mean(c(D0,D)>=D0)
#  }
#  mean(P)

## ----eval=FALSE---------------------------------------------------------------
#  library(RANN)
#  library(boot)
#  library(energy)
#  library(Ball)
#  Tn <- function(z, ix, sizes,k) {
#  n1 <- sizes[1]; n2 <- sizes[2]; n <- n1 + n2
#  if(is.vector(z)) z <- data.frame(z,0);
#  z <- z[ix, ];
#  NN <- nn2(data=z, k=k+1) # what's the first column?
#  block1 <- NN$nn.idx[1:n1,-1]
#  block2 <- NN$nn.idx[(n1+1):n,-1]
#  i1 <- sum(block1 < n1 + .5); i2 <- sum(block2 > n1+.5)
#  (i1 + i2) / (k * n)
#  }
#  set.seed(12345)
#  m<-1e3
#  k<-3
#  p<-2
#  mu<-0.3
#  n1<-50
#  n2<-50
#  R<-999
#  n<-n1+n2
#  N<-c(n1,n2)
#  eqdist.nn<-function(n,sizes,k){
#    boot.obj<-boot(data=z,statistic=Tn,R=R,
#  sim = "permutation", sizes = sizes,k=k)
#    ts <- c(boot.obj$t0,boot.obj$t)
#  p.value <- mean(ts>=ts[1])
#  list(statistic=ts[1],p.value=p.value)
#  }
#  p.values <- matrix(NA,m,3)
#  for(i in 1:m){
#  x <- matrix(rnorm(n1*p,mean=mu,sd=1.5),ncol=p);
#  y <- cbind(rnorm(n2,mean=mu),rnorm(n2,mean=mu));
#  z <- rbind(x,y)
#  p.values[i,1] <- eqdist.nn(z,N,k)$p.value
#  p.values[i,2] <- eqdist.etest(z,sizes=N,R=R)$p.value
#  p.values[i,3] <- bd.test(x=x,y=y,num.permutations=999,seed=i*12345)$p.value
#  }
#  alpha <- 0.1;
#  pow1 <- colMeans(p.values<alpha)
#  pow1

## ----eval=FALSE---------------------------------------------------------------
#  library(RANN)
#  library(boot)
#  library(energy)
#  library(Ball)
#  Tn <- function(z, ix, sizes,k) {
#  n1 <- sizes[1]; n2 <- sizes[2]; n <- n1 + n2
#  if(is.vector(z)) z <- data.frame(z,0);
#  z <- z[ix, ];
#  NN <- nn2(data=z, k=k+1) # what's the first column?
#  block1 <- NN$nn.idx[1:n1,-1]
#  block2 <- NN$nn.idx[(n1+1):n,-1]
#  i1 <- sum(block1 < n1 + .5); i2 <- sum(block2 > n1+.5)
#  (i1 + i2) / (k * n)
#  }
#  set.seed(12345)
#  m<-1e3
#  k<-3
#  p<-2
#  mu<-0.3
#  mu1<-0.5
#  n1<-50
#  n2<-50
#  R<-999
#  n<-n1+n2
#  N<-c(n1,n2)
#  eqdist.nn<-function(n,sizes,k){
#    boot.obj<-boot(data=z,statistic=Tn,R=R,
#  sim = "permutation", sizes = sizes,k=k)
#    ts <- c(boot.obj$t0,boot.obj$t)
#  p.value <- mean(ts>=ts[1])
#  list(statistic=ts[1],p.value=p.value)
#  }
#  p.values <- matrix(NA,m,3)
#  for(i in 1:m){
#  x <- matrix(rnorm(n1*p,mean=mu1,sd=1.5),ncol=p);
#  y <- cbind(rnorm(n2,mean=mu,sd=1),rnorm(n2,mean=mu,sd=1));
#  z <- rbind(x,y)
#  p.values[i,1] <- eqdist.nn(z,N,k)$p.value
#  p.values[i,2] <- eqdist.etest(z,sizes=N,R=R)$p.value
#  p.values[i,3] <- bd.test(x=x,y=y,num.permutations=999,seed=i*12345)$p.value
#  }
#  alpha <- 0.1;
#  pow2 <- colMeans(p.values<alpha)
#  pow2

## ----eval=FALSE---------------------------------------------------------------
#  library(RANN)
#  library(boot)
#  library(energy)
#  library(Ball)
#  Tn <- function(z, ix, sizes,k) {
#  n1 <- sizes[1]; n2 <- sizes[2]; n <- n1 + n2
#  if(is.vector(z)) z <- data.frame(z,0);
#  z <- z[ix, ];
#  NN <- nn2(data=z, k=k+1) # what's the first column?
#  block1 <- NN$nn.idx[1:n1,-1]
#  block2 <- NN$nn.idx[(n1+1):n,-1]
#  i1 <- sum(block1 < n1 + .5); i2 <- sum(block2 > n1+.5)
#  (i1 + i2) / (k * n)
#  }
#  set.seed(12345)
#  m<-1e3
#  k<-3
#  p<-2
#  p1<-0.8
#  n1<-50
#  n2<-50
#  R<-999
#  n<-n1+n2
#  N<-c(n1,n2)
#  eqdist.nn<-function(n,sizes,k){
#    boot.obj<-boot(data=z,statistic=Tn,R=R,
#  sim = "permutation", sizes = sizes,k=k)
#    ts <- c(boot.obj$t0,boot.obj$t)
#  p.value <- mean(ts>=ts[1])
#  list(statistic=ts[1],p.value=p.value)
#  }
#  p.values <- matrix(NA,m,3)
#  for(i in 1:m){
#  x <- matrix(rt(n1*p,df=1),ncol=p);
#  sigma<- sample(c(1, 10), size = n2,
#  replace = TRUE, prob = c(p1, 1-p1))
#  y <- cbind(rnorm(n2,mean=0,sd=sigma),rnorm(n2,mean=0,sd=sigma));
#  z <- rbind(x,y)
#  p.values[i,1] <- eqdist.nn(z,N,k)$p.value
#  p.values[i,2] <- eqdist.etest(z,sizes=N,R=R)$p.value
#  p.values[i,3] <- bd.test(x=x,y=y,num.permutations=999,seed=i*12345)$p.value
#  }
#  alpha <- 0.1;
#  pow3 <- colMeans(p.values<alpha)
#  pow3

## ----eval=FALSE---------------------------------------------------------------
#  library(RANN)
#  library(boot)
#  library(energy)
#  library(Ball)
#  Tn <- function(z, ix, sizes,k) {
#  n1 <- sizes[1]; n2 <- sizes[2]; n <- n1 + n2
#  if(is.vector(z)) z <- data.frame(z,0);
#  z <- z[ix, ];
#  NN <- nn2(data=z, k=k+1) # what's the first column?
#  block1 <- NN$nn.idx[1:n1,-1]
#  block2 <- NN$nn.idx[(n1+1):n,-1]
#  i1 <- sum(block1 < n1 + .5); i2 <- sum(block2 > n1+.5)
#  (i1 + i2) / (k * n)
#  }
#  set.seed(12345)
#  m<-1e3
#  k<-3
#  p<-2
#  p1<-0.8
#  n1<-50
#  n2<-50
#  R<-999
#  n<-n1+n2
#  N<-c(n1,n2)
#  eqdist.nn<-function(n,sizes,k){
#    boot.obj<-boot(data=z,statistic=Tn,R=R,
#  sim = "permutation", sizes = sizes,k=k)
#    ts <- c(boot.obj$t0,boot.obj$t)
#  p.value <- mean(ts>=ts[1])
#  list(statistic=ts[1],p.value=p.value)
#  }
#  p.values <- matrix(NA,m,3)
#  for(i in 1:m){
#  x <- matrix(rt(n1*p,df=1),ncol=p);
#  y <- cbind(rt(n2,df=2),rt(n2,df=2));
#  z <- rbind(x,y)
#  p.values[i,1] <- eqdist.nn(z,N,k)$p.value
#  p.values[i,2] <- eqdist.etest(z,sizes=N,R=R)$p.value
#  p.values[i,3] <- bd.test(x=x,y=y,num.permutations=999,seed=i*12345)$p.value
#  }
#  alpha <- 0.1;
#  pow6 <- colMeans(p.values<alpha)
#  pow6

## ----eval=FALSE---------------------------------------------------------------
#  library(RANN)
#  library(boot)
#  library(energy)
#  library(Ball)
#  Tn <- function(z, ix, sizes,k) {
#  n1 <- sizes[1]; n2 <- sizes[2]; n <- n1 + n2
#  if(is.vector(z)) z <- data.frame(z,0);
#  z <- z[ix, ];
#  NN <- nn2(data=z, k=k+1) # what's the first column?
#  block1 <- NN$nn.idx[1:n1,-1]
#  block2 <- NN$nn.idx[(n1+1):n,-1]
#  i1 <- sum(block1 < n1 + .5); i2 <- sum(block2 > n1+.5)
#  (i1 + i2) / (k * n)
#  }
#  set.seed(12345)
#  m<-1e3
#  k<-3
#  p<-2
#  p1<-0.8
#  p2<-0.9
#  n1<-50
#  n2<-50
#  R<-999
#  n<-n1+n2
#  N<-c(n1,n2)
#  eqdist.nn<-function(n,sizes,k){
#    boot.obj<-boot(data=z,statistic=Tn,R=R,
#  sim = "permutation", sizes = sizes,k=k)
#    ts <- c(boot.obj$t0,boot.obj$t)
#  p.value <- mean(ts>=ts[1])
#  list(statistic=ts[1],p.value=p.value)
#  }
#  p.values <- matrix(NA,m,3)
#  for(i in 1:m){
#  sigma<- sample(c(1, 3), size = n2,
#  replace = TRUE, prob = c(p1, 1-p1))
#  sigma1<- sample(c(1, 2), size = n2,
#  replace = TRUE, prob = c(p2, 1-p2))
#  x <-cbind(rnorm(n1,mean=0,sd=sigma1),rnorm(n1,mean=0,sd=sigma1));
#  y <- cbind(rnorm(n2,mean=0,sd=sigma),rnorm(n2,mean=0,sd=sigma));
#  z <- rbind(x,y)
#  p.values[i,1] <- eqdist.nn(z,N,k)$p.value
#  p.values[i,2] <- eqdist.etest(z,sizes=N,R=R)$p.value
#  p.values[i,3] <- bd.test(x=x,y=y,num.permutations=999,seed=i*12345)$p.value
#  }
#  alpha <- 0.1;
#  pow3 <- colMeans(p.values<alpha)
#  pow3

## ----eval=FALSE---------------------------------------------------------------
#  library(RANN)
#  library(boot)
#  library(energy)
#  library(Ball)
#  Tn <- function(z, ix, sizes,k) {
#  n1 <- sizes[1]; n2 <- sizes[2]; n <- n1 + n2
#  if(is.vector(z)) z <- data.frame(z,0);
#  z <- z[ix, ];
#  NN <- nn2(data=z, k=k+1) # what's the first column?
#  block1 <- NN$nn.idx[1:n1,-1]
#  block2 <- NN$nn.idx[(n1+1):n,-1]
#  i1 <- sum(block1 < n1 + .5); i2 <- sum(block2 > n1+.5)
#  (i1 + i2) / (k * n)
#  }
#  set.seed(12345)
#  m<-1e3
#  k<-3
#  p<-2
#  mu<-0.3
#  n1<-10
#  n2<-100
#  R<-999
#  n<-n1+n2
#  N<-c(n1,n2)
#  eqdist.nn<-function(n,sizes,k){
#    boot.obj<-boot(data=z,statistic=Tn,R=R,
#  sim = "permutation", sizes = sizes,k=k)
#    ts <- c(boot.obj$t0,boot.obj$t)
#  p.value <- mean(ts>=ts[1])
#  list(statistic=ts[1],p.value=p.value)
#  }
#  p.values <- matrix(NA,m,3)
#  for(i in 1:m){
#  x <- matrix(rnorm(n1*p,mean=mu,sd=1),ncol=p);
#  y <- cbind(rnorm(n2,mean=mu,sd=1),rnorm(n2,mean=mu,sd=1));
#  z <- rbind(x,y)
#  p.values[i,1] <- eqdist.nn(z,N,k)$p.value
#  p.values[i,2] <- eqdist.etest(z,sizes=N,R=R)$p.value
#  p.values[i,3] <- bd.test(x=x,y=y,num.permutations=999,seed=i*12345)$p.value
#  }
#  alpha <- 0.1;
#  pow4 <- colMeans(p.values<alpha)
#  pow4

## -----------------------------------------------------------------------------
library(kedd)
n=1000
lambda=matrix(c(0,1,-1,1/10,-1/2,1/10,0,1/10,1/2,1/10,1,1/10),nrow=6,ncol=2,byrow=TRUE)
P=c(.5,.1,.1,.1,.1,.1)
y=numeric(n)
y=sample(1:6,size=n,replace = TRUE,prob=P)
rate=as.matrix(lambda[y,])
mean1=rate[,1]
sd1=rate[,2]
x=numeric(n)
for (i in 1:n) {
 x[i]=rnorm(1,mean=mean1[i],sd=sd1[i]) 
}
hist(x)

## -----------------------------------------------------------------------------
hist(x,breaks = "Scott")

## -----------------------------------------------------------------------------
hist(x,breaks = "FD")

## -----------------------------------------------------------------------------
set.seed(1234)
library(kedd)##先产生100个数据点
n<-100
theta<-matrix(c(0,1,1,0.3),nrow=2,ncol = 2,byrow=TRUE)
y<-sample(1:2,size = n,replace = TRUE,prob = c(0.3,0.7))
z<-as.matrix(theta[y,])
mean2<-z[,1]
sd2<-z[,2]
x<-numeric(n)
for(i in 1:100){
  x[i]<-rnorm(1,mean=mean2[i],sd=sd2[i])  
}
## 带宽选择
h=numeric(7)
h[1]=h.amise(x,deriv.order = 0)$h
h[2]=h.mlcv(x,kernel= "gaussian" )$h
h[3]=h.ucv(x,deriv.order = 0)$h
h[4]=h.bcv(x,whichbcv = 1,deriv.order = 0)$h
h[5]=h.ccv(x,deriv.order = 0,upper=0.5)$h
h[6]=h.mcv(x,deriv.order = 0,upper=0.5)$h
h[7]=h.tcv(x,deriv.order = 0)$h
fx=function(x)0.3*dnorm(x,0,1)+0.7*dnorm(x,1,0.3)
hatf=dkde(x,deriv.order = 0,h=h[1])
hatf1=dkde(x,deriv.order = 0,h=h[2])
hatf2=dkde(x,deriv.order = 0,h=h[3])
hatf3=dkde(x,deriv.order = 0,h=h[4])
hatf4=dkde(x,deriv.order = 0,h=h[5])
hatf5=dkde(x,deriv.order = 0,h=h[6])
hatf6=dkde(x,deriv.order = 0,h=h[7])
##画图
z=seq(-3,3,0.05)
plot(z,fx(z),xlab="",ylab="density function",col="black",main = "kernel density estimation",ylim = c(0,1.2),type = "c")
lines(hatf,col="black",lty=2)
lines(hatf1,col="blue")
lines(hatf2,col="black")
lines(hatf3,col="brown")
lines(hatf4,col="blue4")
lines(hatf5,col="antiquewhite4")
lines(hatf6,col="aquamarine4")
par(new=TRUE)
plot(fx,xlab="",xlim=c(-3,3),ylim=c(0,1.2),col="red",lty=2)
legend("topright",legend=c("f(x)","h.amise","h.mlcv","h.ucv","h.bcv","h.ccv","h.mcv","h.tcv"),col=c("red","black","blue","black","brown","blue4","antiquewhite4","aquamarine4"),inset = 0.02,lty=c(2,2,1,1,1,1,1,1))

## ----eval=FALSE---------------------------------------------------------------
#  set.seed(12345)
#  N=1e4
#  n1=444
#  n2=132
#  n3=361
#  n4=63
#  n=c(n1,n2,n3,n4)
#  L=c(0.5,0.3,0.2)
#  tol= .Machine$double.eps^0.5
#  L.old=L+1
#  for (j in 1:N) {
#    mu1=1/(2*sum(n))*(n4+n1*(1+(L[1])/(2-2*L[2]-L[1])))
#    mu2=1/(2*sum(n))*(n4+n2*(1+(L[2])/(2-2*L[1]-L[2])))
#    mu3=1-mu1-mu2
#    L=c(mu1,mu2,mu3)
#    if (sum(abs(L - L.old)/L.old) < tol) break
#    L.old=L
#  }
#  print(list(p= L[1],q=L[2], iter = j, tol = tol))

## ----eval=FALSE---------------------------------------------------------------
#  set.seed(12345)
#  N=1e4
#  n1=444
#  n2=132
#  n3=361
#  n4=63
#  n=c(n1,n2,n3,n4)
#  L=c(0.5,0.3,0.2)
#  tol= .Machine$double.eps^0.5
#  L.old=L+1
#  p=numeric(0)
#  q=numeric(0)
#  l=numeric(0)
#  for (j in 1:N) {
#    mu1=1/(2*sum(n))*(n4+n1*(1+(L[1])/(2-2*L[2]-L[1])))
#    mu2=1/(2*sum(n))*(n4+n2*(1+(L[2])/(2-2*L[1]-L[2])))
#    mu3=1-mu1-mu2
#    L=c(mu1,mu2,mu3)
#    p[j]=mu1
#    q[j]=mu2
#    l[j]=n1*log(mu1^2+2*mu1*mu3)+n2*log(mu2^2+2*mu2*mu3)+(2*n3)*log(mu3)+n4*log(2*mu1*mu2)
#    if (sum(abs(L - L.old)/L.old) < tol) break
#    L.old=L
#  }
#  print(list(p=p,q=q,l=l))

## -----------------------------------------------------------------------------
mpg=mtcars$mpg
disp=mtcars$disp
wt=mtcars$wt
formulas=list(
mpg ~ disp,
mpg ~ I(1 / disp),
mpg ~ disp + wt,
mpg ~ I(1 / disp) + wt
)
fit.models=function(i){
  lm(formulas[[i]])
}
lapply(seq_along(formulas), fit.models)

## -----------------------------------------------------------------------------
mpg=mtcars$mpg
disp=mtcars$disp
wt=mtcars$wt
formulas=list(
mpg ~ disp,
mpg ~ I(1 / disp),
mpg ~ disp + wt,
mpg ~ I(1 / disp) + wt
)
out=vector("list", length(formulas))
for (i in seq_along(formulas)) {
out[[i]]=lm(formulas[[i]])
}
out

## -----------------------------------------------------------------------------
set.seed(1234)
##the anonymous function is defined
trials <- replicate(
100,
t.test(rpois(10, 10), rpois(7, 10)),
simplify = FALSE
)
sapply(trials,function(f)f$p.value)

## -----------------------------------------------------------------------------
x=list(runif(11), runif(6)) 
y=list(rpois(11, 2) + 1, rpois(6, 2) + 1)
vapply(Map(weighted.mean,x,y),function(x) x,FUN.VALUE = c(a=0))

## ----eval=FALSE---------------------------------------------------------------
#  set.seed(3000)
#  lap_f = function(x) exp(-abs(x))
#  rw.Metropolis = function(sigma, x0, N){
#   x = numeric(N)
#   x[1] = x0
#   u = runif(N)
#   k = 0
#   for (i in 2:N) {
#    y = rnorm(1, x[i-1], sigma)
#    if (u[i] <= (lap_f(y) / lap_f(x[i-1]))) x[i] = y
#    else {
#    x[i] = x[i-1]
#    k = k+1
#    }
#   }
#   return(list(x = x, k = k))
#  }
#  N = 2000
#  sigma = c(.05, .5, 2, 16)
#  x0 = 25
#  rw1 = rw.Metropolis(sigma[1],x0,N)
#  rw2 = rw.Metropolis(sigma[2],x0,N)
#  rw3 = rw.Metropolis(sigma[3],x0,N)
#  rw4 = rw.Metropolis(sigma[4],x0,N)
#  #number of candidate points rejected
#  Rej = cbind(rw1$k, rw2$k, rw3$k, rw4$k)
#  Acc = round((N-Rej)/N,4)
#  rownames(Acc) = "Accept rates"
#  colnames(Acc) = paste("sigma",sigma)
#  knitr::kable(Acc)
#  #plot
#  par(mfrow=c(2,2))  #display 4 graphs together
#      rw = cbind(rw1$x, rw2$x, rw3$x,  rw4$x)
#      for (j in 1:4) {
#          plot(rw[,j], type="l",
#               xlab=bquote(sigma == .(round(sigma[j],3))),
#               ylab="X", ylim=range(rw[,j]))
#      }

## ----eval=FALSE---------------------------------------------------------------
#  o=ppoints(100)
#  zr=numeric(100)
#  for (i in 1:100) {
#    zr[i]=log(2*(o[i]))*as.numeric(o[i]<=0.5)+(-log(2-2*o[i]))*as.numeric(o[i]>0.5)
#  }
#  zx1=quantile(y1,o)
#  zx2=quantile(y2,o)
#  zx3=quantile(y3,o)
#  zx4=quantile(y4,o)
#  zy1=quantile(rw1$x,o)
#  zy2=quantile(rw2$x,o)
#  zy3=quantile(rw3$x,o)
#  zy4=quantile(rw4$x,o)
#  qqplot(zx1,zy1,xlab = "RcPP Function",ylab ="R Function" )
#  qqline(zy1)
#  par(mfrow=c(1,1))
#  qqplot(zx2,zy2,xlab = "RcPP Function",ylab ="R Function")
#  qqline(zy2)
#  par(mfrow=c(1,1))
#  qqplot(zx3,zy3,xlab = "RcPP Function",ylab ="R Function")
#  qqline(zy3)
#  par(mfrow=c(1,1))
#  qqplot(zx4,zy4,xlab = "RcPP Function",ylab ="R Function")
#  qqline(zy4)


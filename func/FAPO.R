###############################################################
#                      3/8/2017   fapo 
# use k_T to generate basis functions; use L2 as the inner 
# product; use simpson's rule to compute L2 inner product
# contain both the balanced case and unbalanced case
# for either the balanced or the unbalanced case, use ttt from 
# the output of generatex-p30-n300.R; in the balanced case
# ttt contains identical rows; in the unbalanced case, ttt contains
# randomly generated rows. 
###############################################################
###############################################################
# function: generate simpson's weights (ns must be odd)
###############################################################
simpson=function(ns){
  wei=rep(0,ns);wei[1]=1;wei[ns]=1
  for(i in 1:((ns-1)/2)) wei[2*i]=4
  for(i in 1:((ns-1)/2-1)) wei[2*i+1]=2
  h=1/(ns-1); wei=wei*(h/3)
  return(wei)
}
############################################################### 
# function: operator norm                                 
############################################################### 
onorm=function(a) return(eigen(round((a+t(a))/2,8))$values[1])
##############################################################
#      function: trace of a matrix
##############################################################
tr=function(a) return(sum(diag(a)))
###############################################################
#          function: Moore-Penrose type power              
#           keep only the first m eigenvalues                  
###############################################################
matpower = function(a,alpha){
  a = (a + t(a))/2
  tmp = eigen(a)
  if(length(tmp$values) == 1){
    m <- tmp$vectors%*%(abs(tmp$values)^alpha)%*%t(tmp$vectors)
  } else{
    m <- tmp$vectors%*%diag(abs(tmp$values)^alpha)%*%t(tmp$vectors)
  }
  return(m)}
#
mppower = function(matrix,power,ignore){
  eig = eigen(matrix)
  eval = eig$values
  evec = eig$vectors
  m = length(eval[abs(eval)>ignore])
  if(m == 1){
    tmp = (eval[1:m]^power)*evec[,1:m]%*%
      t(evec[,1:m])
  } else{
    tmp = evec[,1:m]%*%diag(eval[1:m]^power)%*%
      t(evec[,1:m])
  }
  return(tmp)
}
############################################################### 
# function:  Gram matrix for first-level Hilbert space
#            for a given kernel and time set
# tte: is the time set for evaluation
# tto: the time set for observations
# if you just want to compute gram matrix then set tte=tto
#  kern is either brown or gauss
############################################################### 
gramt=function(tte,tto,kern){
  ltte=length(tte);ltto=length(tto)
  if (kern=="gauss"){
    a1=matrix(tte^2,ltte,ltto);a2=tte%*%t(tto);a3=t(matrix(tto^2,ltto,ltte))
    a=a1-2*a2+a3
    b1=matrix(tto^2,ltto,ltto);b2=tto%*%t(tto);b3=t(matrix(tto^2,ltto,ltto))
    b=b1-2*b2+b3
    sigma=sum(sqrt(b))/(2*choose(ltto,2));gamma=1/(2*sigma^2)
    ktmat=exp(-gamma*(a))}
  if(kern=="brown"){
    arr=array(0,c(ltte,ltto,2))
    arr[,,1]=matrix(tte,ltte,ltto);arr[,,2]=t(matrix(tto,ltto,ltte))
    ktmat=apply(arr,c(1,2),min)}
  return(t(ktmat))
}
##############################################################
#   function:  Generalized Cross Validation for expsilon_T 
#  et is a grid of values
# xxx are the observed random functions
##############################################################
gcvt=function(xxx,ttt,et,kern){
  n=dim(xxx)[1];nt=dim(xxx)[2];p=dim(xxx)[3]
  nuset=numeric();deset=numeric()
  for(i in 1:n){
    kt=gramt(ttt[i,],ttt[i,],kern)
    scale=onorm(kt)
    ktinv=matpower(kt+et*scale*diag(nt),-1)
    for(j in 1:p) {
      nuset=c(nuset,sum((xxx[i,,j]-kt%*%ktinv%*%xxx[i,,j])^2))
      deset=c(deset,(1-tr(kt%*%ktinv)/nt)^2)}}
  out=sum(nuset/deset)
  return(out)
}
############################################################### 
# function: estimates one function
############################################################### 
evalx=function(f,tte,tto,ridge,kern){
  kt=gramt(tto,tto,kern)
  scale=eigen(kt)$values[1]	
  ktinv=matpower(kt+scale*ridge*diag(nrow(kt)),-1)
  kt1=t(gramt(tte,tto,kern))
  out=kt1%*%ktinv%*%f
  return(c(out))}
############################################################### 
# function: estimates a sample of functions
############################################################### 
evalxmat=function(ff,tte,ttt,ridge,kern){
  n=dim(ff)[1];ffcoo=numeric()
  for(i in 1:n) ffcoo=rbind(ffcoo,evalx(ff[i,],tte,ttt[i,],ridge,kern)) 
  return(ffcoo) 
}
############################################################### 
# function: gram matrix for the second RHKS for a sample of x
# x are the estimated functions at the evaluation points
# kt is the diagonal of simpson weights to compute the L2 inner product
############################################################### 
gramx=function(x,kt){
  n=dim(x)[1]
  k2=x%*%kt%*%t(x);k1=t(matrix(diag(k2),n,n));k3=t(k1);k=k1-2*k2+k3
  sigma=sum(sqrt(abs(k)))/(2*choose(n,2));gamma=1/(2*sigma^2)
  return(exp(-gamma*(k1-2*k2+k3)))
}
############################################################### 
# function: Q = I - J/n                                   
############################################################### 
qmat = function(n) return(diag(n)-rep(1,n)%*%t(rep(1,n))/n)
###############################################################
#  function:  K to G=QKQ
###############################################################
cgram=function(K){n=dim(K)[1];Q=qmat(n);return(Q%*%K%*%Q)}
###############################################################
# function: Compute matrix A_i and Lambda_i
###############################################################
aili=function(Gi,ridge){
  n=dim(Gi)[1];scale=eigen(Gi)$values[1];Q=qmat(n)
  mat= Gi+ridge*scale*Q
  Ai=(n^(-1/2))*mppower(Gi,1/2,10^(-7))%*%mppower(mat,-1/2,10^(-7))
  Li=Q-Ai%*%Ai	
  return(list(Ai=Ai,Li=Li))
}
############################################################### 
# fapo: using simpson's rule to compute the L2 inner product
# xvec are the estimated functions
# ridge is the optimal ex
# ns number of evaluation points
############################################################### 
fapo=function(xvec,ridge){
  n=dim(xvec)[1];p=dim(xvec)[3];eps=10^(-7);ns=dim(xvec)[2]
  kt=diag(simpson(ns))
  H=numeric();Linv=numeric()
  for(i in 1:p){
    Gi=cgram(gramx(xvec[,,i],kt));store=aili(Gi,ridge)
    H=rbind(H,store$Ai);Linv=rbind(Linv,mppower(store$Li,-1,eps))}
  M=0;for(i in 1:p){
    Hi=H[((i-1)*n+1):(i*n),];Linvi=Linv[((i-1)*n+1):(i*n),];M=M+Hi%*%Linvi%*%Hi}
  Q=qmat(n);Minv=mppower(M+Q,-1,eps)
  LinvH=numeric()
  for(i in 1:p) LinvH=rbind(LinvH,Linv[((i-1)*n+1):(i*n),]%*%H[((i-1)*n+1):(i*n),])
  big=LinvH%*%Minv%*%t(LinvH);norm1=matrix(0,p,p)
  for (i in 2:p) {
    for (j in 1:(i-1)){
      norm1[i,j] <- norm1[j,i] <- as.numeric(norm(big[((i-1)*n+1):(i*n),((j-1)*n+1):(j*n)],"F"))}}
  return(norm1)
}
##############################################################
#   function:  Generalized Cross Validation for expsilon_X 
#  xxxeva are the estimated functions
#  ex is a grid of points
##############################################################
gcvx=function(xxxeva,ex){
  ns=dim(xxxeva)[2];n=dim(xxxeva)[1]
  kt=diag(simpson(ns))
  nuset=numeric();deset=numeric() 
  for(i in 1:(p-1)) for(j in (i+1):p){
    Gi=gramx(xxxeva[,,i],kt);Gj=gramx(xxxeva[,,j],kt)
    scale=onorm(Gi);Gin2=matpower(Gi+ex*scale*diag(n),-2)
    nuset=c(nuset,sum((Gi-Gj%*%Gi%*%Gin2%*%Gi)^2))
    deset=c(deset,(1-tr(Gi%*%Gin2%*%Gi)/n)^2)}
  return(sum(nuset/deset))
}

#
# density under alternative
#


alt_model = function(X,Jm,Jv,theta,psi,tSigma,tSigma2,tSigma3,dSigma,Adjmu,Adjmu2,Adjmu3,Adjmu4,mu,mu3){
#cat(length(theta),length(psi),'\n')
theta = theta
psi = psi
L = length(X)
vect = rep(0,L)


mean = Jm + exp(theta)*tSigma + Adjmu*mu^2
var = Jv + 4*exp(theta)*tSigma3 + dSigma*(exp(psi)-2*exp(2*theta)) + exp(2*theta)*tSigma2 + Adjmu2*exp(theta)*mu^2+Adjmu3*mu*mu3+4*Adjmu4*mu^2
vect = dnorm(X,mean,sqrt(var))
#cat(theta,psi,var,mean,'\n')

vect[is.na(vect)]=-10000
vect[is.nan(vect)]=-10000
return (vect)
}



#
# density under null
#
#redone!
null_model = function(X,Jm,Jv,tSigma,tSigma2,tSigma3,dSigma,Adjmu,Adjmu2,Adjmu3,Adjmu4,mu,mu3){
vect = alt_model(X,Jm,Jv,theta=-Inf,psi=-Inf,tSigma,tSigma2,tSigma3,dSigma,Adjmu,Adjmu2,Adjmu3,Adjmu4,0,0)
return (vect)
}


#
#  Likelihood mu and musq
#

log_L = function(p,X,Jm,Jv,theta,psi,tSigma,tSigma2,tSigma3,dSigma,Adjmu,Adjmu2,Adjmu3,Adjmu4,mu,mu3){
#cat(X,'\n')
if (is.na(theta) || is.na(psi)){
return (-10000)
}


V  = NULL
L = length(psi)
for (i in 1:L){ 
alt_d = alt_model(X,Jm,Jv,theta,psi,tSigma,tSigma2,tSigma3,dSigma,Adjmu,Adjmu2,Adjmu3,Adjmu4,mu,mu3)
alt_n = null_model(X,Jm,Jv,tSigma,tSigma2,tSigma3,dSigma,Adjmu,Adjmu2,Adjmu3,Adjmu4,mu,mu3)
A = sum(log(p*alt_d + (1-p)*alt_n))
V = c(V,A)
}
#cat(V,'\n')
return (V)
}





log_L00 = function(V){
#cat(X[1],'\n')
p = rep(V[5],length(X))
theta = V[1]
psi = V[2]
mu = V[3]
mu3 = V[4]
#cat(Jm[1],Jv[1],tSigma[1],tSigma2[1],tSigma3[1],dSigma[1],Adjmu[1],Adjmu2[1],Adjmu3[1],Adjmu4[1],'\n')

if (is.na(theta) || is.na(psi)){
return (-10000)
}



V  = NULL
L = length(psi)
for (i in 1:L){ 
alt_d = alt_model(X,Jm,Jv,theta,psi,tSigma,tSigma2,tSigma3,dSigma,Adjmu,Adjmu2,Adjmu3,Adjmu4,mu,mu3)
alt_n = null_model(X,Jm,Jv,tSigma,tSigma2,tSigma3,dSigma,Adjmu,Adjmu2,Adjmu3,Adjmu4,mu,mu3)
A = sum(log(p*alt_d + (1-p)*alt_n))
V = c(V,A)
}
if (is.na(V)){V=-10000}
if (V>10000){V=10000}
if (V < -10000){V = -10000}
#cat(V,'\n')
return (V)
}



find_GlobalO = function(theta,psi,mu,mu3,p,logL){
a=NULL
invisible(capture.output(o<- optimx(c(theta,psi,mu,mu3,p),log_L00,control=list(maximize=TRUE),upper=c(Inf,Inf,Inf,Inf,1),lower=c(-Inf,-Inf,-Inf,-Inf,0),method='L-BFGS-B')))
vm00 = as.vector(o[1:5])
L00 = o[6]

if (L00>=logL){
a=c(a,vm00)
return (list(a=a,logL = L00))
}
cat('Cannot find','\n')
return (list(a=c(theta,psi,mu,mu3,p),logL = logL))
}




find_GlobalOS = function(theta,psi,mu,mu3,p,logL){
a=NULL

thetaX = seq(-2*theta,2*theta,length.out=4)
psiX = seq(-2*psi,2*psi,length.out=4)
pX = seq(0.01,3*p,length.out=4)
muX = seq(-1*mu,1*mu,length.out=4)
mu3X = seq(-1*mu3,1*mu3,length.out=4)


#cat(X,'\n')
vm00X = NULL
L00X = NULL

for (i in 1:4){
for (j in 1:4){
for (k in 1:4){
for (m in 1:4){
for (l in 1:4){

tryCatch(invisible(capture.output(o<- optimx(c(thetaX[i],psiX[j],muX[m],mu3X[l],pX[k]),log_L00,control=list(maximize=TRUE),upper=c(Inf,Inf,Inf,Inf,1),lower=c(-Inf,-Inf,-Inf,-Inf,0),method='L-BFGS-B'))))
vm00X = cbind(vm00X,unlist(o[1:5]))
L00X = c(L00X,unlist(o[6]))
}
}
}
}
}
#cat(L00X,'\n')
kkk = which(L00X == max(L00X,na.rm=T))[1]

vm00 = vm00X[,kkk]
#cat(vm00,'\n')
L00 = L00X[kkk]
if (L00>=logL){
a=c(a,vm00)
return (list(a=a,logL = L00))
}

cat('Cannot find','\n')
return (list(a=c(theta,psi,mu,mu3,p),logL = logL))
}






#
# Get p-value
#


library(MASS)

get_pvalTRUENC = function(p,X,J,Jm,Jv,theta,psi,SigmaList,tSigma,tSigma2,tSigma3,dSigma,Adjmu,Adjmu2,Adjmu3,Adjmu4,mu,mu3,perm){

X<<-X
J<<-J
Jm<<-Jm
Jv<<-Jv
SigmaList<<-SigmaList
tSigma<<-tSigma
tSigma2<<-tSigma2
tSigma3<<-tSigma3
dSigma<<-dSigma
Adjmu<<-Adjmu
Adjmu2<<-Adjmu2
Adjmu3<<-Adjmu3
Adjmu4<<-Adjmu4
mu<<-mu
mu3<<-mu3

logH0 = log_L(rep(0,length(X)),X,Jm,Jv,theta,psi,tSigma,tSigma2,tSigma3,dSigma,Adjmu,Adjmu2,Adjmu3,Adjmu4,mu,mu3)
#cat(m00,p,'\n')

aaaa = rep(p,length(X))
ll = log_L(aaaa*0,X,Jm,Jv,theta,psi,tSigma,tSigma2,tSigma3,dSigma,Adjmu,Adjmu2,Adjmu3,Adjmu4,mu,mu3)
res2 <- try(B <- find_GlobalOS(theta,psi,mu,mu3,p,ll))
Bxx = B
logHa = B$logL

pp =getP(unlist(B$a))

T = unlist(logHa) - unlist(logH0)


if (T <= 0){ 
perm=2
}

ST = NULL
L = length(X)


P=p
for (i in 1:perm){
if (i %% 2 == 1){cat('Permutation',i,'\n')}

X = NULL
mu0 = NULL
mu03 = NULL
deltaV = NULL
for (j in 1:L){
Sigma = SigmaList[[j]]
D = diag(Sigma)
deltaV = c(deltaV,D)
tt = rep(0,length(D))
v = (Sigma+ diag(tt))
s = rmvn(1,rep(0,length(D)),v)
mu0 = c(mu0,mean(s))
mu03 = c(mu03,mean((s-mean(s))^3))
invS = solve(Sigma)
s = sum((s)%*%invS%*%invS%*%t(s))
X = c(X,s)
}
X<<-X
mu<<-mean(mu0)
mu3<<-mean(mu03)

theta = max((sum(X - Jm - Adjmu*mu^2))/sum(tSigma),1/sum(tSigma))
psi = theta
#cat(theta,psi,'\n')
#cat(X,'\n')
theta = log(theta)
psi = theta



aaaa = rep(P,length(X))
ll = log_L(aaaa*0,X,Jm,Jv,theta,psi,tSigma,tSigma2,tSigma3,dSigma,Adjmu,Adjmu2,Adjmu3,Adjmu4,mu,mu3)
res2 <-  try(B <- find_GlobalOS(theta,psi,mu,mu3,p,ll))

if(inherits(res2, "try-error") ){
cat('','\n')
cat('Error',i,'\n')
cat('','\n')
}else{
logH0 = log_L(rep(0,length(X)),X,Jm,Jv,theta,psi,tSigma,tSigma2,tSigma3,dSigma,Adjmu,Adjmu2,Adjmu3,Adjmu4,mu,mu3)
logHa = B$logL
D = unlist(logHa) - unlist(logH0)
ST =  c(ST,unlist(D))
if (i==500){
pval = mean(T<=ST)
if (pval>0.1){break}
}
if (i==50){
pval = mean(T<=ST)
if (pval>0.25){break}
}
if (i==100){
pval = mean(T<=ST)
if (pval>0.15){break}
}
}
}
cat(sum(ST>100),'\n')
if (T<0){T==0}
ST[ST<0] = 0
#cat(ST,'\n')
p = sum(ST>=T)/length(ST)
return (list(p=c(p,T),prob = pp,est = Bxx))
}


getP = function(V){
p = rep(V[5],length(X))
theta = V[1]
psi = V[2]
mu = V[3]
mu3 = V[4]

alt_d = alt_model(X,Jm,Jv,theta,psi,tSigma,tSigma2,tSigma3,dSigma,Adjmu,Adjmu2,Adjmu3,Adjmu4,mu,mu3)
alt_n = null_model(X,Jm,Jv,tSigma,tSigma2,tSigma3,dSigma,Adjmu,Adjmu2,Adjmu3,Adjmu4,mu,mu3)
A =p*alt_d + (1-p)*alt_n

pos = p*alt_d/A 
return (pos)
}



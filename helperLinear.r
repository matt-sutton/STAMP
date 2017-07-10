#
# For Sum of Beta
#


#
# Alt model
#

alt_modelS  = function(X,J,mu,invdelta,Sum,Omega1,Omega2,tau){
L = length(X)
vect = NULL
#cat(mu,tau,'\n')
#cat(invdelta,'\n')
for (i in 1:L){
mu1 = mu*sqrt(invdelta[i])*Omega1[i]
tau1 = (Sum[i] + exp(tau)*Omega2[i])*invdelta[i]
z = dnorm(X[i],mu1,sqrt(tau1))
vect = c(vect,z)
}
return (vect)
}


#
# Null model
#

null_modelS = function(X,J,invdelta,Sum,Omega1,Omega2,tau){
vect = alt_modelS(X,J,0,invdelta,Sum,Omega1,Omega2,-Inf)
return (vect)
}



#
# get conditional p
#

get_ES = function(p,X,J,mu,invdelta,Sum,Omega1,Omega2,tau){

vect_alt = alt_modelS(X,J,mu,invdelta,Sum,Omega1,Omega2,tau)
vect_null = null_modelS(X,J,invdelta,Sum,Omega1,Omega2,tau)
Total = p*vect_alt + (1-p)*vect_null
gamma = p*vect_alt/Total 
return (gamma)
}

solve_mu0_tau0 = function(muV){
mu = muV[1]
tau = muV[2]
####cat('mu',gamma1,'\n')
Smu = 0
Stau = 0

sigmaTT = 1

for (i in 1:length(gamma1)){
muT = X[i] - mu*sqrt(invdelta[i])*Omega1[1]
sigmaT = (Sum[i] + exp(tau)*Omega2[i])*invdelta[i]
Smu = Smu + gamma1[i]*muT/sigmaT*sqrt(invdelta[i])*Omega1[1]
Stau = Stau +  gamma1[i]*((muT)^2/sigmaT^2 - 1/sigmaT)*invdelta[i]*Omega2[i]
}

aa = c(Smu,Stau*exp(muV[2]))
aa[is.na(aa)] = -1

return (aa)
}

log_LS = function(V){
p = rep(V[3],length(X))
mu = V[1]
tau = V[2]


if (is.na(mu) || is.na(tau)){
return (-Inf)
}

V  = NULL
L = length(tau)
for (i in 1:L){ 
alt_d = alt_modelS(X,J,mu,invdelta,Sum,Omega1,Omega2,tau[i])
alt_n = null_modelS(X,J,invdelta,Sum,Omega1,Omega2,tau[i])
A = sum(log(p*alt_d + (1-p)*alt_n))
V = c(V,A)
}
if (is.na(V)){V=-10000}
if (V>10000){V=10000}
if (V < -10000){V = -10000}
return (V)
}

get_pS = function(vect){
L = length(vect)
p = sum(vect[1:(L)])/(L)
p = rep(p,L)
return (p)
}

get_EH = function(p,X,J,mu,invdelta,Sum,Omega1,Omega2,tau){
cat(1,tau,'\n')
#save(p,X,J,mu,delta,deltasq,tau,tausq,file='a')
XXX = NULL
pin = p 
muin = mu 
tauin = tau 

pold = p 
muold = mu + 0.1
tauold = tau + 0.1


vin = c(pin[1],muin,tauin)
vold = c(pold[1],muold,tauold)
#cat('Vin',vin,'\n')
kk = 1
aaa = NULL
rr=0

X <<- X
Sum<<-Sum
Omega1<<-Omega1
Omega2<<-Omega2
mu <<- mu
tau <<- tau
invdelta <<- invdelta
J<-J
XXX = NULL
V = c(mu,tau,p[1])
lin = log_LS(V)
#cat(lin,'\n')
la = lin+1

cat(2,tau,'\n')
kk = 0
while (((abs(lin-la)/abs(max((lin),(la))))>5*10^(-4))){
gamma1 <<- get_ES(pin,X,J,muin,invdelta,Sum,Omega1,Omega2,tauin)

muv = findGlobalS(muin,tauin)

muold = muv$mv[1]
tauold = muv$mv[2]


pold = get_pS(gamma1)

pin = pold
muin = muold
tauin = tauold
V = c(muin,tauin,pin[1])

la = lin
lin = log_LS(V)
#cat('lin',lin,la,'\n')

kk=kk+1
if (kk==10){break}
}

#cat('kk',kk,'\n')
return (c(pin[1],muin,tauin))
}

findGlobalS = function(muin,tauin){
cat(muin,tauin,'\n')
mv = nleqslv(c(muin,tauin),solve_mu0_tau0,control=list(cndtol=10^(-100),btol=.001,ftol=10^(-25),xtol=10^(-25),maxit=1000))$x
if ((sum(abs(solve_mu0_tau0(mv)))<0.0001) & sum(abs(mv))<1000){
C = check_tauS(mv)
if (C){
##cat(-1,'\n')
return (list(mv =mv,est =c(muin,tauin)))
}
}

aam = seq(muin-muin,muin+muin,length.out=20)
bbm = seq(tauin-tauin,tauin+tauin,length.out=20)

L  = length(aam)-1
R = NULL
for (i in 1:(L/2)){
for (j in 1:(L/2)){

mv = nleqslv(c(aam[i],bbm[j]),solve_mu0_tau0,control=list(cndtol=10^(-100),btol=.001,ftol=10^(-25),xtol=10^(-25),maxit=1000))$x
if ((sum(abs(solve_mu0_tau0(mv)))<0.0001) & sum(abs(mv))<1000){
C = check_tauS(mv)
if (C){
##cat(-1,'\n')
return (list(mv =mv,est =c(aam[i],bbm[j])))
}
}

mv = nleqslv(c(aam[L-i],bbm[j]),solve_mu0_tau0,control=list(cndtol=10^(-100),btol=.001,ftol=10^(-25),xtol=10^(-25),maxit=1000))$x
if ((sum(abs(solve_mu0_tau0(mv)))<0.0001) & sum(abs(mv))<1000){
C = check_tauS(mv)
if (C){
##cat(-1,'\n')
return (list(mv =mv,est =c(aam[L-i],bbm[j])))
}
}


mv = nleqslv(c(aam[i],bbm[L-j]),solve_mu0_tau0,control=list(cndtol=10^(-100),btol=.001,ftol=10^(-25),xtol=10^(-25),maxit=1000))$x
if ((sum(abs(solve_mu0_tau0(mv)))<0.0001) & sum(abs(mv))<1000){
C = check_tauS(mv)
if (C){
##cat(-1,'\n')
return (list(mv =mv,est =c(aam[i],bbm[L-j])))
}
}


mv = nleqslv(c(aam[L/2 + i],bbm[L/2 + j]),solve_mu0_tau0,control=list(cndtol=10^(-100),btol=.001,ftol=10^(-25),xtol=10^(-25),maxit=1000))$x
if ((sum(abs(solve_mu0_tau0(mv)))<0.0001) & sum(abs(mv))<1000){
C = check_tauS(mv)
if (C){
##cat(-1,'\n')
return (list(mv =mv,est =c(aam[L/2 + i],bbm[L/2 + j])))
}
}



mv = nleqslv(c(aam[L/2 - i+1],bbm[L/2 - j+1]),solve_mu0_tau0,control=list(cndtol=10^(-100),btol=.001,ftol=10^(-25),xtol=10^(-25),maxit=1000))$x
if ((sum(abs(solve_mu0_tau0(mv)))<0.0001) & sum(abs(mv))<1000){
C = check_tauS(mv)
if (C){
##cat(-1,'\n')
return (list(mv =mv,est =c(aam[L/2 - i+1],bbm[L/2 - j+1])))
}
}


}
}

return (list(mv =c(muin,tauin),est =c(muin,tauin)))

}

check_tauS = function(mv){

alt_d = alt_modelS(X,J,mv[1],invdelta,Sum,Omega1,Omega2,mv[2])
alt_n = null_modelS(X,J,invdelta,Sum,Omega1,Omega2,mv[2])
A1 = sum(gamma1*log(alt_d) + (1-gamma1)*log(alt_n))
if (sum(is.na(alt_d))|| sum(is.na(alt_n))>0){
return (FALSE)
}
alt_d = alt_modelS(X,J,mu,invdelta,Sum,Omega1,Omega2,tau)
alt_n = null_modelS(X,J,invdelta,Sum,Omega1,Omega2,tau)
if (sum(is.na(alt_d))|| sum(is.na(alt_n))>0){
return (FALSE)
}
A2 = sum(gamma1*log(alt_d) + (1-gamma1)*log(alt_n))
##cat(3,A1,A2,'\n')
if (is.na(A1) || is.na(A2) || is.na(A1>=A2)){
return (FALSE)
}
if (A1>=A2){
##cat(33,A1,A2,'\n')
return (TRUE)
}else{
return (FALSE)
}
}



findGlobalSO = function(muin,tauin,pin){
o= optimx(c(muin,tauin,pin),log_LS,control=list(maximize=TRUE),upper=c(Inf,Inf,1),lower=c(-Inf,-Inf,0),method='L-BFGS-B')
vm00 = as.vector(o[1:3])
L00 = o[4]
return (list(est=vm00,L=L00))
}

findGlobalSOO = function(tauin){
o= optimx(c(tauin),log_L0S,control=list(maximize=TRUE),method='L-BFGS-B')
vm00 = as.vector(o[1])
L00 = o[2]
return (list(est=vm00,L=L00))
}

log_L0S = function(V){
#cat(1,'\n')
p = rep(0,length(X))
tau = V[1]
#cat(tau,'\n',X,'\n',J,invdelta,'\n')

if (is.na(tau)){
return (-Inf)
}

V  = NULL
L = length(tau)
#cat(L,'\n')
for (i in 1:L){ 
alt_d = alt_modelS(X,J,mu=1,invdelta,Sum,Omega1,Omega2,tau[i])
alt_n = null_modelS(X,J,invdelta,Sum,Omega1,Omega2,tau[i])
A = sum(log(p*alt_d + (1-p)*alt_n))
V = c(V,A)
}
if (is.na(V)){V=-10000}
if (V>10000){V=10000}
if (V < -10000){V = -10000}
return (V)
}


get_pval = function(p,X,J,mu,invdelta,Sum,Omega1,Omega2,tau,perm){


aa = mu
bb = tau
X<<-X
J<<-J
invdelta<<-invdelta
Sum<<-Sum
Omega1<<-Omega1
cat(11,tau,'\n')
Omega2<<-Omega2

LogH0 = log_L0S(1)
L = length(X)
cat(0,tau,'\n')
a = get_EH(rep(p,L),X,J,mu,invdelta,Sum,Omega1,Omega2,tau)
cat(a,'\n')
est = unlist(findGlobalSO(a[2],a[3],a[1])$est)
cat(est,'\n')
LogHa = unlist(findGlobalSO(a[2],a[3],a[1])$L)

pin = est[3]
muin = est[1]
tauin = est[2]

gamma1p = get_ES(pin,X,J,muin,invdelta,Sum,Omega1,Omega2,tauin)

cat(LogHa,LogH0,'\n')
Tobs = LogHa - LogH0
cat(Tobs,'\n')
vectT = NULL
if (Tobs<=0){
perm=2
}

for (i in 1:perm){

X = rnorm(L)*sqrt(Sum*invdelta)
X<<-X
XL = X
All = max((var(XL)*19 - sum(Sum*invdelta))/sum(invdelta),0.1*var(XL)/mean(invdelta))
muL = 1/2*sqrt(All)*sign(mean(XL))
tauL = log(All)


bb = tauL
aa = muL

LogH0 = log_L0S(1)

a = get_EH(rep(p,L),X,J,aa,invdelta,Sum,Omega1,Omega2,bb)
LogHa = unlist(findGlobalSO(a[2],a[3],a[1])$L)
T = LogHa - LogH0
vectT = c(vectT,T)

if (i==50){
pval = mean(Tobs<=vectT)
if (pval>0.25){break}
}
if (i==100){
pval = mean(Tobs<=vectT)
if (pval>0.15){break}
}
if (i==500){
pval = mean(Tobs<=vectT)
if (pval>0.1){break}
}

}
vectT[vectT<0]=0
pval = mean(Tobs<=vectT)
return (list(pv=pval,pp = gamma1p))
}


get_pvalT = function(X,Var){
T = sum((1/Var*X))/sqrt(sum(1/Var))
p = 2*(1-pnorm(abs(T)))
return (p)
}


get_pvalH =function(X,Total){
T = sum(X)
df1 = Total
p = 1-pchisq(T,df=df1)
return (p)
}




library('psych')
library('pscl')
library('optimx')
library('nleqslv')
library('mvnfast')
library('CompQuadForm')

#
# Conduct asociation analysis with linear test 
# under mixture model
#
# Input parameters: 
#	BetaList - list of genetic effect estimates from each study (S vectors)
#   SeList - list of standart errors from each study (S vectors)
#   CorList - list of  LD pairwise correlation between SNPs (S vectors)
#   WeightList - list of square root of Weights, identity matrix by default (S vectors)
#	ListofMAF - list of MAFs (S vectors)
#	perm - number of permutations

# Output parameters:
#	pvalue - pvalue
#	LogLR - log likelihood ratio
#	Est - list of estimates
#	posterPr - posterior probability


# load helper functions
source('helperLinear.r')

mixtureLinear = function(perm,BetaList,SeList,CorList,ListofMAF,WeightList=NULL){


library('psych')

listSE = list()
listBeta = list()
listCorr = list() 
listOmega = list()


Quad_Mixture = NULL
SigmaList= list()
Sigmas = NULL
Xq = NULL
tSigma=NULL
tSigma2=NULL 
tSigma3=NULL
dSigma = NULL
Adjmu = NULL
Adjmu2 = NULL
Adjmu3 = NULL
Adjmu4 = NULL
Sum = NULL
invdelta = NULL
Jm = NULL
Jv = NULL
mu0 = NULL
mu03 = NULL
JV = NULL
Omega1 = NULL
Omega2 = NULL


L = length(BetaList)


for (st in 1:L){

Beta = BetaList[[st]]
SE =  SeList[[st]]
Covv = CorList[[st]]
V = 2*ListofMAF[[st]]*(1-ListofMAF[[st]])


if (is.null(WeightList)){
A = diag(rep(1,length(Beta)))
}else{
A = WeightList[[st]]
}


listOmega[[st]] = diag(1/sqrt(V))%*%Covv%*%diag(sqrt(V))
listCorr[[st]] = Covv
listBeta[[st]] = Beta
listSE[[st]] = SE
Var = diag(listSE[[st]])%*%listCorr[[st]]%*%diag(listSE[[st]])
JV = c(JV,length(listBeta[[st]]))
Mcor = solve(A)%*%solve(listOmega[[st]])
beta = listBeta[[st]]
Sigma = Var

S = Mcor%*%beta*20
Omega1 = c(Omega1,sum(listOmega[[st]]))
Omega2 = c(Omega2,sum(listOmega[[st]]%*%listOmega[[st]]))
Sigma = Mcor%*%Sigma%*%t(Mcor)*400
Sigmas=c(Sigmas,diag(Sigma))
VarM = Sigma
mu0 = c(mu0,mean(S))
mu03 = c(mu03,mean((S-mean(S))^3))
invS = solve(Sigma)

SigmaList[[st]] = VarM
Jm = c(Jm,tr(invS))
tSigma=c(tSigma,tr(invS%*%invS))
tSigma3=c(tSigma3,tr(invS%*%invS%*%invS))
Adjmu4=c(Adjmu4,sum(invS%*%invS%*%invS))
invS = invS%*%invS
Jv = c(Jv,2*tr(invS))
invS2 = invS%*%invS
tSigma2=c(tSigma2,2*tr(invS2))
dSigma = c(dSigma,tr(diag(diag(invS))%*%invS))
Adjmu = c(Adjmu,sum(invS))
Adjmu2 = c(Adjmu2,4*sum(invS2))
Adjmu3 = c(Adjmu3,4*sum(invS%*%diag(invS)))

}




mu = mean(mu0)*1
mu3 = mean(mu03)*1

Sigmasq = Sigmas^2
#m001 = get_estim(Xq,Sigmas,Sigmasq) 


XL = Linear_Mixture
All = max((var(XL)*19 - sum(Sum*invdelta))/sum(invdelta),0.1*var(XL)/mean(invdelta))
muL = 1/3*sqrt(All)*sign(mean(XL))
tauL = log(All)/2


p=0.25

output = get_pval(p,XL,JV,muL,invdelta,Sum,Omega1,Omega2,tauL,perm)

pvalue = output$pv
LogLR  = output$T
posterPr = output$pp
Est = output$Est
return (list(pvalue=pvalue,LogLR=LogLR,posterPr=posterPr,Est=Est))
}





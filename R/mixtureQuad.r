#'Conduct asociation analysis with quadartic test under mixture model
#'
#' @param perm number of permutations
#' @param BetaList list of genetic effect estimates from each study (S vectors)
#' @param SeList list of standart errors from each study (S vectors)
#' @param CorList list of  LD pairwise correlation between SNPs (S vectors)
#' @param ListofMAF list of MAFs (S vectors)
#' @param WeightList list of square root of Weights, identity matrix by default (S vectors)
#'
#' @export
#'
#' @return \code{mixtureQuad} returns a list object with parameters:
#' \item{pvalue}{p-values}
#' \item{LogLR}{log likelihood ratio}
#' \item{Est}{list of estimates}
#' \item{posterPr}{Posterior probability}
#'
#' @examples
#' data(Example_SC11)
#' \donttest{mixtureQuad(1000,Example_SC11$ListBeta,Example_SC11$ListSE,Example_SC11$ListCor,Example_SC11$MafList) }
#'
#'

mixtureQuad = function(perm,BetaList,SeList,CorList,ListofMAF,WeightList=NULL){
  
  library('psych')
  library('pscl')
  library('optimx')
  library('nleqslv')
  library('mvnfast')
  library('CompQuadForm')
  
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
    Xq =c(Xq,(S)^2)
    Sigma = Mcor%*%Sigma%*%t(Mcor)*400
    Sigmas=c(Sigmas,diag(Sigma))
    VarM = Sigma
    mu0 = c(mu0,mean(S))
    mu03 = c(mu03,mean((S-mean(S))^3))
    invS = solve(Sigma)
    
    Quad_Mixture = c(Quad_Mixture,sum(t(S)%*%(invS)%*%invS%*%S))
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
  
  XQ = Quad_Mixture
  theta = max((sum(XQ - Jm - Adjmu*mu^2))/sum(tSigma),1/sum(tSigma))
  psi = theta
  theta = log(theta)
  psi = log(psi)
  
  
  p=0.25
  output = get_pvalTRUENC(p,Quad_Mixture,JV,Jm,Jv,theta,psi,SigmaList,tSigma,tSigma2,tSigma3,dSigma,Adjmu,Adjmu2,Adjmu3,Adjmu4,mu,mu3,perm)
  
  pvalue = output$p[1]
  LogLR  = output$p[2]
  posterPr = output$prob
  Est = list(pi=output$est$a[5],mu=output$est$a[3]/20,mu3=output$est$a[4]/20^3,logpsi=output$est$a[2]-2*log(20),logtheta=output$est$a[2]-4*log(20))
  return (list(pvalue=pvalue,LogLR=LogLR,posterPr=posterPr,Est=Est))
}

rm(list = ls())


library(EMMIXskew)
library(mclust)
library(clusteval)
#library(ContaminatedMixt)



WD.PATH = paste("D:/Data and code MCNFA/code/simulation2")

source(paste(WD.PATH, '/MCNFA.na.EM.r', sep=""))
source(paste(WD.PATH, '/MCNFA.na.ECM.r', sep=""))
source(paste(WD.PATH, '/MCNFA.na.AECM.r', sep=""))
source(paste(WD.PATH, '/MtFA.na.EM.r', sep=""))
source(paste(WD.PATH, '/MtFA.na.ECM.r', sep=""))
source(paste(WD.PATH, '/MFA.na.EM.2018.r', sep=""))
source(paste(WD.PATH, '/MFA.na.ECM.2018.r', sep=""))
source(paste(WD.PATH, '/MFA.na.AECM.2018.r', sep=""))
source(paste(WD.PATH, '/MtFA.na.AECM.r', sep=""))


################important function######################
########################################################

gener.na = function(X, na.rate)
{
  #set.seed(1)
  n = nrow(X)
  p = ncol(X)
  num.na = floor(n*p*na.rate)
  keep.part = sample(p, n, replace = T) + p * 0:(n-1)
  na.posi = tabulate(sample((1:(n*p))[-keep.part], num.na), nbins = n*p)
  na.posi = matrix(na.posi, ncol = p, byrow = T)
  X[na.posi == 1] = NA
  return(X)
}

mix.N.gener = function(p, mu, S, sample.size = 1000){
  n = sample.size
  g = length(p)
  y = list(g)
  ni = rmultinom(1, size = n, prob = p)
  class = rep(1:g,c(ni))
  for(i in 1:g){
    Sigma = S[,,i]
    y[[i]] = rmvnorm(n=ni[i], mean=mu[,i], sigma=Sigma)
  }
  data = y[[1]]
  for(j in 2:g){
    data = rbind(data,y[[j]])
  }
  return(list(data = data, class = class))
  
}

#############################################


fit.CN.AECM10 = fit.CN.ECM10 = fit.CN.EM10 = list()
fit.t.AECM10 = fit.t.ECM10 = fit.t.EM10 = list()
fit.N.AECM10 = fit.N.ECM10 = fit.N.EM10 = list()

fit.CN.AECM20 = fit.CN.ECM20 = fit.CN.EM20 = list()
fit.t.AECM20 = fit.t.ECM20 = fit.t.EM20 = list()
fit.N.AECM20 = fit.N.ECM20 = fit.N.EM20 = list()

Y.real = list(); Y.10 = list(); Y.20 = list()
tru.class = list()


n=2000; MM = 200
for(l in 1:100){
  
  
  p=10;g=3;q=3;MM=20
  W = c(1/3,1/3,1/3)
  #poor
  mu = cbind(c(-7,-7,7,7,7,7,7,7,7,7),
             c(-5,0,0,0,0,0,0,0,0,0),
             c(7,2,2,2,2,2,2,2,2,2));
  #well
  #mu = cbind(c(-3,3,3,3,3,3,3,3,3,3),
  #           c(-2,0,0,0,0,0,0,0,0,0),
  #           c(0,-1,-1,-1,-1,-1,-1,-1,-1,-1));
  
  A = array(NA,dim = c(10,3,3))
  A[,,1] = cbind(c(rmvnorm(1,rep(0, p), diag(p))),c(rmvnorm(1,rep(0, p), diag(p))),c(rmvnorm(1,rep(0, p), diag(p))))
  A[,,2] = cbind(c(rmvnorm(1,rep(0, p), diag(p))),c(rmvnorm(1,rep(0, p), diag(p))),c(rmvnorm(1,rep(0, p), diag(p))))
  A[,,3] = cbind(c(rmvnorm(1,rep(0, p), diag(p))),c(rmvnorm(1,rep(0, p), diag(p))),c(rmvnorm(1,rep(0, p), diag(p))))
  D = cbind(runif(p,0.5,1),runif(p,0.5,1),runif(p,0.5,1))
  
  S = array(0,c(p,p,g))
  S[,,1]=A[,,1]%*% t(A[,,1])+diag(D[,1])
  S[,,2]=A[,,2]%*% t(A[,,2])+diag(D[,2])
  S[,,3]=A[,,3]%*% t(A[,,3])+diag(D[,3])
  
  DATA = mix.N.gener(W, mu, S, sample.size = n)
  uni.data = cbind(runif(MM,-5,10),runif(MM,-5,10),runif(MM,-5,10),runif(MM,-5,10),runif(MM,-5,10),
                   runif(MM,-5,10),runif(MM,-5,10),runif(MM,-5,10),runif(MM,-5,10),runif(MM,-5,10))
  
  
  
  Y = DATA$data
  tru.class[[l]]=true.clus=DATA$class
  
  Y.new = rbind(Y,uni.data)
  Y.real[[l]] = Y.new
  
  ####################################################################################
  ####################################################################################
  ####################################################################################
  
  Y.na10 = gener.na (Y.new, 0.1)
  
  Y.10[[l]] = Y.na10
  
  ## AECM
  
  fit.CN.AECM10[[l]]=MCNFA.na.AECM (Y.na10, g, q, eta=0.005, tol=0.0001, max.iter=5000, per=5000)
  
  
  fit.t.AECM10[[l]]=MtFA.na.AECM (Y.na10, g, q, eta=0.005, tol=0.0001, max.iter=5000, per=5000)
  
  fit.N.AECM10[[l]]=MFA.na.AECM (Y.na10, g, q, eta=0.005, tol=0.0001, max.iter=5000, per=5000)
  
  ## ECM
  
  
  fit.CN.ECM10[[l]]=MCNFA.na.ECM (Y.na10, g, q, eta=0.005, tol=0.0001, max.iter=5000, per=5000)
  
  
  fit.t.ECM10[[l]]=MtFA.na.ECM (Y.na10, g, q, eta=0.005, tol=0.0001, max.iter=5000, per=5000)
  
  fit.N.ECM10[[l]]=MFA.na.ECM (Y.na10, g, q, eta=0.005, tol=0.0001, max.iter=5000, per=5000)
  
  ## EM
  
  
  fit.CN.EM10[[l]]=MCNFA.na.EM (Y.na10, g, q, eta=0.005, tol=0.0001, max.iter=5000, per=5000)
  
  
  fit.t.EM10[[l]]=MtFA.na.EM (Y.na10, g, q, eta=0.005, tol=0.0001, max.iter=5000, per=5000)
  
  fit.N.EM10[[l]]=MFA.na.EM (Y.na10, g, q, eta=0.005, tol=0.0001, max.iter=5000, per=5000)
  
  
  ####################################################################################
  ####################################################################################
  
  Y.na20 = gener.na (Y.new, 0.2)
  Y.20[[l]] = Y.na20
  
  ## AECM
  fit.CN.AECM20[[l]]=MCNFA.na.AECM (Y.na20, g, q, eta=0.005, tol=0.0001, max.iter=5000, per=5000)
  
  
  fit.t.AECM20[[l]]=MtFA.na.AECM (Y.na20, g, q, eta=0.005, tol=0.0001, max.iter=5000, per=5000)
  
  fit.N.AECM20[[l]]=MFA.na.AECM (Y.na20, g, q, eta=0.005, tol=0.0001, max.iter=5000, per=5000)
  
  ## ECM
  
  
  fit.CN.ECM20[[l]]=MCNFA.na.ECM (Y.na20, g, q, eta=0.005, tol=0.0001, max.iter=5000, per=5000)
  
  
  fit.t.ECM20[[l]]=MtFA.na.ECM (Y.na20, g, q, eta=0.005, tol=0.0001, max.iter=5000, per=5000)
  
  fit.N.ECM20[[l]]=MFA.na.ECM (Y.na20, g, q, eta=0.005, tol=0.0001, max.iter=5000, per=5000)
  
  ## EM
  
  
  fit.CN.EM20[[l]]=MCNFA.na.EM (Y.na20, g, q, eta=0.005, tol=0.0001, max.iter=5000, per=5000)
  
  
  fit.t.EM20[[l]]=MtFA.na.EM (Y.na20, g, q, eta=0.005, tol=0.0001, max.iter=5000, per=5000)
  
  fit.N.EM20[[l]]=MFA.na.EM (Y.na20, g, q, eta=0.005, tol=0.0001, max.iter=5000, per=5000)
  
  print(l)
}



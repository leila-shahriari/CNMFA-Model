  rm(list = ls())
  library(ghyp)
  
  library(sn)
  library(VGAM)
  library(mclust)
  library(EMMIXskew)
  library(clusteval)
  library(ContaminatedMixt)
  library(mnormt)
  
  
  WD.PATH = paste("D:/Data and code CNMFA")
  
  source(paste(WD.PATH, '/MCNFA.na.EM.r', sep=""))
  source(paste(WD.PATH, '/MCNFA.na.ECM.r', sep=""))
  source(paste(WD.PATH, '/MCNFA.na.AECM.r', sep=""))
  source(paste(WD.PATH, '/MCNFA.na.EM.SE.r', sep=""))
  
  ################important function######################
  ########################################################
  
  
  rCN <- function(n, mu = rep(0,p), Sigma, alpha = 0.99, eta = 1.01){
    
    if(missing(Sigma))
      stop("Sigma is missing")
    if(alpha<0 | alpha>1)
      stop("alpha must be in (0,1)")
    if(eta<1)
      stop("eta must be greater than 1")
    
    p <- if(is.matrix(Sigma)) 
      ncol(Sigma)
    else 1
    
    X    <- array(0,c(n,p),dimnames=list(1:n,paste("X.",1:p,sep="")))
    good <- rbinom(n=n,size=1,prob=alpha)
    for(i in 1:n){
      if(good[i]==1)
        X[i,] <- rmnorm(n = 1, mean = mu, varcov = Sigma)
      else
        X[i,] <- rmnorm(n = 1, mean = mu, varcov = eta*Sigma)
    }
    return(X)  
    
  }
  
  
  mix.CN.gener = function(p, alpha, gam, mu, A, D, sample.size = 1000){
    n = sample.size
    g = length(p)
    y = list(g)
    ni = rmultinom(1, size = n, prob = p)
    class = rep(1:g,c(ni))
    for(i in 1:g){
      Sigma = A[,,i]%*%t(A[,,i]) + diag(D[,i])
      y[[i]] = rCN(ni[i], mu = mu[,i], Sigma = Sigma, alpha = alpha[i], eta = gam[i]) 
    }
    data = y[[1]]
    for(j in 2:g){
      data = rbind(data,y[[j]])
    }
    return(list(data = data, class = class))
    
  }
  
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
  
  
  ###############################################
  ###############################################
  
  
  n = 1600;#100,200,400,800,1600
  pi=w = c(0.3,0.7); mu = cbind(c(-1,1.1,1.3),c(1,-1.1,-1.3))
  D = cbind(c(0.7,0.5,0.6),c(0.7,0.6,0.5))
  A = array(NA,dim = c(3,1,2))
  A[,,1] = cbind(c(0.9,1.1,1.2))
  A[,,2] = cbind(c(0.9,1.2,1.1))
  gam = c(1.2,1.2); alpha = c(0.55,0.65)
  

  
  
  #################################
  #################################
  
  
  fit.AECM = fit.AECM10 = fit.AECM20 = list()
  
  rrr = 100
  
  
  for(rr in 1:rrr){
    repeat
    {
      cat('\nr =', rr,'\n')
  
  ##################
  
  DATA = mix.CN.gener(pi, alpha, gam, mu, A, D, sample.size = n)
  
  
  
  Y = DATA$data
  true.clus=DATA$class
  
  ###########################################################
  
  fit.AECM[[rr]]=try(MCNFA.na.AECM (Y, 2, 1, clus=true.clus, eta=0.005, tol=0.0001, max.iter=2000, per=1000))
  
  
  
  ###########################################################
  Y.na10 = gener.na (Y, 0.1)
  
  fit.AECM10[[rr]]=try(MCNFA.na.AECM (Y.na10, 2, 1, clus=true.clus, eta=0.005, tol=0.0001, max.iter=2000, per=1000))
  
  
  
  ###########################################################
  Y.na20 = gener.na (Y, 0.2)
  
  fit.AECM20[[rr]]=try(MCNFA.na.AECM (Y.na20, 2, 1, clus=true.clus, eta=0.005, tol=0.0001, max.iter=2000, per=1000))
  #if (!inherits(fit.AECM10[[rr]], "try-error")) break
  #if (!inherits(fit.AECM[[rr]], "try-error")) break
  if (!inherits(fit.AECM20[[rr]], "try-error")){
    if (!inherits(fit.AECM10[[rr]], "try-error")){
      if (!inherits(fit.AECM[[rr]], "try-error")) break
    }
  }
    
  }
  
  
  }
  
  
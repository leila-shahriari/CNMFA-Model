rm(list = ls())

library(mclust); library(EMMIXskew)

YY = read.table("D:/Data and code MCNFA/Barbara.txt" ,sep = "",header = F)

YY = as.matrix(YY)



WD.PATH = paste("D:/Data and code MCNFA")

source(paste(WD.PATH, '/function/MCNFA.na.AAECM.r', sep=""))
source(paste(WD.PATH, '/function/MCNFA.na.AECM.r', sep=""))
source(paste(WD.PATH, '/function/MFA.na.AAECM.2018.r', sep=""))
source(paste(WD.PATH, '/function/MFA.na.AECM.2018.r', sep=""))
source(paste(WD.PATH, '/function/MtFA.na.AECM.r', sep=""))
source(paste(WD.PATH, '/function/MtFA.na.AAECM.r', sep=""))


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

g = 4
nn = nrow(YY)
clus = true.clus = kmeans(YY,g)$cluster;eta=0.005





nn = nrow(YY)
true.clus = sample(1:g, nn, replace = T)
clus = true.clus;


fit.N.AECMq2 = fit.N.AECMq3 = fit.N.AECMq4 = fit.N.AECMa5 = fit.N.AECMq6 = fit.N.AECMq7 = fit.N.AECMq8 = list()
fit.CN.AECMq2 = fit.CN.AECMq3 = fit.CN.AECMq4 = fit.CN.AECMq5 = fit.CN.AECMq6 = fit.CN.AECMq7 =  fit.CN.AECMq8 = list()
fit.t.AECMq2 = fit.t.AECMq3 = fit.t.AECMq4 = fit.t.AECMq5 = fit.t.AECMq6 = fit.t.AECMq7 =  fit.t.AECMq8 = list()



for(i in 1:100){
  
  
  Y = gener.na (YY, 0.2)
  
  
  fit.N.AECMq2[[i]] = MFA.na.AECM (Y, g = g, q = 2, clus = true.clus, eta = 0.005, tol = 0.000001, max.iter = 5000, per = 5000)
  fit.N.AECMq3[[i]] = MFA.na.AECM (Y, g = g, q = 3, clus = true.clus, eta = 0.005, tol = 0.000001, max.iter = 5000, per = 5000)
  fit.N.AECMq4[[i]] = MFA.na.AECM (Y, g = g, q = 4, clus = true.clus, eta = 0.005, tol = 0.000001, max.iter = 5000, per = 5000)
  fit.N.AECMq5[[i]] = MFA.na.AECM (Y, g = g, q = 5, clus = true.clus, eta = 0.005, tol = 0.000001, max.iter = 5000, per = 5000)
  fit.N.AECMq6[[i]] = MFA.na.AECM (Y, g = g, q = 6, clus = true.clus, eta = 0.005, tol = 0.000001, max.iter = 5000, per = 5000)
  fit.N.AECMq7[[i]] = MFA.na.AECM (Y, g = g, q = 7, clus = true.clus, eta = 0.005, tol = 0.000001, max.iter = 5000, per = 5000)
  fit.N.AECMq8[[i]] = MFA.na.AECM (Y, g = g, q = 8, clus = true.clus, eta = 0.005, tol = 0.000001, max.iter = 5000, per = 5000)
  
  fit.t.AECMq2[[i]] = MtFA.na.AECM (Y, g = g, q = 2, clus = true.clus, eta = 0.005, tol = 0.000001, max.iter = 5000, per = 5000)
  fit.t.AECMq3[[i]] = MtFA.na.AECM (Y, g = g, q = 3, clus = true.clus, eta = 0.005, tol = 0.000001, max.iter = 5000, per = 5000)
  fit.t.AECMq4[[i]] = MtFA.na.AECM (Y, g = g, q = 4, clus = true.clus, eta = 0.005, tol = 0.000001, max.iter = 5000, per = 5000)
  fit.t.AECMq5[[i]] = MtFA.na.AECM (Y, g = g, q = 5, clus = true.clus, eta = 0.005, tol = 0.000001, max.iter = 5000, per = 5000)
  fit.t.AECMq6[[i]] = MtFA.na.AECM (Y, g = g, q = 6, clus = true.clus, eta = 0.005, tol = 0.000001, max.iter = 5000, per = 5000)
  fit.t.AECMq7[[i]] = MtFA.na.AECM (Y, g = g, q = 7, clus = true.clus, eta = 0.005, tol = 0.000001, max.iter = 5000, per = 5000)
  fit.t.AECMq8[[i]] = MtFA.na.AECM (Y, g = g, q = 8, clus = true.clus, eta = 0.005, tol = 0.000001, max.iter = 5000, per = 5000)
  
  
  
  
  
  fit.CN.AECMq2[[i]] = MCNFA.na.AECM (Y, g = g, q = 2, clus = fit.N.AECMq2[[i]]$post.clus, eta = 0.005, tol = 0.000001, max.iter = 5000, per = 5000,
                                    mu = fit.N.AECMq2[[i]]$para$mu,
                                    S = fit.N.AECMq2[[i]]$para$Sig,
                                    D = fit.N.AECMq2[[i]]$para$D,
                                    A = fit.N.AECMq2[[i]]$para$A)
  
  fit.CN.AECMq3[[i]] = MCNFA.na.AECM (Y, g = g, q = 3, clus = fit.N.AECMq3[[i]]$post.clus, eta = 0.005, tol = 0.000001, max.iter = 5000, per = 5000,
                                    mu = fit.N.AECMq3[[i]]$para$mu,
                                    S = fit.N.AECMq3[[i]]$para$Sig,
                                    D = fit.N.AECMq3[[i]]$para$D,
                                    A = fit.N.AECMq3[[i]]$para$A)
  
  fit.CN.AECMq4[[i]] = MCNFA.na.AECM (Y, g = g, q = 4, clus = fit.N.AECMq4[[i]]$post.clus, eta = 0.005, tol = 0.000001, max.iter = 5000, per = 5000,
                                    mu = fit.N.AECMq4[[i]]$para$mu,
                                    S = fit.N.AECMq4[[i]]$para$Sig,
                                    D = fit.N.AECMq4[[i]]$para$D,
                                    A = fit.N.AECMq4[[i]]$para$A)
  
  fit.CN.AECMq5[[i]] = MCNFA.na.AECM (Y, g = g, q = 5, clus = fit.N.AECMq5[[i]]$post.clus, eta = 0.005, tol = 0.000001, max.iter = 5000, per = 5000,
                                    mu = fit.N.AECMq5[[i]]$para$mu,
                                    S = fit.N.AECMq5[[i]]$para$Sig,
                                    D = fit.N.AECMq5[[i]]$para$D,
                                    A = fit.N.AECMq5[[i]]$para$A)
  
  fit.CN.AECMq6[[i]] = MCNFA.na.AECM (Y, g = g, q = 6, clus = fit.N.AECMq6[[i]]$post.clus, eta = 0.005, tol = 0.000001, max.iter = 5000, per = 5000,
                                    mu = fit.N.AECMq6[[i]]$para$mu,
                                    S = fit.N.AECMq6[[i]]$para$Sig,
                                    D = fit.N.AECMq6[[i]]$para$D,
                                    A = fit.N.AECMq6[[i]]$para$A)
  
  fit.CN.AECMq7[[i]] = MCNFA.na.AECM (Y, g = g, q = 7, clus = fit.N.AECMq7[[i]]$post.clus, eta = 0.005, tol = 0.000001, max.iter = 5000, per = 5000,
                                    mu = fit.N.AECMq7[[i]]$para$mu,
                                    S = fit.N.AECMq7[[i]]$para$Sig,
                                    D = fit.N.AECMq7[[i]]$para$D,
                                    A = fit.N.AECMq7[[i]]$para$A)
  
  fit.CN.AECMq8[[i]] = MCNFA.na.AECM (Y, g = g, q = 8, clus = fit.N.AECMq8[[i]]$post.clus, eta = 0.005, tol = 0.000001, max.iter = 5000, per = 5000,
                                    mu = fit.N.AECMq8[[i]]$para$mu,
                                    S = fit.N.AECMq8[[i]]$para$Sig,
                                    D = fit.N.AECMq8[[i]]$para$D,
                                    A = fit.N.AECMq8[[i]]$para$A)
  
}


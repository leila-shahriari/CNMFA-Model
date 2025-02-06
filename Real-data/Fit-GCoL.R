
rm(list = ls())

library(mclust); library(EMMIXskew)



WD.PATH = paste("D:/Data and code MCNFA")

source(paste(WD.PATH, '/function/MCNFA.na.AECM.r', sep=""))
source(paste(WD.PATH, '/function/MCNFA.na.EM.r', sep=""))

DATA=read.csv("D:/Data and code MCNFA/cost-of-living.csv" ,sep = ",")
classss = DATA[,3]
Y1= (DATA[,1:59])#
#Y2= (DATA[classss=="Germany",1:58])#
#Data = rbind(Y1,Y2)
Y = scale(DATA[,4:58])#
#Y = (scale(Y1[,4:58]))
p=ncol(Y)
# O matrices 
na.posi = is.na(Y)
na.class = colSums(t(na.posi) * 2 ^ (0:(p-1)))     #ind.na
full.missing = which(na.class==sum(rep(2,p)^ (0:(p-1))))

Y = scale(as.matrix(Y[ -full.missing,]))

g = 7
nn = nrow(Y)
true.clus = sample(1:g, nn, replace = T)
clus = true.clus; q=4
#true.clus = as.numeric(Y1[,59])
#true.clus[true.clus == 0] = 2


fit.N.EM = fit.N.ECM = fit.N.AECM = list()
fit.CN.EM = fit.CN.ECM = fit.CN.AECM = list()
fit.t.EM = fit.t.ECM = fit.t.AECM = list()



for(i in 2:16){
fit.N.AECM[[i]] = MFA.na.AECM (Y, 
                               g = g, 
                               q = i, 
                               clus = true.clus, 
                               eta = 0.005, 
                               tol = 0.000001, 
                               max.iter = 5000, 
                               per = 5000)



fit.N.EM[[i]] = MFA.na.EM (Y, g = g, 
                          q = i, 
                          clus = true.clus, 
                          eta = 0.005, 
                          tol = 0.000001, 
                          max.iter = 5000, 
                          per = 5000)

fit.t.EM[[i]] = MtFA.na.EM (Y, g = g, q = i, clus=true.clus, eta=0.005, tol=0.000001, max.iter=5000, per=10)
fit.t.AECM[[i]] = MtFA.na.AECM (Y, g = g, q = i, clus=true.clus, eta=0.005, tol=0.000001, max.iter=5000, per=10)

fit.CN.AECM[[i]] = MCNFA.na.AECM (Y, 
                             g = g, 
                             q = i, 
                             clus = fit.N.AECM[[i]]$post.clus, 
                             eta = 0.005, 
                             tol = 0.000001, 
                             max.iter = 5000, 
                             per = 5000,
                             mu = fit.N.AECM[[i]]$para$mu,
                             S = fit.N.AECM[[i]]$para$Sig,
                             D = fit.N.AECM[[i]]$para$D,
                             A = fit.N.AECM[[i]]$para$A)

fit.CN.EM[[i]] = MCNFA.na.EM (Y, 
                           g = g, 
                           q = i, 
                           clus = fit.N.AECM[[i]]$post.clus, 
                           eta = 0.005, 
                           tol = 0.000001, 
                           max.iter = 5000, 
                           per = 5000,
                         mu = fit.N.AECM[[i]]$para$mu,
                         S = fit.N.AECM[[i]]$para$Sig,
                         D = fit.N.AECM[[i]]$para$D,
                         A = fit.N.AECM[[i]]$para$A)




}


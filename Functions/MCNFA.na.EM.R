
####---------------------------------------------------------###
## ---- Mixture of Contaminated normal factor analysis ------###
## ------------- with missing value -------------------------###
####---------------------------------------------------------###

MCNFA.na.EM = function(Y = NULL, # matrix of data
                       g = NULL, # number of groups
                       q = NULL, # number of factors
                       T.clus = NULL, # True groups Labels
                       clus = NULL,
                       init.meth = c("kmeans", "random"),
                       eta = NULL,
                       tol = NULL, # stopping rule
                       max.iter = NULL, # maximum number of iterations
                       Print = T, 
                       per = NULL,
                       mu = NULL,
                       S = NULL,
                       D = NULL,
                       A = NULL)
{
  packages <- c("cli", "mvtnorm")
  
  # Install packages not yet installed
  installed_packages <- packages %in% rownames(installed.packages())
  if (any(installed_packages == FALSE)) {
    install.packages(packages[!installed_packages])
  }
  
  # Packages loading
  invisible(lapply(packages, library, character.only = TRUE))
  
  massage1 = cli::combine_ansi_styles("red", "bold")
  massage2 = cli::combine_ansi_styles("burlywood1", "bold")
  
  begin = proc.time()[1]
  
  if (is.data.frame(Y))
    Y = as.matrix(Y)
  n = nrow(Y)
  p = ncol(Y)
  iter = 1
  
  if (is.null(Y))
    stop(massage1('Hey, we need some data, please! Y is null'))
  if (!is.matrix(Y))
    stop(massage1('Y needs to be in matrix form'))
  if (!is.numeric(Y))
    stop(massage1('Y is required to be numeric'))
  if (n == 1)
    stop(massage1('number of onservation is too small'))
  if (p == 1)
    stop(massage1('This function is for multivariate data'))
  if (is.null(g) | g < 1)
    stop(massage1('g should be specified correctly'))
  
  ####---------------------------------------------------------###
  ####--------------------- (Required Functions) --------------###
  ####---------------------------------------------------------###
  
  Cluster.error.rate = function(clust1, clust2) {
    clust1 <- unclass(as.ordered(clust1))
    clust2 <- unclass(as.ordered(clust2))
    if ((n = length(clust1)) != length(clust2)) {
      stop(
        massage1(
          "Error in clustering evaluation: the length Groups are not equal."
        ),
        call. = FALSE
      )
    }
    if ((g = length(table(clust1))) != length(table(clust2))) {
      stop(
        massage1(
          "Error in clustering evaluation: the number of clusters are not equal."
        ),
        call. = FALSE
      )
    }
    permute <- function(a) {
      n <- length(a)
      if (n == 1)
        f <- a
      else {
        nm <- gamma(n)
        f <- array(0, c(n, n * nm))
        j <- 1
        for (i in a) {
          f[1, (j - 1) * nm + 1:nm] <- i
          f[-1, (j - 1) * nm + 1:nm] <- permute(setdiff(a, i))
          j <- j + 1
        }
      }
      f
    }
    id <- 1:n
    cmb <- permute(1:g)
    nperm <- ncol(cmb)
    rate <- rep(0, nperm)
    for (i in 1:nperm) {
      tmp <- rep(0, g)
      tc <- rep(0, n)
      for (j in 1:g)
        tc[clust2 == j] = cmb[j, i]
      for (j in 1:g) {
        tmp1 <- 0
        for (k in (1:g)[-j])
          tmp1 <- tmp1 + length(intersect(id[clust1 == j], id[tc == k]))
        tmp[j] <- tmp1
      }
      rate[i] <- sum(tmp) / n
    }
    min(rate)
  }
  
  ARI.fun = function(LabelA, LabelB) {
    u <- unclass(as.ordered(LabelA))
    v <- unclass(as.ordered(LabelB))
    if ((N <- length(u)) != length(v))
      stop(massage1("Error in ARI computing: Labels of the groups do not match!"),
           call. = FALSE)
    row <- max(u)
    col <- max(v)
    nvect <- array(0, c(row, col))
    for (i in 1:row) {
      for (j in 1:col) {
        nvect[i, j] <- sum(u == i & v == j)
      }
    }
    SumsA <- rowSums(nvect)
    SumsB <- colSums(nvect)
    a = 0
    for (i in 1:row)
      a = a + choose(SumsA[i], 2)
    b = 0
    for (j in 1:col)
      b = b + choose(SumsB[j], 2)
    c <- a * b / choose(N, 2)
    d = 0
    for (i in 1:row) {
      for (j in 1:col) {
        d = d + choose(nvect[i, j], 2)
      }
    }
    ind <- (d - c) / ((a + b) / 2 - c)
    ind
  }
  
  ####---------------------------------------------------------###
  ####---- Missing values and parameters initialization -------###
  ####---------------------------------------------------------###
  if (is.null(mu)) {
    mu = matrix(NA, p, g)
    S = array(NA, dim = c(p, p, g))
    D = matrix(NA, p, g)
    A = array(NA, dim = c(p, q, g))
  }
  
  Y.hat = array(NA, dim = c(n, p, g))
  alpha = rep(0.9, g)
  gam = rep(3, g)
  
  na.posi = is.na(Y)
  ind.na.posi = matrix(FALSE, n, p)
  Y.hh = Y
  
  if (is.null(clus)){
    if (init.meth == "kmeans") clus = kmeans(Y)$cluster
    if (init.meth == "random") clus = sample(1:g, n, replace = T) 
  }
  
  for (i in 1:g)
  {
    ind.Y = Y[clus == i, ]
    ind.na.posi[clus == i, ] = na.posi[clus == i, ]
    for (j in 1:p)
    {
      Y.hh[ind.na.posi[, j], j] = colMeans(ind.Y, na.rm = T)[j]
    }
  }
  
  if (is.null(mu)) {
    if (g == 1)
    {
      mu[, 1] = colMeans(Y.hh)
      S[, , 1] = var(Y.hh)
      eig.S = eigen(S[, , 1])
      A[, , 1] = eig.S$ve[, 1:q] %*% diag(sqrt(eig.S$va[1:q]), q, q)
      D[, 1] = diag(S[, , 1] - A[, , 1] %*% t(A[, , 1]))
    }
    else
    {
      for (i in 1:g)
      {
        ind.Y = Y.hh[clus == i,]
        mu[, i] = colMeans(ind.Y)
        S[, , i] = var(ind.Y)
        eig.S = eigen(S[, , i])
        A[, , i] = eig.S$ve[, 1:q] %*% diag(sqrt(eig.S$va[1:q]), q, q)
        D[, i] = diag(S[, , i] - A[, , i] %*% t(A[, , i]))
      }
    }
  }
  
  if (g > 1) w = table(clus) / nrow(Y) else w = 1
  
  ####----------------------------------------###
  ####----------------- O matrix --------------###
  ####---------------------------------------###
  na.class = colSums(t(na.posi) * 2 ^ (0:(p - 1)))     
  uni.na.class = unique(na.class)                    
  num.na.class = length(uni.na.class)           
  O.list = ind.list = as.list(num.na.class)
  for (i in 1:num.na.class)
  {
    ind.list[[i]] = which(na.class == uni.na.class[i])
    O.list[[i]] = matrix(diag(p)[!na.posi[ind.list[[i]][1], ], ], ncol = p)
  }
  
  VI = wden = matrix(NA, n, g)
  
  ####---------------------------------------------------------###
  ####----------------- Observed log-likelihood ---------------###
  ####---------------------------------------------------------###
  
  for (i in 1:g)
  {
    S[, , i] = A[, , i] %*% t(A[, , i]) + diag(D[, i])
    Y.hat[, , i] = Y.hh
    for (j in 1:num.na.class)
    {
      O = O.list[[j]]
      ind = ind.list[[j]]
      Y.pat = matrix(Y.hh[ind, ] , ncol = p)
      Y.o = Y.pat %*% t(O)
      OSO = O %*% S[, , i] %*% t(O)
      mu.o = as.vector(O %*% mu[, i])
      wden[ind, i] = w[i] * (alpha[i] * dmvnorm(Y.o, mu.o, OSO) + 
                               (1 - alpha[i]) * dmvnorm(Y.o, mu.o, gam[i] * OSO))
    }
  }
  indv.den = rowSums(wden)
  indv.den[which(indv.den == 0)] <- .Machine$double.xmin
  logli.old = sum(log(indv.den))
  lk = logli.old
  
  if (isTRUE(Print)) {
    cat(cli::rule(line = 1, line_col = "orange"), "\n")
    cat(cli::rule(
      center = cli::col_green(
        "* Fitting Mixture of Contaminated normal factor analysis (EM) *"
      ),
      line_col = "deeppink3",
      line = 2
    ),
    "\n")
    cat(cli::rule(
      center = cli::col_blue("Initial Log-likelihood = ", lk),
      line_col = "deeppink3",
      line = 2
    ),
    "\n")
    cat(cli::rule(line = 1, line_col = "lightcoral"), "\n")
  }
  
  repeat {
    iter = iter + 1
    ####----------------------------------------###
    ####----------------- E-step --------------###
    ####---------------------------------------###
    
    Z = wden / indv.den 
    ni = colSums(Z)  
    w = ni / n
    
    hu= array(NA, dim = c(q,n,g))
    Y.hat=array(NA, dim = c(p,n,g))
    for (i in 1:g)
    {
      Reta=b.hat =DPS =matrix(NA, p, n)
      u1j=matrix(NA,q,n)
      tau=ga=g1j=g2j=matrix(NA,1,n)
      Phi.temp.Sum=matrix(0,p,p)
      RPsi.sum=matrix(0,p,q)
      hu2.sum=matrix(0,q,q)
      zeta1.sum=rep(0,q)
      Ga = solve(S[, , i]) %*% A[, , i]
      
      for (j in 1:num.na.class)
      {
        ind = ind.list[[j]]
        O = O.list[[j]]
        po=nrow(O)
        OMO= O %*% S[,,i] %*% t(O)
        ODO= O %*% diag(D[,i]) %*% t(O)
        Y.pat = matrix(Y[ind,] , ncol = p)
        Y.o=Y.pat%*%t(O)  
        mu.o=as.vector(O%*%mu[,i])
        
        OMoo=t(O) %*% solve(OMO) %*% O  #Soo# 
        Coo=t(O) %*% solve(ODO) %*% O #Coo#
        DCoo=diag(D[,i])%*%Coo
        IAoo=(diag(p)-DCoo)%*%A[,,i] #Roo
        IDoo=(diag(p)-DCoo)%*%diag(D[,i])
        OSO = O %*% S[, , i] %*% t(O)
        
        
        Y.pat = matrix(Y.hh[ind, ] , ncol = p)
        Y.o = Y.pat %*% t(O)
        mu.o = as.vector(O %*% mu[, i])
        
        Correct = dmvnorm(Y.o, mu.o, gam[i] * OSO, log = T) - 
          dmvnorm(Y.o, mu.o, OSO, log = T)
        VI[ind,i] = (alpha[i]) / (alpha[i] + (1 - alpha[i]) * exp(Correct))
        
        tau[,ind] = (VI[ind, i] + (1 - VI[ind, i]) / gam[i])
        
        Woo = solve(diag(q) + t(A[,,i]) %*% Coo %*% A[,,i]) 
        v.o = t(A[,,i]) %*% Coo %*% (t(Y.pat)-mu[,i])
        
        hu[,ind,i] = Woo %*% (t(t(v.o)))
        u1j[,ind] = Woo %*% (t(t(v.o) * tau[,ind]))
        u2j = Woo %*% (t(t(v.o) ) ) 
        
        if(q==1) 
          hu2=(t(Z[ind,i]*(u1j[,ind]))%*%t(v.o)+ diag(sum(Z[ind,i]),q))%*%Woo
        else 
          hu2=(t(Z[ind,i]*t(u1j[,ind]))%*%t(v.o)+ diag(sum(Z[ind,i]),q))%*%Woo
        hu2.sum=hu2.sum+hu2
        RPsi.temp=IAoo%*%hu2     #A*om
        RPsi.sum=RPsi.sum+RPsi.temp
        Reta[,ind]= DCoo%*%A[,,i]%*%(u1j[,ind])
        b.hat[,ind] = mu[,i] + DCoo %*% (t(Y.pat)- mu[,i])
        if(q==1) 
          DPS[,ind]=IAoo%*%t(Z[ind,i]*(u1j[,ind])) 
        else 
          DPS[,ind]=IAoo%*%t(Z[ind,i]*t(u1j[,ind])) 
        if(q==1) 
          ze.hat = mu[,i] + A[,,i] %*% t(hu[,ind,i])
        else 
          ze.hat = mu[,i] + A[,,i] %*% (hu[,ind,i])
        Y.hat[,ind,i] = ze.hat + DCoo %*% (t(Y.pat)- ze.hat)
        Phi.temp =  sum(Z[ind,i]) * IDoo + IAoo%*%hu2%*%t(IAoo)
        Phi.temp.Sum = Phi.temp.Sum + Phi.temp 
        
      }
      
      
      ####----------------------------------------###
      ####----------------- M-step --------------###
      ####---------------------------------------###
      
      
      mu[,i] = colSums(c(Z[,i]*t(tau))*t(b.hat)- Z[,i]*t(Reta))/sum(Z[,i]*t(tau))
      A[,,i] = ((b.hat - mu[,i])%*%(Z[,i]*t(u1j))+ RPsi.sum)%*% solve(hu2.sum)
      
      A.b = b.hat - mu[,i]
      ttmp=Phi.temp.Sum- RPsi.sum%*%t(A[,,i]) -A[,,i]%*% t(RPsi.sum)+
        A[,,i]%*% hu2.sum%*%t(A[,,i])
      Tmp=ttmp+A.b%*%t(DPS) -A.b%*%t(A[,,i]%*%t(Z[,i]*t(u1j)))+DPS%*%t(A.b)-
        A[,,i]%*%t(Z[,i]*t(u1j))%*%t(A.b)+(A.b)%*%(c(Z[,i]*t(tau))*t(A.b))
      D[,i] = diag(Tmp)/(ni[i]) 
      
      u.cent = A.b - A[, , i] %*% hu[,,i]
      
      alpha[i] = max(0.5, sum(Z[, i] * VI[, i]) / ni[i])
      gam[i] = max(1.001, 
                   sum(Z[, i] * (1 - VI[, i]) * diag(
                     t(u.cent) %*% solve(diag(D[, i])) %*% u.cent)) / 
                     sum(p * Z[, i] * (1 - VI[, i])))
    }
    D[D < eta] = eta
    
    ####---------------------------------------------------------###
    ####-------------- Updated Observed log-likelihood ----------###
    ####---------------------------------------------------------###
    for (i in 1:g)
    {
      S[, , i] = A[, , i] %*% t(A[, , i]) + diag(D[, i])
      for (j in 1:num.na.class)
      {
        O = O.list[[j]]
        ind = ind.list[[j]]
        Y.pat = matrix(Y.hh[ind, ] , ncol = p)
        Y.o = Y.pat %*% t(O)
        OSO = O %*% S[, , i] %*% t(O)
        mu.o = as.vector(O %*% mu[, i])
        wden[ind, i] = w[i] * (alpha[i] * dmvnorm(Y.o, mu.o, OSO) + 
                                 (1 - alpha[i]) * dmvnorm(Y.o, mu.o, gam[i] * OSO))
      }
    }
    indv.den = rowSums(wden)
    indv.den[which(indv.den == 0)] <- .Machine$double.xmin
    logli.new = sum(log(indv.den))
    lk = c(lk, logli.new)
    diff = logli.new - logli.old
    
    if (iter %% per == 0 | is.na(diff))
    {
      if (isTRUE(Print)) {
        cat(cli::rule(
          left = cli::col_cyan(
            'iter = ',
            iter,
            ', logli = ',
            logli.new,
            ", log.Like's diff = ",
            diff
          ),
          line_col = "deeppink3",
          line = 1
        ),
        "\n")
        cat(cli::rule(line = 1, line_col = "lightcoral"), "\n")
      }
    }
    
    if (diff < tol | iter == max.iter)
      break
    logli.old = logli.new
  }
  
  end = proc.time()[1]
  if (isTRUE(Print)) {
    cat(cli::rule(
      center = cli::col_green('** Fiting MCN-FA model takes ', end - begin, ' seconds **'),
      line_col = "deeppink3",
      line = 1
    ),
    "\n")
  }
  Z = wden / indv.den
  post.clus = matrix(apply(Z, 1, order), nrow = g)[g, ]
  
  if(!is.null(T.clus)){
    if (length(unique(post.clus)) == length(unique(T.clus)))
      MCR = Cluster.error.rate(T.clus, post.clus)
    else
      MCR = 1 - sum(apply(table(post.clus, T.clus), 1, max)) / n
    
    ARI = ARI.fun(post.clus, T.clus)
  }else{
    MCR = ARI = NULL
  }
  
  hu.sum = matrix(0,n,q)
  for(i in 1:g)
  {
    
    if(q==1)
      hu.sum=hu.sum+Z[,i]* (hu[,,i])
    else
      hu.sum=hu.sum+Z[,i]* t(hu[,,i])
  }
  
  
  para = list(
    w = w,
    mu = mu,
    A = A,
    D = D,
    Sig = S,
    alpha = alpha,
    gamma = gam
  )
  
  no.para = g * (2 + 2 * p + p * q - q * (q - 1) / 2) + g - 1
  AIC = -2 * logli.new + 2 * no.para
  lk.bic = -2 * lk + no.para * log(n)
  bic2 = -2 * logli.new + no.para * log(n)
  BIC = -bic2 / 2
  abic2 = -2 * logli.new + no.para * log(n / 2 + 12)
  ABIC = -abic2 / 2
  ENT = -sum(Z * log(Z + 1e-10))
  ICL = BIC - ENT
  AWE = ICL - no.para * (1.5 + 0.5 * log(n))
  
  model.inf = c(
    no.para = no.para,
    logli = logli.new,
    AIC = AIC,
    BIC = bic2 ,
    ICL = ICL ,
    AWE = AWE ,
    ABIC = ABIC,
    ARI = ARI,
    MCR = MCR,
    iter = iter
  )
  
  return(
    list(
      logli = logli.new,
      para = para,
      model.inf = model.inf,
      post.clus = post.clus,
      lk = lk.bic,
      cpu = end - begin,
      uhat = hu.sum
    )
  )
}


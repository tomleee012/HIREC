library(MCMCpack)
library(tictoc)
library(quadprog)
library(parallel)
library(doParallel)
# library(HIREewas)
library(mvtnorm)

###Simulate data according to the HIREv model
rm(list = ls())

N = 600
M = 1500
K = 4
Q = 3
seed = 524

#Data Generating Process
pi = c(.5,.3,.2)
#Number of cancer subtypes
X1 = rep(1, N*pi[1])
for(q in 2:Q){
  temp = rep(q, N*pi[q])
  X1 = c(X1, temp)
}
X = array(0, dim = c(N,Q))
for(n in 1:N){
  X[n,X1[n]] = 1
}

P = array(NA, dim = c(N,K), dimnames = list(paste('n =', seq(1,N,1)), paste('k =', seq(1,K,1))))
set.seed(seed)
for(n in 1:N){
  temp = rep(1,K)
  temp[X1[n]] = 9
  P[n,] = rdirichlet(1, temp)
}

#Sparsity is common
Beta = array(0, dim=c(M,K,Q))
M_common = 50;    M_seperate = 50;    base = M_common + M_seperate*K;    max_signal = 0.2;    min_signal = 0.15
set.seed(seed)
for(q in 2:Q){
  max_signal = max_signal*(q-1);    min_signal = min_signal*(q-1)
  signs = array(sample(c(-1,1), M_common*K, replace=TRUE), dim=c(M_common,K))
  Beta[base*(q-2) + 1:M_common,,q] = signs*runif(M_common*K, min=min_signal, max=max_signal)
  for(k in 1:K){
    signs = array(sample(c(-1,1), M_seperate, replace=TRUE))
    Beta[base*(q-2) + M_common + M_seperate*(k-1) + (1:M_seperate),k,q] = signs*runif(M_seperate, min=min_signal, max=max_signal)
  }
}

#Baseline Methylation Profile (Between 0 and 1)
set.seed(seed)
Mu = array(rbeta(M*K, 10, 10), dim = c(M,K), dimnames = list(paste('m =', seq(1,M,1)), paste('k =', seq(1,K,1))))

#Sigma_MK = array(runif(M*K,0,.001), dim = c(M,K), dimnames = list(paste('m =', seq(1,M,1)), paste('k =', seq(1,K,1))))
#Sigma_EpsilonM = array(runif(M,0,.001), dim = c(M), dimnames = list(paste('m =', seq(1,M,1))))

U = array(NA, dim = c(N,M,K), dimnames = list(paste('n =', seq(1,N,1)), paste('m =', seq(1,M,1)), paste('k =', seq(1,K,1))))
O = array(NA, dim = c(N,M), dimnames = list(paste('n =', seq(1,N,1)), paste('m =', seq(1,M,1))))
set.seed(seed)
for(n in 1:N){
  for(m in 1:M){
    for(k in 1:K){
      U[n,m,k] = rnorm(1, Mu[m,k] + Beta[m,k,X1[n]], sqrt(.0001))
    }
    O[n,m] = rnorm(1, U[n,m,] %*% P[n,], sqrt(.0001))
  }
}
print(sum(O<0) + sum(O>1))
O[O > 1] = 1
O[O < 0] = 0

#Check whether the DGP is correct
HIRE <- function(Ometh, X, num_celltype, tol = 10^(-5), num_iter=1000, alpha=0.01){
  
  #generate samples from a Dirichlet distribution
  rDirichlet <- function(alpha_vec){
    num <- length(alpha_vec)
    temp <- rgamma(num, shape = alpha_vec, rate = 1)
    return(temp / sum(temp))
  }
  
  #initialize P_t
  
  CorDescent <- function(MethMatr, num_celltype, tol = 0.01, showIter = FALSE){
    err0 <- 0
    err1 <- 1000
    m <- nrow(MethMatr) #number of CpG sites
    n <- ncol(MethMatr) #number of samples
    K <- num_celltype
    if(m < n){
      stop("The CpG site number must be larger than the sample number!")
    }
    #initialize P_matr
    P_matr_t <- vapply(seq_len(n), function(i){ rDirichlet(rep(2,K)) }, FUN.VALUE=rep(-1,K))
    while(abs(err1 - err0) >= tol){
      err0 <- err1
      #update U_matr
      Dmat <- 2*P_matr_t%*%t(P_matr_t)
      Amat <- cbind(diag(rep(1,K)), diag(rep(-1,K)))
      bvec <- c(rep(0,K), rep(-1,K))
      U_matr_t <- t( vapply(seq_len(m), function(j){
        dvec <- 2*P_matr_t %*% as.numeric(MethMatr[j,])
        
        solu <- solve.QP(Dmat, dvec, Amat, bvec)
        solu$solution
      }, FUN.VALUE=rep(-1,K)) )
      
      #update P_matr
      Dmat <- 2*t(U_matr_t) %*% U_matr_t
      Amat <- cbind(matrix(1, K, K), diag(rep(1,K)))
      bvec <- c(rep(1, K), rep(0, K))
      P_matr_t <- vapply(seq_len(n), function(i){
        dvec <- 2 * t(U_matr_t) %*% as.numeric(MethMatr[ ,i])
        solu <- solve.QP(Dmat, dvec, Amat, bvec, meq = K)
        solu$solution 
      }, FUN.VALUE=rep(-1,K))
      
      #calculate errors
      err1 <- sum((MethMatr - U_matr_t %*% P_matr_t)^2)
      if(showIter == TRUE){
        message("  ", err1, "\n")
      }
    }
    
    return(list(U=U_matr_t, P=P_matr_t))	
  }
  
  Initialize <- function(Ometh, num_celltype){
    K <- num_celltype
    sdrow <- apply(Ometh, 1, sd)
    ind <- order(sdrow, decreasing = TRUE)
    m <- nrow(Ometh)
    if(m <= 1000){
      num_cpg_for_init <- m
    }else{
      num_cpg_for_init <- max(3*ncol(Ometh), floor(m/10))
      if(num_cpg_for_init > m){
        num_cpg_for_init <- m
      }
    }
    
    Ometh_part <- Ometh[ind[seq_len(num_cpg_for_init)],] #select CpG sites with the most num_cpg_for_init variant methylation levels 
    
    result <- CorDescent(Ometh_part, num_celltype=K, tol = 0.1, showIter = FALSE)
    P_initial <- result$P
    
    mu_initial <- vapply(seq_len(m), function(j){
      if(K > 2){
        fit <- lm(Ometh[j,]~as.matrix(t(P_initial[-1, ])))
      }else{
        fit <- lm(Ometh[j,]~as.numeric(P_initial[-1, ]))
      }
      tmp <- as.numeric(summary(fit)$coeff[ ,1])
      tmp[-1] <- tmp[1] + tmp[-1]
      tmp
    }, FUN.VALUE = rep(-1,K) )
    return(list(P_initial, mu_initial))
  }
  
  EmEwasRcallC <- function(Ometh, X, P_t, mu_t, beta_t, sig_sqTiss_t, sig_sqErr_t, tol, num_iter){
    args <- list("P_init"=as.numeric(P_t), "mu_init"=as.numeric(mu_t), "beta_init"=as.numeric(beta_t),
                 "beta_init_dim"=as.integer(dim(beta_t)), "Ometh_r"=as.numeric(Ometh),
                 "Ometh_r_dim"=as.integer(dim(Ometh)), "X_r"=as.numeric(X), "X_r_dim"=as.integer(dim(X)),
                 "sig_sqTiss_init"=as.numeric(sig_sqTiss_t), "sig_sqErr_init"=as.numeric(sig_sqErr_t),
                 "tol_r" = as.numeric(tol), "num_iter" = as.integer(num_iter))
    ret_list <- .Call("EmEwasRcallC", args)		
    return(ret_list)
  }
  
  
  m <- nrow(Ometh) #CpG site number
  n <- ncol(Ometh) #sample number
  p <- nrow(X)
  K <- num_celltype
  
  P_t <- matrix(-1, K, n)
  mu_t <- matrix(-1, m, K)
  beta_t <- array(-1, dim=c(m, K, p))
  sig_sqTiss_t <- matrix(-1, m, K)
  sig_sqErr_t <- rep(-1, m)
  
  init <- Initialize(Ometh, K)
  P_t <- init[[1]]
  mu_t <- t(init[[2]])
  message("  Initialization Done.\n")
  message("  Implementing EM algorithm... \n")
  ret_list <- EmEwasRcallC(Ometh, X, P_t, mu_t, beta_t, sig_sqTiss_t, sig_sqErr_t, tol, num_iter)
  message("  Done! \n")
  
  return(ret_list)
}
#hire = HIRE(t(O), t(X), 3)

#Implementation of the EM Algorithm
rm(list= ls()[!(ls() %in% c('M', 'N', 'K', 'Q', 'O', 'X', 'Beta', 'Mu', 'Pi', 'P'))])

#Step 0
pi = Pi; mu = Mu; beta = Beta; p = P
sigma_mk = array(.001, dim = c(M,K));     sigma_epsilonm = array(.001, dim = c(M))
e_sigma = array(0, dim = c(N,M,K,K));     e_mu = array(0, dim = c(N,M,Q,K))

cl = makeCluster(12)
registerDoParallel(cl)
clusterExport(cl, list("O", "p", 'pi', "mu", 'beta', "sigma_mk", "sigma_epsilonm", "e_mu", 'e_sigma'), envir=environment())

z_ = array(0, dim = c(N,M,Q))
for(q in 1:Q){
  z_[,,q] = foreach(m = 1:M, .combine= cbind) %dopar% {
    sapply(1:N, function(n){pi[q] * dnorm(O[n,m], as.numeric(p[n,] %*% (mu[m,] + beta[m,,q])), sqrt(as.numeric(p[n,] %*% diag(sigma_mk[m,]) %*% p[n,] + sigma_epsilonm[m])))})
  }
}
for(n in 1:N){
  z_[n,,] = t(sapply(1:M, function(m){z_[n,m,] / (sum(z_[n,m,]) + 2.225074e-308)}))
}

stopCluster(cl) #close parallel computing

#Step 1: Choosing the optimal initial values 
# Choosing the optimal initial values
CorDescent = function(MethMatr, num_celltype, ini_P, tol = 0.01){
  MethMatr = t(MethMatr)
  err0 = 0
  err1 = 1000
  m = nrow(MethMatr) #number of CpG sites
  n = ncol(MethMatr) #number of samples
  K = num_celltype
  #initialize P_matr
  P_matr_t = ini_P
  while(abs(err1 - err0) >= tol){
    err0 = err1
    #update U_matr
    Dmat = 2*P_matr_t%*%t(P_matr_t)
    Amat = cbind(diag(rep(1,K)), diag(rep(-1,K)))
    bvec = c(rep(0,K), rep(-1,K))
    U_matr_t = t( vapply(seq_len(m), function(j){
      dvec = 2*P_matr_t %*% as.numeric(MethMatr[j,])
      
      solu = solve.QP(Dmat, dvec, Amat, bvec)
      solu$solution
    }, FUN.VALUE=rep(-1,K)) )
    
    #update P_matr
    Dmat = 2*t(U_matr_t) %*% U_matr_t
    Amat = cbind(matrix(1, K, K), diag(rep(1,K)))
    bvec = c(rep(1, K), rep(0, K))
    P_matr_t = vapply(seq_len(n), function(i){
      dvec = 2 * t(U_matr_t) %*% as.numeric(MethMatr[ ,i])
      solu = solve.QP(Dmat, dvec, Amat, bvec, meq = K)
      solu$solution 
    }, FUN.VALUE=rep(-1,K))
    
    #calculate errors
    err1 = sum((MethMatr - U_matr_t %*% P_matr_t)^2)
  }
  
  mu_initial <- vapply(seq_len(m), function(j){
    if(K > 2){
      fit <- lm(MethMatr[j,]~as.matrix(t(P_matr_t[-1, ])))
    }else{
      fit <- lm(MethMatr[j,]~as.numeric(P_matr_t[-1, ]))
    }
    tmp <- as.numeric(summary(fit)$coeff[ ,1])
    tmp[-1] <- tmp[1] + tmp[-1]
    tmp
  }, FUN.VALUE = rep(-1,K) )
  
  return(list(hat.mu=t(mu_initial), hat.p=t(P_matr_t)))
}

Optimize_p = function(O, p, sigma_epsilonm, z, e_mu, e_sigma){
  for(n in 1:N){
    Dmat = array(0, dim = c(K,K)); dvec = array(0, dim = c(K))
    for(m in 1:M){
      for(q in 1:Q){
        Dmat = Dmat + z[n,m,q]* (outer(e_mu[n,m,q,], e_mu[n,m,q,]) + e_sigma[n,m,,]) / sigma_epsilonm[m]
        dvec = dvec + z[n,m,q]* O[n,m] * e_mu[n,m,q,] / sigma_epsilonm[m]
      }
    }
    Amat = cbind(diag(rep(1,K)), diag(rep(-1,K)), rep(1,K), rep(-1,k))
    bvec = c(rep(0,K), rep(-1,K), 1, -1)
    solu = solve.QP(Dmat, dvec, Amat, bvec)
    p[n,] = solu$solution
  }
  return(p)
}

e_sigma = array(0, dim = c(N,M,K,K));     e_mu = array(0, dim = c(N,M,Q,K))
sigma_mk = array(.0001, dim = c(M,K));     sigma_epsilonm = array(.0001, dim = c(M))
cluster.tilde = as.vector(kmeans(O, Q)$cluster) - 1
pi = as.vector(table(cluster.tilde)) / sum(as.vector(table(cluster.tilde)))
initial.val = CorDescent(O, K, t(rdirichlet(N, rep(1, K))))
p = initial.val$hat.p
mu = initial.val$hat.mu

set.seed(9)
beta = array(NA, dim = c(M, K, Q))
beta[,,1] = 0
for(q in 2:Q) beta[,,q] = runif(M*K, -1, 1)

cl = makeCluster(12)
registerDoParallel(cl)
clusterExport(cl, list("O", "p", 'pi', "mu", 'beta', "sigma_mk", "sigma_epsilonm", "e_mu", 'e_sigma'), envir=environment())

tol = 1;     err0 = 0;     err1 = 1000; num_iter = 1000
# while(abs(err1 - err0) >= tol){
for(t in 1:num_iter){
  err0 = err1
  #E_Step
  w = array(0, dim = c(N,Q));     z = array(0, dim = c(N,M,Q))
  for(q in 1:Q){
    w[,q] = foreach(n = 1:N, .combine= rbind) %dopar% {
      temp = log(pi[q])
      for(m in 1:M){
        temp = temp + log(dnorm(O[n,m], as.numeric(p[n,] %*% (mu[m,] + beta[m,,q])), sqrt(as.numeric(p[n,] %*% diag(sigma_mk[m,]) %*% p[n,] + sigma_epsilonm[m]))) + 2.225074e-308)
      }
      temp
    }
    z[,,q] = foreach(m = 1:M, .combine= cbind) %dopar% {
      sapply(1:N, function(n){pi[q] * dnorm(O[n,m], as.numeric(p[n,] %*% (mu[m,] + beta[m,,q])), sqrt(as.numeric(p[n,] %*% diag(sigma_mk[m,]) %*% p[n,] + sigma_epsilonm[m])))})
    }
  }
  for(n in 1:N){
    z[n,,] = t(sapply(1:M, function(m){z[n,m,] / (sum(z[n,m,]) + 2.225074e-308)}))
  }
  y = exp(w - apply(w, 1, max));     Z = y/apply(y, 1, sum)
  
  e_sigma = array(0, dim = c(N,M,K,K))
  result = foreach(n = 1:N, .multicombine = TRUE) %dopar% {
    temp = t(sapply(1:M, function(m) {solve(solve(diag(sigma_mk[m,])) + outer(p[n,], p[n,]) / sigma_epsilonm[m])}))
    dim(temp) = c(M, K, K)
    temp
  }
  for(n in 1:N){
    e_sigma[n,,,] = result[[n]]
  }
  
  e_mu = array(0, dim = c(N,M,Q,K))
  for(q in 1:Q){
    result = foreach(n = 1:N, .multicombine= TRUE) %dopar% {
      t(sapply(1:M, function(m) {as.numeric(((mu[m,] + beta[m,,q]) %*% solve(diag(sigma_mk[m,])) + (p[n,] * O[n,m] / sigma_epsilonm[m])) %*% e_sigma[n,m,,])}))
    }
    for(n in 1:N){
      e_mu[n,,q,] = result[[n]]
    }
  }
  
  #M_Step
  #update mu, beta
  for(m in 1:M){
    for(k in 1:K){
      if(t >= 500){
        mu[m,k] = mean(z[,m,1] * (e_mu[,m,1,k] - beta[m,k,1]) + z[,m,2] * (e_mu[,m,2,k] - beta[m,k,2]) + z[,m,3] * (e_mu[,m,3,k] - beta[m,k,3]))
      }
      for(q in 2:Q){
        beta[m,k,q] = sum(z[,m,q] * (e_mu[,m,q,k] - mu[m,k])) / (sum(z[,m,q]) + 2.225074e-308)
      }
    }
  }
  
    #update p
    p = Optimize_p(O, p, sigma_epsilonm, z, e_mu, e_sigma)
    
  if(t >= 500){
    #update pi
    pi = apply(Z, 2, mean)
  }
  
  w = foreach(n = 1:N, .combine= rbind) %dopar% {
    temp = 0
    for(m in 1:M){
      temp2 = 0
      for(q in 1:Q){
        temp2 = temp2 + pi[q] * dnorm(O[n,m], as.numeric(p[n,] %*% (mu[m,] + beta[m,,q])), sqrt(as.numeric(p[n,] %*% diag(sigma_mk[m,]) %*% p[n,] + sigma_epsilonm[m])))
      }
      temp = temp + log(temp2 + 2.225074e-308)
    }
    temp
  }
  print(paste0("Log-likelihood: ", sum(w)))
  
  print(pi)
  print(cor(beta[,,2], Beta[,,2]))
  print(cor(beta[,,2], Beta[,,3]))
  print(cor(mu, Mu))
  #print(cor(p, P))
  #print(cor(z, z_))
  #print(Z)
  err1 = sum(w)
}

stopCluster(cl) #close parallel computing

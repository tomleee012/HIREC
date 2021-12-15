setwd("C:/Users/tomle/OneDrive - The Chinese University of Hong Kong/Research/HIRE/HIREewasV/src")
Ometh = O; num_celltype = 3; num_subtype = 3; tol = 10^(-50); num_iter = 1000

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
  
  return(list(mu=t(mu_initial), p=t(P_matr_t)))
}

N <- nrow(Ometh) #CpG site number
M <- ncol(Ometh) #sample number
Q <- num_subtype #cancer number
K <- num_celltype

set.seed(seed+1)
X_til = as.vector(kmeans(Ometh, Q)$cluster);     num_q = unique(X_til);     X_tilde = NA
for(q in 1:Q){
  temp = which(X_til == num_q[q])
  for(ind in temp){
    X_tilde[ind] = q
  }
}
pi = as.vector(table(X_tilde)) / sum(as.vector(table(X_tilde)))
sig_sqTiss_t = array(.0001, dim = c(M,K));     sig_sqErr_t = array(.0001, dim = c(M))

system("R CMD SHLIB Initialize.c");     dyn.load('Initialize.dll')
InitializeRcallC <- function(Ometh, pi, P_t, mu_t, beta_t, sig_sqTiss_t, sig_sqErr_t){
  args <- list("P_init"=as.numeric(P_t), "mu_init"=as.numeric(mu_t), "beta_init"=as.numeric(beta_t),
               "beta_init_dim"=as.integer(dim(beta_t)), "Ometh_r"=as.numeric(Ometh),
               "Ometh_r_dim"=as.integer(dim(Ometh)), "pi_init"=as.numeric(pi), 'loglike_init' = c(0,0),
               "sig_sqTiss_init"=as.numeric(sig_sqTiss_t), "sig_sqErr_init"=as.numeric(sig_sqErr_t))
  ret_list <- .Call("InitializeRcallC", args)		
  return(ret_list)
}

T = 20
loglike = rep(0, T)
for(t in 1:T){
  set.seed(seed+t)
  init_p = array(0, dim = c(K,N))
  for(n in 1:N){
    temp = rep(1, K);     temp[X_tilde[n]] = 10;
    init_p[,n] = rdirichlet(1, temp)
  }
  mu = CorDescent(Ometh[X_tilde == 1,], K, init_p[,X_tilde == 1], tol = .001)$mu
  beta = array(0, dim = c(M,K,Q))
  for(q in 2:Q){
    beta[,,q] = CorDescent(Ometh[X_tilde == q,], K, init_p[,X_tilde == q], tol = .001)$mu - mu
  }
  p = CorDescent(Ometh, K, init_p, tol = .001)$p
  loglike[t] = InitializeRcallC(t(Ometh), pi, t(p), mu, beta, sig_sqTiss_t, sig_sqErr_t)$loglike[1]
  if(which.max(loglike) == t){
    PI = pi;     BETA = beta;     MU = mu;     PP = p
  }
}
message("  Initialization Done.\n")

system("R CMD SHLIB HIREV.c");     dyn.load('HIREV.dll')
EmEwasVRcallC <- function(Ometh, X_tilde, pi, P_t, mu_t, beta_t, sig_sqTiss_t, sig_sqErr_t, tol, num_iter){
  args <- list("P_init"=as.numeric(P_t), "mu_init"=as.numeric(mu_t), "beta_init"=as.numeric(beta_t),
               "beta_init_dim"=as.integer(dim(beta_t)), "Ometh_r"=as.numeric(Ometh),
               "Ometh_r_dim"=as.integer(dim(Ometh)), 'X_tilde' = as.numeric(X_tilde), "pi_init"=as.numeric(pi),
               "sig_sqTiss_init"=as.numeric(sig_sqTiss_t), "sig_sqErr_init"=as.numeric(sig_sqErr_t),
               "tol_r" = as.numeric(tol), "num_iter" = as.integer(num_iter))
  ret_list <- .Call("EmEwasVRcallC", args)		
  return(ret_list)
}
message("  Initialization Done.\n")
message("  Implementing EM algorithm... \n")
tic()
ret_list <- EmEwasVRcallC(t(Ometh), X_tilde, PI, t(PP), MU, BETA, sig_sqTiss_t, sig_sqErr_t, tol, num_iter)
toc()
message("  Done! \n")
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
//#include <ppl.h>
#include "HIREV.h"
//using namespace concurrency;

extern SEXP EmEwasVRcallC(SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"EmEwasVRcallC",          (DL_FUNC) &EmEwasVRcallC,          1},
    {NULL, NULL, 0}
};

void R_init_HIREVewas(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}

SEXP getListElement (SEXP list, char *str) {

	SEXP elmt = R_NilValue, names = getAttrib(list, R_NamesSymbol);
	int i;

	for(i = 0; i < length(list); i++) {
		if(strcmp(CHAR(STRING_ELT(names, i)), str) == 0){
			elmt = VECTOR_ELT(list, i);
			break;
		}
	}

	return elmt;
}

SEXP EmEwasVRcallC(SEXP args) {
	int nProtected = 0;
	double *pi, **P_t,  **mu_t, ***beta_t, *X;
	double *sig_sqErr_t, **sig_sqTiss_t;
	double **Ometh;
	int K, N, M, Q;
	double tol, *val_loglikelihood;
	int num_iteration;

	SEXP list, P_init, mu_init, beta_init, beta_init_dim, sig_sqTiss_init, sig_sqErr_init, Ometh_r, \
		Ometh_r_dim, X_tilde, pi_init, tol_r, num_iter, loglikelihood;

	//input variables from R
	list = args;
	P_init = getListElement(list, "P_init");
	//P_init_dim = getListElement(list, "P_init_dim");
	mu_init = getListElement(list, "mu_init");
	//mu_init_dim = getListElement(list, "mu_init_dim");
	beta_init = getListElement(list, "beta_init");
	beta_init_dim = getListElement(list, "beta_init_dim");
	Ometh_r = getListElement(list, "Ometh_r");
	Ometh_r_dim = getListElement(list, "Ometh_r_dim");
	X_tilde = getListElement(list, "X_tilde");
	pi_init = getListElement(list, "pi_init");
	sig_sqTiss_init = getListElement(list, "sig_sqTiss_init");
	sig_sqErr_init = getListElement(list, "sig_sqErr_init");
	tol_r = getListElement(list, "tol_r");
	num_iter = getListElement(list, "num_iter");
	loglikelihood = getListElement(list, "loglikelihood");


	M = INTEGER(Ometh_r_dim)[0];  // CpG site number
	N = INTEGER(Ometh_r_dim)[1];  // sample number
	K = INTEGER(beta_init_dim)[1];// cell number
	Q = INTEGER(beta_init_dim)[2];      // subtype number
	tol = REAL(tol_r)[0]; 		  // tolerance
	num_iteration = INTEGER(num_iter)[0]; // maximum iteration number

	val_loglikelihood = (double*)malloc(sizeof(double)); 
	val_loglikelihood[0] = REAL(loglikelihood)[0]; 		  // loglikelihood

	//declare and initialize variables in C
	Ometh = make2Darray(M, N);
	X = (double*) malloc(N*sizeof(double));
	pi = (double *) malloc(Q*sizeof(double));
	P_t = make2Darray(K, N);
	mu_t = make2Darray(M, K);
	beta_t = make3Darray(M, K, Q);
	sig_sqTiss_t = make2Darray(M, K);
	sig_sqErr_t = (double *) malloc(M*sizeof(double));

	for(int n = 0; n < N; n++){
		for(int m = 0; m < M; m++)
			Ometh[m][n] = REAL(Ometh_r)[m+M*n];
		X[n] = REAL(X_tilde)[n];
		for(int k = 0; k < K; k++)
			P_t[k][n] = REAL(P_init)[k+K*n];
	}

    for(int q = 0; q < Q; q++)
        pi[q]= REAL(pi_init)[q];

	for(int m = 0; m < M; m++){
		for(int k = 0; k < K; k++)
			mu_t[m][k] = REAL(mu_init)[m+M*k];
	}

	for(int m = 0; m < M; m++){
		sig_sqErr_t[m] = REAL(sig_sqErr_init)[m];
		for(int k = 0; k < K; k++){
			sig_sqTiss_t[m][k] = REAL(sig_sqTiss_init)[m+M*k];
			for(int q = 0; q < Q; q++)
				beta_t[m][k][q] = REAL(beta_init)[m+M*k+M*K*q];
		}
	}

	EmEwasV(P_t, mu_t, beta_t, sig_sqErr_t, sig_sqTiss_t, K, N,
					M, Q, Ometh, X, pi, tol, num_iteration, val_loglikelihood);

	//return P_t
	SEXP ret_P_t, ret_P_t_dim;
	PROTECT(ret_P_t = allocVector(REALSXP, K*N));
	++nProtected;
	for(int k = 0; k < K; k++){
		for(int n = 0; n < N; n++)
			REAL(ret_P_t)[k + K*n] = P_t[k][n];
	}
	PROTECT(ret_P_t_dim = allocVector(INTSXP, 2));
	++nProtected;
	INTEGER(ret_P_t_dim)[0] = K;
	INTEGER(ret_P_t_dim)[1] = N;
	setAttrib(ret_P_t, R_DimSymbol, ret_P_t_dim);

	//return mu_t
	SEXP ret_mu_t, ret_mu_t_dim;
	PROTECT(ret_mu_t = allocVector(REALSXP, M*K));
	++nProtected;
	for(int m = 0; m < M; m++){
		for(int k = 0; k < K; k++)
			REAL(ret_mu_t)[m + M*k] = mu_t[m][k];
	}
	PROTECT(ret_mu_t_dim = allocVector(INTSXP, 2));
	++nProtected;
	INTEGER(ret_mu_t_dim)[0] = M;
	INTEGER(ret_mu_t_dim)[1] = K;
	setAttrib(ret_mu_t, R_DimSymbol, ret_mu_t_dim);

	//return beta_t
	SEXP ret_beta_t, ret_beta_t_dim;
	PROTECT(ret_beta_t = allocVector(REALSXP, M*K*Q));
	++nProtected;
	for(int m = 0; m < M; m++){
		for(int k = 0; k < K; k++){
			for(int q = 0; q < Q; q++)
				REAL(ret_beta_t)[m + M*k + M*K*q] = beta_t[m][k][q];
		}
	}
	PROTECT(ret_beta_t_dim = allocVector(INTSXP, 3));
	++nProtected;
	INTEGER(ret_beta_t_dim)[0] = M;
	INTEGER(ret_beta_t_dim)[1] = K;
	INTEGER(ret_beta_t_dim)[2] = Q;
	setAttrib(ret_beta_t, R_DimSymbol, ret_beta_t_dim);

    //return pi
    SEXP ret_pi;
    PROTECT(ret_pi = allocVector(REALSXP, Q));
    ++nProtected;
    for(int q = 0; q < Q; q++)
        REAL(ret_pi)[q] = pi[q];

	//return X
	SEXP ret_X;
	PROTECT(ret_X = allocVector(REALSXP, N));
	++nProtected;
	for (int n = 0; n < N; n++)
		REAL(ret_X)[n] = X[n];

	//return log likelihood
	SEXP ret_loglike;
	PROTECT(ret_loglike = allocVector(REALSXP, 1));
	++nProtected;
	REAL(ret_loglike)[0] = val_loglikelihood[0];

    /*
	//return sig_sqErr_t
	SEXP ret_sig_sqErr_t;
	PROTECT(ret_sig_sqErr_t = allocVector(REALSXP, M));
	++nProtected;
	for(int m = 0; m < M; m++){
		REAL(ret_sig_sqErr_t)[m] = sig_sqErr_t[m];
	}

	//return sig_sqTiss_t
	SEXP ret_sig_sqTiss_t, ret_sig_sqTiss_t_dim;
	PROTECT(ret_sig_sqTiss_t = allocVector(REALSXP, M*K));
	++nProtected;
	for(int m = 0; m < M; m++){
		for(int k = 0; k < K; k++){
			REAL(ret_sig_sqTiss_t)[m+M*k] = sig_sqTiss_t[m][k];
		}
	}
	PROTECT(ret_sig_sqTiss_t_dim = allocVector(INTSXP, 2));
	++nProtected;
	INTEGER(ret_sig_sqTiss_t_dim)[0] = M;
	INTEGER(ret_sig_sqTiss_t_dim)[1] = K;
	setAttrib(ret_sig_sqTiss_t, R_DimSymbol, ret_sig_sqTiss_t_dim);
	*/

	//return list
	SEXP retList;
	PROTECT(retList = allocVector(VECSXP, 6));
	++nProtected;

	SEXP names; // components names in the return list
	PROTECT(names = allocVector(STRSXP, 6));
	++nProtected;

	SET_STRING_ELT(names, 0, mkChar("pi"));
	SET_STRING_ELT(names, 1, mkChar("P_t"));
	SET_STRING_ELT(names, 2, mkChar("mu_t"));
	SET_STRING_ELT(names, 3, mkChar("beta_t"));
	SET_STRING_ELT(names, 4, mkChar("label"));
	SET_STRING_ELT(names, 5, mkChar("loglike"));
	//SET_STRING_ELT(names, 4, mkChar("sig_sqErr_t"));
	//SET_STRING_ELT(names, 5, mkChar("sig_sqTiss_t"));

	//add elements to the return list
	SET_VECTOR_ELT(retList, 0, ret_pi);
	SET_VECTOR_ELT(retList, 1, ret_P_t);
	SET_VECTOR_ELT(retList, 2, ret_mu_t);
	SET_VECTOR_ELT(retList, 3, ret_beta_t);
	SET_VECTOR_ELT(retList, 4, ret_X);
	SET_VECTOR_ELT(retList, 5, ret_loglike);
	//SET_VECTOR_ELT(retList, 4, ret_sig_sqErr_t);
	//SET_VECTOR_ELT(retList, 5, ret_sig_sqTiss_t);
	setAttrib(retList, R_NamesSymbol, names);

	UNPROTECT(nProtected);
	free(pi);
	free(X);
	delet2Darray(P_t, K, N);
	delet2Darray(mu_t, M, K);
	delet2Darray(Ometh, M, N);
	delet2Darray(sig_sqTiss_t, M, K);
	free(sig_sqErr_t);
	delet3Darray(beta_t, M, K, Q);

	return retList;
}
//========================================================================================================================
//Define multiple dimensional arrays (Memory is not continuous here)
//========================================================================================================================

double**** make4Darray(int a1, int a2, int a3, int a4) {
	double**** tmp;
	tmp = (double****)malloc(a1 * sizeof(double***));
	for (int i = 0; i < a1; i++) {
		tmp[i] = (double***)malloc(a2 * sizeof(double**));
		for (int j = 0; j < a2; j++) {
			tmp[i][j] = (double**)malloc(a3 * sizeof(double*));
			for (int k = 0; k < a3; k++) {
				tmp[i][j][k] = (double*)malloc(a4 * sizeof(double));
			}
		}
	}
	return tmp;
}

void delet4Darray(double**** tmp, int a1, int a2, int a3, int a4) {
	for (int i = 0; i < a1; i++) {
		for (int j = 0; j < a2; j++) {
			for (int k = 0; k < a3; k++) {
				free(tmp[i][j][k]);
			}
			free(tmp[i][j]);
		}
		free(tmp[i]);
	}
	free(tmp);
}

double*** make3Darray(int a1, int a2, int a3) {
	double*** tmp;
	tmp = (double***)malloc(a1 * sizeof(double**));
	for (int i = 0; i < a1; i++) {
		tmp[i] = (double**)malloc(a2 * sizeof(double*));
		for (int j = 0; j < a2; j++) {
			tmp[i][j] = (double*)malloc(a3 * sizeof(double));
		}
	}
	return tmp;
}

void delet3Darray(double*** tmp, int a1, int a2, int a3) {
	for (int i = 0; i < a1; i++) {
		for (int j = 0; j < a2; j++) {
			free(tmp[i][j]);
		}
		free(tmp[i]);
	}
	free(tmp);
}

double** make2Darray(int a1, int a2) {
	double** tmp;
	tmp = (double**)malloc(a1 * sizeof(double*));
	for (int i = 0; i < a1; i++) {
		tmp[i] = (double*)malloc(a2 * sizeof(double));
	}
	return tmp;
}

void delet2Darray(double** tmp, int a1, int a2) {
	for (int i = 0; i < a1; i++) {
		free(tmp[i]);
	}
	free(tmp);
}


//========================================================================================================================
//calculate the inverse of a matrix (the code is from http://www.cs.rochester.edu/u/brown/Crypto/assts/projects/adj.html)
//========================================================================================================================

/*Recursive definition of determinate using expansion by minors*/
double Determinant(double** a, int n) {
	int i, j, j1, j2;
	double det = 0;
	double** m = NULL;

	if (n < 1) { /* Error */

	}
	else if (n == 1) { /* Shouldn't get used */
		det = a[0][0];
	}
	else if (n == 2) {
		det = a[0][0] * a[1][1] - a[1][0] * a[0][1];
	}
	else {
		det = 0;
		for (j1 = 0; j1 < n; j1++) {
			m = (double**)malloc((n - 1) * sizeof(double*));
			for (i = 0; i < n - 1; i++)
				m[i] = (double*)malloc((n - 1) * sizeof(double));
			for (i = 1; i < n; i++) {
				j2 = 0;
				for (j = 0; j < n; j++) {
					if (j == j1)
						continue;
					m[i - 1][j2] = a[i][j];
					j2++;
				}
			}
			det += pow(-1.0, j1 + 2.0) * a[0][j1] * Determinant(m, n - 1);
			for (i = 0; i < n - 1; i++)
				free(m[i]);
			free(m);
		}
	}
	return(det);
}

/*Find the cofactor matrix of a square matrix*/

void CoFactor(double** a, int n, double** b) {
	int i, j, ii, jj, i1, j1;
	double det;
	double** c;

	c = (double**)malloc((n - 1) * sizeof(double*));
	for (i = 0; i < n - 1; i++)
		c[i] = (double*)malloc((n - 1) * sizeof(double));

	for (j = 0; j < n; j++) {
		for (i = 0; i < n; i++) {

			/* Form the adjoint a_ij */
			i1 = 0;
			for (ii = 0; ii < n; ii++) {
				if (ii == i)
					continue;
				j1 = 0;
				for (jj = 0; jj < n; jj++) {
					if (jj == j)
						continue;
					c[i1][j1] = a[ii][jj];
					j1++;
				}
				i1++;
			}

			/* Calculate the determinate */
			det = Determinant(c, n - 1);

			/* Fill in the elements of the cofactor */
			b[i][j] = pow(-1.0, i + j + 2.0) * det;
		}
	}
	for (i = 0; i < n - 1; i++)
		free(c[i]);
	free(c);
}

/*Transpose of a square matrix, do it in place*/
void Transpose(double** a, int n) {
	int i, j;
	double tmp;

	for (i = 1; i < n; i++) {
		for (j = 0; j < i; j++) {
			tmp = a[i][j];
			a[i][j] = a[j][i];
			a[j][i] = tmp;
		}
	}
}

/*calculate the inverse*/
void inverse(double** a, int n, double** a_inv) {
	double det;
	double** cofac_a;
	cofac_a = (double**)malloc(n * sizeof(double*));
	for (int i = 0; i < n; i++) {
		cofac_a[i] = (double*)malloc(n * sizeof(double));
	}
	CoFactor(a, n, cofac_a);
	Transpose(cofac_a, n); //turn the cofacotor matrix into the adjoint matrix
	det = Determinant(a, n);
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			a_inv[i][j] = cofac_a[i][j] / det;
		}
	}

	for (int i = 0; i < n; i++) {
		free(cofac_a[i]);
	}
	free(cofac_a);
}

//========================================================================================================================
//quadratic programming
//========================================================================================================================
void quadprog(double** Dmat, double* dvec, double* xvec, int K) { //K is the dimension of xvec
	if (K == 1) {
		xvec[0] = 1;
	}
	else {
		double** Dmat_inv, * x_star, s1, s2;
		double Dmat_inv_sum = 0, x_star_sum = 0, * Dmat_inv_rowsum;
		int num_negatives = 0;
		double* ind_negative;
		Dmat_inv = (double**)malloc(K * sizeof(double*));
		x_star = (double*)malloc(K * sizeof(double));
		Dmat_inv_rowsum = (double*)malloc(K * sizeof(double));
		ind_negative = (double*)malloc(K * sizeof(double));
		for (int k = 0; k < K; k++) {
			Dmat_inv[k] = (double*)malloc(K * sizeof(double));
		}
		inverse(Dmat, K, Dmat_inv);
		for (int k = 0; k < K; k++) {
			s1 = 0;
			s2 = 0;
			for (int k1 = 0; k1 < K; k1++) {
				s1 += Dmat_inv[k][k1] * dvec[k1];
				Dmat_inv_sum += Dmat_inv[k][k1];
				s2 += Dmat_inv[k][k1];
			}
			x_star[k] = s1;
			Dmat_inv_rowsum[k] = s2;
			x_star_sum += s1;
		}
		for (int k = 0; k < K; k++) {
			xvec[k] = x_star[k] + (1 - x_star_sum) / Dmat_inv_sum * Dmat_inv_rowsum[k];
			if (xvec[k] < 0) {
				num_negatives++;
				ind_negative[k] = 1;
			}
			else {
				ind_negative[k] = 0;
			}
		}
		free(x_star);
		free(Dmat_inv_rowsum);
		for (int k = 0; k < K; k++) {
			free(Dmat_inv[k]);
		}
		free(Dmat_inv);

		if (num_negatives == 0) {
			free(ind_negative);
		}
		else {
			int Knew = K - num_negatives, i, j;
			double** Dmat_new, * dvec_new, * xvec_sub;
			Dmat_new = (double**)malloc((Knew) * sizeof(double*));
			for (int k = 0; k < Knew; k++) {
				Dmat_new[k] = (double*)malloc((Knew) * sizeof(double));
			}
			dvec_new = (double*)malloc((Knew) * sizeof(double));
			xvec_sub = (double*)malloc((Knew) * sizeof(double));
			i = 0;

			for (int k1 = 0; k1 < K; k1++) {
				if (ind_negative[k1] == 0) {
					dvec_new[i] = dvec[k1];
					j = 0;
					for (int k2 = 0; k2 < K; k2++) {
						if (ind_negative[k2] == 0) {
							Dmat_new[i][j] = Dmat[k1][k2];
							j++;
						}
					}
					i++;
				}
				else {
					xvec[k1] = 0;
				}
			}

			quadprog(Dmat_new, dvec_new, xvec_sub, Knew);
			i = 0;
			for (int k = 0; k < K; k++) {
				if (ind_negative[k] == 0) {
					xvec[k] = xvec_sub[i];
					i++;
				}
			}
			free(dvec_new);
			free(xvec_sub);
			for (int k = 0; k < Knew; k++) {
				free(Dmat_new[k]);
			}
			free(Dmat_new);
		}

	}

}

//the objective function when updating P_t
double val2(double** P_t, double* sig_sqErr_t, int num_celltypes, int num_CpGsites, double** Ometh,
	double*** E_z, double**** E_Sigma, double**** E_mu, int n, int num_subtypes) {
	int K, M, Q;
	K = num_celltypes;
	M = num_CpGsites;
	Q = num_subtypes;

	double s1, s2, sum = 0;

	for (int m = 0; m < M; m++) {
		s1 = 0;
		s2 = 0;
		for (int k = 0; k < K; k++) {
			for (int k2 = 0; k2 < K; k2++) {
				for (int q = 0; q < Q; q++) {
					s1 += E_z[n][m][q] * P_t[k][n] * E_Sigma[n][m][k][k2] * P_t[k2][n];
				}
			}
			for (int q = 0; q < Q; q++) {
				s2 += E_z[n][m][q] * E_mu[n][m][q][k] * P_t[k][n];
			}
		}
		s2 = Ometh[m][n] - s2;
		s2 = s2 * s2;
		sum += (s1 + s2) / (2 * sig_sqErr_t[m]);
	}

	return sum;
}

double val3(double** P_t, double* sig_sqErr_t, int num_celltypes, int num_CpGsites, double** Ometh,
	double*** E_z, double**** E_Sigma, double**** E_mu, int n, int num_subtypes) {
	int K, M, Q;
	K = num_celltypes;
	M = num_CpGsites;
	Q = num_subtypes;

	double s1, s2, sum = 0;

	for (int m = 0; m < M; m++) {
		s1 = 0;
		s2 = 0;
		for (int k = 0; k < K; k++) {
			for (int k2 = 0; k2 < K; k2++) {
				s1 += P_t[k][n] * E_Sigma[n][m][k][k2] * P_t[k2][n];
			}
			for (int q = 0; q < Q; q++) {
				s2 += E_z[n][m][q] * E_mu[n][m][q][k] * P_t[k][n];
			}
		}
		s2 = Ometh[m][n] - s2;
		s2 = s2 * s2;
		sum += (s1 + s2) / (2 * sig_sqErr_t[m]);
	}

	return sum;
}

//========================================================================================================================
//generate random numbers
//========================================================================================================================

//generate one random number between 0 and 1 with 0 and 1 exclusive
double r_unif() {
	double temp;
	do {
		temp = (double)rand() / (double)RAND_MAX;
	} while (temp >= 1 || temp <= 0);
	return temp;
}

//generate one noraml random number via Box-Muller method
double rnormal(double mu, double sd) {
	double temp[2];
	temp[0] = r_unif();
	temp[1] = r_unif();
	double standardnormal = sqrt(-2 * log(temp[0])) * cos(2 * M_PI * temp[1]);
	return standardnormal * sd + mu;
}

//generate one gamma random number 
double rgamma(double shape, double rate) {
	double scale = 1.0 / rate;
	int shape_int = floor(shape);
	double s = 0;
	for (int i = 0; i < shape_int; i++) {
		s = s - log(r_unif());
	}

	if (shape_int == shape) {
		return scale * s;
	}
	else {
		double U, V, W, xi, eta;
		double delta = shape - shape_int;
		do {
			U = r_unif();
			V = r_unif();
			W = r_unif();
			if (U <= exp(1) / (exp(1) + delta)) {
				xi = pow(V, 1 / delta);
				eta = W * pow(xi, delta - 1);
			}
			else {
				xi = 1 - log(V);
				eta = W * exp(-xi);
			}
		} while (eta > pow(xi, delta - 1) * exp(-xi));
		return scale * (xi + s);
	}

}

//generate one inverse gamma random number
double inv_rgamma(double shape, double rate) {
	double temp = rgamma(shape, rate);
	return 1.0 / temp;
}

//generate a beta random number
double rbeta(double shape1, double shape2) {
	double s1, s2;
	s1 = rgamma(shape1, 1);
	s2 = rgamma(shape2, 1);
	return s1 / (s1 + s2);
}

 //========================================================================================================================
 //functions in the algorithm
 //========================================================================================================================

double normal_density(double x, double mean, double sd) {
	return exp(-pow(x - mean, 2) / (2 * sd * sd)) / sqrt(2 * M_PI * sd * sd);
}

double normal_density_log(double x, double mean, double sd) {
	return -0.5 * log(2 * M_PI) - log(sd) - pow(x - mean, 2) / (2 * sd * sd);
}

// no use without lasso penalty
double positive_part(double x) {
	if (x >= 0) {
		return x;
	}
	else {
		return 0;
	}
}

double absolute(double x) {
	if (x >= 0) {
		return x;
	}
	else {
		return -x;
	}
}

double sign(double x) {
	if (x > 0) {
		return 1.0;
	}
	else if (x == 0) {
		return 0.0;
	}
	else {
		return -1.0;
	}
}

//calculate the observed data likelihood
double observed_log_likelihood(double* pi, double** P_t, double** mu_t, double*** beta_t, double* sig_sqErr_t,
	double** sig_sqTiss_t, int num_celltypes, int num_samples,
	int num_CpGsites, int num_subtypes, double** Ometh) {
	int K = num_celltypes;
	int N = num_samples;
	int M = num_CpGsites;
	int Q = num_subtypes;
	double loglike = 0;
	for (int m = 0; m < M; m++) {
		for (int n = 0; n < N; n++) {
			double s4 = 0;
			for (int q = 0; q < Q; q++) {
				double s1 = 0;
				double s3 = 0;
				for (int k = 0; k < K; k++) {
					s1 += P_t[k][n] * (mu_t[m][k] + beta_t[m][k][q]);
					s3 += P_t[k][n] * P_t[k][n] * sig_sqTiss_t[m][k];
				}
				s4 += pi[q] * normal_density(Ometh[m][n], s1, sqrt(s3 + sig_sqErr_t[m]));
			}
			loglike += log(s4 + 5.56268e-309);
		}
	}
	return loglike;
}

//conduct E step
void Estep(double* pi, double** P_t, double** mu_t, double*** beta_t, double* sig_sqErr_t,
	double** sig_sqTiss_t, int num_celltypes, int num_samples,
	int num_CpGsites, int num_subtypes, double** Ometh,
	double**** E_Sigma, double**** E_mu, double*** E_z, double** E_Z, double* X) {
	int K, N, M, Q;
	K = num_celltypes;
	N = num_samples;
	M = num_CpGsites;
	Q = num_subtypes;
	double g, g_tmp, s1, s2, s3, sum, u, w;
	double* tmp1 = (double*)malloc(K * sizeof(double));
	double* tmp2 = (double*)malloc(K * sizeof(double));

	for (int q = 0; q < Q; q++) {
		for (int n = 0; n < N; n++) {
			s1 = log(pi[q] + 5.56268e-309);
			for (int m = 0; m < M; m++) {
				s2 = 0, s3 = 0;
				for (int k = 0; k < K; k++) {
					s2 += P_t[k][n] * (mu_t[m][k] + beta_t[m][k][q]);
					s3 += P_t[k][n] * P_t[k][n] * sig_sqTiss_t[m][k];
				}
				s1 += normal_density_log(Ometh[m][n], s2, sqrt(s3 + sig_sqErr_t[m]));
				E_z[n][m][q] = pi[q] * normal_density(Ometh[m][n], s2, sqrt(s3 + sig_sqErr_t[m]));
			}
			E_Z[n][q] = s1;
		}
	}

	double** vec_X;
	vec_X = make2Darray(N, Q);

	for (int n = 0; n < N; n++) {
		sum = 0;
		w = -100000;
		for (int q = 0; q < Q; q++) {
			if (E_Z[n][q] > w) {
				w = E_Z[n][q];
			}
		}
		for (int q = 0; q < Q; q++) {
			E_Z[n][q] = exp(E_Z[n][q] - w);
			sum += E_Z[n][q];
		}
		for (int q = 0; q < Q; q++) {
			E_Z[n][q] = E_Z[n][q] / sum;
		}

		for (int m = 0; m < M; m++) {
			u = 5.56268e-309;
			for (int q = 0; q < Q; q++) {
				u += E_z[n][m][q];
			}
			for (int q = 0; q < Q; q++) {
				E_z[n][m][q] = E_z[n][m][q] / u;
			}

			s1 = 0;
			for (int k = 0; k < K; k++)
				s1 += P_t[k][n] * P_t[k][n] * sig_sqTiss_t[m][k];
			g = s1 / sig_sqErr_t[m];
			for (int k1 = 0; k1 < K; k1++) {
				for (int k2 = k1; k2 < K; k2++) {
					s2 = 1.0 / (1.0 + g) * P_t[k1][n] * P_t[k2][n] * sig_sqTiss_t[m][k1] * sig_sqTiss_t[m][k2];
					g_tmp = s2 / sig_sqErr_t[m];
					if (k2 != k1) {
						E_Sigma[n][m][k1][k2] = (-g_tmp);
						E_Sigma[n][m][k2][k1] = E_Sigma[n][m][k1][k2];
					}
					else
						E_Sigma[n][m][k1][k2] = sig_sqTiss_t[m][k2] - g_tmp;
				}
			}

			for (int q = 0; q < Q; q++) {
				for (int k = 0; k < K; k++) {
					tmp1[k] = Ometh[m][n] * P_t[k][n] / sig_sqErr_t[m];
					tmp2[k] = (mu_t[m][k] + beta_t[m][k][q]) / sig_sqTiss_t[m][k];
				}

				for (int k1 = 0; k1 < K; k1++) {
					s1 = 0;
					for (int k2 = 0; k2 < K; k2++)
						s1 += (tmp1[k2] + tmp2[k2]) * E_Sigma[n][m][k2][k1];
					E_mu[n][m][q][k1] = s1;
				}
			}
		}
		//update label X
		/*
		
		for (int q = 0; q < Q; q++) {
			vec_X[n][q] = 0;
			for (int m = 0; m < M; m++) {
				vec_X[n][q] += E_z[n][m][q];
			}
		}
		
		*/

		w = -100000;
		for (int q = 0; q < Q; q++) {
			if (E_Z[n][q] > w) {
				w = E_Z[n][q];
				X[n] = q;
			}
		}
	}

	free(tmp1);
	free(tmp2);
	delet2Darray(vec_X, N, Q);
}

// conduct M step
void Mstep(double* pi, double** P_t, double** mu_t, double*** beta_t, double* sig_sqErr_t,
	double** sig_sqTiss_t, int num_celltypes, int num_samples,
	int num_CpGsites, int num_subtypes, double** Ometh,
	double**** E_Sigma, double**** E_mu, double*** E_z, double** E_Z, int t, int num_iteration) {
	int K, N, M, Q;
	K = num_celltypes;
	N = num_samples;
	M = num_CpGsites;
	Q = num_subtypes;
	double s1, s2, s3, s4, sum;
	//=================================
	// use coordinate descent
	//=================================

	/*
	
	//update pi
	for (int q = 0; q < Q; q++) {
		s1 = 5.56268e-309;
		for (int n = 0; n < N; n++) {
			s1 += E_Z[n][q];
		}
		pi[q] = s1 / N;
	}
	
	*/

	//update pi
	double* vec_pi;
	vec_pi = (double*)malloc(Q * sizeof(double));
	sum = 0;
	for (int q = 0; q < Q; q++) {
		vec_pi[q] = 0;
		for (int n = 0; n < N; n++) {
			for (int m = 0; m < M; m++) {
				vec_pi[q] += E_z[n][m][q];
				sum += E_z[n][m][q];
			}
		}
	}
	for (int q = 0; q < Q; q++) {
		pi[q] = vec_pi[q] / M / N;
	}
	free(vec_pi);

	//use coordinate descent algorithm to update mu_t and beta_t
	for (int m = 0; m < M; m++) {
		for (int k = 0; k < K; k++) {

			for (int q = 1; q < Q; q++) {
				s1 = 0, s3 = 5.56268e-309;
				for (int n = 0; n < N; n++) {
					s1 += E_z[n][m][q] * (E_mu[n][m][q][k] - mu_t[m][k]);
					s3 += E_z[n][m][q];
				}
				beta_t[m][k][q] = s1 / s3;
			}

			s1 = 0;
			for (int n = 0; n < N; n++) {
				for (int q = 0; q < Q; q++)
					s1 += E_z[n][m][q] * (E_mu[n][m][q][k] - beta_t[m][k][q]);
			}
			if (t > num_iteration/2) {
				mu_t[m][k] = s1 / N;
				if (mu_t[m][k] > 1)
					mu_t[m][k] = 1;
				else if (mu_t[m][k] < 0)
					mu_t[m][k] = 0;
			}
		}
	}

	//==============================
	//update P_t
	//==============================
	if (t > num_iteration / 2) {
		double* dvec, ** Dmat, ** Dmat_inv, * xvec, * xvec0, tar2a, tar2b;
		xvec = (double*)malloc(K * sizeof(double));
		xvec0 = (double*)malloc(K * sizeof(double));
		dvec = (double*)malloc(K * sizeof(double));
		Dmat = (double**)malloc(K * sizeof(double*));
		Dmat_inv = (double**)malloc(K * sizeof(double*));
		for (int k = 0; k < K; k++) {
			Dmat[k] = (double*)malloc(K * sizeof(double));
			Dmat_inv[k] = (double*)malloc(K * sizeof(double));
		}
		for (int n = 0; n < N; n++) {
			for (int k = 0; k < K; k++) {
				for (int k2 = k; k2 < K; k2++) {
					s1 = 0;
					for (int m = 0; m < M; m++) {
						for (int q = 0; q < Q; q++) {
							s1 += E_z[n][m][q] * ((E_Sigma[n][m][k][k2] + E_mu[n][m][q][k] * E_mu[n][m][q][k2]) / sig_sqErr_t[m]);
						}
					}
					Dmat[k][k2] = s1;
					if (k2 > k) {
						Dmat[k2][k] = Dmat[k][k2];
					}
				}

				s2 = 0;
				for (int m = 0; m < M; m++) {
					for (int q = 0; q < Q; q++) {
						s2 += E_z[n][m][q] * Ometh[m][n] * E_mu[n][m][q][k] / sig_sqErr_t[m];
					}
				}
				dvec[k] = s2;
			}

			for (int k = 0; k < K; k++) {
				xvec0[k] = P_t[k][n];
			}

			tar2a = val3(P_t, sig_sqErr_t, K, M, Ometh, E_z, E_Sigma, E_mu, n, Q);

			quadprog(Dmat, dvec, xvec, K);


			for (int k = 0; k < K; k++) {
				P_t[k][n] = xvec[k];
			}
			tar2b = val3(P_t, sig_sqErr_t, K, M, Ometh, E_z, E_Sigma, E_mu, n, Q);

			if (tar2b > tar2a) { //avoid offset effects when estimates are stable
				for (int k = 0; k < K; k++) {
					P_t[k][n] = xvec0[k];
				}
			}


		}
		free(xvec);
		free(xvec0);
		free(dvec);
		for (int k = 0; k < K; k++) {
			free(Dmat[k]);
			free(Dmat_inv[k]);
		}
		free(Dmat);
		free(Dmat_inv);
	}

	//==============================
	//update sig_sqErr_t
	//==============================
	for (int m = 0; m < M; m++) {
		sum = 0;
		s3 = 0;
		for (int n = 0; n < N; n++) {
			s1 = 0;
			for (int k1 = 0; k1 < K; k1++) {
				for (int k2 = 0; k2 < K; k2++) {
					s1 += P_t[k1][n] * E_Sigma[n][m][k1][k2] * P_t[k2][n];
				}
			}
			s2 = 0;
			for (int q = 0; q < Q; q++) {
				s3 += E_z[n][m][q];
				for (int k = 0; k < K; k++) {
					s2 += E_mu[n][m][q][k] * P_t[k][n];
				}
				s2 = pow(Ometh[m][n] - s2, 2);
				sum += E_z[n][m][q] * (s1 + s2);
			}
		}
		sig_sqErr_t[m] = sum / s3;
	}

	//==============================
	//update sig_sqTiss_t
	//==============================

	for (int m = 0; m < M; m++) {
		for (int k = 0; k < K; k++) {
			s1 = 0;
			sum = 0;
			for (int n = 0; n < N; n++) {
				for (int q = 0; q < Q; q++) {
					s1 += E_z[n][m][q] * (pow(E_mu[n][m][q][k] - mu_t[m][k] - beta_t[m][k][q], 2) + E_Sigma[n][m][k][k]);
					sum += E_z[n][m][q];
				}
			}
			sig_sqTiss_t[m][k] = s1 / sum;
		}
	}

}

//carry out the EM algorithm
void EmEwasV(double** P_t, double** mu_t, double*** beta_t, double* sig_sqErr_t,
	double** sig_sqTiss_t, int num_celltypes, int num_samples,
	int num_CpGsites, int num_subtypes, double** Ometh, double* X, double* pi,
	double tol, int num_iteration, double* val_loglikelihood) { //tol is the tolerance deciding when the algorithm stops

	int t = 1;
	int K = num_celltypes;
	int N = num_samples;
	int M = num_CpGsites;
	int Q = num_subtypes;

	double**** E_Sigma;
	E_Sigma = make4Darray(N, M, K, K);
	double**** E_mu;
	E_mu = make4Darray(N, M, Q, K);
	double*** E_z;
	E_z = make3Darray(N, M, Q);
	double** E_Z;
	E_Z = make2Darray(N, Q);

	double ploglike0 = -100000;
	// Rprintf("Start implementing EM Algorithm... \n");

	// while ((absolute(*val_loglikelihood - ploglike0) / absolute(ploglike0) > tol) && (t <= num_iteration)) { //relative tolerance
	while ((absolute(*val_loglikelihood - ploglike0) > tol) && (t <= num_iteration)) {
		ploglike0 = *val_loglikelihood;

		Estep(pi, P_t, mu_t, beta_t, sig_sqErr_t, sig_sqTiss_t, K, N, M, Q, Ometh, E_Sigma, E_mu, E_z, E_Z, X);
		Mstep(pi, P_t, mu_t, beta_t, sig_sqErr_t, sig_sqTiss_t, K, N, M, Q, Ometh, E_Sigma, E_mu, E_z, E_Z, t, num_iteration);

		*val_loglikelihood = observed_log_likelihood(pi, P_t, mu_t, beta_t, sig_sqErr_t, sig_sqTiss_t, K, N, M, Q, Ometh);
		Rprintf("Iteration: %d\t observed-data log likelihood: %lf\n", t, *val_loglikelihood);
		if (*val_loglikelihood < ploglike0) {
			printf("\t\t\t\t\t warning!!!\n");
		}
		t++;
	}

	delet4Darray(E_Sigma, N, M, K, K);
	delet4Darray(E_mu, N, M, Q, K);
	delet3Darray(E_z, N, M, Q);
	delet2Darray(E_Z, N, Q);
}
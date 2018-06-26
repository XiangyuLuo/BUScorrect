#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <Rmath.h>

#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <stdlib.h>

extern SEXP fdrDEindicator_c(SEXP);
extern SEXP observed_log_likelihood_c(SEXP);
extern SEXP update_alpha_c(SEXP);
extern SEXP update_gamma_c(SEXP);
extern SEXP update_L_c(SEXP);
extern SEXP update_mu_c(SEXP);
extern SEXP update_sigma_sq_c(SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"fdrDEindicator_c",          (DL_FUNC) &fdrDEindicator_c,          1},
    {"observed_log_likelihood_c", (DL_FUNC) &observed_log_likelihood_c, 1},
    {"update_alpha_c",            (DL_FUNC) &update_alpha_c,            1},
    {"update_gamma_c",            (DL_FUNC) &update_gamma_c,            1},
    {"update_L_c",                (DL_FUNC) &update_L_c,                1},
    {"update_mu_c",               (DL_FUNC) &update_mu_c,               1},
    {"update_sigma_sq_c",         (DL_FUNC) &update_sigma_sq_c,         1},
    {NULL, NULL, 0}
};

void R_init_BUScorrect(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}



//Functions 
double max_c (double*vec, int n_vec){
	double s=vec[0];
	int i;
	for(i=1; i<n_vec; i++){
		if(vec[i]>s){
			s = vec[i];
		}
	}
	return s;
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

SEXP update_L_c (SEXP args) {
	int nProtected = 0;
	int kk, gg, K, G;
	double prob[2],  r;
	double DE_prop_t, tau_mu_zero_t, tau_mu_one_t;
	SEXP ret_L_matr, ret_L_matr_dim, list, mu_t, mu_t_dim;

	list = args;
	mu_t = getListElement(list, "mu_t");
	mu_t_dim = getListElement(list, "mu_t_dim");
	DE_prop_t = REAL(getListElement(list, "DE_prop_t"))[0];
	tau_mu_zero_t = REAL(getListElement(list, "tau_mu_zero_t"))[0];
	tau_mu_one_t = REAL(getListElement(list, "tau_mu_one_t"))[0];
	
	G = INTEGER(mu_t_dim)[0]; 
	K = INTEGER(mu_t_dim)[1]; 

	PROTECT( ret_L_matr = allocVector(REALSXP, G*(K-1)) ); /*ret_L_matr is a G*(K-1) matrix*/
	++nProtected;

	GetRNGstate();
	for(kk = 0; kk < K-1; ++kk) {
		for(gg = 0; gg < G; ++gg) {
			prob[0] = DE_prop_t*dnorm(REAL(mu_t)[gg + G*(kk+1)], 0, tau_mu_one_t, 0);
			prob[1] = (1-DE_prop_t)*dnorm(REAL(mu_t)[gg + G*(kk+1)], 0, tau_mu_zero_t, 0);
			r = runif(0,1);
			if(r < prob[0]/(prob[0] + prob[1])) {
				REAL(ret_L_matr)[gg + G*kk] = 1;
			}else {
				REAL(ret_L_matr)[gg + G*kk] = 0;
			}
		}
	}
  
	PutRNGstate();
	

	PROTECT(ret_L_matr_dim = allocVector(INTSXP, 2));
	++nProtected;
	INTEGER(ret_L_matr_dim)[0] = G;
	INTEGER(ret_L_matr_dim)[1] = K - 1;
	setAttrib(ret_L_matr, R_DimSymbol, ret_L_matr_dim);

	UNPROTECT(nProtected);
	return ret_L_matr;
}





SEXP update_alpha_c (SEXP args) {
	/*return a G dimensional vector*/
	int nProtected = 0, G, B;
	int g, b, j, *n_vec;
	SEXP list, Y, nn_vec, mu_t, mu_t_dim, gamma_t, gamma_t_dim;
	SEXP Z_t, sigma_sq_t, tau_alpha;
	SEXP ret_alpha_vec, eta_alpha;

	list = args;
	Y = getListElement(list, "Y");
	nn_vec = getListElement(list, "n_vec");
	mu_t = getListElement(list, "mu_t");
	mu_t_dim = getListElement(list, "mu_t_dim");
	gamma_t = getListElement(list, "gamma_t");
	gamma_t_dim = getListElement(list, "gamma_t_dim");
	Z_t = getListElement(list, "Z_t");
	sigma_sq_t = getListElement(list, "sigma_sq_t");
	tau_alpha = getListElement(list, "tau_alpha");
	eta_alpha = getListElement(list, "eta_alpha");	

	G = INTEGER(mu_t_dim)[0];
	B = INTEGER(gamma_t_dim)[0];
	

	n_vec = (int *) malloc(B * sizeof(int));

	for(b = 0; b < B; ++b){
		n_vec[b] = INTEGER(nn_vec)[b];
	}
	

	PROTECT(ret_alpha_vec = allocVector(REALSXP, G));
	++nProtected;
    double tau_alpha_dt = REAL(tau_alpha)[0];
	double tau_alpha2 = pow(tau_alpha_dt, 2);


	GetRNGstate();
    double s1, s2, mean, sd;
	double *Y_vec;
	int *Z_vec;
	for(g = 0; g < G; ++g) {
		s1 = 0;
        s2 = 0;
	
		for(b = 0; b < B; ++b) {

			Y_vec = REAL(VECTOR_ELT(Y, b));
			Z_vec = INTEGER(VECTOR_ELT(Z_t, b));
			for(j = 0; j < n_vec[b]; ++j){
				s1 = s1 + (Y_vec[j + n_vec[b]*g] - REAL(mu_t)[g + G*(Z_vec[j]-1)] - REAL(gamma_t)[b + B*g])/\
					REAL(sigma_sq_t)[b + B*g];
				
			}

			s2 = s2 + n_vec[b]/REAL(sigma_sq_t)[b + B*g];
		}
			
			
		mean = (s1*tau_alpha2 + REAL(eta_alpha)[0])/(s2*tau_alpha2 + 1);
		sd = sqrt(tau_alpha2/(1+s2*tau_alpha2));
		REAL(ret_alpha_vec)[g] = rnorm(mean, sd);
	}
 		
	PutRNGstate();	

	free(n_vec);		
			
	UNPROTECT(nProtected);
	return ret_alpha_vec;
}
	



SEXP update_mu_c(SEXP args) {
	int nProtected = 0;
	int B, G, K;
	int b, g, k, j;
	double tau_mu_zero_t, tau_mu_one_t, tau_mu2; 
	SEXP list, Y, L_t, alpha_t, gamma_t, Z_t, sigma_sq_t, n_vec;
	SEXP mu, mu_dim;
	
	list = args;
	Y = getListElement(list, "Y");
	L_t = getListElement(list, "L_t");
	alpha_t = getListElement(list, "alpha_t");
	gamma_t = getListElement(list, "gamma_t");
	Z_t = getListElement(list, "Z_t");
	sigma_sq_t = getListElement(list, "sigma_sq_t");

	n_vec = getListElement(list, "n_vec");

	B = INTEGER(getListElement(list, "B"))[0];	 
	G = INTEGER(getListElement(list, "G"))[0];
	K = INTEGER(getListElement(list, "K"))[0];



	tau_mu_zero_t = REAL(getListElement(list, "tau_mu_zero_t"))[0];
	tau_mu_one_t = REAL(getListElement(list, "tau_mu_one_t"))[0];

	PROTECT(mu = allocVector(REALSXP, G*K));
	++nProtected;


	/*The subtype effect of the first group is kept to be zero*/
	for(g = 0; g < G; ++g) {
		REAL(mu)[g] = 0;
	}	

	GetRNGstate();			
    int num_k;
    int *Z_vec;
	double *Y_vec;
	double s1, s2, mean, sd;
	for(k = 1; k < K; ++k) {
		for(g = 0; g < G; ++g) {
			s1 = 0, s2 = 0; 

			if(INTEGER(L_t)[g + G*(k - 1)] == 1) {
				for(b = 0; b < B; ++b) {
					num_k = 0;
					Y_vec = REAL(VECTOR_ELT(Y, b));
					Z_vec = INTEGER(VECTOR_ELT(Z_t, b));
					
					for(j = 0; j < INTEGER(n_vec)[b]; ++j){
						if(Z_vec[j] == k + 1) {
							s1 = s1 + (Y_vec[j + INTEGER(n_vec)[b]*g] - REAL(alpha_t)[g] -\
								REAL(gamma_t)[b + B*g])/REAL(sigma_sq_t)[b + B*g];
							num_k = num_k + 1;
						}
						
					}
	
					s2 = s2 + num_k/REAL(sigma_sq_t)[b + B*g];

				}	
				
				tau_mu2 = pow(tau_mu_one_t,2);
				mean = tau_mu2*s1/(s2*tau_mu2 + 1);
				sd = sqrt(tau_mu2/(1 + s2*tau_mu2));
				REAL(mu)[g + G*k] = rnorm(mean, sd);

			}else {
				for(b = 0; b < B; ++b) {
					num_k = 0;
					Y_vec = REAL(VECTOR_ELT(Y, b));
					Z_vec = INTEGER(VECTOR_ELT(Z_t, b));
					
					for(j = 0; j < INTEGER(n_vec)[b]; ++j){
						if(Z_vec[j] == k + 1) {
							s1 = s1 + (Y_vec[j + INTEGER(n_vec)[b]*g] - REAL(alpha_t)[g] -\
								REAL(gamma_t)[b + B*g])/REAL(sigma_sq_t)[b + B*g];
							num_k += 1;
						}
						
					}
	
					s2 = s2 + num_k/REAL(sigma_sq_t)[b + B*g];

				}	
				
				tau_mu2 = pow(tau_mu_zero_t,2);
				mean = tau_mu2*s1/(s2*tau_mu2 + 1);
				sd = sqrt(tau_mu2/(1 + s2*tau_mu2));
				REAL(mu)[g + G*k] = rnorm(mean, sd);

			}
		}

	}
				
	
	PutRNGstate();	


	PROTECT(mu_dim = allocVector(INTSXP, 2));
	++nProtected;

	INTEGER(mu_dim)[0] = G;
	INTEGER(mu_dim)[1] = K;
	setAttrib(mu, R_DimSymbol, mu_dim);

	UNPROTECT(nProtected);

	return mu;

}




SEXP update_gamma_c(SEXP args) {
	int nProtected = 0;
	int B, G;
	int b, g, j;
	double tau_gamma2;
	SEXP list, Y, alpha_t, mu_t, Z_t, sigma_sq_t, n_vec;
	SEXP gamma, gamma_dim, tau_gamma;
	
	list = args;
	Y = getListElement(list, "Y");
	alpha_t = getListElement(list, "alpha_t");
	mu_t = getListElement(list, "mu_t");
	Z_t = getListElement(list, "Z_t");
	sigma_sq_t = getListElement(list, "sigma_sq_t");
	tau_gamma = getListElement(list, "tau_gamma");
	n_vec = getListElement(list, "n_vec");

	B = INTEGER(getListElement(list, "B"))[0];	 
	G = INTEGER(getListElement(list, "G"))[0];



	PROTECT(gamma = allocVector(REALSXP, G*B));
	++nProtected;

	

	/*The batch effect of batch one is assumed to be zero*/
	

	for(g = 0; g < G; ++g) {
		REAL(gamma)[0 + B*g] = 0;
	}
		
	GetRNGstate();	
	int *Z_vec;
	double *Y_vec;
	double s1, sigma_sq_t_temp, mean, sd;
	double tau_gamma_dt = REAL(tau_gamma)[0];
	for(b = 1; b < B; ++b) {

		Y_vec = REAL(VECTOR_ELT(Y, b));
		Z_vec = INTEGER(VECTOR_ELT(Z_t, b));

		for(g = 0; g < G; ++g) {
            s1 = 0;

			for(j = 0; j < INTEGER(n_vec)[b]; ++j) {
				s1 = s1 + (Y_vec[j + INTEGER(n_vec)[b]*g] - REAL(alpha_t)[g] - \
					REAL(mu_t)[g + G*(Z_vec[j]-1)])/REAL(sigma_sq_t)[b + B*g];
			}

			tau_gamma2 = pow(tau_gamma_dt, 2);
			sigma_sq_t_temp = REAL(sigma_sq_t)[b + B*g];
			mean = tau_gamma2*s1/(INTEGER(n_vec)[b]*tau_gamma2/sigma_sq_t_temp + 1);
			sd = sqrt( tau_gamma2/(1+ INTEGER(n_vec)[b]*tau_gamma2/sigma_sq_t_temp) );

			REAL(gamma)[b + B*g] = rnorm(mean, sd);
		}
	}
			

	PutRNGstate();	


	PROTECT(gamma_dim = allocVector(INTSXP, 2));
	++nProtected;

	INTEGER(gamma_dim)[0] = B;
	INTEGER(gamma_dim)[1] = G;
	setAttrib(gamma, R_DimSymbol, gamma_dim);

	UNPROTECT(nProtected);

	return gamma;

}

SEXP update_sigma_sq_c(SEXP args) {
	int nProtected = 0;
	int B, G;
	int b, g, j;

	SEXP list, Y, alpha_t, mu_t, Z_t, gamma_t, n_vec;
	SEXP sigma_sq, sigma_sq_dim, a_inv_gamma, b_inv_gamma;
	
	list = args;
	Y = getListElement(list, "Y");
	alpha_t = getListElement(list, "alpha_t");
	mu_t = getListElement(list, "mu_t");
	Z_t = getListElement(list, "Z_t");
	gamma_t = getListElement(list, "gamma_t");
	a_inv_gamma = getListElement(list, "a_inv_gamma");
	b_inv_gamma = getListElement(list, "b_inv_gamma");

	n_vec = getListElement(list, "n_vec");

	B = INTEGER(getListElement(list, "B"))[0];	 
	G = INTEGER(getListElement(list, "G"))[0];



	PROTECT(sigma_sq = allocVector(REALSXP, B*G));
	++nProtected;

	
		
	GetRNGstate();	
	int *Z_vec;
	double *Y_vec;
	double s1, shape, rate, temp;
	for(b = 0; b < B; ++b) {

		Y_vec = REAL(VECTOR_ELT(Y, b));
		Z_vec = INTEGER(VECTOR_ELT(Z_t, b));

		for(g = 0; g < G; ++g) {
			s1 = 0;
			for(j = 0; j < INTEGER(n_vec)[b]; ++j) {
				s1 = s1 + pow(Y_vec[j + INTEGER(n_vec)[b]*g] - REAL(alpha_t)[g] -\
					REAL(mu_t)[g + G*(Z_vec[j] - 1)] - REAL(gamma_t)[b + B*g] ,2);
			}
			
			shape = REAL(a_inv_gamma)[0] + INTEGER(n_vec)[b]/2.0;
			rate = REAL(b_inv_gamma)[0] + 1/2.0*s1;
			temp = rgamma(shape, 1/rate);
			REAL(sigma_sq)[b + B*g] = 1/temp;
		}
	}


	PROTECT(sigma_sq_dim = allocVector(INTSXP, 2));
	++nProtected;

	INTEGER(sigma_sq_dim)[0] = B;
	INTEGER(sigma_sq_dim)[1] = G;
	setAttrib(sigma_sq, R_DimSymbol, sigma_sq_dim);

	UNPROTECT(nProtected);

	return sigma_sq;

}




SEXP observed_log_likelihood_c (SEXP args) {
	/*return the log observed-data likelihood value*/
	int nProtected = 0, G, B, K;
	int g, b, k, j, *n_vec;
	SEXP list, Y, nn_vec, pi_t_post, alpha_t_post, mu_t_post, mu_t_dim; 
	SEXP gamma_t_post, gamma_t_dim, sigma_sq_t_post;
	SEXP ret_value;

	list = args;
	Y = getListElement(list, "Y");
	nn_vec = getListElement(list, "n_vec");
	pi_t_post = getListElement(list, "pi_t_post");
	alpha_t_post = getListElement(list, "alpha_t_post");
	mu_t_post = getListElement(list, "mu_t_post");
	mu_t_dim = getListElement(list, "mu_t_dim");
	gamma_t_post = getListElement(list, "gamma_t_post");
	gamma_t_dim = getListElement(list, "gamma_t_dim");
	sigma_sq_t_post = getListElement(list, "sigma_sq_t_post");	

	G = INTEGER(mu_t_dim)[0];
	K = INTEGER(mu_t_dim)[1];
	B = INTEGER(gamma_t_dim)[0];

	n_vec = (int *) malloc(B * sizeof(int));

	for(b = 0; b < B; ++b){
		n_vec[b] = INTEGER(nn_vec)[b];
	}
	
	PROTECT(ret_value = allocVector(REALSXP, 1));
	++nProtected;

    double s=0, mean, sd;
    double *s2, s3;
    s2 = (double*) malloc(K * sizeof(double));
	double *Y_vec, vmax, s4;
	for(b=0; b < B; b++){
		Y_vec = REAL(VECTOR_ELT(Y, b));
		for(j=0; j < n_vec[b]; j++){
			for(k=0; k < K; k++){
				s3=0;
				for(g=0; g < G; g++){
					mean = REAL(alpha_t_post)[g]+REAL(mu_t_post)[g+G*k]+REAL(gamma_t_post)[b+B*g];
					sd = sqrt(REAL(sigma_sq_t_post)[b+B*g]);
					s3 = s3 + dnorm(Y_vec[j+n_vec[b]*g], mean, sd, 1); //1 means log=TRUE					
				}
				s2[k]=s3;
			}
			vmax = max_c(s2, K);
			for(k=0; k < K; k++){
				s2[k]=s2[k]-vmax;
			}
			s4 = 0;
			for(k=0;k < K; k++){
				s4 = s4+REAL(pi_t_post)[b+B*k]*exp(s2[k]);
			}
			s = s + log(s4) + vmax;
		}
	}
	REAL(ret_value)[0]=s;
	free(n_vec);
	free(s2);		
			
	UNPROTECT(nProtected);
	return ret_value;
}	

SEXP fdrDEindicator_c (SEXP args) {
	/*return the fdr value when the tunning threshold is kappa*/
	int nProtected = 0, G, K, num_rec;
	int g, k, j;
	double kappa;
	SEXP list, L_PosterSamp; 
	SEXP fdr;

	list = args;
	G = INTEGER(getListElement(list, "G"))[0];
	K = INTEGER(getListElement(list, "K"))[0];
	num_rec = INTEGER(getListElement(list, "num_rec"))[0];
	kappa = REAL(getListElement(list, "kappa"))[0];
	L_PosterSamp = getListElement(list, "L_PosterSamp");
	
	PROTECT(fdr = allocVector(REALSXP, 1));
	++nProtected;
	
	double denom=0, numer=0, s;
	double ppi, xi;
	for(k=0; k<K-1; k++){
		for(g=0; g<G; g++){
			s=0.0;
			for(j=0; j<num_rec; j++){
				if(INTEGER(L_PosterSamp)[g+G*k+G*(K-1)*j]==(int) 1){
					s = s + 1;
				}
			}
			ppi = s / num_rec;
			xi = 1-ppi;
			if(xi <= kappa){
				denom = denom + 1;
				numer = numer + xi;
			}
		}
	}	
	if(denom==0){
		REAL(fdr)[0]=0;
	}else{
		REAL(fdr)[0]=numer/denom;
	}
	
	UNPROTECT(nProtected);
	return fdr;
}








 
	
	

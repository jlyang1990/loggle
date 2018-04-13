#include <stdlib.h>
#include <stdio.h>
#include <R.h>
#include <Rmath.h>
#include <math.h>
#include <R_ext/Lapack.h>
#include "loggle_admm.h"


//apply ADMM on lambda grid given h and d
void ADMM_lambda(double *Corr, double *Z, int *P, int *Len, int *Pos, int *Pos_Len, int *member_ind, int *csize_ind, 
	int *No, double *Lambda, int *Lambda_Len, int *fit_type, double *thres, double *Rho, double *Epi_abs, 
	double *Epi_rel, int *Max_step){

	int i, j, k, m, p = *P, L = *Len, Pos_L = *Pos_Len, Lambda_L = *Lambda_Len, no_sum = 0, S_L = 0;

	double *Z_i = (double *) malloc(p*p*L*sizeof(double));
	double *U_i = (double *) malloc(p*p*L*sizeof(double));

	for(i=0; i<L; i++){
		for(j=0; j<p; j++){
			for(k=0; k<p; k++){
				Z_i[p*p*i+p*j+k] = 0;
				U_i[p*p*i+p*j+k] = 0;
			}
		}
	}

	//iteration across lambda
	for(i=(Lambda_L-1); i>=0; i--){

		//meet stopping criterion in grid search across lambda
		if(S_L > ((*thres)*p)){
			Z[p*p*Pos_L*(i+1)-1] = -1;
			continue;
		}

		ADMM_cluster(Corr, Z_i, U_i, P, Len, member_ind+p*i, csize_ind+no_sum, No+i, Lambda+i, fit_type, Rho+i, 
			Epi_abs, Epi_rel, Max_step);

		S_L = 0;
		for(j=0; j<p; j++){
			for(k=0; k<p; k++){
				if(j<k && Z_i[p*p*Pos[0]+p*j+k] != 0){
					S_L++;
				}
				for(m=0; m<Pos_L; m++){
					Z[p*p*Pos_L*i+p*p*m+p*j+k] = Z_i[p*p*Pos[m]+p*j+k];
				}
			} 
		}

		no_sum += (No[i]+1);
	}//end iteration across lambda

	free(Z_i);
	free(U_i);
}


//apply ADMM across block diagonals given h, d and lambda
void ADMM_cluster(double *Corr, double *Z, double *U, int *P, int *Len, int *member_ind, int *csize_ind, int *No, 
	double *Lambda, int *fit_type, double *Rho, double *Epi_abs, double *Epi_rel, int *Max_step){

	int p = *P, no = *No, L = *Len, p_n, n, i, j, k;
	int *member_ind_n;

	//iteration across block diagonals
	for(n=0; n<no; n++){

		p_n = csize_ind[n+1] - csize_ind[n];
		member_ind_n = &member_ind[csize_ind[n]];

		//block diagonal with dimension one
		if(p_n == 1){
			if(*fit_type == 0){
				for(i=0; i<L; i++){
					Z[p*p*i+p*member_ind_n[0]+member_ind_n[0]] = 1/Corr[p*p*i+p*member_ind_n[0]+member_ind_n[0]];
				}
			}
		}
		//block diagonal with dimension larger than one
		else{
			double *Corr_n = (double *) malloc(p_n*p_n*L*sizeof(double));
			double *Z_n = (double *) malloc(p_n*p_n*L*sizeof(double));
			double *U_n = (double *) malloc(p_n*p_n*L*sizeof(double));
		
			for(i=0; i<L; i++){
				for(j=0; j<p_n; j++){
					for(k=0; k<p_n; k++){
						Corr_n[p_n*p_n*i+p_n*j+k] = Corr[p*p*i+p*member_ind_n[j]+member_ind_n[k]];
						Z_n[p_n*p_n*i+p_n*j+k] = Z[p*p*i+p*member_ind_n[j]+member_ind_n[k]];
						U_n[p_n*p_n*i+p_n*j+k] = U[p*p*i+p*member_ind_n[j]+member_ind_n[k]];
					}
				}
			}

			//model fitting
			if(*fit_type == 0){
				ADMM_local_glasso(Corr_n, Z_n, U_n, &p_n, Len, Lambda, Rho, Epi_abs, Epi_rel, Max_step);
			}
			else if(*fit_type == 1){
				ADMM_pseudo_glasso(Corr_n, Z_n, U_n, &p_n, Len, Lambda, Rho, Epi_abs, Epi_rel, Max_step);
			}
			else if(*fit_type == 2){
				double *d_n = (double *) malloc(p_n*L*sizeof(double));
				for(i=0; i<L; i++){
					for(j=0; j<p_n; j++){
						d_n[p_n*i+j] = Corr_n[p_n*p_n*i+p_n*j+j];
					}
				}
				for(i=0; i<3; i++){
					ADMM_SPACE_rho(Corr_n, d_n, Z_n, U_n, &p_n, Len, Lambda, Rho, Epi_abs, Epi_rel, Max_step);
					if(i<2){
						ADMM_SPACE_d(Corr_n, d_n, Z_n, &p_n, Len);
					}
				}
				free(d_n);
			}

			for(i=0; i<L; i++){
				for(j=0; j<p_n; j++){
					for(k=0; k<p_n; k++){
						Z[p*p*i+p*member_ind_n[j]+member_ind_n[k]] = Z_n[p_n*p_n*i+p_n*j+k];
						U[p*p*i+p*member_ind_n[j]+member_ind_n[k]] = U_n[p_n*p_n*i+p_n*j+k];
					}
				}
			}
			
			free(Corr_n);
			free(Z_n);
			free(U_n);		
		}
	}//end iteration across block diagonals
}


//apply ADMM across block diagonals given h, d and lambda (no initialization of Z and U)
void ADMM_simple(double *Corr, double *Z, int *P, int *Len, int *Pos, int *Pos_Len, int *member_ind, int *csize_ind, 
	int *No, double *Lambda, int *fit_type, double *Rho, double *Epi_abs, double *Epi_rel, int *Max_step){
  
	int p = *P, no = *No, L = *Len, Pos_L = *Pos_Len, p_n, n, i, j, k;
	int *member_ind_n;

	//iteration across block diagonals
	for(n=0; n<no; n++){
    
		p_n = csize_ind[n+1] - csize_ind[n];
		member_ind_n = &member_ind[csize_ind[n]];
    
		//block diagonal with dimension one
		if(p_n==1){
			if(*fit_type == 0){
				for(i=0; i<Pos_L; i++){
					Z[p*p*i+p*member_ind_n[0]+member_ind_n[0]] = 1/Corr[p*p*Pos[i]+p*member_ind_n[0]+member_ind_n[0]];
				}
			}
		}
		//block diagonal with dimension larger than one
		else{
			double *Corr_n = (double *) malloc(p_n*p_n*L*sizeof(double));
			double *Z_n = (double *) malloc(p_n*p_n*L*sizeof(double));
			double *U_n = (double *) malloc(p_n*p_n*L*sizeof(double));
      
			for(i=0; i<L; i++){
				for(j=0; j<p_n; j++){
					for(k=0; k<p_n; k++){
						Corr_n[p_n*p_n*i+p_n*j+k] = Corr[p*p*i+p*member_ind_n[j]+member_ind_n[k]];
						Z_n[p_n*p_n*i+p_n*j+k] = 0;
						U_n[p_n*p_n*i+p_n*j+k] = 0;
					}
				}
			}
      
			//model fitting
			if(*fit_type == 0){
				ADMM_local_glasso(Corr_n, Z_n, U_n, &p_n, Len, Lambda, Rho, Epi_abs, Epi_rel, Max_step);
			}
			else if(*fit_type == 1){
				ADMM_pseudo_glasso(Corr_n, Z_n, U_n, &p_n, Len, Lambda, Rho, Epi_abs, Epi_rel, Max_step);
			}
			else if(*fit_type == 2){
				double *d_n = (double *) malloc(p_n*L*sizeof(double));
				for(i=0; i<L; i++){
					for(j=0; j<p_n; j++){
						d_n[p_n*i+j] = Corr_n[p_n*p_n*i+p_n*j+j];
					}
				}
				for(i=0; i<3; i++){
					ADMM_SPACE_rho(Corr_n, d_n, Z_n, U_n, &p_n, Len, Lambda, Rho, Epi_abs, Epi_rel, Max_step);
					if(i<2){
						ADMM_SPACE_d(Corr_n, d_n, Z_n, &p_n, Len);
					}
				}
				free(d_n);
			}
      		
      for(i=0; i<Pos_L; i++){
				for(j=0; j<p_n; j++){
					for(k=0; k<p_n; k++){
						Z[p*p*i+p*member_ind_n[j]+member_ind_n[k]] = Z_n[p_n*p_n*Pos[i]+p_n*j+k];
					}
				}
			}
      
			free(Corr_n);
			free(Z_n);
			free(U_n);	
		}
	}//end iteration across block diagonals
}


//ADMM on likelihood estimation
void ADMM_local_glasso(double *Corr, double *Z, double *U, int *P, int *Len, double *Lambda, double *Rho, 
	double *Epi_abs, double *Epi_rel, int *Max_step){
 
	int p = *P, L = *Len, max_step = *Max_step;
	double lambda = *Lambda, rho = *Rho, epi_abs = *Epi_abs, epi_rel = *Epi_rel, alpha = 1.5;
	double prd = 1, drd = 1, r, s, epi_pri, epi_dual;
	int il = -1, iu = -1, num, lwork = 26*p, liwork = 10*p, info, i, j, k;
	double vl = -1, vu= -1, abstol = pow(10, -6), one = 1, zero = 0, temp_1, epi_pri_1, epi_pri_2;
	double *Omega_ij, *Z_ij, *U_ij, *Corr_ij;
	int *isuppz = (int *) malloc(2*p*sizeof(int));
	double *work = (double *) malloc(lwork*sizeof(double));
	int *iwork = (int *) malloc(liwork*sizeof(int));
	double *Omega = (double *) malloc(p*p*L*sizeof(double));
	double *A = (double *) malloc(p*p*sizeof(double));
	double *eig_vec = (double *) malloc(p*p*sizeof(double));
	double *eig_val = (double *) malloc(p*sizeof(double));
	double *coef = (double *) malloc(p*p*sizeof(double));

	//ADMM iteration
	while((drd>0 || prd>0) && max_step>0){

		r = 0, s = 0, epi_dual = 0, epi_pri_1 = 0, epi_pri_2 = 0;

		//omega-minimization
		for(i=0; i<L; i++){
			for(j=0; j<p; j++){
				Corr_ij = Corr+p*p*i+p*j, Z_ij = Z+p*p*i+p*j, U_ij = U+p*p*i+p*j;
				for(k=j; k<p; k++){
					A[p*j+k] = Corr_ij[k] - rho*(Z_ij[k] - U_ij[k]);
				}
			}
			
			F77_CALL(dsyevr)("V", "A", "L", P, A, P, &vl, &vu, &il, &iu, &abstol, &num, eig_val, eig_vec, P, isuppz, 
				work, &lwork, iwork, &liwork, &info);
			
			for(j=0; j<p; j++){
				eig_val[j] = (-eig_val[j] + sqrt(pow(eig_val[j], 2) + 4*rho))/(2*rho);
				for(k=0; k<p; k++){
					eig_vec[p*j+k] *= sqrt(eig_val[j]);
				}
			}
			F77_CALL(dsyrk)("L", "N", P, P, &one, eig_vec, P, &zero, (Omega+p*p*i), P);
			for(j=0; j<p; j++){
				Omega_ij = Omega+p*p*i+p*j, Z_ij = Z+p*p*i+p*j, U_ij = U+p*p*i+p*j;
				for(k=j; k<p; k++){
					U_ij[k] += alpha*Omega_ij[k] + (1-alpha)*Z_ij[k];
				}
			}
		}

		//z-minimization and dual update
		for(j=0; j<p; j++){
			for(k=j+1; k<p; k++){
				coef[p*j+k] = 0;
			}
		}
		for(i=0; i<L; i++){
			for(j=0; j<p; j++){
				for(k=j+1; k<p; k++){
					coef[p*j+k] += pow(U[p*p*i+p*j+k], 2);
				}
			}
		}
		for(j=0; j<p; j++){
			for(k=j+1; k<p; k++){
				coef[p*j+k] = 1 - lambda/(rho*sqrt(coef[p*j+k]));
			}
		}
		for(i=0; i<L; i++){
			for(j=0; j<p; j++){
				Z_ij = Z+p*p*i+p*j, U_ij = U+p*p*i+p*j;
				s += pow((U_ij[j] - Z_ij[j]), 2)/2;
				Z_ij[j] = U_ij[j];
				U_ij[j] = 0;
				for(k=j+1; k<p; k++){
					if(coef[p*j+k] <= 0){
						s += pow(Z_ij[k], 2);
						Z_ij[k] = 0;
					}
					else{
						temp_1 = Z_ij[k];
						Z_ij[k] = U_ij[k]*coef[p*j+k];
						U_ij[k] -= Z_ij[k];
						s += pow(Z_ij[k] - temp_1, 2);
					}
					epi_dual += pow(U_ij[k], 2);
				}
			}
		}	

		s = rho*sqrt(2*s);
		epi_dual = sqrt(L*p)*epi_abs + epi_rel*rho*sqrt(2*epi_dual);
		drd = s - epi_dual;
		if(drd < 0){
			for(i=0; i<L; i++){
				for(j=0; j<p; j++){
					Omega_ij = Omega+p*p*i+p*j, Z_ij = Z+p*p*i+p*j;
					r += pow(Omega_ij[j]-Z_ij[j], 2)/2;
					epi_pri_1 += pow(Omega_ij[j], 2)/2;
					epi_pri_2 += pow(Z_ij[j], 2)/2;
					for(k=j+1; k<p; k++){
						r += pow(Omega_ij[k]-Z_ij[k], 2);
						epi_pri_1 += pow(Omega_ij[k], 2);
						epi_pri_2 += pow(Z_ij[k], 2);
					}
				}
			}
			r = sqrt(2*r);
			epi_pri = sqrt(L*p)*epi_abs + epi_rel*sqrt(2*((epi_pri_1>epi_pri_2)?epi_pri_1:epi_pri_2));
			prd = r - epi_pri;
		}
		max_step--;	
	}//end ADMM iteration

	if(max_step == 0){
		Rprintf("Warning: algorithm does not converge at max step = %d!\n", *Max_step);
	}

	for(i=0; i<L; i++){
		for(j=0; j<p; j++){
			for(k=0; k<j; k++){
				Z[p*p*i+p*j+k] = Z[p*p*i+p*k+j];
				U[p*p*i+p*j+k] = U[p*p*i+p*k+j];
			}
		}
	}

	free(isuppz);
	free(work);
	free(iwork);
	free(Omega);
	free(A);
	free(eig_vec);
	free(eig_val);
	free(coef);
}


//ADMM on pseudo-likelihood estimation
void ADMM_pseudo_glasso(double *Corr, double *Z, double *U, int *P, int *Len, double *Lambda, double *Rho, 
	double *Epi_abs, double *Epi_rel, int *Max_step){ 

	int p = *P, L = *Len, p_1 = p-1, max_step = *Max_step;
	double lambda = *Lambda, rho = *Rho, epi_abs = *Epi_abs, epi_rel = *Epi_rel, alpha = 1.5;
	double prd = 1, drd = 1, r, s, epi_pri, epi_dual;
	int one = 1, info, i, j, k, index_ijk, index_ikj;
	double temp_1, temp_2, epi_pri_1, epi_pri_2;
	double *Omega_ij, *Z_ij, *U_ij, *Corr_i;
	double *Omega = (double *) malloc(p_1*p*L*sizeof(double));
	double *Chol_ij = (double *) malloc(p_1*p_1*sizeof(double));
	double *C = (double *) malloc(p_1*sizeof(double));
	double *coef = (double *) malloc(p_1*p*sizeof(double));

	for(i=0; i<L; i++){
		for(j=0; j<p; j++){
			for(k=0; k<j; k++){
				Z[p_1*p*i+p_1*j+k] = Z[p*p*i+p*j+k];
				U[p_1*p*i+p_1*j+k] = U[p*p*i+p*j+k];
			}
			for(k=j; k<p_1; k++){
				Z[p_1*p*i+p_1*j+k] = Z[p*p*i+p*j+k+1];
				U[p_1*p*i+p_1*j+k] = U[p*p*i+p*j+k+1];
			}
			Corr[p*p*i+p*j+j] += rho;
		}
		//upper triangular in Corr is used in Givens rotation, lower triangular in Corr stores original information
		F77_CALL(dpotrf)("L", P, (Corr+p*p*i), P, &info); 
	}

	//ADMM iteration
	while((drd>0 || prd>0) && max_step>0){

		r = 0, s = 0, epi_dual = 0, epi_pri_1 = 0, epi_pri_2 = 0;

		//omega-minimization
		for(i=0; i<L; i++){
			Corr_i = Corr+p*p*i;
			for(j=0; j<p; j++){
				Omega_ij = Omega+p_1*p*i+p_1*j, Z_ij = Z+p_1*p*i+p_1*j, U_ij = U+p_1*p*i+p_1*j;
				Givens_rotation(Chol_ij, Corr_i, P, &j);
				for(k=0; k<j; k++){
					C[k] = rho*(Z_ij[k] - U_ij[k]) + Corr_i[p*j+k];
				}
				for(k=j; k<p_1; k++){
					C[k] = rho*(Z_ij[k] - U_ij[k]) + Corr_i[p*(k+1)+j];
				}
				F77_CALL(dtrsv)("L", "N", "N", &p_1, Chol_ij, &p_1, C, &one);
				F77_CALL(dtrsv)("L", "T", "N", &p_1, Chol_ij, &p_1, C, &one);
				for(k=0; k<p_1; k++){
					Omega_ij[k] = C[k];
					U_ij[k] += alpha*C[k] + (1-alpha)*Z_ij[k];
				}
			}
		}

		//z-minimization and dual update
		for(j=1; j<p; j++){
			for(k=0; k<j; k++){
				coef[p_1*j+k] = 0;
			}
		}
		for(i=0; i<L; i++){
			for(j=1; j<p; j++){
				for(k=0; k<j; k++){
					coef[p_1*j+k] += (pow(U[p_1*p*i+p_1*j+k], 2) + pow(U[p_1*p*i+p_1*k+j-1], 2));
				}
			}
		}
		for(j=1; j<p; j++){
			for(k=0; k<j; k++){
				coef[p_1*j+k] = 1 - sqrt(2)*lambda/(rho*sqrt(coef[p_1*j+k]));
			}
		}
		for(i=0; i<L; i++){
			for(j=1; j<p; j++){
				for(k=0; k<j; k++){
					index_ijk = p_1*p*i+p_1*j+k, index_ikj = p_1*p*i+p_1*k+j-1;
					if(coef[p_1*j+k] <= 0){
						s += (pow(Z[index_ijk], 2) + pow(Z[index_ikj], 2));
						Z[index_ijk] = 0;
						Z[index_ikj] = 0;
					}
					else{
						temp_1 = Z[index_ijk], temp_2 = Z[index_ikj];
						Z[index_ijk] = coef[p_1*j+k]*U[index_ijk];
						Z[index_ikj] = coef[p_1*j+k]*U[index_ikj];
						U[index_ijk] -= Z[index_ijk];
						U[index_ikj] -= Z[index_ikj];
						s += (pow(Z[index_ijk] - temp_1, 2) + pow(Z[index_ikj] - temp_2, 2));
					}
					epi_dual += (pow(U[index_ijk], 2) + pow(U[index_ikj], 2));
				}
			}
		}
		
		s = rho*sqrt(s);
		epi_dual = sqrt(L*p)*epi_abs + epi_rel*rho*sqrt(2*epi_dual);
		drd = s - epi_dual;
		if(drd < 0){
			for(i=0; i<L; i++){
				for(j=0; j<p; j++){
					Omega_ij = Omega+p_1*p*i+p_1*j, Z_ij = Z+p_1*p*i+p_1*j;
					for(k=0; k<p_1; k++){
						r += pow(Omega_ij[k]-Z_ij[k], 2);
						epi_pri_1 += pow(Omega_ij[k], 2);
						epi_pri_2 += pow(Z_ij[k], 2);
					}
				}
			}
			r = sqrt(r);
			epi_pri = sqrt(L*p)*epi_abs + epi_rel*sqrt(2*((epi_pri_1>epi_pri_2)?epi_pri_1:epi_pri_2));
			prd = r - epi_pri;
		}
		max_step--;
	}//end ADMM iteration

	if(max_step == 0){
		Rprintf("Warning: algorithm does not converge at max step = %d!\n", *Max_step);
	}

	for(i=L-1; i>=0; i--){
		for(j=p-1; j>=0; j--){
			for(k=p_1-1; k>=j; k--){
				Z[p*p*i+p*j+k+1] = Z[p_1*p*i+p_1*j+k];
				U[p*p*i+p*j+k+1] = U[p_1*p*i+p_1*j+k];
			}
			Z[p*p*i+p*j+j] = 0;
			U[p*p*i+p*j+j] = 0;
			for(k=j-1; k>=0; k--){
				Z[p*p*i+p*j+k] = Z[p_1*p*i+p_1*j+k];
				U[p*p*i+p*j+k] = U[p_1*p*i+p_1*j+k];
			}
		}
	}

	free(Omega);
	free(Chol_ij);
	free(C);
	free(coef);
}


//ADMM on SPACE(sparse partial correlation estimation) (rho step)
void ADMM_SPACE_rho(double *Corr, double *d, double *Z, double *U, int *P, int *Len, double *Lambda, double *Rho, 
	double *Epi_abs, double *Epi_rel, int *Max_step){
 
	int p = *P, L = *Len, p_1 = p-1, max_step = *Max_step;
	double lambda = *Lambda, rho = *Rho, epi_abs = *Epi_abs, epi_rel = *Epi_rel, alpha = 1.5;
	double prd = 1, drd = 1, r, s, epi_pri, epi_dual;
	int one = 1, info, i, j, k, index_ijk, index_ikj;
	double temp_1, temp_2, epi_pri_1, epi_pri_2;
	double *Omega_ij, *Z_ij, *U_ij, *Corr_ij;
	double *Omega = (double *) malloc(p_1*p*L*sizeof(double));
	double *Chol = (double *) malloc(p*p*L*sizeof(double));
	double *Chol_ij = (double *) malloc(p_1*p_1*sizeof(double));
	double *C = (double *) malloc(p_1*sizeof(double));
	double *coef = (double *) malloc(p_1*p*sizeof(double));

	for(i=0; i<L; i++){
		for(j=0; j<p; j++){
			for(k=0; k<j; k++){
				Z[p_1*p*i+p_1*j+k] = Z[p*p*i+p*j+k];
				U[p_1*p*i+p_1*j+k] = U[p*p*i+p*j+k];
			}
			for(k=j; k<p_1; k++){
				Z[p_1*p*i+p_1*j+k] = Z[p*p*i+p*j+k+1];
				U[p_1*p*i+p_1*j+k] = U[p*p*i+p*j+k+1];
			}
			for(k=j; k<p; k++){
				Chol[p*p*i+p*j+k] = Corr[p*p*i+p*j+k]*sqrt(d[p*i+j]*d[p*i+k]);
			}
			Chol[p*p*i+p*j+j] += rho;
		}
		F77_CALL(dpotrf)("L", P, (Chol+p*p*i), P, &info);
	}

	//ADMM iteration
	while((drd>0 || prd>0) && max_step>0){

		r = 0, s = 0, epi_dual = 0, epi_pri_1 = 0, epi_pri_2 = 0;

		//omega-minimization
		for(i=0; i<L; i++){
			for(j=0; j<p; j++){
				Omega_ij = Omega+p_1*p*i+p_1*j, Z_ij = Z+p_1*p*i+p_1*j, U_ij = U+p_1*p*i+p_1*j, Corr_ij = Corr+p*p*i+p*j;
				Givens_rotation(Chol_ij, (Chol+p*p*i), P, &j);
				for(k=0; k<j; k++){
					C[k] = rho*(Z_ij[k] - U_ij[k]) + Corr_ij[k]*sqrt(d[p*i+j]*d[p*i+k]);
				}
				for(k=j; k<p_1; k++){
					C[k] = rho*(Z_ij[k] - U_ij[k]) + Corr_ij[k+1]*sqrt(d[p*i+j]*d[p*i+k+1]);
				}
				F77_CALL(dtrsv)("L", "N", "N", &p_1, Chol_ij, &p_1, C, &one);
				F77_CALL(dtrsv)("L", "T", "N", &p_1, Chol_ij, &p_1, C, &one);
				for(k=0; k<p_1; k++){
					Omega_ij[k] = C[k];
					U_ij[k] += alpha*C[k] + (1-alpha)*Z_ij[k];
				}
			}
		}

		//z-minimization and dual update
		for(j=1; j<p; j++){
			for(k=0; k<j; k++){
				coef[p_1*j+k] = 0;
			}
		}
		for(i=0; i<L; i++){
			for(j=1; j<p; j++){
				for(k=0; k<j; k++){
					coef[p_1*j+k] += pow((U[p_1*p*i+p_1*j+k] + U[p_1*p*i+p_1*k+j-1])/2, 2);
				}
			}
		}
		for(j=1; j<p; j++){
			for(k=0; k<j; k++){
				coef[p_1*j+k] = 1 - lambda/(rho*sqrt(coef[p_1*j+k]));
			}
		}
		for(i=0; i<L; i++){
			for(j=1; j<p; j++){
				for(k=0; k<j; k++){
					index_ijk = p_1*p*i+p_1*j+k, index_ikj = p_1*p*i+p_1*k+j-1;
					if(coef[p_1*j+k] <= 0){
						s += (pow(Z[index_ijk], 2) + pow(Z[index_ikj], 2));
						Z[index_ijk] = 0;
						Z[index_ikj] = 0;
					}
					else{
						temp_1 = Z[index_ijk], temp_2 = Z[index_ikj];
						Z[index_ijk] = coef[p_1*j+k]*(U[index_ijk] + U[index_ikj])/2;
						Z[index_ikj] = Z[index_ijk];
						U[index_ijk] -= Z[index_ijk];
						U[index_ikj] -= Z[index_ikj];
						s += (pow(Z[index_ijk] - temp_1, 2) + pow(Z[index_ikj] - temp_2, 2));
					}
					epi_dual += (pow(U[index_ijk], 2) + pow(U[index_ikj], 2));
				}
			}
		}

		s = rho*sqrt(s);
		epi_dual = sqrt(L*p)*epi_abs + epi_rel*rho*sqrt(2*epi_dual);
		drd = s - epi_dual;
		if(drd < 0){
			for(i=0; i<L; i++){
				for(j=0; j<p; j++){
					Omega_ij = Omega+p_1*p*i+p_1*j, Z_ij = Z+p_1*p*i+p_1*j;
					for(k=0; k<p_1; k++){
						r += pow(Omega_ij[k]-Z_ij[k], 2);
						epi_pri_1 += pow(Omega_ij[k], 2);
						epi_pri_2 += pow(Z_ij[k], 2);
					}
				}
			}
			r = sqrt(r);
			epi_pri = sqrt(L*p)*epi_abs + epi_rel*sqrt(2*((epi_pri_1>epi_pri_2)?epi_pri_1:epi_pri_2));
			prd = r - epi_pri;
		}
		max_step--;
	}//end ADMM iteration

	if(max_step == 0){
		Rprintf("Warning: algorithm does not converge at max step = %d!\n", *Max_step);
	}

	for(i=L-1; i>=0; i--){
		for(j=p-1; j>=0; j--){
			for(k=p_1-1; k>=j; k--){
				Z[p*p*i+p*j+k+1] = Z[p_1*p*i+p_1*j+k];
				U[p*p*i+p*j+k+1] = U[p_1*p*i+p_1*j+k];
			}
			Z[p*p*i+p*j+j] = 0;
			U[p*p*i+p*j+j] = 0;
			for(k=j-1; k>=0; k--){
				Z[p*p*i+p*j+k] = Z[p_1*p*i+p_1*j+k];
				U[p*p*i+p*j+k] = U[p_1*p*i+p_1*j+k];
			}
		}
	}

	free(Omega);
	free(Chol);
	free(Chol_ij);
	free(C);
	free(coef);
}


//ADMM on SPACE(sparse partial correlation estimation) (d step)
void ADMM_SPACE_d(double *Corr, double *d, double *Z, int *P, int *Len){

	int p = *P, L = *Len;
	int i, j, k, m;
	double *rho = (double *) malloc(p*sizeof(double));
	double *d_new = (double *) malloc(p*sizeof(double));

	for(i=0; i<L; i++){
		for(j=0; j<p; j++){
			d_new[j] = 0;
			for(k=0; k<p; k++){
				rho[k] = Z[p*p*i+p*j+k];
			}
			for(k=0; k<p; k++){
				if(rho[k] != 0){
					d_new[j] += pow(rho[k], 2)*d[p*i+k]*Corr[p*p*i+p*k+k];
					for(m=k+1; m<p; m++){
						if(rho[m] != 0){
							d_new[j] += 2*rho[k]*rho[m]*sqrt(d[p*i+k]*d[p*i+m])*Corr[p*p*i+p*k+m];
						}
					}
				}
			}
			d_new[j] = 1/(-d_new[j]/d[p*i+j] + Corr[p*p*i+p*j+j]);
		}
		for(j=0; j<p; j++){
			d[p*i+j] = d_new[j];
		}
	}

	free(rho);
	free(d_new);
}


//Givens rotation
void Givens_rotation(double *U, double *Chol, int *P, int *J){

	int p = *P, j = *J, p_1 = p-1, k, m;
	double r, c, s, temp_1, temp_2;
	double *U_k, *Chol_k;

	//inherit from previous upper triangular matrix with j-1
	//recover (j-1)th column
	for(k=0; k<j; k++){
		U[p_1*k+j-1] = Chol[p*k+j-1];
	}
	//recover (j-1)th row
	if(j >= 1){
		for(m=j; m<p_1; m++){
			U[p_1*(j-1)+m] = Chol[p*(j-1)+m+1];
		}
	}
	//recover jth row
	for(m=j; m<p_1; m++){
		U[p_1*j+m] = Chol[p*j+m+1];
	}

	//implement Givens rotation from jth row to the last row
	for(k=j; k<p_1; k++){
		U_k = U+p_1*k+k, Chol_k = Chol+p*(k+1)+k+1;
		r = sqrt(pow(U_k[0], 2) + pow(Chol_k[0], 2));
		c = U_k[0]/r;
		s = -Chol_k[0]/r;
		U_k[0] = r;	
		for(m=1; m<p_1-k; m++){
			temp_1 = U_k[m], temp_2 = Chol_k[m];
			U_k[m] = c*temp_1 - s*temp_2;
			U_k[m+p_1] = s*temp_1 + c*temp_2;
		}
	}
}
#include <stdlib.h>
#include <stdio.h>
#include <R.h>
#include <Rmath.h>
#include <math.h>
#include <R_ext/Lapack.h>

//apply ADMM across lambda's when h and d are fixed
//member_ind, csize_ind, No: identify block diagonal structure
//LL: number of time points within neighborhood
void ADMM_lambda(int *P, int *member_ind, int *csize_ind, int *No, int *LL, int *Pos, int *Pos_Len, double *Corr, double *Sigma, double *Z, double *Z_pos, int *Lambda_Len, double *Lambda, double *Rho, double *Epi_abs, double *Epi_rel, int *pseudo_fit, int *pseudo_refit, double *thres);

//apply ADMM to fixed h, d and lambda
void ADMM_cluster(int *P, int *member_ind, int *csize_ind, int *No, int *LL, int *Pos, int *Pos_Len, double *Corr, double *sd, double *Z, double *Z_pos, double *U, double *Lambda, double *Rho, double *Epi_abs, double *Epi_rel, int *pseudo_fit, int *pseudo_refit, int *S_Len);

//apply ADMM to fixed h, d and lambda (simple version)
void ADMM_simple(int *P, int *member_ind, int *csize_ind, int *No, int *LL, int *Pos, int *Pos_Len, double *Corr, double *sd, double *Z, double *Z_pos, double *Lambda, double *Rho, double *Epi_abs, double *Epi_rel, int *pseudo_fit, int *pseudo_refit, int *S_Len);

//local group graphical lasso
void ADMM_local_glasso(int *P, int *LL, double *Sigma, double *Z, double *U, double *Lambda, double *Rho, double *Epi_abs, double *Epi_rel);

//pseudo-likelihood group lasso (symmetry)
void ADMM_pseudo_glasso(int *P, int *LL, double *Sigma, double *Z, double *U, double *Lambda, double *Rho, double *Epi_abs, double *Epi_rel);

//SPACE (rho step)
void ADMM_SPACE_rho(int *P, int *LL, double *d, double *Sigma, double *Z, double *U, double *Lambda, double *Rho, double *Epi_abs, double *Epi_rel);

//SPACE (d step)
void ADMM_SPACE_d(int *P, int *LL, double *d, double *Sigma, double *Z);

//likelihood refit
void ADMM_refit(int *P, int *Pos, double *Sigma, double *Z, double *U, int *S, int *S_Len, double *Rho, double *Epi_abs, double *Epi_rel);

//pseudo-likelihood refit
void ADMM_pseudo_refit(int *P, int *Pos, double *Sigma, double *Z, double *Z_pos, int *S_Len);

void Givens_rotation(double *U, double *Chol, int *P, int *J);


/////////////////////////////////////////////////



void ADMM_lambda(int *P, int *member_ind, int *csize_ind, int *No, int *LL, int *Pos, int *Pos_Len, double *Corr, double *Sigma, double *Z, double *Z_pos, int *Lambda_Len, double *Lambda, double *Rho, double *Epi_abs, double *Epi_rel, int *pseudo_fit, int *pseudo_refit, double *thres){

	int i, j, k, m, p = *P, L = *LL, Pos_L = *Pos_Len, Lambda_L = *Lambda_Len, no_sum = 0;

	double *Z_i = (double *) malloc(p*p*L*sizeof(double));
	double *U_i = (double *) malloc(p*p*L*sizeof(double));
	int *S_L = (int *) malloc(sizeof(int));

	for(i=0; i<L; i++){
		for(j=0; j<p; j++){
			for(k=0; k<p; k++){
				Z_i[p*p*i+p*j+k] = 0;
				U_i[p*p*i+p*j+k] = 0;
			}
		}
	}

	*S_L = 0;

	//iteration across lambda
	for(i=(Lambda_L-1); i>=0; i--){

		if((*S_L)>((*thres)*p)){
			break;
		}

		ADMM_cluster(P, (member_ind+p*i), (csize_ind+no_sum), (No+i), LL, Pos, Pos_Len, Corr, Sigma, Z_i, (Z_pos+p*p*Pos_L*i), U_i, (Lambda+i), (Rho+i), Epi_abs, Epi_rel, pseudo_fit, pseudo_refit, S_L);

		for(m=0; m<Pos_L; m++){
			for(j=0; j<p; j++){
				for(k=0; k<p; k++){
					Z[p*p*Pos_L*i+p*p*m+p*j+k] = Z_i[p*p*Pos[m]+p*j+k];
				}
			} 
		}

		no_sum += (No[i]+1);
	}//end iteration across lambda

	free(Z_i);
	free(U_i);
	free(S_L);
}



void ADMM_cluster(int *P, int *member_ind, int *csize_ind, int *No, int *LL, int *Pos, int *Pos_Len, double *Corr, double *sd, double *Z, double *Z_pos, double *U, double *Lambda, double *Rho, double *Epi_abs, double *Epi_rel, int *pseudo_fit, int *pseudo_refit, int *S_Len){

	int p = *P, no = *No, L = *LL, Pos_L = *Pos_Len, p_n, n, i, j, k, pos, S_L;
	int *member_ind_n;
	double rho = (*Rho)/sqrt(L);

	S_Len[0] = 0;

	//iteration across block diagonals
	for(n=0; n<no; n++){

		p_n = csize_ind[n+1]-csize_ind[n];
		member_ind_n = &member_ind[csize_ind[n]];

		//block diagonal with dimension equals one
		if(p_n==1){
			if(*pseudo_fit == 0){
				for(i=0; i<L; i++){
					Z[p*p*i+p*(*member_ind_n)+(*member_ind_n)] = 1/Corr[p*p*i+p*(*member_ind_n)+(*member_ind_n)];
				}
			}
			for(i=0; i<Pos_L; i++){
				Z_pos[p*p*i+p*(*member_ind_n)+(*member_ind_n)] = 1/(Corr[p*p*Pos[i]+p*(*member_ind_n)+(*member_ind_n)]*sd[*member_ind_n]*sd[*member_ind_n]);
			}
		}
		//block diagonal with dimension larger than one
		else{
			double *Sigma_n = (double *) malloc(p_n*p_n*L*sizeof(double));
			double *Z_n = (double *) malloc(p_n*p_n*L*sizeof(double));
			double *U_n = (double *) malloc(p_n*p_n*L*sizeof(double));
			double *Z_pos_n = (double *) malloc(p_n*p_n*sizeof(double));
			int *S = (int *) malloc(p_n*(p_n-1)*sizeof(int));
		
			for(i=0; i<L; i++){
				for(j=0; j<p_n; j++){
					for(k=0; k<p_n; k++){
						Sigma_n[p_n*p_n*i+p_n*j+k] = Corr[p*p*i+p*(*(member_ind_n+j))+(*(member_ind_n+k))];
						Z_n[p_n*p_n*i+p_n*j+k] = Z[p*p*i+p*(*(member_ind_n+j))+(*(member_ind_n+k))];
						U_n[p_n*p_n*i+p_n*j+k] = U[p*p*i+p*(*(member_ind_n+j))+(*(member_ind_n+k))];
					}
				}
			}

			//model fitting
			if(*pseudo_fit == 0){
				ADMM_local_glasso(&p_n, LL, Sigma_n, Z_n, U_n, Lambda, Rho, Epi_abs, Epi_rel);
			}
			else if(*pseudo_fit == 1){
				ADMM_pseudo_glasso(&p_n, LL, Sigma_n, Z_n, U_n, Lambda, Rho, Epi_abs, Epi_rel);
			}
			else if(*pseudo_fit == 2){
				double *d_n = (double *) malloc(p_n*L*sizeof(double));
				for(i=0; i<L; i++){
					for(j=0; j<p_n; j++){
						d_n[p_n*i+j] = Sigma_n[p_n*p_n*i+p_n*j+j];
					}
				}
				for(i=0; i<3; i++){
					ADMM_SPACE_rho(&p_n, LL, d_n, Sigma_n, Z_n, U_n, Lambda, Rho, Epi_abs, Epi_rel);
					if(i<2){
						ADMM_SPACE_d(&p_n, LL, d_n, Sigma_n, Z_n);
					}
				}
				free(d_n);
			}

			for(i=0; i<L; i++){
				for(j=0; j<p_n; j++){
					for(k=0; k<p_n; k++){
						Z[p*p*i+p*(*(member_ind_n+j))+(*(member_ind_n+k))] = Z_n[p_n*p_n*i+p_n*j+k];
						U[p*p*i+p*(*(member_ind_n+j))+(*(member_ind_n+k))] = U_n[p_n*p_n*i+p_n*j+k];
						Sigma_n[p_n*p_n*i+p_n*j+k] = Corr[p*p*i+p*(*(member_ind_n+j))+(*(member_ind_n+k))]*sd[*(member_ind_n+j)]*sd[*(member_ind_n+k)];
					}
				}
			}
	
			//model refitting
			if(*pseudo_refit == 0){

				double *U_pos_n = (double *) malloc(p_n*p_n*sizeof(double));

				for(i=0; i<Pos_L; i++){
	
					pos = Pos[i];
					S_L = 0;
				
					for(j=0; j<p_n; j++){
						for(k=j; k<p_n; k++){
							Z_pos_n[p_n*j+k] = Z_n[p_n*p_n*pos+p_n*j+k];
							U_pos_n[p_n*j+k] = U_n[p_n*p_n*pos+p_n*j+k];
							if(j!=k && Z_n[p_n*p_n*pos+p_n*j+k]!=0 && Z_n[p_n*p_n*pos+p_n*k+j]!=0){
								S[2*S_L] = j;
								S[2*S_L+1] = k;
								S_L++;
								S_Len[0]++;		
							}
						}
					}

					ADMM_refit(&p_n, &pos, Sigma_n, Z_pos_n, U_pos_n, S, &S_L, &rho, Epi_abs, Epi_rel);

					for(j=0; j<p_n; j++){
						for(k=j; k<p_n; k++){
							Z_pos[p*p*i+p*(*(member_ind_n+j))+(*(member_ind_n+k))] = Z_pos_n[p_n*j+k];
							Z_pos[p*p*i+p*(*(member_ind_n+k))+(*(member_ind_n+j))] = Z_pos_n[p_n*j+k];
						}
					}
				}

				free(U_pos_n);
			}
			else if(*pseudo_refit == 1){

				for(i=0; i<Pos_L; i++){
	
					pos = Pos[i];

					ADMM_pseudo_refit(&p_n, &pos, Sigma_n, Z_n, Z_pos_n, S_Len);

					for(j=0; j<p_n; j++){
						for(k=0; k<p_n; k++){
							Z_pos[p*p*i+p*(*(member_ind_n+j))+(*(member_ind_n+k))] = Z_pos_n[p_n*j+k];
						}
					}
				}
			}
			
			free(Sigma_n);
			free(Z_n);
			free(U_n);
			free(Z_pos_n);
			free(S);		
		}
	}//end iteration across block diagonals

	S_Len[0] = S_Len[0]/Pos_L;
}



void ADMM_simple(int *P, int *member_ind, int *csize_ind, int *No, int *LL, int *Pos, int *Pos_Len, double *Corr, double *sd, double *Z, double *Z_pos, double *Lambda, double *Rho, double *Epi_abs, double *Epi_rel, int *pseudo_fit, int *pseudo_refit, int *S_Len){
  
  int p = *P, no = *No, L = *LL, Pos_L = *Pos_Len, p_n, n, i, j, k, pos, S_L;
  int *member_ind_n;
  double rho = (*Rho)/sqrt(L);
  
  S_Len[0] = 0;
  
  //iteration across block diagonals
  for(n=0; n<no; n++){
    
    p_n = csize_ind[n+1]-csize_ind[n];
    member_ind_n = &member_ind[csize_ind[n]];
    
    //block diagonal with dimension equals one
    if(p_n==1){
      if(*pseudo_fit == 0){
        for(i=0; i<Pos_L; i++){
          Z[p*p*i+p*(*member_ind_n)+(*member_ind_n)] = 1/Corr[p*p*Pos[i]+p*(*member_ind_n)+(*member_ind_n)];
        }
      }
      for(i=0; i<Pos_L; i++){
        Z_pos[p*p*i+p*(*member_ind_n)+(*member_ind_n)] = 1/(Corr[p*p*Pos[i]+p*(*member_ind_n)+(*member_ind_n)]*sd[*member_ind_n]*sd[*member_ind_n]);
      }
    }
    //block diagonal with dimension larger than one
    else{
      double *Sigma_n = (double *) malloc(p_n*p_n*L*sizeof(double));
      double *Z_n = (double *) malloc(p_n*p_n*L*sizeof(double));
      double *U_n = (double *) malloc(p_n*p_n*L*sizeof(double));
      double *Z_pos_n = (double *) malloc(p_n*p_n*sizeof(double));
      int *S = (int *) malloc(p_n*(p_n-1)*sizeof(int));
      
      for(i=0; i<L; i++){
        for(j=0; j<p_n; j++){
          for(k=0; k<p_n; k++){
            Sigma_n[p_n*p_n*i+p_n*j+k] = Corr[p*p*i+p*(*(member_ind_n+j))+(*(member_ind_n+k))];
            Z_n[p_n*p_n*i+p_n*j+k] = 0;
            U_n[p_n*p_n*i+p_n*j+k] = 0;
          }
        }
      }
      
      //model fitting
      if(*pseudo_fit == 0){
        ADMM_local_glasso(&p_n, LL, Sigma_n, Z_n, U_n, Lambda, Rho, Epi_abs, Epi_rel);
      }
      else if(*pseudo_fit == 1){
        ADMM_pseudo_glasso(&p_n, LL, Sigma_n, Z_n, U_n, Lambda, Rho, Epi_abs, Epi_rel);
      }
      else if(*pseudo_fit == 2){
        double *d_n = (double *) malloc(p_n*L*sizeof(double));
        for(i=0; i<L; i++){
          for(j=0; j<p_n; j++){
            d_n[p_n*i+j] = Sigma_n[p_n*p_n*i+p_n*j+j];
          }
        }
        for(i=0; i<3; i++){
          ADMM_SPACE_rho(&p_n, LL, d_n, Sigma_n, Z_n, U_n, Lambda, Rho, Epi_abs, Epi_rel);
          if(i<2){
            ADMM_SPACE_d(&p_n, LL, d_n, Sigma_n, Z_n);
          }
        }
        free(d_n);
      }
      
      for(j=0; j<p_n; j++){
        for(k=0; k<p_n; k++){
          for(i=0; i<L; i++){
            Sigma_n[p_n*p_n*i+p_n*j+k] = Corr[p*p*i+p*(*(member_ind_n+j))+(*(member_ind_n+k))]*sd[*(member_ind_n+j)]*sd[*(member_ind_n+k)];
          }
          for(i=0; i<Pos_L; i++){
            Z[p*p*i+p*(*(member_ind_n+j))+(*(member_ind_n+k))] = Z_n[p_n*p_n*Pos[i]+p_n*j+k];
          }
        }
      }
      
      //model refitting
      if(*pseudo_refit == 0){
        
        double *U_pos_n = (double *) malloc(p_n*p_n*sizeof(double));
        
        for(i=0; i<Pos_L; i++){
          
          pos = Pos[i];
          S_L = 0;
          
          for(j=0; j<p_n; j++){
            for(k=j; k<p_n; k++){
              Z_pos_n[p_n*j+k] = Z_n[p_n*p_n*pos+p_n*j+k];
              U_pos_n[p_n*j+k] = U_n[p_n*p_n*pos+p_n*j+k];
              if(j!=k && Z_n[p_n*p_n*pos+p_n*j+k]!=0 && Z_n[p_n*p_n*pos+p_n*k+j]!=0){
                S[2*S_L] = j;
                S[2*S_L+1] = k;
                S_L++;
                S_Len[0]++;		
              }
            }
          }
          
          ADMM_refit(&p_n, &pos, Sigma_n, Z_pos_n, U_pos_n, S, &S_L, &rho, Epi_abs, Epi_rel);
          
          for(j=0; j<p_n; j++){
            for(k=j; k<p_n; k++){
              Z_pos[p*p*i+p*(*(member_ind_n+j))+(*(member_ind_n+k))] = Z_pos_n[p_n*j+k];
              Z_pos[p*p*i+p*(*(member_ind_n+k))+(*(member_ind_n+j))] = Z_pos_n[p_n*j+k];
            }
          }
        }
        
        free(U_pos_n);
      }
      else if(*pseudo_refit == 1){
        
        for(i=0; i<Pos_L; i++){
          
          pos = Pos[i];
          
          ADMM_pseudo_refit(&p_n, &pos, Sigma_n, Z_n, Z_pos_n, S_Len);
          
          for(j=0; j<p_n; j++){
            for(k=0; k<p_n; k++){
              Z_pos[p*p*i+p*(*(member_ind_n+j))+(*(member_ind_n+k))] = Z_pos_n[p_n*j+k];
            }
          }
        }
      }
      
      free(Sigma_n);
      free(Z_n);
      free(U_n);
      free(Z_pos_n);
      free(S);		
    }
  }//end iteration across block diagonals
  
  S_Len[0] = S_Len[0]/Pos_L;
}



void ADMM_local_glasso(int *P, int *LL, double *Sigma, double *Z, double *U, double *Lambda, double *Rho, double *Epi_abs, double *Epi_rel){
 
	int p = *P, L = *LL;
	double lambda = *Lambda, rho = *Rho, epi_abs = *Epi_abs, epi_rel = *Epi_rel, alpha = 1.5;
	double prd = 1, drd = 1, r, s, epi_pri, epi_dual;
	int il = -1, iu = -1, num, lwork = 26*p, liwork = 10*p, info, i, j, k, m, index_ijk;
	double vl = -1, vu= -1, abstol = pow(10, -6), coef, temp, temp_1, temp_2, epi_pri_1, epi_pri_2;
	int *isuppz = (int *) malloc(2*p*sizeof(int));
	double *work = (double *) malloc(lwork*sizeof(double));
	int *iwork = (int *) malloc(liwork*sizeof(int));
	double *Omega = (double *) malloc(p*p*L*sizeof(double));
	double *A = (double *) malloc(p*p*sizeof(double));
	double *eig_vec = (double *) malloc(p*p*sizeof(double));
	double *eig_val = (double *) malloc(p*sizeof(double));

	//ADMM iteration
	while(prd>0 || drd>0){

		r = 0, s = 0, epi_dual = 0, epi_pri_1 = 0, epi_pri_2 = 0;

		for(i=0; i<L; i++){
			for(j=0; j<p; j++){
				for(k=j; k<p; k++){
					A[p*j+k] = Sigma[p*p*i+p*j+k] - rho*(Z[p*p*i+p*j+k] - U[p*p*i+p*j+k]);
				}
			}
			
			F77_CALL(dsyevr)("V", "A", "L", &p, A, &p, &vl, &vu, &il, &iu, &abstol, &num, eig_val, eig_vec, &p, isuppz, work, &lwork, iwork, &liwork, &info);
			
			for(j=0; j<p; j++){
				eig_val[j] = (-eig_val[j] + sqrt(pow(eig_val[j], 2) + 4*rho))/(2*rho);
			}
			for(j=0; j<p; j++){
				for(k=j; k<p; k++){
					index_ijk = p*p*i+p*j+k;
					Omega[index_ijk] = 0;
					for(m=0; m<p; m++){
						Omega[index_ijk] += eig_val[m]*eig_vec[p*m+j]*eig_vec[p*m+k];
					}
					U[index_ijk] += alpha*Omega[index_ijk] +(1-alpha)*Z[index_ijk];
				}
				index_ijk = p*p*i+p*j+j;
				s += pow((U[index_ijk] - Z[index_ijk]), 2)/2;
				Z[index_ijk] = U[index_ijk];
				U[index_ijk] = 0;
				r += pow((Omega[index_ijk] - Z[index_ijk]), 2)/2;
				epi_pri_1 += pow(Omega[index_ijk], 2)/2;
				epi_pri_2 += pow(Z[index_ijk], 2)/2;
			}
		}

		for(j=0; j<p; j++){
			for(k=j+1; k<p; k++){
				temp_1 = 0, temp_2 = 0;
				for(i=0; i<L; i++){
					temp_1 += pow(U[p*p*i+p*j+k], 2);
					temp_2 += pow(Omega[p*p*i+p*j+k], 2);
				}
				coef = 1 - lambda/(rho*sqrt(temp_1));
				if(coef <= 0){
					for(i=0; i<L; i++){
						index_ijk = p*p*i+p*j+k;
						s += pow(Z[index_ijk], 2);
						Z[index_ijk] = 0;
					}
					r += temp_2;
					epi_dual += temp_1;
				}
				else{
					for(i=0; i<L; i++){
						index_ijk = p*p*i+p*j+k;
						temp = Z[index_ijk];
						Z[index_ijk] = U[index_ijk]*coef;
						U[index_ijk] -= Z[index_ijk];
						r += pow((Omega[index_ijk] - Z[index_ijk]), 2);
						s += pow(Z[index_ijk]-temp, 2);
						epi_pri_2 += pow(Z[index_ijk], 2);
						epi_dual += pow(U[index_ijk], 2);
					}
				}
				epi_pri_1 += temp_2;	
			}
		}

		r = sqrt(2*r);
		s = rho*sqrt(2*s);
		epi_pri = sqrt(L*p)*epi_abs + epi_rel*sqrt(2*((epi_pri_1>epi_pri_2)?epi_pri_1:epi_pri_2));
		epi_dual = sqrt(L*p)*epi_abs + epi_rel*rho*sqrt(2*epi_dual);
		
		prd = (r>epi_pri);
		drd = (s>epi_dual);
	}//end ADMM iteration

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
}



void ADMM_pseudo_glasso(int *P, int *LL, double *Sigma, double *Z, double *U, double *Lambda, double *Rho, double *Epi_abs, double *Epi_rel){ 

	int p = *P, L = *LL, p_1 = p-1;
	double lambda = *Lambda, rho = *Rho, epi_abs = *Epi_abs, epi_rel = *Epi_rel, alpha = 1.5;
	double prd = 1, drd = 1, r, s, epi_pri, epi_dual;
	int nrhs = 1, info, i, j, k, m, index_ijk, index_ikj;
	double coef, temp, temp_1, temp_2, epi_pri_1, epi_pri_2, alp = 1;
	double *Omega = (double *) malloc(p_1*p*L*sizeof(double));
	double *Chol = (double *) malloc(p*p*L*sizeof(double));
	double *Chol_ij = (double *) malloc(p_1*p_1*sizeof(double));
	double *C = (double *) malloc(p_1*sizeof(double));

	for(i=0; i<L; i++){
		for(j=0; j<p; j++){
			for(k=0; k<p_1; k++){
				Z[p_1*p*i+p_1*j+k] = Z[p*p*i+p*j+((k>=j)?(k+1):k)];
				U[p_1*p*i+p_1*j+k] = U[p*p*i+p*j+((k>=j)?(k+1):k)];
			}
			for(k=j; k<p; k++){
				Chol[p*p*i+p*j+k] = Sigma[p*p*i+p*j+k];
			}
			Chol[p*p*i+p*j+j]+=rho;
		}
		F77_CALL(dpotrf)("L", P, (Chol+p*p*i), P, &info);
	}

	//ADMM iteration
	while(prd>0 || drd>0){

		r = 0, s = 0, epi_dual = 0, epi_pri_1 = 0, epi_pri_2 = 0;

		for(i=0; i<L; i++){
			for(j=0; j<p; j++){
				Givens_rotation(Chol_ij, (Chol+p*p*i), P, &j);
				for(k=0; k<j; k++){
					C[k] = rho*(Z[p_1*p*i+p_1*j+k]-U[p_1*p*i+p_1*j+k])+Sigma[p*p*i+p*j+k];
				}
				for(k=j; k<p_1; k++){
					C[k] = rho*(Z[p_1*p*i+p_1*j+k]-U[p_1*p*i+p_1*j+k])+Sigma[p*p*i+p*j+k+1];
				}
				F77_CALL(dtrsm)("L", "L", "N", "N", &p_1, &nrhs, &alp, Chol_ij, &p_1, C, &p_1);
				F77_CALL(dtrsm)("L", "L", "T", "N", &p_1, &nrhs, &alp, Chol_ij, &p_1, C, &p_1);
				for(k=0; k<p_1; k++){
					index_ijk = p_1*p*i+p_1*j+k;
					Omega[index_ijk] = C[k];
					U[index_ijk] += alpha*C[k]+(1-alpha)*Z[index_ijk];
				}
			}
		}

		for(j=1; j<p; j++){
			for(k=0; k<j; k++){
				temp_1 = 0, temp_2 = 0;
				for(i=0; i<L; i++){
					temp_1 += (pow(U[p_1*p*i+p_1*j+k], 2) + pow(U[p_1*p*i+p_1*k+j-1], 2));
					temp_2 += (pow(Omega[p_1*p*i+p_1*j+k], 2) + pow(Omega[p_1*p*i+p_1*k+j-1], 2));
				}
				coef = 1 - sqrt(2)*lambda/(rho*sqrt(temp_1));
				if(coef <= 0){
					for(i=0; i<L; i++){
						index_ijk = p_1*p*i+p_1*j+k, index_ikj = p_1*p*i+p_1*k+j-1;
						s += (pow(Z[index_ijk], 2) + pow(Z[index_ikj], 2));
						Z[index_ijk] = 0;
						Z[index_ikj] = 0;
					}
					r += temp_2;
					epi_dual += temp_1;
				}
				else{
					for(i=0; i<L; i++){
						index_ijk = p_1*p*i+p_1*j+k, index_ikj = p_1*p*i+p_1*k+j-1;
						temp = Z[index_ijk], temp_1 = Z[index_ikj];
						Z[index_ijk] = U[index_ijk]*coef;
						Z[index_ikj] = U[index_ikj]*coef;
						U[index_ijk] -= Z[index_ijk];
						U[index_ikj] -= Z[index_ikj];
						r += (pow((Omega[index_ijk] - Z[index_ijk]), 2) + pow((Omega[index_ikj] - Z[index_ikj]), 2));
						s += (pow(Z[index_ijk]-temp, 2) + pow(Z[index_ikj]-temp_1, 2));
						epi_pri_2 += (pow(Z[index_ijk], 2) + pow(Z[index_ikj], 2));
						epi_dual += (pow(U[index_ijk], 2) + pow(U[index_ikj], 2));
					}
				}
				epi_pri_1 += temp_2;	
			}
		}
		
		r = sqrt(r);
		s = rho*sqrt(s);
		epi_pri = sqrt(L*p)*epi_abs + epi_rel*sqrt(2*((epi_pri_1>epi_pri_2)?epi_pri_1:epi_pri_2));
		epi_dual = sqrt(L*p)*epi_abs + epi_rel*rho*sqrt(2*epi_dual);
		
		prd = (r>epi_pri);
		drd = (s>epi_dual);
	}//end ADMM iteration

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
}



void ADMM_SPACE_rho(int *P, int *LL, double *d, double *Sigma, double *Z, double *U, double *Lambda, double *Rho, double *Epi_abs, double *Epi_rel){
 
	int p = *P, L = *LL, p_1 = p-1;
	double lambda = *Lambda, rho = *Rho, epi_abs = *Epi_abs, epi_rel = *Epi_rel, alpha = 1.5;
	double prd = 1, drd = 1, r, s, epi_pri, epi_dual;
	int nrhs = 1, info, i, j, k, m, index_ijk, index_ikj;
	double coef, temp, temp_1, temp_2, epi_pri_1, epi_pri_2, alp = 1;
	double *Omega = (double *) malloc(p_1*p*L*sizeof(double));
	double *Chol = (double *) malloc(p*p*L*sizeof(double));
	double *Chol_ij = (double *) malloc(p_1*p_1*sizeof(double));
	double *C = (double *) malloc(p_1*sizeof(double));

	for(i=0; i<L; i++){
		for(j=0; j<p; j++){
			for(k=0; k<p_1; k++){
				Z[p_1*p*i+p_1*j+k] = Z[p*p*i+p*j+((k>=j)?(k+1):k)];
				U[p_1*p*i+p_1*j+k] = U[p*p*i+p*j+((k>=j)?(k+1):k)];
			}
			for(k=j; k<p; k++){
				Chol[p*p*i+p*j+k] = Sigma[p*p*i+p*j+k]*sqrt(d[p*i+j]*d[p*i+k]);
			}
			Chol[p*p*i+p*j+j]+=rho;
		}
		F77_CALL(dpotrf)("L", P, (Chol+p*p*i), P, &info);
	}

	//ADMM iteration
	while(prd>0 || drd>0){

		r = 0, s = 0, epi_dual = 0, epi_pri_1 = 0, epi_pri_2 = 0;

		for(i=0; i<L; i++){
			for(j=0; j<p; j++){
				Givens_rotation(Chol_ij, (Chol+p*p*i), P, &j);
				for(k=0; k<j; k++){
					C[k] = rho*(Z[p_1*p*i+p_1*j+k]-U[p_1*p*i+p_1*j+k])+Sigma[p*p*i+p*j+k]*sqrt(d[p*i+j]*d[p*i+k]);
				}
				for(k=j; k<p_1; k++){
					C[k] = rho*(Z[p_1*p*i+p_1*j+k]-U[p_1*p*i+p_1*j+k])+Sigma[p*p*i+p*j+k+1]*sqrt(d[p*i+j]*d[p*i+k+1]);
				}
				F77_CALL(dtrsm)("L", "L", "N", "N", &p_1, &nrhs, &alp, Chol_ij, &p_1, C, &p_1);
				F77_CALL(dtrsm)("L", "L", "T", "N", &p_1, &nrhs, &alp, Chol_ij, &p_1, C, &p_1);
				for(k=0; k<p_1; k++){
					index_ijk = p_1*p*i+p_1*j+k;
					Omega[index_ijk] = C[k];
					U[index_ijk] += alpha*C[k]+(1-alpha)*Z[index_ijk];
				}
			}
		}

		for(j=1; j<p; j++){
			for(k=0; k<j; k++){
				coef = 0, temp_1 = 0, temp_2 = 0;
				for(i=0; i<L; i++){
					index_ijk = p_1*p*i+p_1*j+k, index_ikj = p_1*p*i+p_1*k+j-1;
					coef += pow((U[index_ijk] + U[index_ikj])/2, 2);
					temp_1 += (pow(U[index_ijk], 2) + pow(U[index_ikj], 2));
					temp_2 += (pow(Omega[index_ijk], 2) + pow(Omega[index_ikj], 2));
				}
				coef = 1 - lambda/(rho*sqrt(coef));
				if(coef <= 0){
					for(i=0; i<L; i++){
						index_ijk = p_1*p*i+p_1*j+k, index_ikj = p_1*p*i+p_1*k+j-1;
						s += (pow(Z[index_ijk], 2) + pow(Z[index_ikj], 2));
						Z[index_ijk] = 0;
						Z[index_ikj] = 0;
					}
					r += temp_2;
					epi_dual += temp_1;
				}
				else{
					for(i=0; i<L; i++){
						index_ijk = p_1*p*i+p_1*j+k, index_ikj = p_1*p*i+p_1*k+j-1;
						temp = Z[index_ijk], temp_1 = Z[index_ikj];
						Z[index_ijk] = coef*(U[index_ijk] + U[index_ikj])/2;
						Z[index_ikj] = Z[index_ijk];
						U[index_ijk] -= Z[index_ijk];
						U[index_ikj] -= Z[index_ikj];
						r += (pow((Omega[index_ijk] - Z[index_ijk]), 2) + pow((Omega[index_ikj] - Z[index_ikj]), 2));
						s += (pow(Z[index_ijk]-temp, 2) + pow(Z[index_ikj]-temp_1, 2));
						epi_pri_2 += 2*pow(Z[index_ijk], 2);
						epi_dual += (pow(U[index_ijk], 2) + pow(U[index_ikj], 2));
					}
				}
				epi_pri_1 += temp_2;	
			}
		}
		
		r = sqrt(r);
		s = rho*sqrt(s);
		epi_pri = sqrt(L*p)*epi_abs + epi_rel*sqrt(2*((epi_pri_1>epi_pri_2)?epi_pri_1:epi_pri_2));
		epi_dual = sqrt(L*p)*epi_abs + epi_rel*rho*sqrt(2*epi_dual);
		
		prd = (r>epi_pri);
		drd = (s>epi_dual);
	}//end ADMM iteration

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
}



void ADMM_SPACE_d(int *P, int *LL, double *d, double *Sigma, double *Z){

	int p = *P, L = *LL;
	int i, j, k, m, index_ijk, index_ijm;
	double *rho = (double *) malloc(p*sizeof(double));
	double *d_new = (double *) malloc(p*sizeof(double));

	for(i=0; i<L; i++){
		for(j=0; j<p; j++){
			d_new[j] = 0;
			for(k=0; k<p; k++){
				rho[k] = Z[p*p*i+p*j+k];
			}
			for(k=0; k<p; k++){
				if(rho[k]!=0){
					d_new[j] += pow(rho[k], 2)*d[p*i+k]*Sigma[p*p*i+p*k+k];
					for(m=k+1; m<p; m++){
						if(rho[m]!=0){
							d_new[j] += 2*rho[k]*rho[m]*sqrt(d[p*i+k]*d[p*i+m])*Sigma[p*p*i+p*k+m];
						}
					}
				}
			}
			d_new[j] = 1/(-d_new[j]/d[p*i+j]+Sigma[p*p*i+p*j+j]);
		}
		for(j=0; j<p; j++){
			d[p*i+j] = d_new[j];
		}
	}

	free(rho);
	free(d_new);
}



void ADMM_refit(int *P, int *Pos, double *Sigma, double *Z, double *U, int *S, int *S_Len, double *Rho, double *Epi_abs, double *Epi_rel){
 
	int p = *P, pos = *Pos, S_L = *S_Len;
	double rho = *Rho, epi_abs = *Epi_abs, epi_rel = *Epi_rel;
	double prd = 1, drd = 1, r, s, epi_pri, epi_dual;
	int il = -1, iu = -1, num, lwork = 26*p, liwork = 10*p, info, j, k, m, index_jk;
	double vl = -1, vu= -1, abstol = pow(10, -6), coef, epi_pri_1, epi_pri_2;
	int *isuppz = (int *) malloc(2*p*sizeof(int));
	double *work = (double *) malloc(lwork*sizeof(double));
	int *iwork = (int *) malloc(liwork*sizeof(int));
	double *Omega = (double *) malloc(p*p*sizeof(double));
	double *A = (double *) malloc(p*p*sizeof(double));
	double *eig_vec = (double *) malloc(p*p*sizeof(double));
	double *eig_val = (double *) malloc(p*sizeof(double));

	//ADMM iteration
	while(prd>0 || drd>0){

		r = 0, s = 0, epi_dual = 0, epi_pri_1 = 0, epi_pri_2 = 0;

		for(j=0; j<p; j++){
			for(k=j; k<p; k++){
				A[p*j+k] = Sigma[p*p*pos+p*j+k] - rho*(Z[p*j+k] - U[p*j+k]);
			}
		}
			
		F77_CALL(dsyevr)("V", "A", "L", &p, A, &p, &vl, &vu, &il, &iu, &abstol, &num, eig_val, eig_vec, &p, isuppz, work, &lwork, iwork, &liwork, &info);
			
		for(j=0; j<p; j++){
			eig_val[j] = (-eig_val[j] + sqrt(pow(eig_val[j], 2) + 4*rho))/(2*rho);
		}
		for(j=0; j<p; j++){
			for(k=j; k<p; k++){
				index_jk = p*j+k;
				Omega[index_jk] = 0;
				for(m=0; m<p; m++){
					Omega[index_jk] += eig_val[m]*eig_vec[p*m+j]*eig_vec[p*m+k];
				}
				U[index_jk] += Omega[index_jk]; 
			}
			index_jk = p*j+j;
			s += pow((U[index_jk] - Z[index_jk]), 2)/2;
			Z[index_jk] = U[index_jk];
			U[index_jk] = 0;
			r += pow((Omega[index_jk] - Z[index_jk]), 2)/2;
			epi_pri_1 += pow(Omega[index_jk], 2)/2;
			epi_pri_2 += pow(Z[index_jk], 2)/2;
		}

		if(S_L!=0){
			for(m=0; m<S_L; m++){
				index_jk = p*S[2*m] + S[2*m+1];
				s += pow((U[index_jk] - Z[index_jk]), 2);
				Z[index_jk] = U[index_jk];
			}
		}

		for(j=0; j<p; j++){
			for(k=j+1; k<p; k++){
				index_jk = p*j+k;
				U[index_jk] -= Z[index_jk];
				r += pow((Omega[index_jk] - Z[index_jk]), 2);
				epi_pri_1 += pow(Omega[index_jk], 2);
				epi_pri_2 += pow(Z[index_jk], 2);
				epi_dual += pow(U[index_jk], 2);
			}
		}
		r = sqrt(2*r);
		s = rho*sqrt(2*s);
		epi_pri = sqrt(p)*epi_abs + epi_rel*sqrt(2*((epi_pri_1>epi_pri_2)?epi_pri_1:epi_pri_2));
		epi_dual = sqrt(p)*epi_abs + epi_rel*rho*sqrt(2*epi_dual);
		
		prd = (r>epi_pri);
		drd = (s>epi_dual);
	}//end ADMM iteration

	free(isuppz);
	free(work);
	free(iwork);
	free(Omega);
	free(A);
	free(eig_vec);
	free(eig_val);
}



void ADMM_pseudo_refit(int *P, int *Pos, double *Sigma, double *Z, double *Z_pos, int *S_Len){

	int p = *P, pos = *Pos, S_L;
	int nrhs = 1, info, j, k, m;
	int *S = (int *) malloc(p*sizeof(int));
	double temp;
	double *Sigma_j = (double *) malloc(p*p*sizeof(double));
	double *C = (double *) malloc(p*sizeof(double));

	for(j=0; j<p; j++){
		S_L = 0;
		for(k=0; k<p; k++){
			Z_pos[p*j+k] = 0;
			if(j!=k && Z[p*p*pos+p*j+k]!=0 && Z[p*p*pos+p*k+j]!=0){
				S[S_L] = k;
				S_L++;
				if(k>j){
					S_Len[0]++;
				}		
			}
		}
		if(S_L==0){
			Z_pos[p*j+j] = 1/Sigma[p*p*pos+p*j+j];
		}
		else{
			for(k=0; k<S_L; k++){
				for(m=k; m<S_L; m++){
					Sigma_j[S_L*k+m] = Sigma[p*p*pos+p*S[k]+S[m]];
				}
				C[k] = Sigma[p*p*pos+p*j+S[k]];
			}
			F77_CALL(dpotrf)("L", &S_L, Sigma_j, &S_L, &info);
			F77_CALL(dpotrs)("L", &S_L, &nrhs, Sigma_j, &S_L, C, &S_L, &info);
			for(k=0; k<S_L; k++){
				Z_pos[p*j+j] += C[k]*Sigma[p*p*pos+p*j+S[k]];
			}
			Z_pos[p*j+j] = 1/(-Z_pos[p*j+j]+Sigma[p*p*pos+p*j+j]);
			for(k=0; k<S_L; k++){
				if(S[k]>j){
					Z_pos[p*j+S[k]] = C[k];
				}
				else{
					temp = Z_pos[p*S[k]+j]*C[k];
					if(temp>0){
						if(Z[p*S[k]+j]>0){
							Z_pos[p*S[k]+j] = -sqrt(temp);
						}
						else{
							Z_pos[p*S[k]+j] = sqrt(temp);
						}
					}
					else{
						if(Z[p*S[k]+j]>0){
							Z_pos[p*S[k]+j] = -sqrt(-temp);
						}
						else{
							Z_pos[p*S[k]+j] = sqrt(-temp);
						}
					}
					Z_pos[p*S[k]+j] = Z_pos[p*S[k]+j]*sqrt(Z_pos[p*S[k]+S[k]]*Z_pos[p*j+j]);
					Z_pos[p*j+S[k]] = Z_pos[p*S[k]+j];
				}
			}
		}
	}

	free(S);
	free(Sigma_j);
	free(C);
}



void Givens_rotation(double *U, double *Chol, int *P, int *J){

	int p = *P, j = *J, p_1 = p-1, k, m;
	double r, c, s;

	for(k=0; k<j; k++){
		U[p_1*k+j-1] = Chol[p*k+j-1];
		for(m=j; m<p_1; m++){
			U[p_1*k+m] = Chol[p*k+m+1];
		}
	}
	for(m=j; m<p_1; m++){
		U[p_1*j+m] = Chol[p*j+m+1];
	}

	for(k=j; k<p_1; k++){
		r = sqrt(pow(U[p_1*k+k], 2)+pow(Chol[p*(k+1)+k+1], 2));
		c = U[p_1*k+k]/r;
		s = -Chol[p*(k+1)+k+1]/r;
		U[p_1*k+k] = r;	
		for(m=k+1; m<p_1; m++){
			U[p_1*(k+1)+m] = s*U[p_1*k+m]+c*Chol[p*(k+1)+m+1];
			U[p_1*k+m] = c*U[p_1*k+m]-s*Chol[p*(k+1)+m+1];
		}
	}
}
//apply ADMM on lambda grid given h and d
void ADMM_lambda(double *Corr, double *Z, int *P, int *Len, int *Pos, int *Pos_Len, int *member_ind, int *csize_ind, 
	int *No, double *Lambda, int *Lambda_Len, int *fit_type, double *thres, double *Rho, double *Epi_abs, 
	double *Epi_rel, int *Max_step);

//apply ADMM across block diagonals given h, d and lambda
void ADMM_cluster(double *Corr, double *Z, double *U, int *P, int *Len, int *member_ind, int *csize_ind, int *No, 
	double *Lambda, int *fit_type, double *Rho, double *Epi_abs, double *Epi_rel, int *Max_step);

//apply ADMM across block diagonals given h, d and lambda (no initialization of Z and U)
void ADMM_simple(double *Corr, double *Z, int *P, int *Len, int *Pos, int *Pos_Len, int *member_ind, int *csize_ind, 
	int *No, double *Lambda, int *fit_type, double *Rho, double *Epi_abs, double *Epi_rel, int *Max_step);

//ADMM on likelihood estimation
void ADMM_local_glasso(double *Sigma, double *Z, double *U, int *P, int *Len, double *Lambda, double *Rho, 
	double *Epi_abs, double *Epi_rel, int *Max_step);

//ADMM on pseudo-likelihood estimation
void ADMM_pseudo_glasso(double *Sigma, double *Z, double *U, int *P, int *Len, double *Lambda, double *Rho, 
	double *Epi_abs, double *Epi_rel, int *Max_step);

//ADMM on SPACE(sparse partial correlation estimation) (rho step)
void ADMM_SPACE_rho(double *Sigma, double *d, double *Z, double *U, int *P, int *Len, double *Lambda, double *Rho, 
	double *Epi_abs, double *Epi_rel, int *Max_step);

//ADMM on SPACE(sparse partial correlation estimation) (d step)
void ADMM_SPACE_d(double *Sigma, double *d, double *Z, int *P, int *Len);

//Givens rotation
void Givens_rotation(double *U, double *Chol, int *P, int *J);
//apply ADMM across lambda's when h and d are fixed
void ADMM_lambda(double *Corr, double *Z, int *P, int *LL, int *Pos, int *Pos_Len, int *member_ind, int *csize_ind, int *No, 
	double *Lambda, int *Lambda_Len, int *pseudo_fit, double *thres, double *Rho, double *Epi_abs, double *Epi_rel, int *Max_step);

//apply ADMM to fixed h, d and lambda
void ADMM_cluster(double *Corr, double *Z, double *U, int *P, int *LL, int *member_ind, int *csize_ind, int *No, double *Lambda, 
	int *pseudo_fit, double *Rho, double *Epi_abs, double *Epi_rel, int *Max_step);

//apply ADMM to fixed h, d and lambda (simple version)
void ADMM_simple(double *Corr, double *Z, int *P, int *LL, int *Pos, int *Pos_Len, int *member_ind, int *csize_ind, int *No, 
	double *Lambda, int *pseudo_fit, double *Rho, double *Epi_abs, double *Epi_rel, int *Max_step);

//local group graphical lasso
void ADMM_local_glasso(double *Sigma, double *Z, double *U, int *P, int *LL, double *Lambda, double *Rho, double *Epi_abs, 
	double *Epi_rel, int *Max_step);

//pseudo-likelihood group lasso
void ADMM_pseudo_glasso(double *Sigma, double *Z, double *U, int *P, int *LL, double *Lambda, double *Rho, double *Epi_abs, 
	double *Epi_rel, int *Max_step);

//SPACE (rho step)
void ADMM_SPACE_rho(double *Sigma, double *d, double *Z, double *U, int *P, int *LL, double *Lambda, double *Rho, double *Epi_abs, 
	double *Epi_rel, int *Max_step);

//SPACE (d step)
void ADMM_SPACE_d(double *Sigma, double *d, double *Z, int *P, int *LL);

//Givens rotation
void Givens_rotation(double *U, double *Chol, int *P, int *J);
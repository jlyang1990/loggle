//apply ADMM across lambda's when h and d are fixed
void ADMM_lambda(double *Corr, double *Z, int *P, int *Len, int *Pos, int *Pos_Len, int *member_ind, int *csize_ind, 
	int *No, double *Lambda, int *Lambda_Len, int *fit_type, double *thres, double *Rho, double *Epi_abs, 
	double *Epi_rel, int *Max_step, double *record);

//apply ADMM across block diagonals when h, d and lambda are fixed
void ADMM_cluster(double *Corr, double *Z, double *U, int *P, int *Len, int *member_ind, int *csize_ind, int *No, 
	double *Lambda, int *fit_type, double *Rho, double *Epi_abs, double *Epi_rel, int *Max_step, double *record);

//apply ADMM across block diagonals when h, d and lambda are fixed (simple version)
void ADMM_simple(double *Corr, double *Z, int *P, int *Len, int *Pos, int *Pos_Len, int *member_ind, int *csize_ind, 
	int *No, double *Lambda, int *fit_type, double *Rho, double *Epi_abs, double *Epi_rel, int *Max_step, double *record);

//local group graphical lasso (likelihood estimation)
void ADMM_local_glasso(double *Sigma, double *Z, double *U, int *P, int *Len, double *Lambda, double *Rho, 
	double *Epi_abs, double *Epi_rel, int *Max_step, double *record);

//pseudo-likelihood group lasso (pseudo-likelihood estimation)
void ADMM_pseudo_glasso(double *Sigma, double *Z, double *U, int *P, int *Len, double *Lambda, double *Rho, 
	double *Epi_abs, double *Epi_rel, int *Max_step, double *record);

//SPACE (sparse partial correlation estimation) (rho step)
void ADMM_SPACE_rho(double *Sigma, double *d, double *Z, double *U, int *P, int *Len, double *Lambda, double *Rho, 
	double *Epi_abs, double *Epi_rel, int *Max_step);

//SPACE (sparse partial correlation estimation) (d step)
void ADMM_SPACE_d(double *Sigma, double *d, double *Z, int *P, int *Len);

//Givens rotation
void Givens_rotation(double *U, double *Chol, int *P, int *J);
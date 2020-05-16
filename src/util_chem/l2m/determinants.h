
/** 
  *
  *  Set of routines for the computation of determinants 
  *
  *  February 2020
  *
  */

#ifndef _DETERMINANTS
#define _DETERMINANTS

//void naive_det( int n, double *A, int Alda, double *Ac, int Aclda, int *ipiv );
//void opti_det( int n, double *A, int Alda, double *B, int Blda, double *C, int Clda );
//void opti_det_wk( int n, double *A, int Alda, double *B, int Blda, double *workspace );

//void reduce_Hessenberg3( int m, int n, double *C, int Clda, double *d, double *workspace );

void getlvl2minors_opt (int no, const int nv2, double *csc, int ld_csc, double *l2minors, int ld_l2minors);

void reduce_Hessenberg( int m, int n, double *B, int Blda, double *s );
void reduce_Hessenberg_ip( int m, int n, double *C, int Clda, double *d );
void reduce_Hessenberg_ip_wk( int m, int n, double *C, int Clda, double *d, double *workspace );
void reduce_Hessenberg_wk_Gfast( int m, int n, double *C, int Clda, double *d, double *workspace );
void det_naive( int n, double *A, int Alda, double *Ac, int Aclda, int *ipiv, double *sum, double *fsum, double *table );
void det_base( int n, double *A, int Alda, double *B, int Blda, int transforms, double *sum, double *fsum );
void det_opt( int n, double *A, int Alda, double *B, int Blda, double *C, int Clda, int transforms, double *sum, double *fsum );
void det_opt_ip( int n, double *A, int Alda, double *B, int Blda, double *C, int Clda, int transforms, double *sum, double *fsum );
void det_opt_wk( int n, double *A, int Alda, double *B, int Blda, double *workspace, int transforms, double *sum, double *fsum, double *table );
void det_opt_rc( int n, double *A, int Alda, double *B, int Blda, double *workspace, int transforms, int r1, int r2, double s1, double *table );
void det_opt_rc_row( int n, double *A, int Alda, double *B, int Blda, double *workspace, int transforms, int r1, int r2, double s1, double *table );
void reduce_Hessenberg_G( int m, int n, double *B, int Blda, double *d );
void reduce_Hessenberg_Gfast( int m, int n, double *B, int Blda, double *d );
void reduce_Hessenberg2_G( int m, int n, double *B, int Blda, double *d );
void reduce_Hessenberg2( int m, int n, double *B, int Blda, double *d );
void reduce_Hessenberg_rc_Gfast_row( int m, int n, double *Ct, int Ctlda, double *d, double *workspace );
void reduce_Hessenberg_rc_Gfast( int m, int n, double *C, int Clda, double *d, double *workspace );
void reduce_Hessenberg_row( int m, int n, double *Bt, int Btlda, double *s );
void reduce_Hessenberg_ip_rc_row( int m, int n, double *Ct, int Ctlda, double *d, double *workspace );
void reduce_Hessenberg_ip_rc( int m, int n, double *Ct, int Ctlda, double *d, double *workspace );
void det_opt_wk_rows( int n, double *At, int Atlda, double *Bt, int Btlda, double *workspace, int transforms, double *sum, double *fsum, double *table );
void reduce_Hessenberg_ip_wk_row( int m, int n, double *Ct, int Ctlda, double *d, double *workspace );

#endif

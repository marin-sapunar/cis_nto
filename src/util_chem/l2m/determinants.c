

/** 
  *
  *  Set of routines for the computation of determinants 
  *
  *  February 2020
  *
  */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "determinants.h"

#define min(a,b)     ( (a) > (b) ? (b) : (a) )

#define Aref(a1,a2)  A[ (a2-1)*(Alda)+(a1-1) ]
#define Acref(a1,a2) Ac[ (a2-1)*(Aclda)+(a1-1) ]
#define Bref(a1,a2)  B[ (a2-1)*(Blda)+(a1-1) ]
#define Cref(a1,a2)  C[ (a2-1)*(Clda)+(a1-1) ]

#define Atref(a1,a2)  At[ (a1-1)*(Atlda)+(a2-1) ]
#define Actref(a1,a2) Act[ (a1-1)*(Actlda)+(a2-1) ]
#define Btref(a1,a2)  Bt[ (a1-1)*(Btlda)+(a2-1) ]
#define Ctref(a1,a2)  Ct[ (a1-1)*(Ctlda)+(a2-1) ]

#define Apref(a1)    workspace[ (a1-1) ]
#define Anpref(a1)   workspace[ (n+a1-1) ]

extern void dswap_( int *, double *, int *, double *, int * );
extern void daxpy_( int *, double *, double *, int *, double *, int * );
extern void dcopy_( int *, double *, int *, double *, int * );
extern void dlacpy_( char *, int *, int *, double *, int *, double *, int * );
extern void dgetrf_( int *, int *, double *, int *, int *, int * );
extern void drotg_( double *, double *, double *, double *);                  
/** 
Naive solution
*/
void det_naive( int n, double *A, int Alda, double *Ac, int Aclda, int *ipiv, double *sum, double *fsum, double *table ) {
  double d; 
  int    c1, c2, k, m = n-2, info;

  *sum = 0.0, *fsum = 0.0;
  for (c1 = 1; c1<=n-1; c1++){
  // for (c1 = 1; c1<=1; c1++){
    k = c1-1;
    dlacpy_("All", &m, &k, &Aref(1,1), &Alda, &Acref(1,1), &Aclda);
    for (c2 = c1+1; c2<=n; c2++){
    // for (c2 = 2; c2<=2; c2++){
      k = c2-c1-1;
      dlacpy_("All", &m, &k, &Aref(1,c1+1), &Alda, &Acref(1,c1),   &Aclda);
      k = n-c2;
      dlacpy_("All", &m, &k, &Aref(1,c2+1), &Alda, &Acref(1,c2-1), &Aclda);
      // print_matrix( "Aci", m, m, &Acref(1,1), Aclda );
      dgetrf_( &m, &m, &Acref(1,1), &Aclda, ipiv, &info);
      // print_matrix( "Acf", m, m, &Acref(1,1), Aclda );
      // exit(-1);
      // dgeqrf_( &m, &m, &Acref(1,1), &Aclda,... &info);
      d = 1.0;
      for (k = 1; k<=m; k++)
        d = d*Acref(k,k);
      // printf("d(%d,%d) = %22.15e;\n", c1, c2, d); fflush(stdout);
      *sum = (*sum) + d; *fsum = (*fsum)+fabs(d); table[(c1-1)*n+(c2-1)] = d;
    }
  }
}

void reduce_Hessenberg_ip_rc( int m, int n, double *C, int Clda, double *d, double *workspace ) {
// Reduces an mxn Hessenberg matrix to upper triangular form using Gauss transforms. 
// Performs implicit pivoting and employs reduced workspace (does not modify A)

  int k, p, mn, pivoting, ione = 1;
  double pivot, *ptrD, *ptrSD, *ptrTmp;
  double s;

  *d = 1.0;
  s = 1.0;
  mn = min(m,n)-1;

  if (mn<0)
    return;

  ptrD   = &Apref(1);
  ptrSD  = &Anpref(1);
  dcopy_(&n, &Cref(1,1), &Clda, &ptrD[0], &ione);

  for (k = 1; k <= mn; k++) {
    p = n-k+1;
    pivoting = fabs(Cref(k+1,k))>fabs(ptrD[k-1]);
    ptrSD[k-1] = Cref(k+1,k);
    
    if (pivoting) {
      pivot = ptrD[k-1]/ptrSD[k-1];
      *d = (*d)*ptrSD[k-1];
      s = -s;

      #pragma ivdep
      for (p = k+1; p<=n; p++)
         ptrSD[p-1] = ptrD[p-1] - pivot*Cref(k+1,p);
    }
    else {
      pivot = ptrSD[k-1]/ptrD[k-1];
      *d = (*d)*ptrD[k-1];

      #pragma ivdep
      for (p = k+1; p<=n; p++)
         ptrSD[p-1] = Cref(k+1,p) - pivot*ptrD[p-1];
    }

    ptrTmp = ptrD;
    ptrD = ptrSD;
    ptrSD = ptrTmp;
  }

  if (mn>=0)
    *d = (*d)*ptrD[mn];
  *d = s * (*d);
}

/** 
  Optimized solution, with implicit pivoting and workspace reduced and row storage
*/
void det_opt_rc_row( int n, double *At, int Atlda, double *Bt, int Btlda, double *workspace, int transforms, double *sum, double *fsum, double *table ) {
  double d, d1, d2, d3, s2;
  int    c1, c2, k, r, c, m = n-2;

  *sum = 0.0, *fsum = 0.0;
  if( n == 0 ) return;
  for (c1 = 1; c1 <= n-1; c1++){
    d1 = 1.0;
    for (k = 1; k <= c1-1; k++) {
      d1 = d1*Atref(k,k);
    }
    // print_matrixt("A1", m-c1+1, c1-1, Atref(1,1), Atlda);

    r    = m-c1+1; c = n-c1;
    Btlda = c;
    //dlacpy_("All", &r, &c, &Atref(c1,c1+1), &Atlda, &Btref(1,1), &Btlda);
    dlacpy_("All", &c, &r, &Atref(c1,c1+1), &Atlda, &Btref(1,1), &Btlda);
    // print_matrixt("A2", r, c, &Btref(1,1), Btlda);
    reduce_Hessenberg_row( r, c, &Btref(1,1), Btlda, &s2 );

    for (c2 = 1; c2 <= n-c1; c2++){
      d2 = 1.0;
      for (k = 1; k <= c2-1; k++)
        d2 = d2*Btref(k,k);
      d2 = d2*s2;

      r    = m-c1-c2+2;
      c    = n-c1-c2;
      // print_matrixt("A3", r, c, &Btref(c2,c2+1), Btlda);
      if (transforms==0)
        reduce_Hessenberg_ip_rc_row( r, c, &Btref(c2,c2+1), Btlda, &d3, workspace );
      // else if (transforms==2)
        // reduce_Hessenberg_rc_Gfast_row( r, c, &Btref(c2,c2+1), Btlda, &d3, workspace );
      else {
         printf("Unknown type of transforms\n");
         exit(-1);
      }

      d = d1*d2*d3;
      // printf("det(%d,%d) = %22.15e %22.15e %22.15e %22.15e;\n", c1, c1+c2, d1, d2, d3, d); fflush(stdout);
      *sum = (*sum) + d; *fsum = (*fsum)+fabs(d); table[(c1-1)*n+(c1+c2-1)] = d;

      // if (transforms==0) {
        // printf("det_opt_ip_rc(%d,%d) = %22.15e;\n", c1, c1+c2, d); fflush(stdout);
      // }
      // else if (transforms==2) {
        // printf("det_opt_rc_Gfast(%d,%d) = %22.15e;\n", c1, c1+c2, d); fflush(stdout);
      // }
    }
  }
}

void reduce_Hessenberg_rc_Gfast( int m, int n, double *C, int Clda, double *d, double *workspace ) {
// Reduces an mxn Hessenberg matrix to upper triangular form using Givens rotations. 
// Only updates bottom row and employes reduced workspace (does not modify A)
//
    int k, p, mn, pivoting, ione = 1;
    double gs, gc, *ptrD, *ptrSD, *ptrTmp;

    *d = 1.0;
    mn = min(m,n)-1;
    if (mn<0)
        return;

    ptrD   = &Apref(1);
    ptrSD  = &Anpref(1);
    dcopy_(&n, &Cref(1,1), &Clda, &ptrD[0], &ione);

    for (k = 1; k <= mn; k++) {
        ptrSD[k-1] = Cref(k+1,k);
                                         
       drotg_(&ptrD[k-1], &ptrSD[k-1], &gc, &gs); // drotg already applies the Givens transform to the arguments
       *d = (*d)*ptrD[k-1];

       #pragma ivdep
       for (p = k+1; p<=n; p++) {
           ptrSD[p-1] = gc*Cref(k+1,p) - gs*ptrD[p-1];
       }
                                                                             
       ptrTmp = ptrD;
       ptrD = ptrSD;
       ptrSD = ptrTmp;
   }
                                                                                             if (mn>=0)
       *d = (*d)*ptrD[mn];
}

void reduce_Hessenberg_rc_Gfast_row( int m, int n, double *Ct, int Ctlda, double *d, double *workspace ) {
// Reduces an mxn Hessenberg matrix to upper triangular form using Givens rotations, assuming row storage. 
// // Only updates bottom row and employes reduced workspace (does not modify A)
//
    int k, p, mn, pivoting, ione = 1;
    double gs, gc, *ptrD, *ptrSD, *ptrTmp;

    *d = 1.0;
    mn = min(m,n)-1;
    if (mn<0)
        return;

    //double *Ap  = malloc(n*sizeof(double));
    //double *Anp = malloc(n*sizeof(double));

    // Ap  = A(1,:);
    ptrD   = &Apref(1);
    ptrSD  = &Anpref(1);
    dcopy_(&n, &Ctref(1,1), &ione, &ptrD[0], &ione);

    for (k = 1; k <= mn; k++) {
        // dcopy_(&p, &Cref(k+1,k), &Clda, &ptrSD[k-1], &ione);
        ptrSD[k-1] = Ctref(k+1,k);
                                         
        drotg_(&ptrD[k-1], &ptrSD[k-1], &gc, &gs); // drotg already applies the Givens transform to the arguments
        *d = (*d)*ptrD[k-1];

        #pragma ivdep
        for (p = k+1; p<=n; p++) {
            // ptrSD[p-1] = gc*ptrSD[p-1] - gs*ptrD[p-1];
            ptrSD[p-1] = gc*Ctref(k+1,p) - gs*ptrD[p-1];
        }
                                                                             
/*
        p = n-k; dcopy_(&p, &ptrD[k], &ione, &ptrSD[k], &ione); double mgs = -gs; dscal_(&p, &mgs, &ptrSD[k], &ione); daxpy_(&p, &gc, &Cref(k+1,k+1), &Clda, &ptrSD[k], &ione);
*/
/*
        dcopy_(&p, &Cref(k+1,k+1), &Clda, &ptrSD[k], &ione); dscal_(&p, &gc, &ptrSD[k], &ione); double mgs = -gs; daxpy_(&p, &mgs, &ptrD[k], &ione, &ptrSD[k], &ione);
*/
        ptrTmp = ptrD;
        ptrD = ptrSD;
        ptrSD = ptrTmp;
    }

    if (mn>=0)
        *d = (*d)*ptrD[mn];
}

void reduce_Hessenberg_row( int m, int n, double *Bt, int Btlda, double *s ) {

// Reduces an mxn Hessenberg matrix to upper triangular form using Gauss transforms, assuming row storage

  int k, p, ione = 1;
  double pivot, gc, gs;

  *s = 1.0;
  for (k = 1; k<=min(m,n)-1; k++){
    // Gauss. If pivoting is enabled, the determinants are scaled by +/-1 due to permutations. 
    // Keep track of permutations using s
    //*s = Bref(k+1,k);
    if (fabs((Btref(k+1,k)))>fabs((Btref(k,k)))){
      *s = -(*s);
      p  = n-k+1;
      dswap_(&p, &Btref(k,k), &ione, &Btref(k+1,k), &ione);
    }
    pivot        = -Btref(k+1,k)/Btref(k,k);
    Btref(k+1,k)  = 0.0;
    p = n-k;
    daxpy_(&p, &pivot, &Btref(k,k+1), &ione, &Btref(k+1,k+1), &ione);
  }
}

/** 
 *   Optimized solution, with implicit pivoting and workspace reduced 
 *   */
void det_opt_rc( int n, double *A, int Alda, double *B, int Blda, double *workspace, int transforms, int r1, int r2, double s1, double *table ) {
  
	double d, d1, d2, d3, s2, s3;
	int    c1, c2, k, r, c, m = n-2;
	int pos1, pos2, pos3, pos4;
	unsigned int n2, n3;
	int c2_ext;

	n2 = n * n;
	n3 = n2 * n;

	unsigned long int table_dim = n*n*n*n;

	if( n == 0 ) return;
		
	for (c1 = 1; c1 <= n-1; c1++){
		d1 = 1.0;
		for (k = 1; k <= c1-1; k++)
			d1 = d1*Aref(k,k);

		r = m-c1+1; 
		c = n-c1;
		Blda = r;
		dlacpy_("All", &r, &c, &Aref(c1,c1+1), &Alda, &Bref(1,1), &Blda);
		reduce_Hessenberg( r, c, &Bref(1,1), Blda, &s2 );

		for (c2 = 1; c2 <= n-c1; c2++){
			d2 = 1.0;
			for (k = 1; k <= c2-1; k++)
				d2 = d2*Bref(k,k);
			d2 = d2*s2;

			r = m-c1-c2+2;
			c = n-c1-c2;
					
			if (transforms==0)
				reduce_Hessenberg_ip_rc( r, c, &Bref(c2,c2+1), Blda, &d3, workspace );
			else if (transforms==2)
				reduce_Hessenberg_rc_Gfast( r, c, &Bref(c2,c2+1), Blda, &d3, workspace );
			else {
				printf("Unknown type of transforms\n");
				exit(-1);
			}

			// The final determinant with included signs from Gauss rotations
			d = s1 * d1*d2*d3;

			// Extend c2 index to a 1:n range
			c2_ext = c2 + c1;

			/* Compute positions in the determinant table */
			pos1 = (r1-1)*n3 + (r2-1)*n2 + (c1-1)*n + c2_ext - 1;
			pos2 = (r1-1)*n3 + (r2-1)*n2 + (c2_ext-1)*n + c1 - 1;
			pos3 = (r2-1)*n3 + (r1-1)*n2 + (c1-1)*n + c2_ext - 1;
			pos4 = (r2-1)*n3 + (r1-1)*n2 + (c2_ext-1)*n + c1 - 1;

			/* Fill-in the table with the determinants with the correct sign */ 
			/*if (c1+r1 % 2 == 0) {
				table[pos1] = d;
			} else {
				table[pos1] = -d;
			}
			if (c2_ext+r1+1 % 2 == 0) {
				table[pos2] = d;
			} else {
				table[pos2] = -d;
			}
			if(c1+r2+1 % 2 == 0) {
				table[pos3] = d;
			} else {
				table[pos3] = -d;
			}
			if (c2_ext+r2 % 2 == 0) {
				table[pos4] = d;
			} else {
				table[pos4] = -d;
			}
			*/
				
			table[pos1] = pow(-1.0, c1+r1) * d;
			table[pos2] = pow(-1.0, c2_ext+r1+1) * d;
			table[pos3] = pow(-1.0, c1+r2+1) * d;
			table[pos4] = pow(-1.0, c2_ext+r2) * d;
		}
	}
}

void reduce_Hessenberg( int m, int n, double *B, int Blda, double *s ) {
// Reduces an mxn Hessenberg matrix to upper triangular form using Gauss transforms

	int k, p;
	double pivot, gc, gs;

	*s = 1.0;

	for (k = 1; k<=min(m,n)-1; k++){
	// Gauss. If pivoting is enabled, the determinants are scaled by +/-1 due to permutations. 
	// Keep track of permutations using s
		if (fabs((Bref(k+1,k)))>fabs((Bref(k,k)))){
			*s = -(*s);
			p  = n-k+1;
			dswap_(&p, &Bref(k,k), &Blda, &Bref(k+1,k), &Blda);
		}
		pivot = -Bref(k+1,k)/Bref(k,k);
		Bref(k+1,k)  = 0.0;
		p = n-k;
		daxpy_(&p, &pivot, &Bref(k,k+1), &Blda, &Bref(k+1,k+1), &Blda);
	}
}

void reduce_Hessenberg_ip_rc_row( int m, int n, double *Ct, int Ctlda, double *d, double *workspace ) {
// Reduces an mxn Hessenberg matrix to upper triangular form using Gauss transforms, assuming row storage. 
// Performs implicit pivoting and employs reduced workspace (does not modify A)

  int k, p, mn, pivoting, ione = 1;
  double pivot, *ptrD, *ptrSD, *ptrTmp;
  double s;

  *d = 1.0;
  s = 1.0;
  mn = min(m,n)-1;
  if (mn<0)
    return;

  // Ap  = A(1,:);
  ptrD   = &Apref(1);
  ptrSD  = &Anpref(1);
  dcopy_(&n, &Ctref(1,1), &ione, &ptrD[0], &ione);

  for (k = 1; k <= mn; k++) {
    p = n-k+1;
    pivoting = fabs(Ctref(k+1,k))>fabs(ptrD[k-1]);
    // dcopy_(&p, &Cref(k+1,k), &Clda, &ptrSD[k-1], &ione);
    ptrSD[k-1] = Ctref(k+1,k);
    
    // drotg_(&ptrD[k-1], &ptrSD[k-1], &gc, &gs); // drotg already applies the Givens transform to the arguments
    if (pivoting) {
      pivot = ptrD[k-1]/ptrSD[k-1];
      *d = (*d)*ptrSD[k-1];
      s = -s;

      #pragma ivdep
      for (p = k+1; p<=n; p++)
         ptrSD[p-1] = ptrD[p-1] - pivot*Ctref(k+1,p);
    }
    else {
      pivot = ptrSD[k-1]/ptrD[k-1];
      *d = (*d)*ptrD[k-1];

      #pragma ivdep
      for (p = k+1; p<=n; p++)
         ptrSD[p-1] = Ctref(k+1,p) - pivot*ptrD[p-1];
    }

    ptrTmp = ptrD;
    ptrD = ptrSD;
    ptrSD = ptrTmp;
  }

  if (mn>=0)
    *d = (*d)*ptrD[mn];

  *d = s * (*d);
}


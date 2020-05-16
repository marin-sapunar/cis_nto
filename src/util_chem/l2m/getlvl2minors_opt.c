#define _GNU_SOURCE
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<sched.h>
#include<omp.h>

#include "determinants.h"

#ifdef OPENBLAS
#include <cblas.h>
#elif MKL
#include <mkl.h>
#endif

#define TEST 1
#define TASKPAR

#define Aref(a1,a2)    A[ (a2-1)*(Alda)+(a1-1) ]
#define Apref(a1,a2)   Ap[ my_thread_id*n*n+(a2-1)*(Aplda)+(a1-1) ]
#define Bpref(a1,a2)   Bp[ my_thread_id*n*n+(a2-1)*(Bplda)+(a1-1) ]

#define Atpref(a1,a2)  Ap[ my_thread_id*n*n+(a1-1)*(Atplda)+(a2-1) ]

#define ipivref(a1)    ipivp[ my_thread_id*n+(a1-1) ]
#define dabs(a)        ( (a) > 0.0 ? (a) : -(a) )
#define min(a,b)       ( (a) > (b) ? (b) : (a) )

void getlvl2minors_opt (int n, const int nv2, double *A, int Alda, double *l2minors, int ld_l2minors) {

	int my_thread_id;
	cpu_set_t my_set;

	double *Ap = NULL;
	double *Bp = NULL;
	double *wk = NULL;
	int *ipivp = NULL;
    
	int Aplda = n;
	int Atplda = n;
	int Bplda = n;
	int nm2 = n - 2;
	int info = 0;

	double zero = 0.0;
	int r1, r2, c1, c2;
	int k;
	int i, j;
	double s1;
	int thmax = omp_get_max_threads();
    int blas_threads;

	/* Allocate space for data */
	Ap = (double *) malloc(n*n*thmax*sizeof(double));
	Bp = (double *) malloc(n*n*thmax*sizeof(double));
	wk = (double *) malloc(2*n*thmax*sizeof(double));
	ipivp = (int *) malloc(n*thmax*sizeof(int));
	
	/* Start parallel region */
	if(TEST) {
		for (i=0; i<n*n*n*n; i++) {
			l2minors[i] = 0.0;
		}
	}

#ifdef OPENBLAS
    blas_threads = openblas_get_num_threads();
    openblas_set_num_threads(1);
#elif MKL
    blas_threads = mkl_get_max_threads();
    mkl_set_num_threads(1);
#endif    

	#pragma omp parallel private(my_thread_id, my_set)
	{
		my_thread_id = omp_get_thread_num();
		CPU_ZERO(&my_set);	/* initialize it all to 0, i.e. no CPUs selected. */
		CPU_SET(my_thread_id, &my_set);	/* Set the bit that represents the core with the same id */
		sched_setaffinity(0, sizeof(cpu_set_t), &my_set); /* Set affinity of this process to */

#ifdef TASKPAR
		#pragma omp for collapse(2) schedule(dynamic) private (r2, s1, k)
#else
		#pragma omp single
		{
#endif
		for (r1 = 1; r1 <= n-1; r1++) {
			for (r2=1; r2 <= n; r2++) {
				
				if( r2 >= r1+1 ) {

#ifdef TASKPAR
					#pragma omp task firstprivate (r1, r2) private (k, my_thread_id, s1, i)
#endif
					{
						my_thread_id = omp_get_thread_num();

						/* Copy first r1-1 rows */
						k = r1 - 1;
						dlacpy_("All", &k, &n, &Aref(1,1), &Alda, &Apref(1,1), &Aplda);

						/* Copy rows from r1+1 till r2 */
						k = r2 - r1 - 1;
						dlacpy_("All", &k, &n, &Aref(r1+1,1), &Alda, &Apref(r1,1), &Aplda);

						/* Copy rows from r2+1 till end */
						k = n - r2;
						dlacpy_("All", &k, &n, &Aref(r2+1,1), &Alda, &Apref(r2-1,1), &Aplda);

						/* LU of the matrix A without rows r1 and r2 */
						dgetrf_(&nm2, &n, &Apref(1,1), &Aplda, &ipivref(1), &info);

						/* Compute Pivoting LU factorization sign */
						s1 = 1.0;
						for (i = 0; i < nm2 ; i++) {
							if(ipivref(1+i) != i+1) {
								s1 = -s1;
							}
						}

						k = nm2 - 1;
						dlaset_("Lower", &k, &n, &zero, &zero, &Apref(2,1), &Aplda);

                        for( i = 1; i <= nm2; i++ ){
                            for( j = 1; j <= n; j++ ){
                                Atpref(i,j) = Apref(i,j);
                            }
        
                        }
	
#ifdef GIVENS
                        /* Column wise */
						//det_opt_rc( n, &Apref(1,1), Aplda, &Bpref(1,1), Bplda, &wk[my_thread_id*2*n], 2,
						//	r1, r2, s1, l2minors);

                        /* Row wise */
						det_opt_rc_row( n, &Atpref(1,1), Atplda, &Bpref(1,1), Bplda, &wk[my_thread_id*2*n], 2,
							r1, r2, s1, l2minors);
#endif
#ifdef GAUSS
                        /* Column wise */
						//det_opt_rc( n, &Apref(1,1), Aplda, &Bpref(1,1), Bplda, &wk[my_thread_id*2*n], 0,
						//	r1, r2, s1, l2minors);

                        /* Row wise */
						det_opt_rc_row( n, &Atpref(1,1), Atplda, &Bpref(1,1), Bplda, &wk[my_thread_id*2*n], 0,
							r1, r2, s1, l2minors);
#endif
					}
				}
			}
		}
#ifndef TASKPAR
		}
#endif
	}	

#ifdef OPENBLAS
    openblas_set_num_threads(blas_threads);
#elif MKL
    mkl_set_num_threads(blas_threads);
#endif    

	/* Free data */
	free(Ap);
	free(Bp);
	free(wk);
	free(ipivp);

	return;
}

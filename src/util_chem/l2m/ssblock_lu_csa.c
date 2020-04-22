/*--------------------------------------------------------------------------------------------------
  FUNCTION: ssblock_lu

  DESCRIPTION: 
--------------------------------------------------------------------------------------------------
*/

#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sched.h>
#include <omp.h>

#ifdef OPENBLAS
#include <cblas.h>
#elif MKL
#include <mkl.h>
#endif

#define max(a, b) ( a < b ? b : a )

#define Tmpref(a1)      tmp[ my_thread_id*ld_tmp*tmp_cols*ns2 +a1 ]
#define Msumref(a1)     msum[ my_thread_id*ld_msum*nv2 + a1 ]

extern void getlvl2minors_lu_csa_(int*, const int*, double*, int*, double*, int*);

void ssblock_lu(double *csc, int no, const int nv1, const int nv2, const int ns1, const int ns2, double *wf1, double *wf2, double *l1minor, double *ss)
{
	int ld_csc;
    int ld_msum, ld_tmp;
    int ld_l2minors;
	int tmp_cols;
    int i, o1, o2;
    int pos_tmp, pos_msum;
    int ss_size;
	int num_threads;
    int thmax;
	double *msum = NULL;  //Auxiliary matrix to store computed determinants
	double *tmp = NULL;  // Auxiliary matrix for intermediate results
	double *l2minors = NULL;
	double zero = 0.0, one = 1.0, neg_one = -1.0;
	int    one_i = 1;
	double coef;
	double ss_local;
	int v1;
	
    /* OpenMP variables */
    int my_thread_id;
    cpu_set_t my_set;

    /* Timing routines */
    int n_timers = 6;
    double start, finish;
    double times[n_timers]; 

    /* Set leading dimensions */
    ld_csc =  no + nv1;
    ld_l2minors = no*no*no;
	ld_tmp = max(no, nv1);
	tmp_cols = nv2;
	ld_msum = nv1;
    ss_size = ns1 * ns2;

	/* Set timings to zero */
	for( i = 0; i < n_timers; i++ ){
		times[i] = zero;
	}

    thmax = omp_get_max_threads();

    /* Allocate arrays */
    msum = (double*) malloc (thmax * ld_msum * nv2 * sizeof(double));
    tmp = (double*) malloc(thmax * ld_tmp * tmp_cols * ns2 * sizeof(double));
	l2minors = (double*) malloc(ld_l2minors * no * sizeof(double));

	/* Compute the determinants of all possible l2minors of the reference matrix. Pass through all combinations of rows and columns (r1, c1, r2, c2) to removed in order to form a l2minor */
	printf("[ssblock_lu] Computing l2minors...\n");
	start = omp_get_wtime();
	getlvl2minors_lu_csa_(&no, &nv2, csc, &ld_csc, l2minors, &ld_l2minors);
	finish = omp_get_wtime();
	times[4] = finish - start;

    /* Scale lower csc matrix block with alternating array (1,-1,1,-1...) i.e. starting from the second 
    * column, scale every second column with -1.
    * Arrays wf1 and wf2 are used as auxiliary arrays and are not referenced (unchaged on exit).
	*/
    start = omp_get_wtime();
	cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, nv1, no/2, nv2, zero, wf1, nv1, wf2, nv2, neg_one, &csc[no+nv1 + no], 2*ld_csc);

	/* Scale ss matrix to zero */
    cblas_dscal(ns1*ns2, zero, ss, 1);

    finish = omp_get_wtime();
    times[5] += finish - start;

#ifdef OPENBLAS
	num_threads = openblas_get_num_threads();
    openblas_set_num_threads(1);
#elif MKL
	num_threads = mkl_get_max_threads();
    mkl_set_num_threads(1);
#endif    

    #pragma omp parallel default(shared) private(my_thread_id, my_set)
    {

        my_thread_id = omp_get_thread_num();
        CPU_ZERO(&my_set);
        CPU_SET(my_thread_id, &my_set);
		sched_setaffinity(0, sizeof(cpu_set_t), &my_set); /* Set affinity of this process to */

	    /* Start iterate through the blocks */
        #pragma omp for collapse(2) schedule(dynamic) private(o2, ld_tmp, pos_tmp, pos_msum, coef, i, v1, start, finish) reduction(+:ss[:ss_size]) reduction(+:times[:4]) 
	    for( o1 = 0; o1 < no; o1++ ){
	        for( o2 = 0; o2 < no; o2++ ){

		        printf("[ssblock_lu] Computing block (%d, %d)\n", o1, o2);

		        /* Apply rows from the right of the ref matrix, cscmat positions (1:no, no+1:np+nv2)
		        * Compute first product l2minors * csc(1:no, no+1:no+nv2) */
        	    start = omp_get_wtime(); 
		        pos_tmp = o1*no*no + o2*no;
		        ld_tmp = no;

        	    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, no, nv2, no, one, &l2minors[pos_tmp], ld_l2minors, &csc[no*ld_csc], ld_csc, zero, &Tmpref(0), ld_tmp);
        	
		        finish = omp_get_wtime();
        	    times[0] += finish - start;

        	    start = omp_get_wtime();
                pos_msum = o1*no*nv1 + o2*nv1;

		        cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, nv1, nv2, no, pow(neg_one, o1+2), &csc[no], ld_csc, &Tmpref(0), ld_tmp, zero, &Msumref(0), ld_msum);

        		/* Get l1minor determinants and apply permutation sign */
		        coef = pow(neg_one, o1+o2) * l1minor[o2*no + o1];

        		/* Get coefficients (no+1:no+nv1, no+1:no+nv2) into auxiliary storage msum */
		        for( i = 0; i < nv2; i++ ){
			        cblas_daxpy(nv1, coef, &csc[(no+i)*ld_csc+no], 1, &Msumref(i*ld_msum), 1);
		        }
        	    finish = omp_get_wtime();
        	    times[1] += finish - start;

        		/* Compute ss matrix  */
		        for( v1 = 0; v1 < nv1; v1++ ){ 

			        /* Compute product wf2 * wf1 */
			        start = omp_get_wtime();

			        cblas_dgemv( CblasColMajor, CblasTrans, nv2, ns2, one, &wf2[o2*nv2], no*nv2, &Msumref(v1), nv1, zero, &Tmpref(0), 1 );

        		    finish = omp_get_wtime();
			        times[2] += finish - start;

        			/* Update ss = msum * wf_prod */
		        	start = omp_get_wtime();

			        cblas_dger( CblasColMajor, ns1, ns2, one, &wf1[o1*nv1+v1], no*nv1, &Tmpref(0), one, ss, ns1 );
        			finish = omp_get_wtime();
		    	    times[3] +=  finish - start;
		        }
            }    
	    }
	}

    /* Average time */
    for( i = 0; i < 4; i++) {
        times[i] = times[i] / thmax; 
    }

	/* Free allocated spaces */
	free(msum);
	free(tmp);
	free(l2minors);

	/* Print timings */
    printf("[SSBLOCK] Execution times:\n");
    printf("\tscale:         %.7lf sec\n", times[5]);
    printf("\tl2minors:      %.7lf sec\n", times[4]),
    printf("\tl2minor * wf2: %.7lf sec\n", times[0]),
    printf("\twf1 * tmp:     %.7lf sec\n", times[1]);
    printf("\twf1 * wf2:     %.7lf sec\n", times[2]);
    printf("\tmsum * tmp:    %.7lf sec\n", times[3]);
} 

/*--------------------------------------------------------------------------------------------------
  FUNCTION: ssblock_lu

  DESCRIPTION: 
--------------------------------------------------------------------------------------------------
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

#ifdef OPENBLAS
#include <cblas.h>
#elif MKL
#include <mkl.h>
#endif

#define max(a, b) ( a < b ? b : a )

extern void getlvl2minors_lu_(int*, const int*, double*, int*, double*, int*);

void ssblock_lu(double *csc, int no, const int nv1, const int nv2, const int ns1, const int ns2, double *wf1, double *wf2, double *l1minor, double *ss, const int print_level)
{
	int ld_csc;
        int ld_msum, ld_tmp;
        int ld_l2minors;
	int tmp_cols;
        int i, o1, o2;
        int pos_tmp, pos_msum;
        int st1, st2;
	int num_threads;
	double *msum = NULL;  //Auxiliary matrix to store computed determinants
	double *tmp = NULL;  // Auxiliary matrix for intermediate results
	double *l2minors = NULL;
	double zero = 0.0, one = 1.0, neg_one = -1.0;
	int    one_i = 1;
	double coef;
	double ss_local;
	int v1;
	
        /* Timing routines */
        double start, finish;
        double time00, time0;
        double time_l2m, time_det, time_sum, time_tot;

	/* Auxiliary variables */
	int ndevices;
	int devices[10];

	if (print_level >= 2) {
	    time00 = omp_get_wtime();
            printf("\n         ---- start ssblock subroutine ----\n");
	}
        /* Set leading dimensions */
        ld_csc =  no + nv1;
        ld_l2minors = no*no*no;
	ld_tmp = max(no, nv1);
	tmp_cols = nv2;
	ld_msum = nv1;

        /* Allocate arrays */
        msum = (double*) malloc (ld_msum * nv2 * sizeof(double));
        tmp = (double*) malloc(ld_tmp * tmp_cols * ns2 * sizeof(double));
	l2minors = (double*) malloc(ld_l2minors * no * sizeof(double));

#ifdef OPENBLAS
	num_threads = openblas_get_num_threads();
#elif MKL
	num_threads = mkl_get_max_threads();
#endif    

	/* Compute the determinants of all possible l2minors of the reference matrix. Pass through all combinations of rows and columns (r1, c1, r2, c2) to removed in order to form a l2minor */
	time0 = omp_get_wtime();
	getlvl2minors_lu_(&no, &nv2, csc, &ld_csc, l2minors, &ld_l2minors);
	time_l2m = omp_get_wtime() - time0;

        /* Scale lower csc matrix block with alternating array (1,-1,1,-1...) i.e. starting from the second 
         * column, scale every second column with -1.
         * Arrays wf1 and wf2 are used as auxiliary arrays and are not referenced (unchaged on exit).
	 */
	cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, nv1, no/2, nv2, zero, wf1, nv1, wf2, nv2, neg_one, &csc[no+nv1 + no], 2*ld_csc);

	/* Scale ss matrix to zero */
        cblas_dscal(ns1*ns2, zero, ss, 1);


	/* Start iterate through the blocks */
	for( o1 = 0; o1 < no; o1++ ){
	    for( o2 = 0; o2 < no; o2++ ){

	/*	printf("[ssblock_lu] Computing block (%d, %d)\n", o1, o2); */

		/* Apply rows from the right of the ref matrix, cscmat positions (1:no, no+1:np+nv2)
		 * Compute first product l2minors * csc(1:no, no+1:no+nv2) */
        	time0 = omp_get_wtime(); 
		pos_tmp = o1*no*no + o2*no;
		ld_tmp = no;

        	cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, no, nv2, no, one, &l2minors[pos_tmp], ld_l2minors, &csc[no*ld_csc], ld_csc, zero, tmp, ld_tmp);
        	

                pos_msum = o1*no*nv1 + o2*nv1;

		cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, nv1, nv2, no, pow(neg_one, o1+2), &csc[no], ld_csc, tmp, ld_tmp, zero, msum, ld_msum);

		/* Get l1minor determinants and apply permutation sign */
		coef = pow(neg_one, o1+o2) * l1minor[o2*no + o1];

		/* Get coefficients (no+1:no+nv1, no+1:no+nv2) into auxiliary storage msum */
		for( i = 0; i < nv2; i++ ){
			cblas_daxpy(nv1, coef, &csc[(no+i)*ld_csc+no], 1, &msum[i*ld_msum], 1);
		}
        	time_det = omp_get_wtime() - time0;


	        time0 = omp_get_wtime();
		/* Compute ss matrix  */
		for( v1 = 0; v1 < nv1; v1++ ){ 
			cblas_dgemv( CblasColMajor, CblasTrans, nv2, ns2, one, &wf2[o2*nv2], no*nv2, &msum[v1], nv1, zero, tmp, 1 );
			cblas_dger( CblasColMajor, ns1, ns2, one, &wf1[o1*nv1+v1], no*nv1, tmp, one, ss, ns1 );
		}	
        	finish = omp_get_wtime();
		time_sum =  omp_get_wtime() - time0;
	    }
	}

	/* Free allocated spaces */
	free(msum);
	free(tmp);
	free(l2minors);

	/* Print timings */
	if (print_level >= 2) {
            printf("\n         ssblock time:\n");
	    printf("             L2 minors                  - time (sec): %.4lf \n",  time_l2m);
	    printf("             Determinants from minors   - time (sec): %.4lf \n",  time_det);
	    printf("             Final sum                  - time (sec): %.4lf \n",  time_sum);
	    printf("                                                     --------------\n");
	    printf("             Total                      - time (sec): %.4lf \n",  time_tot);
            printf("         ---- end ssblock subroutine ----\n");
	}

} 

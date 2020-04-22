/*--------------------------------------------------------------------------------------------------
  FUNCTION: ssblock_blas

  DESCRIPTION: 
--------------------------------------------------------------------------------------------------
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

#ifdef OPENBLAS
#include <cblas.h>
#endif
#ifdef MKL
#include <mkl.h>
#endif

extern void getlvl2minors_blk_(int*, int*, int*, int*, double*, double*, int*);


void ssblock_blas(double *csc, int no, const int nv1, const int nv2, const int ns1, const int ns2, double *wf1, double *wf2, double *l1minor, double *ss)
{
	int ld_csc;
        int ld_msum, ld_tmp;
        int ld_l2minors;
	int tmp_cols;
        int i, o1, o2, c, r;
	int dimmat;
        int pos_tmp, pos_msum;
        int st1, st2;
        int num_threads;
	double *msum = NULL;  //Auxiliary matrix to store computed determinants
	double *tmp = NULL;  // Auxiliary matrix for intermediate results
	double *msum_sorted = NULL; // Auxiliary matrix to store sorted msum values
	double *l2minors = NULL;
	double *wf1_col, *wf2_col = NULL; // Vectors holding current colums of wf1/wf2 
	double zero = 0.0, one = 1.0, neg_one = -1.0;
	double coef;
	double ss_local;

        /* Timing routines */
        double start, finish;
        double times[10]; 

	/* Auxiliary variables */
	int ndevices;
	int devices[10];

        /* Set leading dimensions */
        ld_csc =  no + nv1;
        ld_l2minors = no;
	ld_tmp = (no < nv2) ? nv2 : no;
	tmp_cols = (nv2 < nv1) ? nv1 : nv2;
	ld_msum = nv1;

	/* Set timings to zero */
	for( i = 0; i < 10; i++ ){
		times[i] = zero;
	}

        /* Allocate arrays */
	start = omp_get_wtime();
        msum = (double*) malloc (ld_msum * nv2 * sizeof(double));
        tmp = (double*) malloc(ld_tmp * tmp_cols * sizeof(double));
	l2minors = (double*) malloc(ld_l2minors * no * sizeof(double));
	wf1_col = (double*) malloc(nv1 * sizeof(double));
	wf2_col = (double*) malloc(nv2 * sizeof(double));
	times[0] += omp_get_wtime() - start;

#ifdef OPENBLAS
        num_threads = openblas_get_num_threads();
        printf("Number of openblas threads for LAPACK calls: %d\n", num_threads);
#elif MKL
	num_threads = mkl_get_max_threads();
        printf("Number of intel mkl threads for LAPACK calls: %d\n", num_threads);
#endif    

        /* Scale lower csc matrix block with alternating array (1,-1,1,-1...) i.e. starting from the second 
         * column, scale every second column with -1.
         * Arrays wf1 and wf2 are used as auxiliary arrays and are not referenced (unchaged on exit).
	 */
        start = omp_get_wtime();
	cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, nv1, no/2, nv2, zero, wf1, nv1, wf2, nv2, neg_one, &csc[no+nv1 + no], 2*ld_csc);

	/* Scale ss matrix to zero */
        cblas_dscal(ns1*ns2, zero, ss, 1);

        finish = omp_get_wtime();
        times[1] += finish - start;

	/* Start iterate through the blocks */
	for( o1 = 0; o1 < no; o1++ ){
	    for( o2 = 0; o2 < no; o2++ ){

		printf("[ssblock_blas] Computing block (%d, %d)\n", o1, o2);
		/* OpenMP parallel block - getlvl2minors
		 * set to single-threaded OpenBLAS/MKL kernels
		 */
#ifdef OPENBLAS
		openblas_set_num_threads(1);
#elif MKL
		mkl_set_num_threads(1);
#endif
		
		/* Compute block (o1,o2) with l2minors */
		r = o1 + 1;
		c = o2 + 1;
		dimmat = no+nv1;
		start = omp_get_wtime();
		getlvl2minors_blk_(&r, &c, &no, &dimmat, csc, l2minors, &num_threads);
		finish = omp_get_wtime();
		times[9] += finish - start;

		/* Restore the original number of threads */
#ifdef OPENBLAS
		openblas_set_num_threads(num_threads);
#elif MKL
		mkl_set_num_threads(num_threads);
#endif

		/* Apply rows from the right of the ref matrix, cscmat positions (1:no, no+1:np+nv2)
		 * Compute first product l2minors * csc(1:no, no+1:no+nv2) */
        	start = omp_get_wtime(); 
		pos_tmp = o1*no*no + o2*no;

        	cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, no, nv2, no, one, l2minors, ld_l2minors, &csc[no*ld_csc], ld_csc, zero, tmp, ld_tmp);
        	
		finish = omp_get_wtime();
        	times[2] += finish - start;

		/* Apply columns from below of the ref matrix, cscmat positions (no+1:no+nv1, 1:no)*/
        	start = omp_get_wtime();
                pos_msum = o1*no*nv1 + o2*nv1;

		cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, nv1, nv2, no, pow(neg_one, o1+2), &csc[no], ld_csc, tmp, ld_tmp, zero, msum, ld_msum);

		/* Get l1minor determinants and apply permutation sign */
		coef = pow(neg_one, o1+o2) * l1minor[o2*no + o1];

		/* Get coefficients (no+1:no+nv1, no+1:no+nv2) into auxiliary storage msum */
		for( i = 0; i < nv2; i++ ){
			cblas_daxpy(nv1, coef, &csc[(no+i)*ld_csc+no], 1, &msum[i*ld_msum], 1);
		}
        	finish = omp_get_wtime();
        	times[3] += finish - start;

		/* Compute ss matrix  */
	        for( st1 = 0; st1 < ns1; st1++ ){

			/* Copy columns wf1(o1*nv1, st1) to wf1_col, i.e. to a 1-dimensional array */
			start = omp_get_wtime();
			cblas_dcopy(nv1, &wf1[st1*no*nv1 + o1*nv1], 1, wf1_col, 1);
        		finish = omp_get_wtime();
			times[8] += finish - start;

			for( st2 = 0; st2 < ns2; st2++ ){

				/* Copy columns wf2(n2*nv2, st2) to wf2_col */
			        start = omp_get_wtime();
				cblas_dcopy(nv2, &wf2[st2*no*nv2 + o2*nv2], 1, wf2_col, 1);
        			finish = omp_get_wtime();
				times[8] += finish - start;

				/* Set auxiliary matrix tmp to zero
				 * Tmp holds the outer product of vectors wf2(o1*nv1, *) and wf1(o2*nv2, *)
				 * */
			        start = omp_get_wtime();
				cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, ld_tmp, tmp_cols, nv2, zero, wf1, ld_tmp, wf2, nv2, zero, tmp, ld_tmp);
        			finish = omp_get_wtime();
			        times[4] += finish - start;

				/* Compute wf2 * wf1 */
        			start = omp_get_wtime();
				cblas_dger(CblasColMajor, nv2, nv1, one, wf2_col, 1, wf1_col, 1, tmp, ld_tmp);
        			finish = omp_get_wtime();
			        times[5] += finish - start;

				/* Local variable to store partial ss sum in omp */
				ss_local = zero;

				/* Set openblas/mkl to single-threaded */
#ifdef OPENBLAS
				openblas_set_num_threads(1);
#elif MKL
				mkl_set_num_threads(1);
#endif
			        /* Update ss = msum * wf_prod (ddot) */
        			start = omp_get_wtime();

				#pragma omp parallel default(shared) num_threads(num_threads)
				#pragma omp for reduction(+:ss_local)
				for(i=0; i<nv1; i++){
					ss_local += cblas_ddot(nv2, &msum[i*ld_msum], 1, &tmp[i], ld_tmp);
				}
        			finish = omp_get_wtime();
			        times[7] += finish - start;
#ifdef OPENBLAS
				openblas_set_num_threads(num_threads);
#elif MKL
				mkl_set_num_threads(num_threads);
#endif
				/* Update ss with local update */
				ss[st2*ns1 + st1] += ss_local;
			}
		}
	    }
	}

	/* Free allocated spaces */
	free(msum);
	free(tmp);
	free(l2minors);
	free(wf1_col);
	free(wf2_col);

	/* Print timings */
        printf("[SSBLOCK] Execution times:\n");
        printf("\tscale:         %.7lf sec\n", times[1]);
        printf("\tl2minors:      %.7lf sec\n", times[9]),
        printf("\tl2minor * wf2: %.7lf sec\n", times[2]),
        printf("\twf1 * tmp:     %.7lf sec\n", times[3]);
        printf("\ttmp to zero    %.7lf sec\n", times[4]);
        printf("\twf1 * wf2:     %.7lf sec\n", times[5]);
        printf("\tmsum * tmp:    %.7lf sec\n", times[7]);
        printf("\tdcopy:         %.7lf sec\n", times[8]);
        printf("\tallocation:    %.7lf sec\n", times[0]);

}

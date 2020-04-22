/*--------------------------------------------------------------------------------------------------
  FUNCTION: ssblock_gpu

  DESCRIPTION: 
--------------------------------------------------------------------------------------------------
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include "cuda.h"
#include "cuda_runtime.h"
#include "cublas_v2.h"
#include <cblas.h>

void printMatrix(int, int, double*);
void sortMsum(double *, int, int, int, double *);
extern void getlvl2minors_blk_( int*, int*, int*, int*, double*, double*);


void ssblock_gpu(double *csc, int no, const int nv1, const int nv2, const int ns1, const int ns2, double *wf1, double *wf2, double *l1minor, double *ss)
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

        /* Timing routines */
        double start, finish;
        double times[10]; 
	//t_gemm1, t_gemm2, t_gemm3, t_dger, t_dgemv, t_sort, t_zero;
	//double alloc;

	/* Auxiliary variables */
	int ndevices;
	int devices[10];

	/* CUBLAS variables */
	cublasStatus_t status;
	cublasHandle_t handle;
	struct cudaDeviceProp prop;

/*	status = cublasCreate(&handle);
	if (status != CUBLAS_STATUS_SUCCESS){
		printf("CUBLAS initialization failed!");
	}

	cudaGetDeviceCount(&ndevices);
	for( o1 = 0; o1 < ndevices; o1++){
		cudaGetDeviceProperties(&prop, o1);
		printf("DEVICE INFO\n");
		printf("Device name: %s\n", prop.name);
		printf("CUDA version: %d.%d\n", prop.major, prop.minor);
		printf("Total memory (GB): %lf\n", (double)prop.totalGlobalMem/pow(1024,3));
		printf("\n");
	}
*/
        /* Set leading dimensions */
        ld_csc =  no + nv1;
        ld_l2minors = no;
        //ld_l2minors = no * no * no;
        //ld_tmp = ld_l2minors;
	// Two tmp are required.
	// The first is the auxiliary matrix for computing msum, dim (no x nv2)
	// The second is the auxiliary matrix to store product wf2 * wf1, dim (nv2 x nv1)
	ld_tmp = (no < nv2) ? nv2 : no;
	tmp_cols = (nv2 < nv1) ? nv1 : nv2;
//        ld_msum = no * no * nv1;
	ld_msum = nv1;

	/* Set timings to zero */
	for( i = 0; i < 10; i++ ){
		times[i] = zero;
	}

        /* Allocate msum matrix to keep computed determinants */
	start = omp_get_wtime();
        msum = (double*) malloc (ld_msum * nv2 * sizeof(double));
        tmp = (double*) malloc(ld_tmp * tmp_cols * sizeof(double));
	l2minors = (double*) malloc(ld_l2minors * no * sizeof(double));
        //msum_sorted = (double*) malloc(ld_msum * nv2 * sizeof(double));
	wf1_col = (double*) malloc(nv1 * sizeof(double));
	wf2_col = (double*) malloc(nv2 * sizeof(double));
	times[0] += omp_get_wtime() - start;

#ifdef OPENBLAS
        printf("Number of openblas threads for LAPACK calls: %d\n", openblas_get_num_threads());
#elif MKL
        printf("Number of intel mkl threads for LAPACK calls: %d\n", mkl_get_max_threads());
#endif        

        /* Scale lower csc matrix block with alternating array (1,-1,1,-1...) i.e. starting from the second 
         * column, scale every second column with -1.
         * Arrays wf1 and wf2 are used as auxiliary arrays and are not referenced (unchaged on exit).
	 */
        start = omp_get_wtime();
	cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, nv1, no/2, nv2, zero, wf1, nv1, wf2, nv2, neg_one, &csc[no+nv1 + no], 2*ld_csc);
        finish = omp_get_wtime();
        times[1] += finish - start;

	/* Set ss to zero-matrix */
/*	for( st1 = 0; st1 < ns1; st1++ ){
		for( st2 = 0; st2 < ns2; st2++ ){
			ss[st2*ns1 + st1] = zero;
		}
	}
*/
        cblas_dscal(ns1*ns2, zero, ss, 1);

	// POČETAK FOR PETLJE PO BLOKOVIMA //
	for( o1 = 0; o1 < no; o1++ ){
	    for( o2 = 0; o2 < no; o2++ ){

		/* Compute block (o1,o2) with l2minors */
		r = o1 + 1;
		c = o2 + 1;
		dimmat = no+nv1;
		start = omp_get_wtime();
		getlvl2minors_blk_(&r, &c, &no, &dimmat, csc, l2minors);
		//printf("l2minors:\n");
		//printMatrix(no, no, l2minors);
		finish = omp_get_wtime();
		times[9] += finish - start;
	        
		/* Compute first product l2minors * csc(1:no, no+1:no+nv2) */
        	start = omp_get_wtime(); 
		pos_tmp = o1*no*no + o2*no;
        	//cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, pow(no,3), nv2, no, one, l2minor, ld_l2minors, &csc[no*ld_csc], ld_csc, zero, tmp, ld_tmp);
        	cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, no, nv2, no, one, l2minors, ld_l2minors, &csc[no*ld_csc], ld_csc, zero, tmp, ld_tmp);
        	finish = omp_get_wtime();
        	times[2] += finish - start;

        	start = omp_get_wtime();
		//for( o1 = 0; o1 < no; o1++ ){
	    	//	for( o2 = 0; o2 < no; o2++ ){
                //pos_tmp = o1*no*no + o2*no;
                pos_msum = o1*no*nv1 + o2*nv1;
		//cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, nv1, nv2, no, pow(neg_one, o1+2), &csc[no], ld_csc, tmp, ld_tmp, zero, &msum[pos_msum], ld_msum);
		cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, nv1, nv2, no, pow(neg_one, o1+2), &csc[no], ld_csc, tmp, ld_tmp, zero, msum, ld_msum);
		//printf("Pass o1/o2: %d/%d\n", o1+1, o2+1);
		//printMatrix(ld_msum, nv2, msum);
//                call dgemm('N', 'N', nv1, nv2, no, (-1.0d0)**(o1+1), csc(no+1, 1), ld_csc, tmp(pos_tmp, 1), ld_tmp, 0.0d0, msum(pos_msum, 1), ld_msum)
		coef = pow(neg_one, o1+o2) * l1minor[o2*no + o1];
		for( i = 0; i < nv2; i++ ){
			//cblas_daxpy(nv1, coef, &csc[(no+i)*ld_csc+no], 1, &msum[i*ld_msum+pos_msum], 1);
			cblas_daxpy(nv1, coef, &csc[(no+i)*ld_csc+no], 1, &msum[i*ld_msum], 1);
		}
//                msum(pos_msum:pos_msum+nv1-1, 1:nv2) = msum(pos_msum:pos_msum+nv1-1, 1:nv2) + ((-1.0d0)**(o1+o2)*l1minor(o1,o2))*csc(no+1:no+nv1, no+1:no+nv2)
	    //}
	//}
        	finish = omp_get_wtime();
        	times[3] += finish - start;

//		printf("\n");
//		printf("MSUM\n");
//		for(i=0;i< ld_msum*nv2;i++)
//			printf("%.14e ", msum[i]);
//		printf("\n");

		/* Compute ss matrix  */
	        for( st1 = 0; st1 < ns1; st1++ ){

			/* Copy columns wf1(o1*nv1, st1) to wf1_col */
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

//				printf("TMP before\n");
//				for(i=0;i< ld_tmp*tmp_cols;i++)
//					printf("%.14e ", tmp[i]);
//				printf("\n");

				//printMatrix(nv1*nv2,1, tmp);

				/* Compute wf2 * wf1 */
        			start = omp_get_wtime();
				cblas_dger(CblasColMajor, nv2, nv1, one, wf2_col, 1, wf1_col, 1, tmp, ld_tmp);
        			finish = omp_get_wtime();
			        times[5] += finish - start;

//				printf("wf1/wf2\n");
				//printMatrix(1, nv1, wf1_col);
				//printMatrix(1, nv2, wf2_col);
				//printMatrix(nv1*nv2,1, tmp);
//				printf("TMP after\n");
//				for(i=0;i< ld_tmp*tmp_cols;i++)
//					printf("%.14e ", tmp[i]);
//				printf("\n");

			        /* Update ss = msum * wf_prod (ddot) */
        			start = omp_get_wtime();
				for(i=0; i<nv1; i++)
					ss[st2*ns1 + st1] += cblas_ddot(nv2, &msum[i*ld_msum], 1, &tmp[i], ld_tmp);
				//ss[st2*ns1 + st1] += cblas_ddot(nv1*nv2, msum, 1, tmp, 1);
        			finish = omp_get_wtime();
			        times[7] += finish - start;
			}
		}
	    }
	}

//	printf("Printing msum...\n");
//	printMatrix(ld_msum, nv2, msum);

//	free(tmp);

//      ld_tmp = no*no*nv1*nv2;
//	start = omp_get_wtime();
//      tmp = (double*) malloc(ld_tmp*ns1*ns2 * sizeof(double));
//	times[0] += omp_get_wtime() - start;

	/* Set auxiliary matrix tmp to zero */
//        start = omp_get_wtime();
        //cblas_dscal(ld_tmp*ns1*ns2, zero, tmp, 1);
//        cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, ld_tmp, ns1*ns2, nv2, zero, wf1, ld_tmp, wf2, nv2, zero, tmp, ld_tmp); 
//        finish = omp_get_wtime();
//        times[4] += finish - start;

        /* Generate products wf_prod = wf1(st1) * wf2(st2) (dger) */
/*        start = omp_get_wtime();
        for(st2 = 0; st2 < ns2; st2++){
		for(st1 = 0; st1 < ns1; st1++){
                	cblas_dger(CblasColMajor, no*nv2, no*nv1, one, &wf2[st2*no*nv2], 1, &wf1[st1*no*nv1], 1, &tmp[(st2*ns1+st1)*ld_tmp], no*nv2);
            	}
	}
        finish = omp_get_wtime();
        times[5] += finish - start;
*/
//	printf("Printing tmp...\n");
//	printMatrix(ld_tmp, ns1*ns2, tmp);

        /* Reshape/sort msum */
/*        start = omp_get_wtime();
        sortMsum(msum, no, nv1, nv2, msum_sorted);
        finish = omp_get_wtime();
        times[6] += finish - start;
*/
        /* Compute ss = msum * wf_prod (dgemv) */
/*        start = omp_get_wtime();
        cblas_dgemv(CblasColMajor, CblasTrans, ld_tmp, ns1*ns2, one, tmp, ld_tmp, msum_sorted, 1, zero, ss, 1);
        finish = omp_get_wtime();
        times[7] += finish - start;
*/
	/* Free allocated spaces */
	free(msum);
	free(tmp);
	free(l2minors);
//	free(msum_sorted);
	free(wf1_col);
	free(wf2_col);

	/* Print timings */
        printf("[SSBLOCK] Execution times:\n");
        printf("\tscale 1 -1:   %.7lf sec\n", times[1]);
        printf("\tl2minors:     %.7lf sec\n", times[9]),
        printf("\tl2minor * wf2:%.7lf sec\n", times[2]),
        printf("\twf1 * tmp:    %.7lf sec\n", times[3]);
        printf("\ttmp to zero   %.7lf sec\n", times[4]);
        printf("\twf1 * wf2:    %.7lf sec\n", times[5]);
        printf("\tmsum sort:    %.7lf sec\n", times[6]);
        printf("\tmsum * tmp:   %.7lf sec\n", times[7]);
        printf("\tdcopy:        %.7lf sec\n", times[8]);
        printf("\tallocation:   %.7lf sec\n", times[0]);

/*	istatus = cublasDestroy(handle);
	if( status != CUBLAS_STATUS_SUCCESS ){
		printf("Cannot destroy CUBLAS - not initialized!");
	}
*/
}
 
void printMatrix(int m, int n, double *mat)
{
        int i, j;

	printf("Printing matrix...\n");
	for (i = 0 ; i < m; i++){
		for (j = 0; j < n; j++)
			printf("%.15e ", mat[j*m+i]);
		printf("\n");
	}
}

void sortMsum(double *msum, int no, int nv1, int nv2, double *msum_sorted)
{
        int blk;
        int blk_start_pos;
        int row_in_blk;
        int sorted_array_start_pos;
        int offset1, offset2;
        int ld_msum;

        ld_msum = no*no*nv1;

	for( blk = 0; blk < no; blk++ ){
            blk_start_pos = blk*no*nv1;

	    for( row_in_blk = 0; row_in_blk < no*nv1; row_in_blk++ ){
                offset1 = row_in_blk % nv1;
                offset2 = row_in_blk/nv1;
                sorted_array_start_pos = blk*no*nv1*nv2 + offset1*nv2*no + offset2*nv2;
                cblas_dcopy(nv2, &msum[blk_start_pos+row_in_blk], ld_msum, &msum_sorted[sorted_array_start_pos], 1);
	    }
	}
}


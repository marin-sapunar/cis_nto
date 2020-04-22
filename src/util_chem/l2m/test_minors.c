/*--------------------------------------------------------------------------------------------------
  PROGRAM: test_minors

  DESCRIPTION: 

  The program tests the computation of the determinants of all minors of order n-2 of a given 
  input matrix. The matrix is randomly generated at the runtime.

  Input parameters:
 
	NMIN:STEP:NMAX - input matrix size, range from NMIN to NMAX is steps STEP. Default 10:10:100
	NMIN:NMAX      - input matrix size, range from NMIN to NMAX, step is computed as (NMAX-NMIX)/10

--------------------------------------------------------------------------------------------------
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <omp.h>

#ifdef OPENBLAS
#include <cblas.h>
#elif MKL
#include <mkl.h>
#endif

#define max(a, b) ( a < b ? b : a )

extern void getlvl2minors_lu_(int*, const int*, double*, int*, double*, int*);
extern void getlvl2minors_lu_csa_(int*, const int*, double*, int*, double*, int*);
extern void getlvl2minors(const int*, const double*, double*);

#ifdef OPENBLAS
	extern void dlarnv_(int *, int*, int*, double*);
#endif

/*void print_mat(int m, int n, double *A, int lda){
	
	int i, j;

	for(i = 0; i < m; i++){
		for(j=0; j<n;j++){
			printf("%e ", A[i*lda + j]);
		}
		printf("\n");
	}
	printf("\n\n");
}
*/
int readMat_file(char *filename, int m, int n, double *A, int lda)
{
	FILE *file;
	int m_max, n_max;
	int status;
	int i, j;
	double skip_num;

	file = fopen(filename, "r");

	/* Read number of rows/columns */
	status = fscanf(file, "%d", &m_max); //skip the first number
	status = fscanf(file, "%d", &m_max);
	status = fscanf(file, "%d", &n_max);

	if(n_max < n){
		printf("Given max number of columns is larger than the matrix size!\nCannot read from file. Breaking...\n");
		return -1;
	}
	if(m_max < n){
		printf("Given max number of rows is larger than the matrix size!\nCannot read from file. Breaking...\n");
		return -1;
	}

	printf("Max number of rows = %d and columns = %d\n", m_max, n_max);

	for (i = 0; i < m; i++) {
		for (j = 0; j < n; j++) {
			status = fscanf(file, "%lf", &A[i*lda+j]);
		}
		/* Read rest of the line */
		for (j = n; j < n_max; j++)
			status = fscanf(file, "%lf", &skip_num);
	}

	fclose(file);

	return 0;
}

int main(int argc, char **argv)
{
	int nmin, nmax, step;
	int num_step;
	int num_det;
	int lda, ldadet;
	int n;
	int nv2 = 0;
	int iseed[] = {1,13,23,7};
	double gflops_v1, gflops_v2;
	double error, max_error;
	double *A = NULL;
	double *Adet = NULL;
	double *Adet_ref = NULL;
	char *filename = NULL;

	//char input_str[100];
	const char delim[2] = ":";
	char *input_str = NULL;
	int inputs[3];

	/* Timing routines */
	double start, finish;
	double time_v1, time_v2;
	
    /* Set default values */
	nmin = 10;
	nmax = 100;
	step = 10;

	/* Read input parameters */
	if(argc < 2){
		printf("%s: Missing input parameter\n", argv[0]);
		printf("Usage:\n");
		printf("\t%s NMIN:STEP:NMAX <path-to-input-matrix>\n", argv[0]);
		printf("\t%s NMIN:NMAX <path-to-input-matrix>\n", argv[0]);
		exit(-1);
	}

	/* Parse first input parameter which is in for NMIN:STEP:NMAX or NMIN:NMAX */
	input_str = strtok(argv[1], delim);
	for( int i = 0; i < 3; i++ ) {
		inputs[i] = atoi(input_str);
		input_str = strtok(NULL, delim);
		if(input_str == NULL){
			inputs[i+1] = -1;
			break;
		}
	}

	nmin = inputs[0];
	if(inputs[2] == -1){
		nmax = inputs[1];
	}
	else{
		step = inputs[1];
		nmax = inputs[2];
	}
	n = nmax;

	/* Read input file name */
	if( argc == 3 )
		filename = argv[2];

	printf("nmin = %d, nmax = %d, step = %d\n", nmin, nmax, step);
	if(filename != NULL){
		printf("Read matrix from file %s\n", filename);
	}else
		printf("Generating random matrix (normal distribution [-1, 1]\n");


	printf("Number of omp threads =      %d\n", omp_get_max_threads());
#ifdef OPENBLAS
	printf("Number of openblas threads = %d\n", openblas_get_num_threads());
#endif
#ifdef MKL
	printf("Number of mkl threads = %d\n", mkl_get_max_threads());
#endif

    /* Allocate reference matrix */
	A = (double*) malloc( n * n * sizeof(double) );
	Adet = (double*) malloc ( n * n * n * n * sizeof(double) );
	Adet_ref = (double*) malloc ( n * n * n * n * sizeof(double) );

	int num_elem = n*n;
	int idist = 2;

	/* Read matrix from cscmat file */
	if(filename){
		readMat_file(filename, nmax, nmax, A, n);
	}else{
	
		/* Generate random matrices */
		dlarnv_(&idist, iseed, &num_elem, A);
	}

	printf("%%                          |      Version 1     |     Version 2\n");
	printf("%%   STEP     N     no.Det  | TIME (s)   Gflop/s | TIME (s)    Gflops/s  Error\n");
	printf("-----------------------------------------------------------------------------\n");

	num_step = 0;

	for( n = nmin; n <= nmax; n += step ) {

		lda = n;
		ldadet = n*n*n;

		/* Version 1 - partially exploting structure, LU inside c2 loop */
		start = omp_get_wtime();

		getlvl2minors_lu_(&n, &nv2, A, &lda, Adet_ref, &ldadet);
		//getlvl2minors(&n, A, Adet_ref);
		
		finish = omp_get_wtime();
		time_v1 = finish - start;
	
		/* Gflops = 1/6 n^7 */
		gflops_v1 = (1.0/6.0 * (double)ldadet * ldadet * n)/(pow(10,9)*time_v1);
		
		/* Version 2 - exploting structure, LU outside of the c1/c2 loops */
		start = omp_get_wtime();

		getlvl2minors_lu_csa_(&n, &nv2, A, &lda, Adet, &ldadet);
		
		finish = omp_get_wtime();
		time_v2 = finish - start;
	
		/* Gflops = 1/6 n^7 */
		gflops_v2 = (1.0/6.0 * (double)ldadet * ldadet * n)/(pow(10,9)*time_v2);

		/* Number of determinants computed */
		num_det = 1.0/4.0 * pow(n,4);

		/* Check errors */
		/* Compute max error between v1 and v2 */
		max_error = 0.0;
/*		for(int i = 0; i < n; i++)
			for(int j = 0; j < n*n*n; j++){
				error = abs(Adet_ref[i*ldadet + j] - Adet[i*ldadet + j]);
				if(error > max_error)
					max_error = error;
			}*/

		/* Print outputs */
		printf("%7lld %7lld %7lld   %8.3f   %8.3f   %8.3f   %8.3f   %8.3e\n",
				(long long) num_step, (long long) n, (long long) num_det,
				time_v1,    gflops_v1, time_v2, gflops_v2, max_error );
		fflush(stdout);
		num_step++;
	}

	/* Free allocated spaces */
	free(A);
	free(Adet);
	free(Adet_ref);

	return 0;
} 

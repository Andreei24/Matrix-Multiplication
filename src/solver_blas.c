/*
 * Tema 2 ASC
 * 2022 Spring
 */
#include "utils.h"
#include "cblas.h"
#include <string.h>


/* 
 * Add your BLAS implementation here
 */
double* my_solver(int N, double *A, double *B) {
	printf("BLAS SOLVER\n");
	
	//B * A matrix to be overwritten by dtrmm
	double* BA = calloc(N * N,sizeof(double));
	memcpy(BA,B,N * N * sizeof(double));

	//multiply B and A, from the right of B, where A is an
	//Upper triangular, not transposed, non unit matrix
	cblas_dtrmm(CblasRowMajor, CblasRight, CblasUpper, CblasNoTrans,
				CblasNonUnit, N, N, 1.0, A, N, BA, N);
	
	// B * A * (A**T) matrix
	double* BAAt  = calloc(N * N, sizeof(double));
	memcpy(BAAt,BA,N*N*sizeof(double));
	
	//muliply B * A and (A**T) where the params are the same
	//as the last time, except now A is transposed
	cblas_dtrmm(CblasRowMajor, CblasRight, CblasUpper, CblasTrans,
				CblasNonUnit, N, N, 1.0, A, N, BAAt, N);

	//sum matirx, has the value of B * A * (A**T) operation
	double* sum = calloc(N * N,sizeof(double));
	memcpy(sum,BAAt,N*N*sizeof(double));

	//multiply (B**T) and B, where first matrix is transposed
	//and the second is not, add the result to sum matrix 
	cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, N, N, N,
				1.0, B, N, B, N, 1.0, sum, N);

	free(BA);
	free(BAAt);

	return sum;
}

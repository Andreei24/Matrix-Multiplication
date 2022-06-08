/*
 * Tema 2 ASC
 * 2022 Spring
 */
#include "utils.h"

/*
 * Add your unoptimized implementation here
 */
double* my_solver(int N, double *A, double* B) {
	printf("NEOPT SOLVER\n");

	//B * A matrix
	double* BA = calloc(N * N, sizeof(double));

	//only the first j elements are of interest 
	//because A is Upper Triangular
	for(int i = 0; i < N; i++)
		for(int j = 0; j < N; j++)
			for(int k = 0; k <= j; k++){
				BA[i * N + j] += B[i * N + k] * A[k * N + j];
			}
	
	//B * A * (A**T) matrix
	double* BAAt = calloc(N * N, sizeof(double));

	//only the elements from j through N are of interest
	//because now A is transposed
	for(int i = 0; i < N; i++)
		for(int j = 0; j < N; j++)
			for(int k = j; k < N; k++){
				BAAt[i * N + j] += BA[i * N + k] * A[j * N + k];
			}
	
	//(B**T) * B
	double* BtB = calloc(N * N, sizeof(double));

	//go through all elements
	for(int i = 0; i < N; i++)
		for(int j = 0; j < N; j++)
			for(int k = 0; k < N; k++){
				BtB[i * N + j] += B[k * N + i] * B[k * N + j];
			}

	// add B * A * (A**T) and (B**T)*B
	double* sum = calloc(N * N, sizeof(double));

	for(int i = 0; i < N; i++)
		for(int j = 0; j < N; j++){
			sum[i * N + j] = BAAt[i * N + j] + BtB[i * N + j];
		}

	free(BA);
	free(BAAt);
	free(BtB);

	return sum;
}

/*
 * Tema 2 ASC
 * 2022 Spring
 */
#include "utils.h"
#include <string.h>

/*
 * Add your optimized implementation here
 */

//function to compute the transpose of a matrix
double* transpose(double* M, int N){

	double* Mt = calloc(N * N, sizeof(double));
	if(!Mt)
		return NULL;

	for(int i = 0; i < N; i++){
		
		register double* aux_M = M + i * N;
		register double* aux_Mt =  Mt + i;

		for(int j = 0; j < N; j++){
			*aux_Mt = *aux_M;
			aux_Mt += N;
			aux_M++;
		}
	}

	return Mt;
}

double* my_solver(int N, double *A, double* B) {
	printf("OPT SOLVER\n");

	int bi,bj,bk;
	int i,j,k;
	int blockSize = 40;

	// B * A matrix
	double* BA = calloc(N * N, sizeof(double));
	if(!BA)
		return NULL;

	//Optimizations:
	//Use registers, block calculation 40 x 40
	//Vector access with pointers
	//Take into account A is Upper Triangluar and (A**T) Lower

	for(bi=0; bi<N; bi+=blockSize)
        for(bj=0; bj<N; bj+=blockSize)
            for(bk=0; bk<=bj; bk+=blockSize) {

                double *aux_pb = B +  bi * N + bk; // B[bi][bk]
                double *aux_pa = A + bk * N + bj; // A[bk][bj]
                double *aux_pc = BA + bi * N + bj; // C[bi][bj]

                for(i=0; i<blockSize; i++) {
                    for(j=0; j<blockSize; j++) {

                        double *pb = aux_pb;
                        double *pa = aux_pa + j;
                        register double suma = 0.0;

                        for(k=0; k<blockSize; k++) {
                            suma += *pb * *pa;
                            pb++;
                            pa += N;
                        }
                        *(aux_pc+j) += suma;
                    }
                    aux_pb += N;
                    aux_pc += N;
                }
            }
	
	//B * A * (A**T) matrix
	double* BAAt = calloc(N * N, sizeof(double));
	if(!BAAt)
		return NULL;

	double* At = transpose(A,N);

	for(bi = 0; bi < N; bi += blockSize)
		for(bj = 0; bj < N; bj += blockSize)
			for(bk = bj; bk < N; bk += blockSize){

				double* aux_pba = BA + bi * N + bk; // BA[bi][bk]
				double* aux_pat = At + bk * N + bj; // At[bk][bj]
				double* aux_pc = BAAt + bi * N + bj; // BAAt[bi][bj]

				for(i = 0; i < blockSize; i++){
					for(j = 0; j < blockSize; j++){

						double* pba = aux_pba;
						double* pat = aux_pat + j;
						register double suma = 0.0;

						for(k = 0; k < blockSize; k++){
							suma += *pba * *pat;
							pba++;
							pat += N;
						}
						*(aux_pc + j) += suma;
					}
					aux_pba += N;
					aux_pc += N;
				}
			}

	// B * (B**T) matrix
	double* BtB = calloc(N * N, sizeof(double));
	if(!BtB)
		return NULL;

	double* Bt = transpose(B,N);

	for(bi = 0; bi < N; bi += blockSize)
		for(bj = 0; bj < N; bj += blockSize)
			for(bk = 0; bk < N; bk += blockSize){

				double* aux_pbt = Bt + bi * N + bk; // Bt[bi][bk]
				double* aux_pb = B + bk * N + bj; // B[bk][bj]
				double* aux_pc = BtB + bi * N + bj; // BtB[bi][bj]

				for(i = 0; i < blockSize; i++){
					for(j = 0; j < blockSize; j++){

						double* pbt = aux_pbt;
						double* pb = aux_pb + j;
						register double suma = 0.0;

						for(k = 0; k < blockSize; k++){
							suma += *pbt * *pb;
							pbt++;
							pb += N;
						}
						*(aux_pc + j) += suma;
					}
					aux_pbt += N;
					aux_pc += N;
				}
			}

	//final sum matrix
	double* sum = calloc(N * N, sizeof(double));
	if(!sum)
		return NULL;

	register double* aux_BAAt = BAAt;
	register double* aux_BtB = BtB;
	register double* aux_sum = sum;

	for(i = 0; i < N*N; i++){
		*aux_sum  = *aux_BAAt + *aux_BtB;
		aux_BAAt++;
		aux_BtB++;
		aux_sum++;
	}

	free(BA);
	free(BAAt);
	free(At);
	free(Bt);
	free(BtB);

	
	return sum;	
}

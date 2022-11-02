#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <cblas.h>
#include <cuda.h>
#include <cublas.h>

/* Perform A@x = b */
/* Device Code */
__global__ void vector_matrix_multiplication(
        const float** A,
        const float* x,
        float* b,
        const float n) {
    int my_index = blockDim.x * blockIdx.x +threadIdx.x;
    if (my_index < n){
        b[my_index] = 0;
        for (int i=0; i<n; i++){
            b[my_index] += A[my_index][i]*x[i];
        }
    }
}

void allocate_vectors(
        float *** hA_p,
        float ** hx_p,
        float ** hb_p,
        float *** dA_p,
        float ** dx_p,
        float ** db_p,
        int n) {
    cudaMalloc(dA_p, n*sizeof (float));
    cudaMalloc(dx_p, n*sizeof (float));
    cudaCalloc(db_p, n*sizeof (float));

    *hA_p = (float **) malloc(n*sizeof (float));
    *hx_p = (float *) malloc(n*sizeof (float));
    *hb_p = (float *) calloc(n*sizeof (float));
}

/* Host Code */
int main(){
    int n = 100, thread_per_block=10, block_counter=10;
    float *hx, *dx, *hb, *db, **hA, **dA;

    allocate_vectors(&hA, &hx, &hb, &dA, &dx, &db, n);
    for (int i=0; i < n; i++){
        hx[i] = 1;
        for (int j=0; j < n; j++){
            hA[i][j] = 1;
        }
    }

    cudaMemcpy(dx, hx, n*sizeof (float), cudaMemcpyHostToDevice);
    cudaMemcpy(dA, hA, n*sizeof (float), cudaMemcpyHostToDevice);

    vector_matrix_multiplication <<<block_counter, thread_per_block>>>(dA, dx, db, n);
    cudaDeviceSyncronize();

    cudaMemcpy(hb, db, n*sizeof (float), cudaMemcpyDeviceToHost);

    printf("[");
    for (int i=0; i<n; i++){
        printf("%f\n", hb[i]);
    }
    printf("]");

    return 0;
}
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <cblas.h>
#include <cuda.h>

/* Perform x@y = r */
/* Device Code */
__global__ void vector_vector_multiplication(
        const float* x,
        const float* y,
        float* a) {
    int my_index = blockDim.x * blockIdx.x +threadIdx.x;
    if (my_index * sizeof (float *) < sizeof x){
        a[my_index] = x[my_index]*y[my_index];
    }
}

__global__ void vector_matrix_multiplication(
        const float** A,
        const float* x,
        float* b,
        const float n) {
    int my_index = blockDim.x * blockIdx.x +threadIdx.x;
    if (my_index * sizeof (float *) < sizeof x){
        b[my_index] = 0;
        for (int i=0; i<n; i++){
            b[my_index] += A[my_index][i]*x[i];
        }
    }
}

void allocate_vectors(
        float *** hA_p,
        float *** dA_p,
        float** hx_p,
        float** hy_p,
        float** ha_p,
        float** dx_p,
        float** dy_p,
        float** da_p,
        int n) {
    cudaMalloc(dx_p, n*sizeof (float));
    cudaMalloc(dy_p, n*sizeof (float));
    cudaMalloc(da_p, n*sizeof (float));
    cudaMalloc(dA_p, n*sizeof (float));

    *hA_p = (float **) malloc(n*sizeof (float));
    *hx_p = (float *) malloc(n*sizeof (float));
    *hy_p = (float *) malloc(n*sizeof (float));
    *ha_p = (float *) malloc(n*sizeof (float));
}

/* Host Code */
int main(){
    int n = 100, thread_per_block=10, block_counter=10, flag=0;
    float *hx, *dx, *hy, *dy, *ha, *da, r, **hA, **dA;

    allocate_vectors(&hA, &dA, &hx, &hy, &ha, &dx, &dy, &da, n);
    for (int i=0; i < n; i++){
        hx[i] = 1;
        hy[i] = 1;
        for (int j=0; j < n; j++){
            hA[i][j] = 1;
        }
    }

    cudaMemcpy(dx, hx, n*sizeof (float), cudaMemcpyHostToDevice);
    cudaMemcpy(dy, hy, n*sizeof (float), cudaMemcpyHostToDevice);
    cudaMemcpy(dA, hA, n*sizeof (float), cudaMemcpyHostToDevice);

    if (flag==0){
        vector_vector_multiplication <<<block_counter, thread_per_block>>>(dx, dy, da);
    }
    else{
        vector_matrix_multiplication <<<block_counter, thread_per_block>>>(dA, dx, da, n);

    }
    cudaDeviceSyncronize();

    cudaMemcpy(ha, da, n*sizeof (float), cudaMemcpyDeviceToHost);

    for (int i=0; i < n; i++) {
        r += ha[i];
    }
    if ()
    printf("result = %f\n", r);
    return 0;
}
#include <cuda.h>
#include <stdio.h>
#include <stdlib.h>

__global__ void mandelKernel(float lowerX, float lowerY, float stepX, float stepY, int *result, int maxIterations) {
    // To avoid error caused by the floating number, use the following pseudo code
    //
    // float x = lowerX + thisX * stepX;
    // float y = lowerY + thisY * stepY;
}

// Host front-end function that allocates the memory and launches the GPU kernel
void hostFE (float upperX, float upperY, float lowerX, float lowerY, int* img, int resX, int resY, int maxIterations)
{
    float stepX = (upperX - lowerX) / resX;
    float stepY = (upperY - lowerY) / resY;

    // Allocate result array on host memory
    int *result = (int *)malloc(resX * resY * sizeof(int));

    // Allocate result array on device memory
    int *cudaResult = NULL;
    cudaMalloc(&cudaResult, resX * resY * sizeof(int));

    int gridSize;
    int blockSize;

    mandelKernel<<<1, 1>>>(lowerX, lowerY, stepX, stepY, cudaResult, maxIterations);

    // Copy result array from device to host memory
    cudaMemcpy(result, cudaResult, resX * resY * sizeof(int), cudaMemcpyKind::cudaMemcpyDeviceToHost);
    cudaFree(cudaResult);

    memcpy(img, result, resX * resY * sizeof(int));
    free(result);
}

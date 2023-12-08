#include <cuda.h>
#include <stdio.h>
#include <stdlib.h>

#define BLOCK_WIDTH 16
#define BLOCK_HEIGHT 16

// limited version of checkCudaErrors from helper_cuda.h in CUDA examples
#define checkCudaErrors(val) check_cuda((val), #val, __FILE__, __LINE__)

void check_cuda(cudaError_t result, char const *const func,
                const char *const file, int const line)
{
    if (result)
    {
        printf("CUDA error %s at %s:%d '%s' \n", cudaGetErrorString(result), file, line, func);

        // Make sure we call CUDA Device Reset before exiting
        cudaDeviceReset();
        exit(99);
    }
}

__device__ int mandel(float c_re, float c_im, int count)
{
    float z_re = c_re, z_im = c_im;
    int i;
    for (i = 0; i < count; ++i)
    {
        if (z_re * z_re + z_im * z_im > 4.f) break;

        float new_re = z_re * z_re - z_im * z_im;
        float new_im = 2.f * z_re * z_im;
        z_re = c_re + new_re;
        z_im = c_im + new_im;
    }

    return i;
}

__global__ void mandelKernel(float lowerX, float lowerY, float stepX, float stepY, int width, int *result, int maxIterations) {
    // To avoid error caused by the floating number, use the following pseudo code
    //
    int thisX = blockIdx.x * blockDim.x + threadIdx.x;
    int thisY = blockIdx.y * blockDim.y + threadIdx.y;
    float x = lowerX + thisX * stepX;
    float y = lowerY + thisY * stepY;

    int index = (thisY * width + thisX);
    result[index] = mandel(x, y, maxIterations);
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
    checkCudaErrors(cudaMalloc(&cudaResult, resX * resY * sizeof(int)));

    // 1600 x 1200 = 1920000
    dim3 blockSize(BLOCK_WIDTH, BLOCK_HEIGHT);
    dim3 gridSize(resX / blockSize.x, resY / blockSize.y);

    mandelKernel<<<gridSize, blockSize>>>(lowerX, lowerY, stepX, stepY, resX, cudaResult, maxIterations);

    // Copy result array from device to host memory
    checkCudaErrors(cudaMemcpy(result, cudaResult, resX * resY * sizeof(int), cudaMemcpyKind::cudaMemcpyDeviceToHost));
    checkCudaErrors(cudaFree(cudaResult));

    memcpy(img, result, resX * resY * sizeof(int));
    free(result);
}

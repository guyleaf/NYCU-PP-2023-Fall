#include <cuda.h>
#include <stdio.h>
#include <stdlib.h>

#define BLOCK_WIDTH 16
#define BLOCK_HEIGHT 16

// comment out if you need faster implementation or define macro from outside
// #define USE_FASTER

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
        if (z_re * z_re + z_im * z_im > 4.f)
        {
            break;
        }

        float new_re = z_re * z_re - z_im * z_im;
        float new_im = 2.f * z_re * z_im;
        z_re = c_re + new_re;
        z_im = c_im + new_im;
    }

    return i;
}

__global__ void mandelKernel(float lowerX, float lowerY, float stepX, float stepY, size_t width, int *result, int maxIterations) {
    // To avoid error caused by the floating number, use the following pseudo code
    //
    int thisX = blockIdx.x * blockDim.x + threadIdx.x;
    int thisY = blockIdx.y * blockDim.y + threadIdx.y;
    float x = lowerX + thisX * stepX;
    float y = lowerY + thisY * stepY;

    int result_ = mandel(x, y, maxIterations);

#ifdef USE_FASTER
    size_t index = thisY * width + thisX;
    result[index] = result_;
#else
    ((int *)((char *)result + thisY * width))[thisX] = result_;
#endif
}

// Host front-end function that allocates the memory and launches the GPU kernel
void hostFE (float upperX, float upperY, float lowerX, float lowerY, int* img, int resX, int resY, int maxIterations)
{
    int *result = nullptr, *cudaResult = nullptr;
    float stepX = (upperX - lowerX) / resX;
    float stepY = (upperY - lowerY) / resY;

    // Plan 1: cudaHostAlloc + mapped (faster)
    // Plan 2: cudaHostAlloc + cudaMallocPitch + cudaMemcpy2D (slower) (homework required)

#ifdef USE_FASTER
    // Allocate result array on host memory
    checkCudaErrors(cudaHostAlloc(&result, resX * resY * sizeof(int), cudaHostAllocMapped));

    // Get the pointer to mapped memory on device
    checkCudaErrors(cudaHostGetDevicePointer(&cudaResult, result, 0));
#else
    // Allocate result array on host memory
    checkCudaErrors(cudaHostAlloc(&result, resX * resY * sizeof(int), cudaHostAllocDefault));

    size_t pitch;
    checkCudaErrors(cudaMallocPitch(&cudaResult, &pitch, resX * sizeof(int), resY));
#endif

    // 1600 x 1200 = 1920000
    int block_width = BLOCK_WIDTH;
    int block_height = BLOCK_HEIGHT;
    dim3 blockSize(block_width, block_height);
    dim3 gridSize(resX / block_width, resY / block_height);

#ifdef USE_FASTER
    mandelKernel<<<gridSize, blockSize>>>(lowerX, lowerY, stepX, stepY, resX, cudaResult, maxIterations);

    // Use cudaMemcpy with cudaMemcpyHostToHost flag without calling another synchronization
    checkCudaErrors(cudaMemcpy(img, result, resX * resY * sizeof(int), cudaMemcpyHostToHost));
#else
    mandelKernel<<<gridSize, blockSize>>>(lowerX, lowerY, stepX, stepY, pitch, cudaResult, maxIterations);

    // Copy result array from device to host memory
    checkCudaErrors(cudaMemcpy2D(result, resX * sizeof(int), cudaResult, pitch, resX * sizeof(int), resY, cudaMemcpyDeviceToHost));
    checkCudaErrors(cudaFree(cudaResult));

    memcpy(img, result, resX * resY * sizeof(int));
#endif

    checkCudaErrors(cudaFreeHost(result));
}

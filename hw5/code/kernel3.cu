#include <cuda.h>
#include <stdio.h>
#include <stdlib.h>

#define BLOCK_WIDTH 16
#define BLOCK_HEIGHT 16

#define GROUP_WIDTH 8
#define GROUP_HEIGHT 8

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

__global__ void mandelKernel(float lowerX, float lowerY, float stepX, float stepY, int width, int height, size_t pitch, int *result, int maxIterations) {
    // To avoid error caused by the floating number, use the following pseudo code
    //
    int thisX = (blockIdx.x * blockDim.x + threadIdx.x) * GROUP_WIDTH;
    int thisY = (blockIdx.y * blockDim.y + threadIdx.y) * GROUP_HEIGHT;

    width = min(thisX + GROUP_WIDTH, width);
    height = min(thisY + GROUP_HEIGHT, height);

    float x, y;
    for (size_t localY = thisY; localY < height; localY++)
    {
        y = lowerY + localY * stepY;
        for (size_t localX = thisX; localX < width; localX++)
        {
            x = lowerX + localX * stepX;
            ((int *)((char *)result + localY * pitch))[localX] = mandel(x, y, maxIterations);
        }
    }
}

// Host front-end function that allocates the memory and launches the GPU kernel
void hostFE (float upperX, float upperY, float lowerX, float lowerY, int* img, int resX, int resY, int maxIterations)
{
    float stepX = (upperX - lowerX) / resX;
    float stepY = (upperY - lowerY) / resY;

    int *result = nullptr;
    int *cudaResult = nullptr;
    size_t pitch;

    // Allocate result array on host memory
    checkCudaErrors(cudaHostAlloc(&result, resX * resY * sizeof(int), cudaHostAllocDefault));

    // Allocate padded array on device memory
    checkCudaErrors(cudaMallocPitch(&cudaResult, &pitch, resX * sizeof(int), resY));

    // 1600 x 1200 = 1920000
    int block_width = BLOCK_WIDTH;
    int block_height = BLOCK_HEIGHT;
    dim3 blockSize(block_width, block_height);

    block_width *= GROUP_WIDTH;
    block_height *= GROUP_HEIGHT;
    dim3 gridSize((int)std::ceil((float)resX / block_width), (int)std::ceil((float)resY / block_height));

    mandelKernel<<<gridSize, blockSize>>>(lowerX, lowerY, stepX, stepY, resX, resY, pitch, cudaResult, maxIterations);

    // Copy result array from device to host memory
    checkCudaErrors(cudaMemcpy2D(result, resX * sizeof(int), cudaResult, pitch, resX * sizeof(int), resY, cudaMemcpyDeviceToHost));
    checkCudaErrors(cudaFree(cudaResult));

    memcpy(img, result, resX * resY * sizeof(int));
    checkCudaErrors(cudaFreeHost(result));
}

#include "hostFE.h"

#include <stdio.h>
#include <stdlib.h>

#include "helper.h"

#define KERNEL_NAME "convolution"
#define WORK_GROUP_SIZE 32

void hostFE(int filterWidth, float *filter, int imageHeight, int imageWidth,
            float *inputImage, float *outputImage, cl_device_id *device,
            cl_context *context, cl_program *program)
{
    cl_int status;
    size_t imageSize = imageWidth * imageHeight * sizeof(int);
    size_t filterSize = filterWidth * filterWidth * sizeof(float);

    // load image and filter data to the input buffers on device memory
    cl_mem imageInput =
        clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
                       imageSize, (void *)inputImage, &status);
    CHECK(status, "clCreateBuffer(imageInput)");
    cl_mem filterInput =
        clCreateBuffer(context, CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR,
                       filterSize, (void *)filter, &status);
    CHECK(status, "clCreateBuffer(filterInput)");

    // create output buffer on device memory
    cl_mem imageOutput =
        clCreateBuffer(context, CL_MEM_WRITE_ONLY, imageSize, NULL, &status);
    CHECK(status, "clCreateBuffer(imageOutput)");

    // setup conv kernel
    cl_kernel convKernel = clCreateKernel(*program, KERNEL_NAME, &status);
    clSetKernelArg(convKernel, 0, sizeof(imageInput), (void *)&imageInput);
    clSetKernelArg(convKernel, 1, sizeof(filterInput), (void *)&filterInput);
    clSetKernelArg(convKernel, 2, sizeof(imageOutput), (void *)&imageOutput);
    clSetKernelArg(convKernel, 3, sizeof(imageWidth), (void *)&imageWidth);
    clSetKernelArg(convKernel, 4, sizeof(imageHeight), (void *)&imageHeight);
    clSetKernelArg(convKernel, 5, sizeof(filterWidth), (void *)&filterWidth);

    // create command queue
    cl_command_queue command_queue =
        clCreateCommandQueue(context, *device, 0, &status);
    CHECK(status, "clCreateCommandQueue");

    // execute kernel function
    size_t globalWorkSize[] = {imageWidth, imageHeight};
    size_t localWorkSize[] = {WORK_GROUP_SIZE, WORK_GROUP_SIZE};
    status =
        clEnqueueNDRangeKernel(command_queue, convKernel, 2, NULL,
                               globalWorkSize, localWorkSize, 0, NULL, NULL);
    CHECK(status, "clEnqueueNDRangeKernel");

    // wait for kernel
    // copy image from device to host memory
    status = clEnqueueReadBuffer(command_queue, imageOutput, CL_BLOCKING, 0,
                                 imageSize, (void *)outputImage, 0, NULL, NULL);
    CHECK(status, "clEnqueueReadBuffer");

    status = clFlush(command_queue);
    status = clFinish(command_queue);

    // clean up all objects created in this function
    clReleaseCommandQueue(command_queue);
    clReleaseKernel(convKernel);
    clReleaseMemObject(imageInput);
    clReleaseMemObject(filterInput);
    clReleaseMemObject(imageOutput);
}

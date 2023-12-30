__kernel void convolution(__global const float *imageInput,
                          __global const float *filterInput,
                          __global float *imageOutput, int imageWidth,
                          int imageHeight, int filterWidth)
{
    size_t x = get_global_id(0);
    size_t y = get_global_id(1);

    float sum = 0;
    imageOutput[y * imageWidth + x] = sum;
}

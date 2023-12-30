__kernel void convolution(__global const float *imageInput,
                          __global const float *filterInput,
                          __global float *imageOutput, int imageWidth,
                          int imageHeight, int filterWidth)
{
    size_t j = get_global_id(0);
    size_t i = get_global_id(1);
    int halfFilterSize = filterWidth / 2;

    float sum = 0;
    // Apply the filter to the neighborhood
    for (int k = -halfFilterSize; k <= halfFilterSize; k++)
    {
        for (int l = -halfFilterSize; l <= halfFilterSize; l++)
        {
            if (i + k >= 0 && i + k < imageHeight && j + l >= 0 &&
                j + l < imageWidth)
            {
                sum += imageInput[(i + k) * imageWidth + j + l] *
                       filterInput[(k + halfFilterSize) * filterWidth + l +
                                   halfFilterSize];
            }
        }
    }

    imageOutput[i * imageWidth + j] = sum;
}

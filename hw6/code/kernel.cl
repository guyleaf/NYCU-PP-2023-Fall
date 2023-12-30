// __kernel void convolution(__global const float *imageInput,
//                           __global const float *filterInput,
//                           __global float *imageOutput, int imageWidth,
//                           int imageHeight, int filterWidth)
// {
//     size_t j = get_global_id(0);
//     size_t i = get_global_id(1);

//     int halfFilterSize = filterWidth / 2;
//     float sum = 0;
//     // Apply the filter to the neighborhood
//     for (int k = -halfFilterSize; k <= halfFilterSize; k++)
//     {
//         for (int l = -halfFilterSize; l <= halfFilterSize; l++)
//         {
//             if (i + k >= 0 && i + k < imageHeight && j + l >= 0 &&
//                 j + l < imageWidth)
//             {
//                 sum += imageInput[(i + k) * imageWidth + j + l] *
//                        filterInput[(k + halfFilterSize) * filterWidth + l +
//                                    halfFilterSize];
//             }
//         }
//     }

//     imageOutput[i * imageWidth + j] = sum;
// }

__kernel void convolution(__global const float *imageInput,
                          __global const float *filterInput,
                          __global float *imageOutput, int imageWidth,
                          int imageHeight, int filterWidth,
                          __local float *filterLocal)
{
    size_t localJ = get_local_id(0);
    size_t localI = get_local_id(1);
    size_t j = get_global_id(0);
    size_t i = get_global_id(1);

    if (localJ < filterWidth && localI < filterWidth)
    {
        filterLocal[localI * filterWidth + localJ] =
            filterInput[localI * filterWidth + localJ];
    }

    // Synchronize the read into LMEM
    barrier(CLK_LOCAL_MEM_FENCE);

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
                       filterLocal[(k + halfFilterSize) * filterWidth + l +
                                   halfFilterSize];
            }
        }
    }

    imageOutput[i * imageWidth + j] = sum;
}

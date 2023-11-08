#include "utility.h"

#include <stdint.h>
#include <stdlib.h>

double random_uniform_r(unsigned int *seedp)
{
    return (double)rand_r(seedp) / RAND_MAX;
}

// const __m256i full = _mm256_set1_epi32(INT32_MAX);

// __m256i random_avx_xorshift128plus_r(avx_xorshift128plus_key_t *key)
// {
//     return avx_xorshift128plus(key);
// }
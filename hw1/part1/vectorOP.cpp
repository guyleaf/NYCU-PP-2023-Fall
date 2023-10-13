#include "PPintrin.h"

// implementation of absSerial(), but it is vectorized using PP intrinsics
void absVector(float *values, float *output, int N)
{
    __pp_vec_float x;
    __pp_vec_float result;
    __pp_vec_float zero = _pp_vset_float(0.f);
    __pp_mask maskAll, maskIsNegative, maskIsNotNegative;

    //  Note: Take a careful look at this loop indexing.  This example
    //  code is not guaranteed to work when (N % VECTOR_WIDTH) != 0.
    //  Why is that the case?
    for (int i = 0; i < N; i += VECTOR_WIDTH)
    {
        // All ones
        maskAll = _pp_init_ones(N - i);

        // All zeros
        maskIsNegative = _pp_init_ones(0);

        // Load vector of values from contiguous memory addresses
        _pp_vload_float(x, values + i, maskAll);  // x = values[i];

        // Set mask according to predicate
        _pp_vlt_float(maskIsNegative, x, zero, maskAll);  // if (x < 0) {

        // Execute instruction using mask ("if" clause)
        _pp_vsub_float(result, zero, x, maskIsNegative);  //   output[i] = -x;

        // Inverse maskIsNegative to generate "else" mask
        maskIsNotNegative = _pp_mask_not(maskIsNegative);  // } else {

        // Prevent accessing invalid memory when VECTOR_WIDTH not divides N
        maskIsNotNegative = _pp_mask_and(maskIsNotNegative, maskAll);

        // Execute instruction ("else" clause)
        _pp_vload_float(result, values + i,
                        maskIsNotNegative);  //   output[i] = x; }

        // Write results back to memory
        _pp_vstore_float(output + i, result, maskAll);
    }
}

void clampedExpVector(float *values, int *exponents, float *output, int N)
{
    __pp_vec_float x;
    __pp_vec_int y;
    __pp_vec_float result;
    __pp_vec_int zero = _pp_vset_int(0);
    __pp_vec_int one = _pp_vset_int(1);
    __pp_vec_float maxValue = _pp_vset_float(9.999999f);
    __pp_mask inputMask, maskIsGreaterThanZero, maskIsNotClamped;
    __pp_mask allMask = _pp_init_ones();

    //
    // PP STUDENTS TODO: Implement your vectorized version of
    // clampedExpSerial() here.
    //
    // Your solution should work for any value of
    // N and VECTOR_WIDTH, not just when VECTOR_WIDTH divides N
    //
    for (int i = 0; i < N; i += VECTOR_WIDTH)
    {
        // Load vector of values from contiguous memory addresses
        _pp_vload_float(x, values + i, allMask);
        _pp_vload_int(y, exponents + i, allMask);

        // Set results to 1.0f as default
        result = _pp_vset_float(1.f);

        // while
        while (true)
        {
            // (count > 0)
            _pp_vgt_int(maskIsGreaterThanZero, y, zero, allMask);
            if (_pp_cntbits(maskIsGreaterThanZero) == 0)
            {
                break;
            }

            // result *= x;
            _pp_vmult_float(result, result, x, maskIsGreaterThanZero);

            // count--;
            _pp_vsub_int(y, y, one, allMask);
        }

        // Clampe results if > 9.999999f
        _pp_vgt_float(maskIsNotClamped, result, maxValue, allMask);
        _pp_vmove_float(result, maxValue, maskIsNotClamped);

        // mask invalid part of vector
        inputMask = _pp_init_ones(N - i);

        // Write results back to memory
        _pp_vstore_float(output + i, result, inputMask);
    }
}

// returns the sum of all elements in values
// You can assume N is a multiple of VECTOR_WIDTH
// You can assume VECTOR_WIDTH is a power of 2
float arraySumVector(float *values, int N)
{
    //
    // PP STUDENTS TODO: Implement your vectorized version of arraySumSerial
    // here
    //

    for (int i = 0; i < N; i += VECTOR_WIDTH)
    {
    }

    return 0.0;
}
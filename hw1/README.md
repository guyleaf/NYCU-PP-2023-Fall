# Programming Assignment I: SIMD Programming
**Parallel Programming by Prof. Yi-Ping You**

Assignment description: https://pp-f23.github.io/assignments/HW1

## Introduction
The purpose of this assignment is to familiarize yourself with SIMD (single instruction, multiple data) programming. Most modern processors include some form of vector operations (i.e., SIMD instructions), and some applications may take advantage of SIMD instructions to improve performance through vectorization. Although modern compilers support automatic vectorization optimizations, the capabilities of compilers to fully auto-vectorize a given piece of code are often limited. Fortunately, almost all compilers (targeted to processors with SIMD extensions) provide SIMD intrinsics to allow programmers to vectorize their code explicitly.

## Q1
### Q1-1
#### Question
Does the vector utilization increase, decrease or stay the same as `VECTOR_WIDTH` changes? Why?

#### Answer
1. The vector utilization is decreasing as `VECTOR_WIDTH` changes.

2. Due we assume `N` % `VECTOR_WIDTH` is always zero, **raising the value to the power of exponent** (while loop) is the only part to make it different.

    Because of the stop conditon of loop, exponent is not greater than zero, it has to wait until each exponent in the vector is less or equal to zero.

    Therefore, in the later loops, the operations which mask will contain more zero or negative values as `VECTOR_WIDTH` changes will decrease the utilization. It tells us the higher `VECTOR_WIDTH` isn't always having no drawback.

    As for the total of **vector lanes** or **utilized vector lanes** in other operations, e.g., clamping, are the same among different `VECTOR_WIDTH`.


#### Results
##### `#define VECTOR_WIDTH 2`
```
CLAMPED EXPONENT (required) 
Results matched with answer!
****************** Printing Vector Unit Statistics *******************
Vector Width:              2
Total Vector Instructions: 162515
Vector Utilization:        92.6%
Utilized Vector Lanes:     301043
Total Vector Lanes:        325030
************************ Result Verification *************************
ClampedExp Passed!!!
```

##### `#define VECTOR_WIDTH 4`
```
CLAMPED EXPONENT (required) 
Results matched with answer!
****************** Printing Vector Unit Statistics *******************
Vector Width:              4
Total Vector Instructions: 94571
Vector Utilization:        90.1%
Utilized Vector Lanes:     340985
Total Vector Lanes:        378284
************************ Result Verification *************************
ClampedExp Passed!!!
```

##### `#define VECTOR_WIDTH 8`
```
CLAMPED EXPONENT (required) 
Results matched with answer!
****************** Printing Vector Unit Statistics *******************
Vector Width:              8
Total Vector Instructions: 51627
Vector Utilization:        88.9%
Utilized Vector Lanes:     367037
Total Vector Lanes:        413016
************************ Result Verification *************************
ClampedExp Passed!!!
```

##### `#define VECTOR_WIDTH 16`
```
CLAMPED EXPONENT (required) 
Results matched with answer!
****************** Printing Vector Unit Statistics *******************
Vector Width:              16
Total Vector Instructions: 26967
Vector Utilization:        88.3%
Utilized Vector Lanes:     380885
Total Vector Lanes:        431472
************************ Result Verification *************************
ClampedExp Passed!!!
```

## Q2
### Q2-1
#### Question
Fix the code to make sure it uses aligned moves for the best performance.

Hint: we want to see `vmovaps` rather than `vmovups`.

#### Answer
According to the [Reference](#References), the AVX2 registers (YMM) are 32-byte registers.  
So, we should declare the data is aligned by 32 bytes instead of 16 bytes.

We could solve it by adjusting `__builtin_assume_aligned` function from `16` to `32`.

```cpp
void test1(float *__restrict a, float *__restrict b, float *__restrict c, int N)
{
    __builtin_assume(N == 1024);
    a = (float *)__builtin_assume_aligned(a, 32);
    b = (float *)__builtin_assume_aligned(b, 32);
    c = (float *)__builtin_assume_aligned(c, 32);

    fasttime_t time1 = gettime();
    for (int i = 0; i < I; i++)
    {
        for (int j = 0; j < N; j++)
        {
            c[j] = a[j] + b[j];
        }
    }
    fasttime_t time2 = gettime();

    double elapsedf = tdiff(time1, time2);
    std::cout << "Elapsed execution time of the loop in test1():\n"
              << elapsedf << "sec (N: " << N << ", I: " << I << ")\n";
}
```

#### Result
##### Difference
```
50,53c50,53
<       vmovups (%rbx,%rcx,4), %ymm0
<       vmovups 32(%rbx,%rcx,4), %ymm1
<       vmovups 64(%rbx,%rcx,4), %ymm2
<       vmovups 96(%rbx,%rcx,4), %ymm3
---
>       vmovaps (%rbx,%rcx,4), %ymm0
>       vmovaps 32(%rbx,%rcx,4), %ymm1
>       vmovaps 64(%rbx,%rcx,4), %ymm2
>       vmovaps 96(%rbx,%rcx,4), %ymm3
58,61c58,61
<       vmovups %ymm0, (%r14,%rcx,4)
<       vmovups %ymm1, 32(%r14,%rcx,4)
<       vmovups %ymm2, 64(%r14,%rcx,4)
<       vmovups %ymm3, 96(%r14,%rcx,4)
---
>       vmovaps %ymm0, (%r14,%rcx,4)
>       vmovaps %ymm1, 32(%r14,%rcx,4)
>       vmovaps %ymm2, 64(%r14,%rcx,4)
>       vmovaps %ymm3, 96(%r14,%rcx,4)
```
##### Unaligned ([Godbolt](https://godbolt.org/z/Pz4dzYqW7))
```
_Z4testPfS_S_i: # @_Z4testPfS_S_i
  xorl %eax, %eax
.LBB0_1: # =>This Loop Header: Depth=1
  xorl %ecx, %ecx
.LBB0_2: # Parent Loop BB0_1 Depth=1
  vmovups (%rdi,%rcx,4), %ymm0
  vmovups 32(%rdi,%rcx,4), %ymm1
  vmovups 64(%rdi,%rcx,4), %ymm2
  vmovups 96(%rdi,%rcx,4), %ymm3
  vaddps (%rsi,%rcx,4), %ymm0, %ymm0
  vaddps 32(%rsi,%rcx,4), %ymm1, %ymm1
  vaddps 64(%rsi,%rcx,4), %ymm2, %ymm2
  vaddps 96(%rsi,%rcx,4), %ymm3, %ymm3
  vmovups %ymm0, (%rdx,%rcx,4)
  vmovups %ymm1, 32(%rdx,%rcx,4)
  vmovups %ymm2, 64(%rdx,%rcx,4)
  vmovups %ymm3, 96(%rdx,%rcx,4)
  addq $32, %rcx
  cmpq $1024, %rcx # imm = 0x400
  jne .LBB0_2
  addl $1, %eax
  cmpl $20000000, %eax # imm = 0x1312D00
  jne .LBB0_1
  vzeroupper
  retq
```

##### Aligned ([Godbolt](https://godbolt.org/z/489xWc9ob))
```
_Z4testPfS_S_i: # @_Z4testPfS_S_i
  xorl %eax, %eax
.LBB0_1: # =>This Loop Header: Depth=1
  xorl %ecx, %ecx
.LBB0_2: # Parent Loop BB0_1 Depth=1
  vmovaps (%rdi,%rcx,4), %ymm0
  vmovaps 32(%rdi,%rcx,4), %ymm1
  vmovaps 64(%rdi,%rcx,4), %ymm2
  vmovaps 96(%rdi,%rcx,4), %ymm3
  vaddps (%rsi,%rcx,4), %ymm0, %ymm0
  vaddps 32(%rsi,%rcx,4), %ymm1, %ymm1
  vaddps 64(%rsi,%rcx,4), %ymm2, %ymm2
  vaddps 96(%rsi,%rcx,4), %ymm3, %ymm3
  vmovaps %ymm0, (%rdx,%rcx,4)
  vmovaps %ymm1, 32(%rdx,%rcx,4)
  vmovaps %ymm2, 64(%rdx,%rcx,4)
  vmovaps %ymm3, 96(%rdx,%rcx,4)
  addq $32, %rcx
  cmpq $1024, %rcx # imm = 0x400
  jne .LBB0_2
  addl $1, %eax
  cmpl $20000000, %eax # imm = 0x1312D00
  jne .LBB0_1
  vzeroupper
  retq
```

---

### Q2-2
#### Question
What speedup does the vectorized code achieve over the unvectorized code? What additional speedup does using `-mavx2` give (`AVX2=1` in the `Makefile`)? You may wish to run this experiment several times and take median elapsed times; you can report answers to the nearest 100% (e.g., 2×, 3×, etc). What can you infer about the bit width of the default vector registers on the PP machines? What about the bit width of the AVX2 vector registers.

Hint: Aside from speedup and the vectorization report, the most relevant information is that the data type for each array is `float`.

#### Answer
##### Unvectorized vs. Vectorized $$SpeedUp = \frac{8.28740}{2.62903} = 3.15226528 \approx 3.15× $$

##### Unvectorized vs. Vectorized with AVX2 $$SpeedUp = \frac{8.28740}{1.40664} = 5.89162828 \approx 5.89× $$

##### Vectorized vs. Vectorized with AVX2 $$SpeedUp = \frac{2.62903}{1.40664} = 1.8690141 \approx 1.87× $$

The size of `float` type is `32 bits`.

Although the speedup of vectorized is only close to 3x, because of the size of registers is usually power of 2, I infer the size of the default vector registers is `32 bits * 4 = 128 bits`.

The speedup of vectorized with AVX2 is close to 2x.  
So, I infer the size of the AVX2 vector registers is `128 bits * 2 = 256 bits`.

#### Experiments
##### Case 1
`make clean && make && ./test_auto_vectorize -t 1`
* 8.28264 sec
* 8.28274 sec
* 8.28431 sec
* 8.28616 sec
* 8.28620 sec
* <font color="#f00">8.28740 sec</font>
* 8.28750 sec
* 8.28782 sec
* 8.29040 sec
* 8.29208 sec
* 8.29663 sec

##### Case 2
`make clean && make VECTORIZE=1 && ./test_auto_vectorize -t 1`
* 2.62432 sec
* 2.62493 sec
* 2.62772 sec
* 2.62799 sec
* 2.62870 sec
* <font color="#f00">2.62903 sec</font>
* 2.62918 sec
* 2.63050 sec
* 2.63114 sec
* 2.63313 sec
* 2.63462 sec

##### Case 3
`make clean && make VECTORIZE=1 AVX2=1 && ./test_auto_vectorize -t 1`
* 1.40373 sec
* 1.40383 sec
* 1.40543 sec
* 1.40571 sec
* 1.40638 sec
* <font color="#f00">1.40664 sec</font>
* 1.40710 sec
* 1.40871 sec
* 1.40958 sec
* 1.42159 sec
* 1.83462 sec

---

### Q2-3
#### Question
Provide a theory for why the compiler is generating dramatically different assembly.

#### Answer

1. Because the `maxps` instruction acts like either A or B conditonal assignment, the compiler needs to be hint by the branch condition.

2. In the original version, the vertorized assembly code (need additional masking) costs more computation than unvectorized's. So, the compiler chooses to use the latter.


## References
* [Wikipedia: Advanced Vector Extensions](https://en.wikipedia.org/wiki/Advanced_Vector_Extensions)

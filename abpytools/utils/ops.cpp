#include <cstdio>
#include <xmmintrin.h>
#include "ops.h"

void subtract_op(double **A, double **B, double **C, int size) {

#pragma omp simd
    for (int i=0; i<size; ++i) {
//        printf("%f = %f - %f\n", *(C[i]), *(A[i]), *(B[i]));
        *(C[i]) = *(A[i]) - *(B[i]);
    }

}


void subtract_op_sse(double *A, double *B, double *C, int N) {

    __m128d x1, x2, result;

    for(int i =0; i<N; ++i) {

        x1 = _mm_load_pd((A+(i*2)));
        x2 = _mm_load_pd((B+(i*2)));

        result = _mm_sub_pd(x1, x2);

        _mm_store_pd(&C[i*2], result);
    }

}


void subtract_op_sequential(double *A, double *B, double *C, int N) {

    for (int i=0; i<N; ++i) {
        C[i] = A[i] - B[i];
    }

}


void subtract_op(double *A, double *B, double *C, int size) {

//#pragma omp simd
//    for (int i=0; i<size; ++i) {
////        printf("%f = %f - %f\n", C[i], A[i], B[i]);
//        C[i] = A[i] - B[i];
//    }

    int extra = size % 2 == 0? 0:1;
    int N = size/2;

#if __SSE4_2__
    subtract_op_sse(A, B, C, N);
    // more calculations if necessary -> sequential
    // for doubles this could be replaced with one hardcoded condition
    while (extra > 0) {
        C[size-extra] = A[size-extra] - B[size-extra];
        extra--;
    }
#else
    subtract_op_sequential(A, B, C, size);
#endif
}
#include <cstdio>
#include <cmath>
#include <malloc.h>
#if __SSE4_2__
#include <nmmintrin.h>
#endif
#include "ops.h"

}


void subtract_op_sse(double *A, double *B, double *C, int N) {

    __m128d x1, x2;

    for(int i =0; i<N; ++i) {

        x1 = _mm_load_pd(A+(i*2));
        x2 = _mm_load_pd(B+(i*2));

        _mm_store_pd(C+(i*2), _mm_sub_pd(x1, x2));
    }
}


void subtract_op_sequential(double *A, double *B, double *C, int N) {
    for (int i=0; i<N; ++i) {
        C[i] = A[i] - B[i];
    }
}


void subtract_op(double *A, double *B, double *C, int size) {

#if __SSE4_2__
    int extra = size % 2 == 0? 0:1;
    int N = size/2;
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
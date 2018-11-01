#include <cstdio>
#include <cmath>
#ifdef __IS_DARWIN__
#include <malloc.h>
#endif
#if __SSE4_2__
#include <nmmintrin.h>
#endif
#include "omp.h"
#include "ops.h"

static const int N_THREADS = omp_get_max_threads();

static const long int S64 = (1L << 52);
static const __m128d cAA = _mm_set1_pd((double) S64);
static const __m128d cBB = _mm_set1_pd((double) S64 * 1023);

inline __m128d abs_pd(__m128d x) {
    // taken from https://stackoverflow.com/questions/5508628/how-to-absolute-2-double-or-4-floats-using-sse-instruction-set-up-to-sse4
    static const __m128d sign_mask = _mm_set1_pd(-0.); // -0. = 1 << 63
    return _mm_andnot_pd(sign_mask, x); // !sign_mask & x
}


#define MM_EXTRACT_DOUBLE(v,i) _mm_cvtsd_f64(_mm_shuffle_pd(v, v, _MM_SHUFFLE2(0, i)))

__m128d exp_op_sse_helper_fast(const __m128d a) {

    // based on https://gitlab.mpcdf.mpg.de/bbramas/inastemp/blob/master/Src/SSE3/InaVecSSE3Double.hpp
    // and https://wapco.e-ce.uth.gr/2015/papers/SESSION3/WAPCO_3_5.pdf

    __m128d x, factor;

    static const __m128d l2e = _mm_set1_pd(1.442695040888963407359924681001892137426645954153);  // log2(e)
    static const __m128d cA  = _mm_set1_pd(-1.33266405715993271723e-03);
    static const __m128d cB  = _mm_set1_pd(-9.61837182960864275905e-03);
    static const __m128d cC  = _mm_set1_pd(-5.55040609720754696266e-02);
    static const __m128d cD  = _mm_set1_pd(-2.40226511645233870018e-01);
    static const __m128d cE  = _mm_set1_pd( 3.06852819617161765020e-01);
    static const __m128d cF  = _mm_set1_pd(-1.10150186041739869460e-12);
    static const __m128d cX  = _mm_set1_pd(-1.87582286605066256753e-06);
    static const __m128d cY  = _mm_set1_pd(-1.41484352491262699514e-05);
    static const __m128d cZ  = _mm_set1_pd(-1.55186852765468104613e-04);

    x = _mm_mul_pd(a, l2e);
    const __m128d fractional_part = _mm_sub_pd(x, _mm_floor_pd(x));

    factor = _mm_add_pd(_mm_mul_pd(_mm_add_pd(
                        _mm_mul_pd(_mm_add_pd( _mm_mul_pd(_mm_add_pd(
                        _mm_mul_pd(_mm_add_pd( _mm_mul_pd(_mm_add_pd(
                        _mm_mul_pd(_mm_add_pd( _mm_mul_pd(_mm_add_pd(_mm_mul_pd(
                        cX,  fractional_part), cY), fractional_part),
                        cZ), fractional_part), cA), fractional_part),
                        cB), fractional_part), cC), fractional_part),
                        cD), fractional_part), cE), fractional_part),
                        cF);

    x = _mm_sub_pd(x, factor);
    x = _mm_add_pd(_mm_mul_pd(cAA, x), cBB);


    alignas(64) long int allvalint[2] = {_mm_cvtsd_si64(x), _mm_cvtsd_si64(_mm_shuffle_pd(x, x, 1))};

    return _mm_castsi128_pd(_mm_set_epi64x(allvalint[1], allvalint[0]));
}


__m128d exp_op_sse_helper_faster(const __m128d a) {

    // based on https://gitlab.mpcdf.mpg.de/bbramas/inastemp/blob/master/Src/SSE3/InaVecSSE3Double.hpp
    // and https://wapco.e-ce.uth.gr/2015/papers/SESSION3/WAPCO_3_5.pdf

    __m128d x, factor;
    static const __m128d l2e = _mm_set1_pd(1.442695040888963407359924681001892137426645954153);  // log2(e)
    static const __m128d cC  = _mm_set1_pd(-7.92041454535668681958e-02);
    static const __m128d cD  = _mm_set1_pd(-2.24339532327269441936e-01);
    static const __m128d cE  = _mm_set1_pd( 3.03543677780836240743e-01);
    static const __m128d cF  = _mm_set1_pd( 1.06906116358144185133e-04);

    x = _mm_mul_pd(a, l2e);
    const __m128d fractional_part = _mm_sub_pd(x, _mm_floor_pd(x));

    factor = _mm_add_pd(_mm_mul_pd(_mm_add_pd(
                        _mm_mul_pd(_mm_add_pd( _mm_mul_pd(
                        cC,  fractional_part),
                        cD), fractional_part),
                        cE), fractional_part),
                        cF);


    x = _mm_sub_pd(x, factor);
    x = _mm_add_pd(_mm_mul_pd(cAA, x), cBB);


    alignas(64) long int allvalint[2] = {_mm_cvtsd_si64(x), _mm_cvtsd_si64(_mm_shuffle_pd(x, x, 1))};

    return _mm_castsi128_pd(_mm_set_epi64x(allvalint[1], allvalint[0]));
}


void subtract_op_sse(const double *A, const double *B, double *C, int N) {

    __m128d x1, x2;

    for(int i =0; i<N; i+=2) {

        x1 = _mm_load_pd(&A[i]);
        x2 = _mm_load_pd(&B[i]);

        _mm_store_pd(&C[i], _mm_sub_pd(x1, x2));
    }
}


void subtract_op_sequential(const double *A, const double *B, double *C, int N) {
    for (int i=0; i<N; ++i) {
        C[i] = A[i] - B[i];
    }
}


void subtract_op(const double *A, const double *B, double *C, int size) {

#if __SSE4_2__
    int extra = size % 2 == 0? 0:1;
    subtract_op_sse(A, B, C, size);
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


void multiply_op_sse(const double *A, const double *B, double *C, int N) {

    __m128d x1, x2;

    for(int i =0; i<N; i+=2) {

        x1 = _mm_load_pd(&A[i]);
        x2 = _mm_load_pd(&B[i]);

        _mm_store_pd(&C[i], _mm_mul_pd(x1, x2));
    }
}


void multiply_op_sequential(const double *A, const double *B, double *C, int N) {
    for (int i=0; i<N; ++i) {
        C[i] = A[i] * B[i];
    }
}


void multiply_op(const double *A, const double *B, double *C, int size) {

    #if __SSE4_2__
    int extra = size % 2 == 0? 0:1;
    multiply_op_sse(A, B, C, size);
    // more calculations if necessary -> sequential
    // for doubles this could be replaced with one hardcoded condition
    while (extra > 0) {
        C[size-extra] = A[size-extra] * B[size-extra];
        extra--;
    }
    #else
    multiply_op_sequential(A, B, C, size);
    #endif

}


void add_op_sse(const double *A, const double *B, double *C, int N) {

    __m128d x1, x2;

    for(int i =0; i<N; i+=2) {

        x1 = _mm_load_pd(&A[i]);
        x2 = _mm_load_pd(&B[i]);

        _mm_store_pd(&C[i], _mm_add_pd(x1, x2));
    }
}


void add_op_sequential(const double *A, const double *B, double *C, int N) {
    for (int i=0; i<N; ++i) {
        C[i] = A[i] + B[i];
    }
}


void add_op(const double *A, const double *B, double *C, int size) {

    #if __SSE4_2__
    int extra = size % 2 == 0? 0:1;
    add_op_sse(A, B, C, size);
    // more calculations if necessary -> sequential
    // for doubles this could be replaced with one hardcoded condition
    while (extra > 0) {
        C[size-extra] = A[size-extra] + B[size-extra];
        extra--;
    }
    #else
    add_op_sequential(A, B, C, size);
    #endif

}


double norm_op_sse(const double *A, int p, int size) {

//    int extra = size % 2 == 0? 0:1;
//    double log2_p = std::log2(p);
//    double reciprocal_power = static_cast<double>(1.0/p);
//
//    __m128 intermediate_f;
//    __m128d a, intermediate_d, result;
//
//    result = _mm_set1_pd(0.0);
//
//    for (int i = 0; i< size / 2; i++) {
//
//        a = _mm_load_pd(A+(i*2));
//
//        intermediate_d = _mm_mul_pd(a, a);
//
//        intermediate_d = exp_op_sse_helper(intermediate_d, log2_p);
//
//        a = _mm_add_pd(a, intermediate_d);
//
//    }
//
//    if (extra) {
//        __m128d extra_ = _mm_set_pd(0.0, A[size-1]);
//        a = _mm_add_pd(a, extra_);
//    }
//
//    // horizontal sum of result
//    a = _mm_hadd_pd(a, a);
//    a = exp_op_sse_helper(a, reciprocal_power);
//
//    return MM_EXTRACT_DOUBLE(a, 0);
}


void exp_op_sse(const double *A, double *B, int size) {

    __m128d a, result;
    int extra = size % 2 == 0? 0:1;
    int N = size-extra;

//    #pragma omp parallel for
    for (int i=0; i< N; i+=2) {

        a = _mm_load_pd(&A[i]);

        result = exp_op_sse_helper_faster(a);

        _mm_store_pd(&B[i], result);

    }

    if (extra) {
        __m128d extra_ = _mm_load_pd(&A[size-2]);
        result = exp_op_sse_helper_faster(extra_);
        B[size-1] = MM_EXTRACT_DOUBLE(result, 0);
    }

}


void exp_op(const double *A, double *B, int size) {

    exp_op_sse(A, B, size);
//    #pragma omp parallel for
//    for (int i = 0; i < size; i++)
//        B[i] = exp(A[i]);
}


double norm_op(const double *A, int p, int size) {

    #if __SSE4_2__
        norm_op_sse(A, p, size);
//    #else
//        subtract_op_sequential(A, B, C, size);
    #endif

}


double l2_norm_op_sse_helper(const double *A, int size) {

    __m128d result = _mm_setzero_pd();
    __m128d a;


    for (int i=0; i<size; i+=2) {

        a = _mm_load_pd(&A[i]);

        result = _mm_add_pd(_mm_mul_pd(a, a), result); // square a and add to result

    }

    result = _mm_hadd_pd(result, result);

    return MM_EXTRACT_DOUBLE(result, 0);
}


double l2_norm_op_sse(const double *A, int size) {

    int offset = size / N_THREADS;
    int extra = size % N_THREADS;
    double sum = 0;
    int i;

    #pragma omp parallel for num_threads(N_THREADS)
    for (i=0; i<N_THREADS; i++) {
        sum+=l2_norm_op_sse_helper(&(A[i*offset]), offset);
    }

    for (int j = size - extra; j<size; i++) {
        sum+=A[j];
    }

    return sqrt(sum);

}

double l2_norm_op_sequential(const double *A, int size) {

    double sum = 0.0;

    for (int i = 0; i<size; i++) {
        sum += pow(A[i], 2);
    }

    return sqrt(sum);
}

double l2_norm_op(const double *A, int size) {

    double result;

    // for small values SSE seems to get stuck..
    if (size < 100)
        goto SEQUENTIAL;

    result = l2_norm_op_sse(A, size);
    goto END;

    SEQUENTIAL:
    result = l2_norm_op_sequential(A, size);

    END:
    return result;

}


double l1_norm_sse_helper(const double *A, int size) {

    __m128d result = _mm_setzero_pd();
    __m128d a;

    for (int i=0; i<size; i+=2) {

        a = _mm_load_pd(&A[i]);

        result = _mm_add_pd(abs_pd(a), result); // add a to result

    }

    result = _mm_hadd_pd(result, result);

    return MM_EXTRACT_DOUBLE(result, 0);
}


double l1_norm_op_sequential(const double *A, int size) {

    double sum = 0.0;

    for (int i=0; i<size; i++) {
        sum+=fabs(A[i]);
    }

    return sum;

}


double l1_norm_op_sse(const double *A, int size) {

    int offset = size / N_THREADS;
    int extra = size % N_THREADS;
    double sum = 0;

    if (size < 100)
        goto SEQUENTIAL;

    #pragma omp parallel for num_threads(N_THREADS)
    for (int i=0; i<N_THREADS; i++) {
        sum+=l1_norm_sse_helper(&(A[i*offset]), offset);
    }

    for (int j = size - extra; j<size; j++) {
        sum+=fabs(A[j]);
    }
    goto END;

    SEQUENTIAL:

    sum = l1_norm_op_sequential(A, size);

    END:
    return sum;

}


double l1_norm_op(const double *A, int size) {

    double result;

    #if __SSE4_2__
        result = l1_norm_op_sse(A, size);
    #else
        result = l1_norm_op_sequential(A, size);
    #endif

    return result;
}
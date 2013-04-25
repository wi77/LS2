/*
  This file is part of LS² - the Localization Simulation Engine of FU Berlin.

  Copyright 2011-2013   Heiko Will, Marcel Kyas, Thomas Hillebrandt,
  Stefan Adler, Malte Rohde, Jonathan Gunthermann

  LS² is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  LS² is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with LS².  If not, see <http://www.gnu.org/licenses/>.

 */

#ifndef VECTOR_SHOOTER_H
#define VECTOR_SHOOTER_H 1

#define NUM_VERSION(major, minor, patchlevel) (major * 10000 + minor * 100 + patchlevel)

#undef GCC_VERSION
#define GCC_VERSION NUM_VERSION(__GNUC__, __GNUC_MINOR__, __GNUC_PATCHLEVEL__)

/*!
 * The cache line alignement.
 */
#ifndef CACHE_ALIGNMENT_BITS
#  define CACHE_ALIGNMENT_BITS 64
#endif

#undef CACHE_ALIGNMENT
#define CACHE_ALIGNMENT ((CACHE_ALIGNMENT_BITS)/sizeof(unsigned char))

/*! The default number of threads can be defined in the makefile - ignored by
 * library version, configurable.
 */
#ifndef NUM_THREADS
#  define NUM_THREADS 8
#endif

/*! The number of runs per pixel. RUNS has to be dividable by the size of
 * the vectors! We suggest 8. */
#ifndef RUNS
#  define RUNS (8*50)
#endif

/*! The length of the image square SIZE has to be dividable by NUM_THREADS
 *  or some pixels may be left out.
 */
#ifndef SIZE
#  define SIZE 1000
#endif

#if defined(STAND_ALONE)
#  if !defined(ALGORITHM)
#    error "You must define a valid algorithm."
#  endif
#  if !defined(ERRORMODEL)
#    error "You must define a valid error model."
#  endif
#endif


// This value indicates which error is considerable for your simulation
// (e.g. the estimated ranging error or better) and which is not. Values
// below this error are marked with colors.
#define VERY_GOOD_VALUES_COLOR ((int)(50)) 

#  include <math.h>
#  include <float.h>

// defines for SSE or AVX usage. Try to minimize ifdefs in c code files.
#ifdef __AVX__

#  ifdef __cplusplus
extern "C" {
#  endif
#  include "avx_mathfun.h"
#  ifdef __cplusplus
}
#  endif

#  define _mm_printf(x) (printf("%f %f %f %f %f %f %f %f\n",x[0],x[1],x[2],x[3],x[4],x[5],x[6],x[7]))
#  define VECTOR __m256
#  define VECTOR_OPS 8
#  define VECTOR_CONST_BROADCAST(V) { V, V, V, V, V, V, V, V }

#  define VECTOR_BROADCAST(x)     _mm256_broadcast_ss(x)
#  define VECTOR_BROADCASTF(x)    _mm256_set1_ps(x)
#  define VECTOR_ZEROUPPER()      _mm256_zeroupper();
#  define VECTOR_SUM(x)           (x[0]+x[1]+x[2]+x[3]+x[4]+x[5]+x[6]+x[7])
#  define VECTOR_SQRT(x)          _mm256_sqrt_ps(x)
#  define VECTOR_MIN(x, y)        _mm256_min_ps(x, y)
#  define VECTOR_MAX(x, y)        _mm256_max_ps(x, y)
#  define VECTOR_AND(x, y)        _mm256_and_ps(x, y)
#  define VECTOR_ANDNOT(x, y)     _mm256_andnot_ps(x, y)
#  define VECTOR_OR(x, y)         _mm256_or_ps(x, y)
#  define VECTOR_XOR(x, y)        _mm256_xor_ps(x, y)
#  define VECTOR_LE(x, y)         _mm256_cmp_ps(x, y, _CMP_LE_OQ)
#  define VECTOR_LT(x, y)         _mm256_cmp_ps(x, y, _CMP_LT_OQ)
#  define VECTOR_GE(x, y)         _mm256_cmp_ps(x, y, _CMP_GE_OQ)
#  define VECTOR_GT(x, y)         _mm256_cmp_ps(x, y, _CMP_GT_OQ)
#  define VECTOR_NE(x, y)         _mm256_cmp_ps(x, y, _CMP_NEQ_OQ)
#  define VECTOR_EQ(x, y)         _mm256_cmp_ps(x, y, _CMP_EQ_OQ)
#  define VECTOR_HADD(x, y)       _mm256_hadd_ps(x, y)
#  define VECTOR_BLENDV(x, y, z)  _mm256_blendv_ps(x, y, z)
#  define VECTOR_TEST_ALL_ONES(x) _mm256_testc_ps(x, _mm256_cmp_ps(x, x,_CMP_EQ_US))
#  define VECTOR_CMPLT(x,y)       _mm256_cmp_ps(x,y,_CMP_NGE_US)
#  define VECTOR_CMPGE(x,y)       _mm256_cmp_ps(x,y,_CMP_NLT_US)
#  define VECTOR_CMPLE(x,y)       _mm256_cmp_ps(x,y,_CMP_NGT_US)
#  define VECTOR_CMPGT(x,y)       _mm256_cmp_ps(x,y,_CMP_NLE_US)
#  define VECTOR_CMPEQ(x,y)       _mm256_cmp_ps(x,y,_CMP_EQ_US)
#  define VECTOR_ZERO()           _mm256_setzero_ps()
#  define VECTOR_LOG(x)            log256_ps(x)
#  define VECTOR_SIN(x)            sin256_ps(x)
#  define VECTOR_COS(x)            cos256_ps(x)

#else

#  ifdef __cplusplus
extern "C" {
#  endif
#  include "sse_mathfun.h"
#  ifdef __cplusplus
}
#  endif

#  define _mm_printf(x) (printf("%f %f %f %f\n",x[0],x[1],x[2],x[3]))
#  define VECTOR __m128
#  define VECTOR_OPS 4
#  define VECTOR_CONST_BROADCAST(V) { V, V, V, V }

#  define VECTOR_BROADCAST(x)     _mm_load1_ps(x)
#  define VECTOR_BROADCASTF(x)    _mm_set1_ps(x)
#  define VECTOR_ZEROUPPER()      do { } while (0)
#  define VECTOR_SUM(x)           (x[0]+x[1]+x[2]+x[3])
#  define VECTOR_SQRT(x)          _mm_sqrt_ps(x)
#  define VECTOR_MIN(x, y)        _mm_min_ps(x, y)
#  define VECTOR_MAX(x, y)        _mm_max_ps(x, y)
#  define VECTOR_LE(x, y)         _mm_cmple_ps(x, y)
#  define VECTOR_LT(x, y)         _mm_cmplt_ps(x, y)
#  define VECTOR_GE(x, y)         _mm_cmpge_ps(x, y)
#  define VECTOR_GT(x, y)         _mm_cmpgt_ps(x, y)
#  define VECTOR_NE(x, y)         _mm_cmpneq_ps(x, y)
#  define VECTOR_EQ(x, y)         _mm_cmpeq_ps(x, y)
#  define VECTOR_AND(x, y)        _mm_and_ps(x, y)
#  define VECTOR_ANDNOT(x, y)     _mm_andnot_ps(x, y)

#  define VECTOR_OR(x, y)         _mm_or_ps(x, y)
#  define VECTOR_XOR(x, y)        _mm_xor_ps(x, y)
#  define VECTOR_HADD(x, y)       _mm_hadd_ps(x, y)
#  define VECTOR_LOG(x)           log_ps(x)
#  define VECTOR_SIN(x)            sin_ps(x)
#  define VECTOR_COS(x)            cos_ps(x)

#ifdef __SSE4_1__
#  define VECTOR_BLENDV(x, y, z)  _mm_blendv_ps(x, y, z)
#  define VECTOR_TEST_ALL_ONES(x) _mm_test_all_ones(_mm_castps_si128(x))
#  define VECTOR_CEIL(x)          _mm_round_ps(x, _MM_FROUND_TO_POS_INF)
#  define VECTOR_TRUNCATE(x)      _mm_round_ps(x, _MM_FROUND_TO_ZERO)
#else
#  define VECTOR_BLENDV(x, y, z)  emul_mm_blendv_ps(x, y, z)
#  define _mm_blendv_ps(x, y, z)  emul_mm_blendv_ps(x, y, z)
static inline __m128 __attribute__((always_inline,pure,artificial))
emul_mm_blendv_ps(__m128 a, __m128 b, __m128 mask)
{
    __m128 result;
    __m128i _mask = _mm_castps_si128(mask);
    for (int i = 0; i < 4; i++)
        result[i] = (_mask[i] & 0x80000000) ? b[i] : a[i];
    return result;
}

#  define VECTOR_TEST_ALL_ONES(x) emul_mm_test_all_ones(_mm_castps_si128(x))
#  define _mm_test_all_ones(x)    emul_mm_test_all_ones(x)
static inline int __attribute__((__always_inline__))
emul_mm_test_all_ones(__m128i bits)
{
    int result = 1;
    for (int i = 0; i < 4; i++)
        result = result && (bits[i] == 0xFFFFFFFF);
    return result;
}

#  define VECTOR_CEIL(x)          emul_mm_ceil(x)
static inline __m128  __attribute__((__always_inline__))
emul_mm_ceil(__m128 a)
{
    __m128 result;
    for (int i = 0; i < 4; i++)
        result[i] = (float) ceil((double) a[i]);
    return result;
}

#  define VECTOR_FLOOR(x)         emul_mm_floor(x)
static inline __m128  __attribute__((__always_inline__))
emul_mm_floor(__m128 a)
{
    __m128 result;
    for (int i = 0; i < 4; i++)
        result[i] = (float) floor((double) a[i]);
    return result;
}

#endif

#  define VECTOR_ZERO()           _mm_setzero_ps()

#  define IVECTOR __m128i
typedef union {
    __m128i v;
    int m[4];
} ivector_u;
#endif /* ! __AVX__ */

#define VECTOR_ONES()                   VECTOR_EQ(zero, zero)
#define VECTOR_ABS(x)                   VECTOR_MAX(x, -(x))
#define VECTOR_NOT(x)                   VECTOR_XOR(x, VECTOR_ONES())
#define VECTOR_CLAMP(val, low, high)    VECTOR_MIN(VECTOR_MAX(val, low), high)

#if defined(__cplusplus) && (GCC_VERSION < 40701)
#  define VECTOR_GET(v, i) (v)[i]
#else
#  define VECTOR_GET(v, i) v[i]
#endif

static const VECTOR zero = VECTOR_CONST_BROADCAST(0.0f);
static const VECTOR half = VECTOR_CONST_BROADCAST(0.5f);
static const VECTOR one =  VECTOR_CONST_BROADCAST(1.0f);
static const VECTOR two =  VECTOR_CONST_BROADCAST(2.0f);
static const VECTOR rmax = VECTOR_CONST_BROADCAST(((float) RAND_MAX) * 2.0f);
static const VECTOR eps =  VECTOR_CONST_BROADCAST(FLT_EPSILON);





/* ***************************************
 *  Internal structures an variables
 * ***************************************/

// Dirty helpers to generate function names
#define EMF(func,name) func##_##name
#define EMFU(func,name) EMF(func,name)
#define EMFUNCTION(name) EMFU(ERRORMODEL,name)

#define LF(func) func##_##run
#define LFU(func) LF(func)
#define ALGORITHM_RUN LFU(ALGORITHM)

// Include helper for errormodel
#define QUOTEME(M) #M
#define IEM(M) QUOTEME(error_model/M##_em.c)
#define INCLUDE_EM(M) IEM(M)

#define IEMH(M) QUOTEME(error_model/M##_em.h)
#define INCLUDE_EM_H(M) IEMH(M)

#ifndef ALG_DIR
#  define ALG_DIR algorithm
#endif

// Include helper for algorithms
#define IALG(D,M) QUOTEME(D/M##_algorithm.c)
#define INCLUDE_ALG(M) IALG(ALG_DIR,M)

#define IALGH(D,M) QUOTEME(D/M##_algorithm.h)
#define INCLUDE_ALG_H(M) IALGH(ALG_DIR,M)


#include "ls2/ls2.h"

#ifndef UNITTEST
#  define UNITTEST 0
#endif

#if (UNITTEST == 1)
int test_triangle_area();  
int test_triangle_perimeter();
int test_circle_get_intersection();
int test_circle_get_approx_intersection();
int test_point_geometric_median();
int test_center_of_mass();
#endif

#endif

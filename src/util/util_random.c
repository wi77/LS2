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

/********************************************************************
 **
 **  This file is made only for including in the LS² project
 **  and not desired for stand alone usage!
 **
 ********************************************************************/

/*******************************************************************
 ***
 *** Random Numbers
 ***
 *******************************************************************/

#ifndef INCLUDED_UTIL_RANDOM_H
#define INCLUDED_UTIL_RANDOM_H

#include <limits.h>
#include <immintrin.h>

#include "vector_shooter.h"

// calculate four 32 bit random integer between -RAN_DMAX and RAN_DMAX in a vector
// based on http://software.intel.com/en-us/articles/fast-random-number-generator-on-the-intel-pentiumr-4-processor/
static inline __m128 __attribute__((__always_inline__,__gnu_inline__,__nonnull__,__artificial__))
rand_sse(__m128i* cur_seed)
{
    __m128i cur_seed_split;
    __m128i multiplier;
    __m128i adder;
    __m128i mod_mask;

    static const unsigned int __attribute__((__aligned__(16),__may_alias__))
        gadd[4] = { 2531011, 10395331, 13737667, 1 };
    static const unsigned int __attribute__((__aligned__(16),__may_alias__))
        mult[4] = { 214013, 17405, 214013, 69069 };
    static const unsigned int __attribute__((__aligned__(16),__may_alias__))
        mask[4] = { 0xFFFFFFFF, 0, 0xFFFFFFFF, 0 };

    adder = _mm_load_si128( (__m128i*) gadd);
    multiplier = _mm_load_si128( (__m128i*) mult);
    mod_mask = _mm_load_si128( (__m128i*) mask);

    cur_seed_split = _mm_shuffle_epi32(*cur_seed, _MM_SHUFFLE(2, 3, 0, 1));
    *cur_seed = _mm_mul_epu32(*cur_seed, multiplier);
    multiplier = _mm_shuffle_epi32(multiplier, _MM_SHUFFLE(2, 3, 0, 1));
    cur_seed_split = _mm_mul_epu32(cur_seed_split, multiplier);
    *cur_seed = _mm_and_si128(*cur_seed, mod_mask);
    cur_seed_split = _mm_and_si128(cur_seed_split, mod_mask);
    cur_seed_split = _mm_shuffle_epi32(cur_seed_split, _MM_SHUFFLE(2, 3, 0, 1));
    *cur_seed = _mm_or_si128(*cur_seed, cur_seed_split );
    *cur_seed = _mm_add_epi32(*cur_seed, adder);
    return _mm_cvtepi32_ps(*cur_seed);
}


// Returns a vector of random numbers between 0..1
#if !defined(__RDRND__)
#  if defined(__AVX__)
static inline VECTOR __attribute__((__always_inline__,__gnu_inline__,__nonnull__,__artificial__))
rnd(__m128i *seed)
{
    __m128 tmp[2];
    tmp[0] = rand_sse(seed);
    tmp[1] = rand_sse(seed);
    VECTOR ret = _mm256_load_ps((float *) tmp);
    ret /= rmax;
    ret += half;
    return ret;
}
#  else
static inline VECTOR __attribute__((__always_inline__,__gnu_inline__,__nonnull__,__artificial__))
rnd(__m128i *seed)
{
    VECTOR ret = rand_sse(seed);
    ret = ret / rmax;
    ret += half;
    return ret;
}
#  endif
#else /* defined(__RDRAND__) */

static inline void
__attribute__((__always_inline__,__gnu_inline__,__artificial__,__nonnull__))
// __attribute__((optimize(0)))
rand64(unsigned long long *result)
{
    int c;
    do {
        c = _rdrand64_step(result);
    } while (c == 0);
}

static inline VECTOR
__attribute__((__always_inline__,__gnu_inline__,__artificial__,__nonnull__))
rnd(__m128i *seed __attribute__((__unused__)))
{
    static const VECTOR rnd_divisor = VECTOR_CONST_BROADCAST((float) UINT_MAX);
    union {
        uint32_t ui[VECTOR_OPS];
        uint64_t ull[VECTOR_OPS/2];
    } t;

    rand64(&(t.ull[0]));
    rand64(&(t.ull[1]));
#if VECTOR_OPS == 8
    rand64(&(t.ull[2]));
    rand64(&(t.ull[3]));
#endif

    const VECTOR v =
        { (float) t.ui[0], (float) t.ui[1], (float) t.ui[2], (float) t.ui[3]
#if VECTOR_OPS == 8
        , (float) t.ui[4], (float) t.ui[5], (float) t.ui[6], (float) t.ui[7]
#endif
        }; 

    return v / rnd_divisor;
}
#endif




#ifndef NSUM
#  define NSUM 25
#endif

static inline VECTOR
__attribute__((__always_inline__,__gnu_inline__,__artificial__,__nonnull__))
gaussrand(__m128i* seed, float mean, float sdev)
{
    VECTOR gauss = VECTOR_ZERO();
    int i;

    for(i = 0; i < NSUM; i++){
        gauss += rnd(seed);
    }
    gauss -= VECTOR_BROADCASTF(NSUM / 2.0f);
    gauss /= VECTOR_SQRT(VECTOR_BROADCASTF(NSUM / 12.0f));
    // Multiply with desired standard deviation
    gauss *= VECTOR_BROADCASTF(sdev);
    // Add mean
    gauss += VECTOR_BROADCASTF(mean);
        
    return gauss;
}

// Creates a exponential distribution with rate,
// also shifts the rate, see parameter help
static inline VECTOR
__attribute__((__always_inline__,__gnu_inline__,__artificial__,__nonnull__))
exp_rand(__m128i* seed, VECTOR rate)
{
    VECTOR result = rnd(seed);
    result = - VECTOR_LOG(result) / rate;
    return result;
}

#endif

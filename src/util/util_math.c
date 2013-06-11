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
 **  This file is made only for including in the lib_lat project
 **  and not desired for stand alone usage!
 **
 ********************************************************************/

#ifndef UTIL_MATH_C_INCLUDED
#define UTIL_MATH_C_INCLUDED 1

/** All long-representable factorials */
static const long long factorials[] =  {
    1LL, 1LL, 2LL,
    6LL, 24LL, 120LL,
    720LL, 5040LL, 40320LL,
    362880LL, 3628800LL, 39916800LL,
    479001600LL, 6227020800LL, 87178291200LL,
    1307674368000LL, 20922789888000LL, 355687428096000LL,
    6402373705728000LL, 121645100408832000LL, 2432902008176640000LL
};

 

/*! Calculate binomial coefficient */
static inline int
__attribute__((__always_inline__,__gnu_inline__,__const__,__artificial__))
binom(int n, int k)
{
    if (n < k)
        return 0;
    long long ret = n;
    for (int i = n - 1; i > n - k; i--) {
        ret *= i;
    }

    return (int) (ret / factorials[k]);
}


/*! Round a value N up to the next largest value divisible by K */
static inline long
__attribute__((__always_inline__,__gnu_inline__,__const__,__artificial__))
iceil(const long N, const long K)
{
    const long r = N % K;
    if (r != 0)
        return N + K - r;
    else
        return N;
}

/** Utility method for k-permutation building */
static inline int incCounter(int v[], int i, int n, int k) {
    v[i]++;
    return (v[i] == n - ((k - 1) - i)); // return overflow information
}

#endif

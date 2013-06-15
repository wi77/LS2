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

#ifndef UTIL_MATRIX_C_INCLUDED
#define UTIL_MATRIX_C_INCLUDED

/*******************************************************************
 ***
 ***   Parallel Matrix operations. 
 ***   A matrix is defined as an array of vectors
 ***   having row*cols entries. So a single vector m(a,b) is accesable
 ***   with m[a*cols+b]
 ***
 *******************************************************************/

// Transpose a Matrix
static inline void
__attribute__((__always_inline__,__gnu_inline__,__nonnull__,__artificial__))
transpose(const VECTOR *restrict a, VECTOR *restrict ret,
          const size_t rows, const size_t cols)
{
    for (size_t i = 0; i < rows; i++) {
        for (size_t j = 0; j < cols; j++) {
            ret[j*rows + i] = a[i*cols + j];
        }
    }
}

// Invert a 2x2 Matrix
static inline void
__attribute__((__always_inline__,__gnu_inline__,__nonnull__,__artificial__))
invert_2x2(const VECTOR * restrict a, VECTOR * restrict ret)
{
    VECTOR m_one = VECTOR_BROADCASTF(-1.0F); 
    VECTOR det = one / ((a[0] * a[3])-(a[1] * a[2]));
    ret[0] = det * a[3];
    ret[1] = det * m_one * a[1];
    ret[2] = det * m_one * a[2];
    ret[3] = det * a[0];   
}

// Multiply a Matrix
static inline void
__attribute__((__always_inline__,__gnu_inline__,__nonnull__,__artificial__))
times(const VECTOR * a, const size_t arows, const size_t acols,
      const VECTOR * b, const size_t bcols, VECTOR * restrict ret)
{
    // Ways faster than memset!
    for (size_t i=0; i < arows * bcols; i++)
        ret[i] = zero;

    for (size_t i = 0; i < arows; i++) {
        for (size_t j = 0; j < bcols; j++) {
            for (size_t k = 0; k < acols; k++) {
                ret[i*bcols+j] += (a[i*acols+k] * b[k*bcols+j]);
            }
        }
    }
}

#endif

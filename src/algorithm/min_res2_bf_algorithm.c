/*
  This file is part of LS² - the Localization Simulation Engine of FU Berlin.

  Copyright 2011-2013  Heiko Will, Marcel Kyas, Thomas Hillebrandt,
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
 **  and not intended for stand alone usage!
 **
 ********************************************************************/

/*******************************************************************
 ***
 ***   Min-Max aka Bounding Box algorithm.
 ***
 *******************************************************************/

/* @algorithm_name: Minimise residuals (norm 2, brute force) */

static inline void __attribute__((__always_inline__,__gnu_inline__,__nonnull__,__artificial__)) __attribute__((__always_inline__,__artificial__,__nonnull__))
min_res2_bf_run(const VECTOR* ax, const VECTOR* ay, const VECTOR *restrict r,
                const size_t num_anchors, const size_t width, const size_t height,
                VECTOR *restrict resx, VECTOR *restrict resy)
{
    VECTOR min_res = VECTOR_BROADCASTF(FLT_MAX),
           rx = VECTOR_BROADCASTF(0.0f),
           ry = VECTOR_BROADCASTF(0.0f);
    for (size_t y = 0; y < height; y++) {
        for (size_t x = 0; x < width; x++) {
            VECTOR res = VECTOR_ZERO();
            const VECTOR vx = VECTOR_BROADCASTF((float) x),
                         vy = VECTOR_BROADCASTF((float) y);
            for (size_t a = 0; a < num_anchors; a++) {
                const VECTOR d =  distance(vx, vy, ax[a], ay[a]) - r[a];
                res += d * d;
            }
            const VECTOR m = VECTOR_GT(min_res, res);
            rx = VECTOR_BLENDV(rx, vx, m);
            ry = VECTOR_BLENDV(ry, vy, m);
            min_res = VECTOR_MIN(min_res, res);
        }
    }
    *resx = rx;
    *resy = ry;
}

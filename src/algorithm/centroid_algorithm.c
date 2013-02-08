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
 ***   Centroid algorithm
 ***
 *******************************************************************/

/* @algorithm_name: Centroid */

static inline void __attribute__((__always_inline__,__gnu_inline__,__nonnull__,__artificial__))
centroid_run (const VECTOR* vx, const VECTOR* vy,
              const VECTOR *restrict r, size_t num_anchors, int width __attribute__((__unused__)), int height __attribute__((__unused__)), VECTOR *restrict resx, VECTOR *restrict resy)
{
    // px and py hold the intersections, sx, sy the sum of the intersection
    // coordinates to accumulate for the final result.
    // n is the number of intersections returned by vcircle_intersections;
    VECTOR px[2], py[2], sx = zero, sy = zero, n = zero;

    // Step 1: Calculate all intersections.
    for (size_t i = 0; i < num_anchors; i++) {
        for (size_t j = i + 1; j < num_anchors; j++) {
            n += vcircle_intersections(vx[i], vy[i], vx[j], vy[j], r[i], r[j],
                                       px, py);
            sx += px[0] + px[1];
            sy += py[0] + py[1];
        }
    }

    // Step 2: Compute the center of mass.
    *resx = sx / VECTOR_MAX(n, one);
    *resy = sy / VECTOR_MAX(n, one);
}

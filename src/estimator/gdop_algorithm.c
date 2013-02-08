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

/*******************************************************************
 ***
 ***   Constant value localization (debug only)
 ***
 *******************************************************************/

/* @algorithm_name: Gdop estimation */

static inline float
__attribute__((__always_inline__,__gnu_inline__,__pure__,__nonnull__,__artificial__))
gdop_run(const vector2 *restrict anchor, const size_t num_anchors,
         const vector2 *restrict location)
{
    float A[num_anchors][2];
    float At[2][num_anchors];
    float P[2][2];
    float AmulAt[2][2];

    for (size_t i = 0; i < num_anchors; i++) {
         A[i][0] = (anchor[i].x - location->x) / distance_v(location, &(anchor[i]));
         A[i][1] = (anchor[i].y - location->y) / distance_v(location, &(anchor[i]));
    }

    // transpose
    for (size_t i = 0; i< num_anchors; i++) {
        for (int j = 0; j < 2; j++) {
            At[j][i] = A[i][j];
        }
    }

    // multiply
    AmulAt[0][0] = At[0][0] * A[0][0] + At[0][1] * A[1][0] + At[0][2] * A[2][0];
    AmulAt[0][1] = At[0][0] * A[0][1] + At[0][1] * A[1][1] + At[0][2] * A[2][1];
    AmulAt[1][0] = At[1][0] * A[0][0] + At[1][1] * A[1][0] + At[1][2] * A[2][0];
    AmulAt[1][1] = At[1][0] * A[0][1] + At[1][1] * A[1][1] + At[1][2] * A[2][1];
    
    // invert
    float det = AmulAt[0][0] * AmulAt[1][1] - AmulAt[0][1] * AmulAt[1][0];
    
    P[0][0] = 1.0f/det * AmulAt[1][1];
    P[0][1] = -1.0f/det * AmulAt[0][1];
    P[1][0] = -1.0f/det * AmulAt[1][0];
    P[1][1] = 1.0f/det * AmulAt[0][0];

    float gdop = sqrtf(P[0][0] + P[1][1]);
    return gdop;
}

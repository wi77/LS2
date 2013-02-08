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
 **  and not desired for stand alone usage!
 **
 ********************************************************************/

/*******************************************************************
 ***
 ***   Voting Based Location Estimation
 ***
 *******************************************************************/
 
 /* @algorithm_name: Voting Based Location Estimation */

#define MAX(a, b)  (((a) > (b)) ? (a) : (b))
#define MIN(a, b)  (((a) < (b)) ? (a) : (b))

static inline void __attribute__((__always_inline__,__gnu_inline__,__nonnull__,__artificial__))
vble_run (const VECTOR* vx, const VECTOR* vy, const VECTOR *restrict r, size_t num_anchors, int width __attribute__((__unused__)), int height __attribute__((__unused__)), VECTOR *restrict resx, VECTOR *restrict resy) {
    int ii;
    float l;
    float error_threshold = 100.0F; // set to max error of error model
    if (num_anchors<3) return;
    for (ii = 0; ii < VECTOR_OPS; ii++) {
        // step 1: find minimum rectangle that covers all anchors
        float maxRanging = FLT_MIN;
        float minX = FLT_MAX, maxX = 0;
        float minY = FLT_MAX, maxY = 0;
        for (size_t i = 0; i < num_anchors; i++) {
            if (vx[i][ii] > maxX) {
                maxX = vx[i][ii];
            }
            if (vx[i][ii] < minX) {
                minX = vx[i][ii];
            }
            if (vy[i][ii] > maxY) {
                maxY = vy[i][ii];
            }
            if (vy[i][ii] < minY) {
                minY = vy[i][ii];
            }
            if (r[i][ii] > maxRanging) {
                maxRanging = r[i][ii];
            }
        }

        // step 2: extend rectangle by maximum transmission range of a beacon
        //         signal, here: extend rectangle by maximum ranging value
        minX -= maxRanging;
        maxX += maxRanging;
        minY -= maxRanging;
        maxY += maxRanging;

        // for the sake of comparison, set L - the side length of a cell in
        // meters (grid step size) to the end value of our optimized version
        float min_side = MIN(maxX - minX, maxY - minY);
        l = 0.05F * min_side;

        // step 3: devide rectangle into M small squares (cells) with the same
        //         side length L. No iterative refinement or other optimizations
        //         are used in this version, because we don't run on resource
        //         constrained sensor nodes. NOTE: Use appropriate side length
        //         value conversion if value is not given in meters!
        int score;
        int maxScore = 0;
        int maxScoreIndex = 0;
        float xMaxScore = 0;
        float yMaxScore = 0;
        float x, y;
        for (x = minX; x < maxX; x += l) {
            for (y = minY; y < maxY; y += l) {
                score = 0;
                for (size_t i = 0; i < num_anchors; i++) {
                    // calculate candidate ring with given error threshold in
                    // meters. NOTE: Use appropriate error threshold value
                    // conversion if value is not given in meters!
                    // Test if cell overlaps with current candidate ring,
                    // code not optimized!
                    float dMin, dMax;
                    float ri = ((r[i][ii] - error_threshold) > 0) ? (r[i][ii] - error_threshold) : 0;
                    float ro = r[i][ii] + error_threshold;

                    // precalculate some distances
                    __m128 vxi = _mm_setr_ps(vx[i][ii],vx[i][ii],vx[i][ii],vx[i][ii]);
                    __m128 vyi = _mm_setr_ps(vy[i][ii],vy[i][ii],vy[i][ii],vy[i][ii]);
                    __m128 xl = _mm_setr_ps(x,x+l,x,x+l);
                    __m128 yl = _mm_setr_ps(y,y,y+l,y+l);
                    __m128 dxy = _mm_sqrt_ps((vxi-xl)*(vxi-xl)+(vyi-yl)*(vyi-yl));
                    if (vx[i][ii] < x && vy[i][ii] < y) {
                        // sector 1
                        dMin = dxy[0];
                        dMax = dxy[3];
                    } else if (vx[i][ii] > x+l && vy[i][ii] < y) {
                        // sector 3
                        dMin = dxy[1];
                        dMax = dxy[2];
                    } else if (vx[i][ii] < x && vy[i][ii] > y+l) {
                        // sector 7
                        dMin = dxy[2];
                        dMax = dxy[1];
                    } else if (vx[i][ii] > x+l && vy[i][ii] > y+l) {
                        // sector 9
                        dMin = dxy[3];
                        dMax = dxy[0];
                    } else if (vy[i][ii] < y &&vx[i][ii] >= x && vx[i][ii] <= x+l) {
                        // sector 2
                        dMin = y - vy[i][ii];
                        if (vx[i][ii] - x > (x+l) - vx[i][ii]) {
                            dMax = dxy[2];
                        } else {
                            dMax = dxy[3];
                        }
                    } else if (vy[i][ii] > y+l && vx[i][ii] >= x && vx[i][ii] <= x+l) {
                        // sector 8
                        dMin = vy[i][ii] - (y+l);
                        if (vx[i][ii] - x > (x+l) - vx[i][ii]) {
                            dMax = dxy[0];
                        } else {
                            dMax = dxy[1];
                        }
                    } else if (vx[i][ii] < x && vy[i][ii] >= y && vy[i][ii] <= y+l) {
                        // sector 4
                        dMin = x - vx[i][ii];
                        if (vy[i][ii] - y > (y+l) - vy[i][ii]) {
                            dMax = dxy[1];
                        } else {
                            dMax = dxy[3];
                        }
                    } else if (vx[i][ii] > x+l && vy[i][ii] >= y && vy[i][ii] <= y+l) {
                        // sector 6
                        dMin = vx[i][ii] - (x+l);
                        if (vy[i][ii] - y > (y+l) - vy[i][ii]) {
                            dMax = dxy[0];
                        } else {
                            dMax = dxy[2];
                        }
                    } else {
                        // sector 5
                        dMin = 0;
                        float distTopL = dxy[0];
                        float distTopR = dxy[1];
                        float distBottomL = dxy[2];
                        float distBottomR = dxy[3];
                        dMax = MAX(MAX(distTopL, distTopR), MAX(distBottomL, distBottomR));
                    }

                    // test if candidate ring overlaps with cell
                    if (!(dMin > ro || dMax < ri)) {
                        score++;
                    }
                }
                if (score >= maxScore) {
                    if (score > maxScore) {
                        maxScore = score;
                        xMaxScore = 0.0L;
                        yMaxScore = 0.0L;
                        maxScoreIndex = 0.0L;
                    }
                    if (maxScore > 0) {
                        xMaxScore += x + l/2.0F;
                        yMaxScore += y + l/2.0F;
                        maxScoreIndex++;
                    }
                }
            }
        }
        
        // step 4: return geometric centroid from the cells with the highest
        //         vote as the estimated location
        (*resx)[ii] = xMaxScore / (float)maxScoreIndex;
        (*resy)[ii] = yMaxScore / (float)maxScoreIndex;
    }
}


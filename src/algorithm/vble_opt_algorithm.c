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
 ***   Optimized Voting Based Location Estimation
 ***
 *******************************************************************/

 /* @algorithm_name: Optimized Voting Based Location Estimation */

/*
 * ivector_u is currently only available to SSE4.1 builds.
 */
#if !__AVX__ && __SSE4_1__
static inline int 
ivector_hmax(ivector_u a) {
    int max = a.m[0];
    if(a.m[1] > max)
        max = a.m[1];
    if(a.m[2] > max)
        max = a.m[2];
    if(a.m[3] > max)
        max = a.m[3];
    return max;
}

typedef struct {
    VECTOR x;
    VECTOR y;
    VECTOR ri;
    VECTOR ro;
    VECTOR iBoxMinX;
    VECTOR iBoxMinY;
    VECTOR iBoxMaxX;
    VECTOR iBoxMaxY;
} anchor_info_t;

static inline void __attribute__((__always_inline__,__gnu_inline__,__nonnull__,__artificial__))
multilaterate(anchor_info_t *anchors, size_t num_anchors,
              VECTOR L, VECTOR last_L, VECTOR minX, VECTOR maxX,
              VECTOR minY, VECTOR maxY, VECTOR *restrict resx, VECTOR *restrict resy)
{
    // do iterative/recursive solution
    ivector_u maxScore;
    ivector_u maxScoreIndex = { _mm_set1_epi32(1) };
    VECTOR finalX = zero;
    VECTOR finalY = zero;

    VECTOR iterMinX = minX;
    VECTOR iterMaxX = maxX;
    VECTOR iterMinY = minY;
    VECTOR iterMaxY = maxY;
    
    // Note: The vectorized implementation can result in more iterations for some
    // values.
    while (! VECTOR_TEST_ALL_ONES(VECTOR_LT(L, last_L))) {
        finalX = zero;
        finalY = zero;
        maxScore.v = _mm_set1_epi32(0);
        maxScoreIndex.v = _mm_set1_epi32(0);
        
        ivector_u xLength = { _mm_cvttps_epi32((maxX - minX) / L + one) };
        ivector_u yLength = { _mm_cvttps_epi32((maxY - minY) / L + one) };            
        ivector_u xMin, xMax, yMin, yMax;
        
        int xLengthMax = ivector_hmax(xLength) + 1;
        int yLengthMax = ivector_hmax(yLength) + 1;
        char scores[xLengthMax*yLengthMax][VECTOR_OPS];
        memset(scores, 0, sizeof(char*) * (size_t)(VECTOR_OPS*xLengthMax*yLengthMax));
        
        VECTOR minXNew = VECTOR_BROADCASTF(FLT_MAX); 
        VECTOR maxXNew = VECTOR_BROADCASTF(FLT_MIN);
        VECTOR minYNew = VECTOR_BROADCASTF(FLT_MAX);
        VECTOR maxYNew = VECTOR_BROADCASTF(FLT_MIN);
        
        // TODO comment the following section
        ivector_u a = { _mm_cvttps_epi32((iterMinX - minX) / L) };
        ivector_u b = { _mm_cvttps_epi32((iterMinY - minY) / L) };
        ivector_u c = { _mm_cvttps_epi32(VECTOR_CEIL((iterMaxX - minX) / L)) };
        ivector_u d = { _mm_cvttps_epi32(VECTOR_CEIL((iterMaxY - minY) / L)) };
        VECTOR ones = VECTOR_ONES();
        VECTOR abcdmsk = _mm_or_ps(_mm_or_ps(
                        _mm_xor_ps(ones, (_mm_castsi128_ps(_mm_cmpeq_epi32(a.v, _mm_setzero_si128())))), // a != 0
                        _mm_xor_ps(ones, (_mm_castsi128_ps(_mm_cmpeq_epi32(b.v, _mm_setzero_si128())))) // b != 0
                        ), _mm_or_ps(
                        _mm_xor_ps(ones, (VECTOR_EQ(maxX - iterMaxX, zero))), // maxX - iterMaxX != 0
                        _mm_xor_ps(ones, (VECTOR_EQ(maxY - iterMaxY, zero)))  // maxY - iterMaxY != 0
                        ));
        
        for (size_t i = 0; i < num_anchors; i++) {
            
            // optimize outer test region
            xMin.v = _mm_cvttps_epi32((anchors[i].x - anchors[i].ro - minX) / L);
            yMin.v = _mm_cvttps_epi32((anchors[i].y - anchors[i].ro - minY) / L);
            ivector_u itmp = { _mm_cvttps_epi32(VECTOR_CEIL(two * anchors[i].ro / L + one)) };
            xMax.v = xMin.v + itmp.v;
            yMax.v = yMin.v + itmp.v;

            // for iteration: find overlaping rectangles
            VECTOR anchorRectMinX = minX + _mm_cvtepi32_ps(xMin.v) * L;
            VECTOR anchorRectMinY = minY + _mm_cvtepi32_ps(yMin.v) * L;
            VECTOR anchorRectMaxX = minX + _mm_cvtepi32_ps(xMax.v) * L;
            VECTOR anchorRectMaxY = minY + _mm_cvtepi32_ps(yMax.v) * L;
            
            // test for intersection
            VECTOR iblacklist = VECTOR_OR(
                VECTOR_OR(VECTOR_LE(iterMaxX, anchorRectMinX), VECTOR_GE(iterMinX, anchorRectMaxX)), 
                VECTOR_OR(VECTOR_LE(iterMaxY, anchorRectMinY), VECTOR_GE(iterMinY, anchorRectMaxY)));

            // TODO comment
            xMin.v = _mm_castps_si128(VECTOR_BLENDV(_mm_castsi128_ps(xMin.v), _mm_castsi128_ps(_mm_max_epi32(xMin.v, a.v)), abcdmsk));
            yMin.v = _mm_castps_si128(VECTOR_BLENDV(_mm_castsi128_ps(yMin.v), _mm_castsi128_ps(_mm_max_epi32(yMin.v, b.v)), abcdmsk));
            xMax.v = _mm_castps_si128(VECTOR_BLENDV(_mm_castsi128_ps(xMax.v), _mm_castsi128_ps(_mm_min_epi32(xMax.v, c.v)), abcdmsk));
            yMax.v = _mm_castps_si128(VECTOR_BLENDV(_mm_castsi128_ps(yMax.v), _mm_castsi128_ps(_mm_min_epi32(yMax.v, d.v)), abcdmsk));

            // possible NaN workaround
            xMax.v += _mm_castps_si128(VECTOR_AND(_mm_castsi128_ps(_mm_set1_epi32(1)), VECTOR_AND(abcdmsk, _mm_castsi128_ps(_mm_cmpeq_epi32(xMin.v, xMax.v)))));
            yMax.v += _mm_castps_si128(VECTOR_AND(_mm_castsi128_ps(_mm_set1_epi32(1)), VECTOR_AND(abcdmsk, _mm_castsi128_ps(_mm_cmpeq_epi32(yMin.v, yMax.v)))));             
                                
            // test each cell and increase score counter, use outer
            // test optimization
            ivector_u j = { xMin.v };
            while (1) {
                
                // Break condition: All j elements have passed their corresponding
                // xMax values.
                VECTOR jgexMax = VECTOR_GE(_mm_cvtepi32_ps(j.v), _mm_cvtepi32_ps(xMax.v));
                if (VECTOR_TEST_ALL_ONES(jgexMax))
                    break;
                
                // We blacklist those vector elements where j has passed the
                // xMax value for the current iteration.
                VECTOR xblacklist = jgexMax;
                
                // Calculate cell's x-coordinate.
                VECTOR x = minX + _mm_cvtepi32_ps(j.v) * L;
                                        
                ivector_u k = { yMin.v };
                while (1) {
                    
                    // Same here.
                    VECTOR kgeyMax = VECTOR_GE(_mm_cvtepi32_ps(k.v), _mm_cvtepi32_ps(yMax.v));
                    if (VECTOR_TEST_ALL_ONES(kgeyMax))
                        break;
                    
                    VECTOR yblacklist = kgeyMax;
                    
                    // Calculate cell's y-coordinate.
                    VECTOR y = minY + _mm_cvtepi32_ps(k.v) * L;
                                                                    
                    // Check for inner cell optimization.
                    // Blacklist optimized elements and remember that we 
                    // already incremented k.
                    VECTOR ico = VECTOR_AND(VECTOR_AND(VECTOR_GT(x, anchors[i].iBoxMinX), VECTOR_LT(x + L, anchors[i].iBoxMaxX)),
                                            VECTOR_AND(VECTOR_GT(y, anchors[i].iBoxMinY), VECTOR_LT(y + L, anchors[i].iBoxMaxY)));
                    VECTOR skip = VECTOR_TRUNCATE(((anchors[i].iBoxMaxY - y) - L) / L);
                    skip = VECTOR_BLENDV(one, skip, VECTOR_GT(skip, one));
                    k.v += _mm_cvttps_epi32(VECTOR_BLENDV(zero, skip, ico));                        
                    ivector_u kincr = { _mm_cvttps_epi32(VECTOR_BLENDV(one, zero, ico)) };
                    yblacklist = VECTOR_OR(yblacklist, ico);
                                        
                    // Optimization.
                    if (VECTOR_TEST_ALL_ONES(yblacklist))
                        continue;                       
                                      
                    // determine minimum and maximum distance to cell boundary
                    VECTOR closestX = VECTOR_MIN(VECTOR_MAX(anchors[i].x, x), x+L);
                    VECTOR closestY = VECTOR_MIN(VECTOR_MAX(anchors[i].y, y), y+L);
                    VECTOR farXoffset = VECTOR_MAX(closestX - x, (x+L) - closestX);
                    VECTOR farYoffset = VECTOR_MAX(closestY - y, (y+L) - closestY);
                    VECTOR minDistX = VECTOR_MAX(closestX - anchors[i].x, anchors[i].x - closestX);
                    VECTOR minDistY = VECTOR_MAX(closestY - anchors[i].y, anchors[i].y - closestY);
                    VECTOR maxDistX = minDistX + farXoffset;
                    VECTOR maxDistY = minDistY + farYoffset;
                    VECTOR dMin = minDistX * minDistX + minDistY * minDistY;
                    VECTOR dMax = maxDistX * maxDistX + maxDistY * maxDistY;                   

                    // xyblacklists can contain NaN. Make sure it only contains true or false
                    // so that the test below can work.
                    VECTOR blacklist = VECTOR_AND(one, VECTOR_OR(iblacklist, VECTOR_OR(xblacklist, yblacklist)));
                    
                    // Calculate offsets in scores array.
                    // TODO I'd love to know why gcc stops me from initializing scoresIdx directly (like I did with every ivector_u above!).
                    ivector_u scoresIdx;
                    scoresIdx.v = _mm_mullo_epi32(j.v, yLength.v + _mm_set1_epi32(1)) + k.v;
                    
                    // Test if candidate ring overlaps with cell.
                    VECTOR sx = VECTOR_AND(one, VECTOR_ANDNOT(blacklist, VECTOR_AND(
                        VECTOR_LE(dMin, anchors[i].ro * anchors[i].ro), VECTOR_GE(dMax, anchors[i].ri))));
                    
                    // Increase cell's scores if candidate ring overlaps.
                    ivector_u cell_scores;
                    cell_scores.m[0] = scores[scoresIdx.m[0]][0] = (char) (scores[scoresIdx.m[0]][0] + sx[0]);
                    cell_scores.m[1] = scores[scoresIdx.m[1]][1] = (char) (scores[scoresIdx.m[1]][1] + sx[1]);
                    cell_scores.m[2] = scores[scoresIdx.m[2]][2] = (char) (scores[scoresIdx.m[2]][2] + sx[2]);
                    cell_scores.m[3] = scores[scoresIdx.m[3]][3] = (char) (scores[scoresIdx.m[3]][3] + sx[3]);
                    
                    // Update next iteration's test region in case the highest ranked cells have changed.
                    VECTOR maxScoreExceeded = VECTOR_ANDNOT(blacklist, _mm_castsi128_ps(_mm_cmpgt_epi32(cell_scores.v, maxScore.v)));
                    maxScore.v = _mm_castps_si128(VECTOR_BLENDV(_mm_castsi128_ps(maxScore.v), _mm_castsi128_ps(cell_scores.v), maxScoreExceeded));

                    VECTOR maxScoreMatched = VECTOR_ANDNOT(blacklist, VECTOR_AND(_mm_castsi128_ps(_mm_cmpgt_epi32(cell_scores.v, _mm_setzero_si128())), 
                                                                                 _mm_castsi128_ps(_mm_cmpeq_epi32(cell_scores.v, maxScore.v))));
                    
                    VECTOR updateMinXNew = VECTOR_OR(maxScoreExceeded, VECTOR_AND(maxScoreMatched, VECTOR_LT(x, minXNew)));
                    minXNew = VECTOR_BLENDV(minXNew, x, updateMinXNew);
                    
                    VECTOR updateMaxXNew = VECTOR_OR(maxScoreExceeded, VECTOR_AND(maxScoreMatched, VECTOR_GT(x, maxXNew)));
                    maxXNew = VECTOR_BLENDV(maxXNew, x, updateMaxXNew);
                    
                    VECTOR updateMinYNew = VECTOR_OR(maxScoreExceeded, VECTOR_AND(maxScoreMatched, VECTOR_LT(y, minYNew)));
                    minYNew = VECTOR_BLENDV(minYNew, y, updateMinYNew);
                    
                    VECTOR updateMaxYNew = VECTOR_OR(maxScoreExceeded, VECTOR_AND(maxScoreMatched, VECTOR_GT(y, maxYNew)));
                    maxYNew = VECTOR_BLENDV(maxYNew, y, updateMaxYNew);
                    
                    VECTOR updateFinalXYAndIndex = VECTOR_OR(maxScoreExceeded, maxScoreMatched);
                    VECTOR newFinalX = VECTOR_BLENDV(x + finalX, x, maxScoreExceeded);
                    VECTOR newFinalY = VECTOR_BLENDV(y + finalY, y, maxScoreExceeded);
                    VECTOR newMaxScoreIndex = VECTOR_BLENDV(_mm_castsi128_ps(maxScoreIndex.v + _mm_set1_epi32(1)), _mm_castsi128_ps(_mm_set1_epi32(1)), maxScoreExceeded);
                    finalX = VECTOR_BLENDV(finalX, newFinalX, updateFinalXYAndIndex);
                    finalY = VECTOR_BLENDV(finalY, newFinalY, updateFinalXYAndIndex);
                    maxScoreIndex.v = _mm_castps_si128(VECTOR_BLENDV(_mm_castsi128_ps(maxScoreIndex.v), newMaxScoreIndex, updateFinalXYAndIndex));

                    
                    k.v += kincr.v;
                }

                j.v += _mm_set1_epi32(1);
            }
        }
        
        iterMinX = minXNew;
        iterMaxX = maxXNew+L;
        iterMinY = minYNew;
        iterMaxY = maxYNew+L;
        L /= two;
    } 
        
    // final position is the centroid of the highest ranked cells
    // Note: we add L (which was L/2 in the last iteration) here to indicate
    //       the centers of the cells.         
    *resx = finalX / _mm_cvtepi32_ps(maxScoreIndex.v) + L;
    *resy = finalY / _mm_cvtepi32_ps(maxScoreIndex.v) + L;   
}


static inline void __attribute__((__always_inline__,__gnu_inline__,__nonnull__,__artificial__))
vble_opt_run (const VECTOR* vx, const VECTOR* vy, const VECTOR *restrict r, size_t num_anchors, int width __attribute__((__unused__)), int height __attribute__((__unused__)), VECTOR *restrict resx, VECTOR *restrict resy) {
    if (num_anchors<3) {
        (*resx) = VECTOR_BROADCASTF(NAN);
        (*resy) = VECTOR_BROADCASTF(NAN);
        return;
    }
    
    const VECTOR error_threshold = VECTOR_BROADCASTF(85.0f);
    const VECTOR sqrtTwo = VECTOR_BROADCASTF(sqrtf(2.0f));
    
    // concentrate anchor information in a struct to increase data locality
    anchor_info_t anchors[num_anchors] __attribute__ ((aligned(ALIGNMENT)));
    
    // step 1: find minimum rectangle that covers all anchors
    VECTOR maxRanging = VECTOR_BROADCASTF(FLT_MIN);
    VECTOR minX = VECTOR_BROADCASTF(FLT_MAX), maxX = VECTOR_BROADCASTF(0);
    VECTOR minY = VECTOR_BROADCASTF(FLT_MAX), maxY = VECTOR_BROADCASTF(0);        
    for(size_t i = 0; i < num_anchors; i++) {
        maxX = VECTOR_BLENDV(maxX, vx[i], VECTOR_GT(vx[i], maxX));
        minX = VECTOR_BLENDV(minX, vx[i], VECTOR_LT(vx[i], minX));
        maxY = VECTOR_BLENDV(maxY, vy[i], VECTOR_GT(vy[i], maxY));
        minY = VECTOR_BLENDV(minY, vy[i], VECTOR_LT(vy[i], minY));
        maxRanging = VECTOR_BLENDV(maxRanging, r[i], VECTOR_GT(r[i], maxRanging));
        
        // meanwhile, precalculate candidate rings + inner test regions
        
        // calculate candidate ring with given error threshold in
        // meters. NOTE: Use appropriate error threshold value
        // conversion if value is not given in meters!
        // Test if cell overlaps with current candidate ring, in
        // contrast to the authors we dont't use negative distance
        // measurement errors!
        anchors[i].ri = VECTOR_BLENDV(zero, r[i] - error_threshold, 
                        VECTOR_GT(r[i] - error_threshold, zero));
        anchors[i].ro = r[i]; // we don't measure too short, so
                              // don't add error threshold
                                
        // optimize inner test region
        VECTOR tmp = (sqrtTwo * anchors[i].ri) / two;
        anchors[i].iBoxMinX = vx[i] - tmp;
        anchors[i].iBoxMinY = vy[i] - tmp;
        anchors[i].iBoxMaxX = vx[i] + tmp;
        anchors[i].iBoxMaxY = vy[i] + tmp;
        
        // calculate squared inner radius as we test squared distances
        anchors[i].ri = anchors[i].ri * anchors[i].ri;
        
        // store anchor position
        anchors[i].x = vx[i];
        anchors[i].y = vy[i];
    }
            
    // step 2: extend rectangle by maximum transmission range of a beacon
    //         signal, here: extend rectangle by maximum ranging value
    minX -= maxRanging;
    maxX += maxRanging;
    minY -= maxRanging;
    maxY += maxRanging;

    // step 3: calculate L - the side length of a cell in meters
    //         (grid step size). Use 40 percent of the rectangles
    //         shorter side.
    VECTOR min_side = VECTOR_MIN(maxX - minX, maxY - minY);
    VECTOR l = VECTOR_BROADCASTF(0.4f) * min_side;
    VECTOR last_l = VECTOR_BROADCASTF(0.05f) * min_side;
    
    multilaterate(anchors, num_anchors, l, last_l, minX, maxX, minY, maxY, resx, resy);   
}
#else
static inline void __attribute__((__always_inline__,__gnu_inline__,__nonnull__,__artificial__))
multilaterate(const float* anchorsx, const float* anchorsy, const float* ranges,
              size_t num_anchors, float L, float last_L,
              float errorThreshold, float minX, float maxX,
              float minY, float maxY, float iterMinX, float iterMaxX,
              float iterMinY, float iterMaxY, float *resx, float *resy)
{
        // do iterative/recursive solution
        int maxScore;
        int maxScoreIndex = 1;
        float finalX = 0, finalY = 0;
        const float sqrtTwo = sqrtf(2.0F);
        while (L >= last_L) {
            finalX = 0;
            finalY = 0;
            maxScore = 0;
            maxScoreIndex = 0;
            int xLength = (int) ((maxX - minX)/L + 1.0F);
	    int yLength = (int) ((maxY - minY)/L + 1.0F);
            float x, y;
            int xMin, xMax, yMin, yMax;
            char *scores = calloc((size_t)((yLength+1)*(xLength+1)), sizeof(char));
            float minXNew = FLT_MAX, maxXNew = FLT_MIN, minYNew = FLT_MAX, maxYNew = FLT_MIN;
            for (size_t i = 0; i < num_anchors; i++) {
                // calculate candidate ring with given error threshold in
                // meters. NOTE: Use appropriate error threshold value
                // conversion if value is not given in meters!
                // Test if cell overlaps with current candidate ring, in
                // contrast to the authors we dont't use negative distance
                // measurement errors!
                float ri = ((ranges[i] - errorThreshold) > 0)
                        ? (ranges[i] - errorThreshold) : 0;
                float ro = ranges[i]; // we don't measure too short, so
                                      // don't add error threshold

                // optimize outer test region
                xMin = (int) ((anchorsx[i] - ro - minX) / L);
                yMin = (int) ((anchorsy[i] - ro - minY) / L);
                int itmp = (int) ceil(2.0F * ro / L + 1.0F);
                xMax = xMin + itmp;
                yMax = yMin + itmp;

                // for iteration: find overlaping rectangles
		float anchorRectMinX = minX + (float) xMin * L;
            	float anchorRectMinY = minY + (float) yMin * L;
           	float anchorRectMaxX = minX + (float) xMax * L;
            	float anchorRectMaxY = minY + (float) yMax * L;

                // test for intersection
            	if (iterMaxX <= anchorRectMinX || iterMinX >= anchorRectMaxX || iterMaxY <= anchorRectMinY || iterMinY >= anchorRectMaxY) {
                    continue; // no intersection
            	}

                int a = (int) ((iterMinX - minX) / L);
                int b = (int) ((iterMinY - minY) / L);
                int c = (int) ceil((iterMaxX - minX) / L);
                int d = (int) ceil((iterMaxY - minY) / L);
                if (a != 0 || b != 0 || maxX-iterMaxX != 0 || maxY-iterMaxY != 0) {
                    xMin = MAX(xMin, a);
                    yMin = MAX(yMin, b);
                    xMax = MIN(xMax, c);
                    yMax = MIN(yMax, d);
		        if (xMin == xMax) {
                        xMax++; // possible NaN workaround
                    }
                    if (yMin == yMax) {
                        yMax++; // possible NaN workaround
                    }
                }

                // optimize inner test region
                float tmp = (sqrtTwo * ri) / 2.0F;
                float iboxMinX = anchorsx[i] - tmp;
                float iBoxMinY = anchorsy[i] - tmp;
                float iboxMaxX = anchorsx[i] + tmp;
                float iBoxMaxY = anchorsy[i] + tmp;

                // test each cell and increase score counter, use outer
                // test optimization
                for (int j = xMin; j < xMax; j++) {
                    x = minX + (float) j * L;
                    int k = yMin;
                    while (k < yMax) {
                        y = minY + (float) k * L;
                        float dMin, dMax;
                        // check for inner cell optimization
                        if (x > iboxMinX && x + L < iboxMaxX && y > iBoxMinY && y + L < iBoxMaxY) {
                            int skip = (int) (((iBoxMaxY - y) - L) / L);
                            k += skip > 1 ? skip : 1;
                            continue;
                        }
                        // precalculate some distances
                        __m128 vxi = _mm_setr_ps(anchorsx[i],anchorsx[i],anchorsx[i],anchorsx[i]);
                        __m128 vyi = _mm_setr_ps(anchorsy[i],anchorsy[i],anchorsy[i],anchorsy[i]);
                        __m128 xl = _mm_setr_ps(x,x+L,x,x+L);
                        __m128 yl = _mm_setr_ps(y,y,y+L,y+L);
                        __m128 dxy = _mm_sqrt_ps((vxi-xl)*(vxi-xl)+(vyi-yl)*(vyi-yl));
                        if (anchorsx[i] < x && anchorsy[i] < y) {
                            // sector 1
                            dMin = dxy[0];
                            dMax = dxy[3];
                        } else if (anchorsx[i] > x+L && anchorsy[i] < y) {
                            // sector 3
                            dMin = dxy[1];
                            dMax = dxy[2];
                        } else if (anchorsx[i] < x && anchorsy[i] > y+L) {
                            // sector 7
                            dMin = dxy[2];
                            dMax = dxy[1];
                        } else if (anchorsx[i] > x+L && anchorsy[i] > y+L) {
                            // sector 9
                            dMin = dxy[3];
                            dMax = dxy[0];
                        } else if (anchorsy[i] < y && anchorsx[i] >= x && anchorsx[i] <= x+L) {
                            // sector 2
                            dMin = y - anchorsy[i];
                            if (anchorsx[i] - x > (x+L) - anchorsx[i]) {
                                dMax = dxy[2];
                            } else {
                                dMax = dxy[3];
                            }
                        } else if (anchorsy[i] > y+L && anchorsx[i] >= x && anchorsx[i] <= x+L) {
                            // sector 8
                            dMin = anchorsy[i] - (y+L);
                            if (anchorsx[i] - x > (x+L) - anchorsx[i]) {
                                dMax = dxy[0];
                            } else {
                                dMax = dxy[1];
                            }
                        } else if (anchorsx[i] < x && anchorsy[i] >= y && anchorsy[i] <= y+L) {
                            // sector 4
                            dMin = x - anchorsx[i];
                            if (anchorsy[i] - y > (y+L) - anchorsy[i]) {
                                dMax = dxy[1];
                            } else {
                                dMax = dxy[3];
                            }
                        } else if (anchorsx[i] > x+L && anchorsy[i] >= y && anchorsy[i] <= y+L) {
                            // sector 6
                            dMin = anchorsx[i] - (x+L);
                            if (anchorsy[i] - y > (y+L) - anchorsy[i]) {
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
                            scores[k*(xLength+1)+j]++;
                        }

                        if (scores[k*(xLength+1)+j] >= maxScore) {
                            if (scores[k*(xLength+1)+j] > maxScore) {
                                maxScore = scores[k*(xLength+1)+j];
                                maxScoreIndex = 0;
                                minXNew = FLT_MAX;
                                maxXNew = FLT_MIN;
                                minYNew = FLT_MAX;
                                maxYNew = FLT_MIN;
                                finalX = 0;
                                finalY = 0;
                            }
                            if (maxScore > 0) {
                                float tx = minX + (float)j*L;
                                float ty = minY + (float)k*L;
                                if (tx < minXNew) {
                                    minXNew = tx;
                                }
                                if (tx > maxXNew) {
                                    maxXNew = tx;
                                }
                                if (ty < minYNew) {
                                    minYNew = ty;
                                }
                                if (ty > maxYNew) {
                                    maxYNew = ty;
                                }
                                finalX += tx + L/2.0F;
                                finalY += ty + L/2.0F;
                                maxScoreIndex++;
                            }
                        }

                        k++;
                    }
                }
            }
            iterMinX = minXNew;
            iterMaxX = maxXNew+L;
            iterMinY = minYNew;
            iterMaxY = maxYNew+L;
            L /= 2;
            free(scores);
         }
    *resx = finalX / (float) maxScoreIndex;
    *resy = finalY / (float) maxScoreIndex;   
}


static inline void __attribute__((__always_inline__,__gnu_inline__,__nonnull__,__artificial__))
vble_opt_run (const VECTOR* vx, const VECTOR* vy, const VECTOR *restrict r, size_t num_anchors, int width __attribute__((__unused__)), int height __attribute__((__unused__)), VECTOR *restrict resx, VECTOR *restrict resy) {
    int ii;
    float l, last_l;
    const float error_threshold = 85.0F;
    if (num_anchors<3) return;
    for (ii = 0; ii < VECTOR_OPS; ii++) {
        float anchorsx[num_anchors];
        float anchorsy[num_anchors];
        float ranges[num_anchors];
        // step 1: find minimum rectangle that covers all anchors
        float maxRanging = FLT_MIN;
        float minX = FLT_MAX, maxX = 0;
        float minY = FLT_MAX, maxY = 0;
        for (size_t i = 0; i < num_anchors; i++) {
            anchorsx[i] = vx[i][ii];
            anchorsy[i] = vy[i][ii];
            ranges[i] = r[i][ii];
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

        // step 3: calculate L - the side length of a cell in meters
        //         (grid step size). Use 40 percent of the rectangles
        //         shorter side.
        float min_side = MIN(maxX - minX, maxY - minY);
        l = 0.4F * min_side;
	last_l = 0.05F * min_side;

        multilaterate(anchorsx, anchorsy, ranges, num_anchors, l, last_l, error_threshold, minX, maxX, minY, maxY, minX, maxX, minY, maxY, &((*resx)[ii]),&((*resy)[ii]));
    }
}

#endif

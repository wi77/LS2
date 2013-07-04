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
 ***   Geolateration
 ***
 *******************************************************************/

/* @algorithm_name: Geolateration */

/*!
 * \arg \c mcount 
 * \arg \c mx     
 * \arg \c my     
 * \arg \c r
 * \arg \c ptscount Number of circle intersections.
 * \arg \c ptsy     X coordinate of circle intersections.
 * \arg \c ptsy     Y coordinate of circle intersections.
 * \arg \c min
 */
static inline size_t
circle_minimum_circle_containment(size_t mcount, float *mx, float *my,
                                  float *r, size_t ptscount, float *ptsx,
				  float *ptsy, float *ptsw, size_t min)
{
    float vx [ptscount];
    float vy [ptscount];
    float vw [ptscount];
    size_t icount = 0;
    for (size_t i = 0; i < ptscount; i++) {
        size_t min_c = 0;
        for (size_t j = 0; j < mcount; j++) {
            if (distance_s(mx[j],my[j],ptsx[i],ptsy[i]) <= r[j] + 0.01F) {
                min_c ++;
            }
        }
        if (min_c >= min || ptsw[i] == 0.5F) {
            vx[icount] = ptsx[i];
            vy[icount] = ptsy[i];
            vw[icount] = ptsw[i];
            icount ++;
        }
    }
    memcpy (ptsx, vx, sizeof(float) * icount);
    memcpy (ptsy, vy, sizeof(float) * icount);
    memcpy (ptsw, vw, sizeof(float) * icount);
    return icount;
}


static inline void __attribute__((__always_inline__,__gnu_inline__,__nonnull__,__artificial__))
geon_run (const VECTOR* vx, const VECTOR* vy,
          const VECTOR *restrict r, size_t num_anchors, int width __attribute__((__unused__)), int height __attribute__((__unused__)), VECTOR *restrict resx, VECTOR *restrict resy)
{
    if (num_anchors < 3) return;

    // step 1: calculate circle intersections
    size_t bin = num_anchors * (num_anchors - 1);

    for (int ii = 0 ; ii < VECTOR_OPS; ii++) {
        float intersectionsx[bin];
        float intersectionsy[bin];
        float intersectionsw[bin];
        size_t icount = 0;
        size_t tmp = 0;
        
        // build all k-permutations
        for (size_t i = 0; i < num_anchors-1; i++) {
            for (size_t j = i+1; j < num_anchors; j++) {
                // set weight for intersection
                float weight = 1;
                // Berechne Schnittpunkte mit aktueller Permutation
                tmp = circle_get_intersection(vx[i][ii], vy[i][ii], vx[j][ii],
                                              vy[j][ii], r[i][ii], r[j][ii],
					      &intersectionsx[icount],
                                              &intersectionsy[icount]);
                intersectionsw[icount] = weight;
                intersectionsw[icount+1] = weight;
                
                if (tmp == 0) {
                    // no intersection => try to get approximated intersection
                    tmp = circle_get_approx_intersection(
                            vx[i][ii], vy[i][ii], vx[j][ii], vy[j][ii],
                            r[i][ii], r[j][ii], &intersectionsx[icount],
                            &intersectionsy[icount]);
                    weight = 0.5F;
                    intersectionsw[icount] = weight;
                }
                icount += tmp;                
            }
        }

        // step 2: copy intersections into one array

        // done in step 1

        // step 3: filter intersections points: only keep points which are
        //         contained in anchor length - 2 circles.
        float ax[num_anchors];
        float ay[num_anchors];
        float ar[num_anchors];
        for (size_t i = 0; i <num_anchors; i++) {
            ax[i] = vx[i][ii];   
            ay[i] = vy[i][ii];   
            ar[i] = r[i][ii];        
        }
        icount = circle_minimum_circle_containment(num_anchors, ax, ay, ar, icount, intersectionsx, intersectionsy, intersectionsw, num_anchors - 2);
       
        // step 4: if there are n*(n-1)/2 points which are very close together
        //          => no ranging error, take one of them as result
        size_t close_num_anchors = (num_anchors * (num_anchors - 1)) / 2;
        int wasbreak = 0;
        for (size_t i = 0; i < icount; i++) {
            size_t current_close_num_anchors = 1;
            for (size_t j = 0; j < icount; j++) {
                if (i != j && distance_s(intersectionsx[i],intersectionsy[i],intersectionsx[j],intersectionsy[j]) < 0.1F) {
                    current_close_num_anchors++;
                }
            }
            if (current_close_num_anchors >= close_num_anchors) {
                (*resx)[ii] = intersectionsx[i];
                (*resy)[ii] = intersectionsy[i];
                wasbreak = 1;
                break;
            }
        }
        if (wasbreak) continue;

        // step 5: apply median filter on remaining points
        float distances[icount];
        memset(distances, 0, icount * sizeof(float));
        for (size_t i = 0; i < icount; i++) {
            for (size_t j = 0; j < icount; j++) {
                if (i != j) {
                    distances[i] += distance_s(intersectionsx[i],intersectionsy[i],intersectionsx[j],intersectionsy[j]);
                }
            }
        }

        const float median = fselect_s(distances, icount, icount / 2);
        
        float vvx[icount];
        float vvy[icount];
        float vvw[icount];
        int vvcount = 0;
        
        for (size_t i = 0; i < icount; i++) {
	    if (icount >=3) {
		    if (distances[i] <= median * 1.0F) { // Median factor
		        vvx[vvcount] = intersectionsx[i];
		        vvy[vvcount] = intersectionsy[i];
		        vvw[vvcount] = intersectionsw[i];
		        vvcount++;
		    } 
	    } else {
		    vvx[vvcount] = intersectionsx[i];
		    vvy[vvcount] = intersectionsy[i];
		    vvw[vvcount] = intersectionsw[i];
		    vvcount++;
	    }
		
        }

        
        // step 6: calculate final position estimation with given algorithm
        float masses[vvcount];
        for (int i = 0; i < vvcount; i++) {
            if (vvw[i] == 1.0F) {
                masses[i] = 3.0F * vvw[i];
            } else {
                masses[i] = 1.0F;
            }
        }
        center_of_mass(vvcount, vvx, vvy, masses, &((*resx)[ii]), &((*resy)[ii]))    ;    
     }
}

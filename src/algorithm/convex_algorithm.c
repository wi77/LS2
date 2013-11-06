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
**  This file is made only for including in the LS² project
**  and not desired for stand alone usage!
**
********************************************************************/




#ifndef CONVEX_ALGORITHM_C_INCLUDED
#define CONVEX_ALGORITHM_C_INCLUDED 1

#if HAVE_CONFIG_H
#  include "ls2/ls2-config.h"
#endif

#if HAVE_POPT_H
#  include <popt.h>
#endif

#define CONVEX_EPS 1.0f

#include "util/util_circle.c"
#include "util/util_vector.c"

#ifndef CONVEX_DEFAULT_OFFSET
#define CONVEX_DEFAULT_OFFSET ((3.0f - 1.0f) / (3.0f / 50.0f))
#endif

static float convex_offset = CONVEX_DEFAULT_OFFSET;

struct poptOption convex_arguments[] = {
        { "convex-offset", 0, POPT_ARG_FLOAT | POPT_ARGFLAG_SHOW_DEFAULT,
          &convex_offset, 0,
          "increase all radii by this offset", NULL },
        POPT_TABLEEND
};



/*
* This function finds a point in the intersection of
* the no_anchors given open disks with center points
* (vx[i][ii],vy[i][ii]) and radii r[i][ii].
* If no such point exists it the result is (NAN,NAN).
* Runtime is in O(no_anchors^3). Additional used space is in O(1).
*/
static inline void
__attribute__((__always_inline__,__gnu_inline__,__nonnull__,__artificial__))
convex_run(const VECTOR* vx, const VECTOR* vy, const VECTOR *restrict r,
           size_t no_anchors,
           int width __attribute__((__unused__)),
           int height __attribute__((__unused__)),
           VECTOR *restrict resx, VECTOR *restrict resy)
{
	// First check if two disks are disjoint,
	// or if one centerpoint is inside of all disks
	float dist = 0.0f;
	char centerinside = 1;
        VECTOR offset = VECTOR_BROADCAST(&convex_offset);
        VECTOR rr[no_anchors];

	for (size_t i = 0; i < no_anchors; i++) {
		rr[i] = r[i] + offset;
	}

    for (int ii = 0; ii < VECTOR_OPS; ii++) {
	for(size_t i = 0; i < no_anchors; i++){
		centerinside = 1;
		for(size_t j = 0; j < no_anchors; j++){
			if(i == j)
			continue;
			dist = distance_squared_s(vx[i][ii],vy[i][ii],vx[j][ii],vy[j][ii]);
			// If and only if  
			if(dist >= rr[i][ii]*rr[i][ii]+ 2*rr[i][ii]*rr[j][ii]  + rr[j][ii]*rr[j][ii]){
#ifndef NDEBUG
				const int size = 4096;
				char buffer[size];
				int p = 0;
				p += snprintf(buffer, (size_t) size, "Failure for parameters:\n");
				for (size_t k = 0; k < no_anchors; k++) {
					p += snprintf(buffer + p, (size_t) (size - p), "\tA[%zu] = (%f,%f) range = %f\n", k, vx[k][ii], vy[k][ii], rr[k][ii]);
				}
				buffer[(p < size) ? p : size - 1] = '\0';
				fprintf(stderr, buffer);
#endif
				(*resx)[ii] = NAN;
				(*resy)[ii] = NAN;
				goto next;
			}
			// If d(A_i,A_j) >= r_j, then
			// A_i is not in B_j.
			if(dist >= rr[j][ii]*rr[j][ii]){
				centerinside = 0;
			}

		}
		// return the centerpoint if it
		// already is in I
		if(centerinside){
			(*resx)[ii] = vx[i][ii];
			(*resy)[ii] = vy[i][ii];
			goto next;
		}
	}

	// From here on the two conditions before the upper
	// loop should be false. Then also no disk can be completely
	// inside all others.
	// Therefore 2 intersection points of the circles
	// are in all other closed disks iff and only if
	// the intersection of all open disks is not empty.
	
	size_t no_res = 0;
	float resultx[2], resulty[2];
	resultx[0] = NAN;
	resultx[1] = NAN;
	resulty[0] = NAN;
	resulty[1] = NAN;
	float resx1[2], resy1[2];
	resx1[0] = NAN;
	resy1[0] = NAN;
	resx1[1] = NAN;
	resy1[1] = NAN;
	size_t found = 0;
	char insideflag;
	
	
	

	// Loop over each pair of circles.
	for (size_t i = 0; i < no_anchors - 1; i++){
		for(size_t j = i+1; j < no_anchors; j++){
			// find the intersection point/points
			no_res = circle_get_intersection(vx[i][ii],vy[i][ii],vx[j][ii],vy[j][ii],rr[i][ii],rr[j][ii],resx1,resy1);
			if(no_res == 0){
				// happens if one circle is inside of the other
				continue;
			}
			if(no_res == 1){
				// should never happen
				fprintf(stderr,"1 intersection point. Unlikely.\r\n");
				// if only one intersection point is found
				// then go to the code for the first intersection point
				goto int1;
			}
			
			insideflag = 1;
			// check if the second intersection point
			// is inside of all other closed disks
			for(size_t k = 0; k < no_anchors; k++){
				// skip the circles from the outer loop
				if(k == i || k == j)
				continue;
				dist = distance_squared_s(vx[k][ii],vy[k][ii],resx1[1],resy1[1]);
				if(dist > rr[k][ii]*rr[k][ii]){
					insideflag = 0;
					break;
				}
			}
			
			if(insideflag){
				// Arriving here, we have found a point.
				// Write the result to the appropriate index.
				resultx[found] = resx1[1];
				resulty[found] = resy1[1];
				if(found == 1){
					if(distance_squared_s(resultx[0],resulty[0],resultx[1],resulty[1]) > CONVEX_EPS)
					found++;
				}
				else
				found++;
				if(found == 2)
				goto end;
			}
			
			// do the same for the first intersection point
int1: 
			// this label is just here for the case of one intersection point
			// and one circle is inside the other

			insideflag = 1;
			for(size_t k = 0; k < no_anchors; k++){
				// skip the circles from the outer loop
				if(k == i || k == j)
				continue;
				dist = distance_squared_s(vx[k][ii],vy[k][ii],resx1[0],resy1[0]);
				if(dist > rr[k][ii]*rr[k][ii]){
					insideflag = 0;
					break;
				}
			}
			if(insideflag){
				// Arriving here, we have found a point.
				// Write the result to the appropriate index.
				resultx[found] = resx1[0];
				resulty[found] = resy1[0];
				if(found == 1){
					if(distance_squared_s(resultx[0],resulty[0],resultx[1],resulty[1]) > CONVEX_EPS)
					found++;
				}
				else
				found++;
				if(found == 2)
				goto end;
			}
		}
	}
	// Arriving here means no two points
	// on the boundary of the intersection
	// of the open disks have been found.
	// That means the intersection is empty.
#ifndef NDEBUG
	(*resx)[ii] = NAN;
	(*resy)[ii] = NAN;	
	do {
		const int size = 4096;
		char buffer[size];
		int p = 0;
		p += snprintf(buffer, (size_t) size, "Empty intersection for parameters:\n");
		for (size_t k = 0; k < no_anchors; k++) {
			p += snprintf(buffer + p, (size_t) (size - p), "\tA[%zu] = (%f,%f), range = %f\n", k, vx[k][ii], vy[k][ii], rr[k][ii]);
		}
		buffer[(p < size) ? p : size - 1] = '\0';
		fprintf(stderr, buffer);
	} while (0);
#endif
	continue;

	// Here we found two according points.
	// Then we construct a convex combination
	// and return.
end:
	(*resx)[ii] = (float)(0.5f*resultx[0] + 0.5f*resultx[1]);
	(*resy)[ii] = (float)(0.5f*resulty[0] + 0.5f*resulty[1]);
next:
        continue;
    }
	return;
}


// This algorithm finds a discrete point (a grid point) inside
// of the given disks, if it exists.
// It's very slow, because it every time iterates over all possible
// grid points and checks if they are results.
// Can be used if other approaches fail.
static inline void __attribute__((__always_inline__,__gnu_inline__,__nonnull__,__artificial__))
convex_pixel(const VECTOR* vx, const VECTOR* vy, const VECTOR *restrict r,
size_t no_anchors,
int width,
int height,
VECTOR *restrict resx,
VECTOR *restrict resy)
{
	char inside = 1;
	float dist = 0.0f;
    for (int ii = 0; ii < VECTOR_OPS; ii++) {
	for(float i = 0.0f; i < width; i++){
		for(float j = 0.0f; j < height; j++){
			inside = 1;
			for(size_t k = 0; k < no_anchors; k++){
				dist = (i - vx[k][ii])*(i - vx[k][ii]) + (j - vy[k][ii])*(j + vy[k][ii]);
				if(dist >= r[k][ii]*r[k][ii]){
					inside = 0;
					break;
				}
			}
			if(inside == 1){
				(*resx)[ii] = i;
				(*resy)[ii] = j;
				goto next;
			}
		}
	}
next:
        continue;
    }
    return;
}

#endif

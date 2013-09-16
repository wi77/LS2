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

#define CONVEX_EPS 1.0f

#include "util/util_circle.c"
#include "util/util_vector.c"

/*
* This function finds a point in the intersection of
* the no_anchors given open disks with center points
* (vx[i][0],vy[i][0]) and radii r[i][0].
* If no such point exists it the result is (NAN,NAN).
* Runtime is in O(no_anchors^3). Additional used space is in O(1).
*/
static inline void __attribute__((__always_inline__,__gnu_inline__,__nonnull__,__artificial__))
convex_run(const VECTOR* vx, const VECTOR* vy, const VECTOR *restrict r,
size_t no_anchors,
int width __attribute__((__unused__)),
int height __attribute__((__unused__)),
VECTOR *restrict resx,
VECTOR *restrict resy)
{
	// First check if two disks are disjoint,
	// or if one centerpoint is inside of all disks
	double dist = 0.0f;
	char centerinside = 1;
	
	for(size_t i = 0; i < no_anchors; i++){
		centerinside = 1;
		for(size_t j = 0; j < no_anchors; j++){
			if(i == j)
			continue;
			dist = distance_squared_sf(vx[i][0],vy[i][0],vx[j][0],vy[j][0]);
			// If and only if  
			if(dist >= r[i][0]*r[i][0]+ 2*r[i][0]*r[j][0]  + r[j][0]*r[j][0]){
				(*resx)[0] = NAN;
				(*resy)[0] = NAN;
				return;
			}
			// If d(A_i,A_j) >= r_j, then
			// A_i is not in B_j.
			if(dist >= r[j][0]*r[j][0]){
				centerinside = 0;
			}

		}
		// return the centerpoint if it
		// already is in I
		if(centerinside){
			(*resx)[0] = vx[i][0];
			(*resy)[0] = vy[i][0];
			return;
		}
	}

	// From here on the two conditions before the upper
	// loop should be false. Then also no disk can be completely
	// inside all others.
	// Therefore 2 intersection points of the circles
	// are in all other closed disks iff and only if
	// the intersection of all open disks is not empty.
	
	size_t no_res = 0;
	double resultx[2], resulty[2];
	resultx[0] = NAN;
	resultx[1] = NAN;
	resulty[0] = NAN;
	resulty[1] = NAN;
	double resx1[2], resy1[2];
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
			no_res = circle_get_intersectionf(vx[i][0],vy[i][0],vx[j][0],vy[j][0],r[i][0],r[j][0],resx1,resy1);
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
				dist = distance_squared_sf(vx[k][0],vy[k][0],resx1[1],resy1[1]);
				if(dist > r[k][0]*r[k][0]){
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
						if(distance_squared_sf(resultx[0],resulty[0],resultx[1],resulty[1]) > CONVEX_EPS)
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
				dist = distance_squared_sf(vx[k][0],vy[k][0],resx1[0],resy1[0]);
				if(dist > r[k][0]*r[k][0]){
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
						if(distance_squared_sf(resultx[0],resulty[0],resultx[1],resulty[1]) > CONVEX_EPS)
							found++;
					}
					else
						found++;
					if(found == 2)
						goto end;
				}
			}
		}
	fprintf(stderr, "error\r\n");	
	// Arriving here means no two points
	// on the boundary of the intersection
	// of the open disks have been found.
	// That means the intersection is empty.
	(*resx)[0] = NAN;
	(*resy)[0] = NAN;	
	return;

	// Here we found two according points.
	// Then we construct a convex combination
	// and return.
end:
	(*resx)[0] = (float)(0.5f*resultx[0] + 0.5f*resultx[1]);
	(*resy)[0] = (float)(0.5f*resulty[0] + 0.5f*resulty[1]);
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
	for(float i = 0.0f; i < width; i++){
		for(float j = 0.0f; j < height; j++){
			inside = 1;
			for(size_t k = 0; k < no_anchors; k++){
				dist = (i - vx[k][0])*(i - vx[k][0]) + (j - vy[k][0])*(j + vy[k][0]);
				if(dist >= r[k][0]*r[k][0]){
					inside = 0;
					break;
				}
			}
			if(inside == 1){
				(*resx)[0] = (float)i;
				(*resy)[0] = (float)j;
				return;
			}
		}
	}
	return;
}

#endif

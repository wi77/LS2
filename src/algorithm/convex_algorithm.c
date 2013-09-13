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
#define CONVEX_EPS 0.01f

#include "util/util_circle.c"
#include "util/util_vector.c"

/*
* This function finds a point, if existing, in the intersection of
* the no_anchors given open disks with center points
* (vx[i][0],vy[i][0]) and radii r[i][0].
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
	// first check if two disks are disjoint,
	// or if one centerpoint is inside of all disks
	double dist = 0.0f;
	char onecircle = 1;
	
	for(size_t i = 0; i < no_anchors; i++){
		onecircle = 1;
		for(size_t j = 0; j < no_anchors; j++){
			if(i == j)
				continue;
			dist = distance_squared_sf(vx[i][0],vy[i][0],vx[j][0],vy[j][0]);
			// If and only if  
			if(dist >= r[i][0]*r[i][0]+ 2*r[i][0]*r[j][0]  + r[j][0])*r[j][0]){
				(*resx)[0] = NAN;
				(*resy)[0] = NAN;
				return;
			}
			// If d(A_i,A_j) >= r_j, then
			// A_i is not in B_j.
			if(dist >= r[j][0]*r[j][0]){
				onecircle = 0;
			}

		}
		// return the centerpoint if it
		// already is in I
		if(onecircle){
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
	char nonefound = 1;
	
	
	

	// Loop over each pair of circles.
	for (size_t i = 0; i < no_anchors; i++){
		for(size_t j = 0; j < no_anchors; j++){
			if(i == j)
				continue;
			// find the intersection point/points
			no_res = circle_get_intersectionf(vx[i][0],vy[i][0],r[i][0],vx[j][0],vy[j][0],r[j][0],resx1,resy1);
			if(no_res == 1){
				// if only one intersection point is found
				// then go to the code for the first intersection point
				goto int1;
			}
			
			// check if the second intersection point
			// is inside of all other closed disks
			for(size_t k = 0; k < no_anchors; k++){
				// skip the circles from the outer loop
				if(k == i || k == j)
					continue;
				dist = distance_squared_sf(vx[k][0],vy[k][0],resx1[1],resy1[1]);
				if(dist > r[k][0]*r[k][0] + CONVEX_EPS){
					goto int1;
				}
			}
			
			// Arriving here, we have found a point
			if(nonefound){
				// if its the first one, write it to index 0
				resultx[0] = resx1[1];
				resulty[0] = resy1[1];
				nonefound = 0;
			}
			else{
				// if its the second one, write it to index 1
				// and go to the end
				resultx[1] = resx1[1];
				resulty[1] = resy1[1];
				goto end;
			}

			// do the same for the first intersection point
int1: 
			for(size_t k = 0; k < no_anchors; k++){
				// skip the circles from the outer loop
				if(k == i || k == j)
					continue;
				dist = distance_squared_sf(vx[k][0],vy[k][0],resx1[0],resy1[0]);
				if(dist > r[k][0]*r[k][0] + CONVEX_EPS){
					goto cont;
				}
			}
			// Arriving here, we have found a point
			if(nonefound){
				// if its the first one, write it to index 0
				resultx[0] = resx1[0];
				resulty[0] = resy1[0];
			}
			else{
				// if its the second one, write it to index 1
				// and go to the end
				resultx[1] = resx1[0];
				resulty[1] = resy1[0];
				goto end;
			}
		}
cont: ;
	}
	// Arriving here means no two points
	// on the boundary of the intersection
	// of the open disks have been found.
	// That means the intersection is empty.
	(*resx)[0] = NAN;
	(*resy)[0] = NAN;	
	fprintf(convex_file,"((x-%f)^2 + (y - %f)^2 < %f) && ((x-%f)^2 + (y - %f)^2 < %f) && ((x-%f)^2 + (y - %f)^2 < %f) && ((x-%f)^2 + (y - %f)^2 < %f)", vx[0][0], vy[0][0], r[0][0], vx[1][0], vy[1][0], r[1][0], vx[2][0], vy[2][0], r[2][0], vx[3][0], vy[3][0], r[3][0]);
	return;

	// Here we found two according points.
	// Then we construct a convex combination
	// and return.
end:
	(*resx)[0] = (float)(0.5f*resultx[0] + 0.5f*resultx[1]);
	(*resy)[0] = (float)(0.5f*resulty[0] + 0.5f*resulty[1]);
	return;
}


// pixel iterator algorithm
static inline void __attribute__((__always_inline__,__gnu_inline__,__nonnull__,__artificial__))
pixel_run(const VECTOR* vx, const VECTOR* vy, const VECTOR *restrict r,
size_t no_anchors,
int width __attribute__((__unused__)),
int height __attribute__((__unused__)),
VECTOR *restrict resx,
VECTOR *restrict resy)
{
	for(int i = 0; i < 1000; i ++){
		for(int j = 0; j < 1000; j++){
			char inside = 1;
			for(size_t k = 0; k < no_anchors; k++){
				if(distance_s(i,j,vx[k][0],vy[k][0]) > r[k][0]){
					inside = 0;
				}
			}
			if(inside){
				(*resx)[0] = (float)i;
				(*resy)[0] = (float)j;
			}
		}
	}
	fprintf(convex_file,"((x-%f)^2 + (y - %f)^2 < %f) && ((x-%f)^2 + (y - %f)^2 < %f) && ((x-%f)^2 + (y - %f)^2 < %f) && ((x-%f)^2 + (y - %f)^2 < %f)", vx[0][0], vy[0][0], r[0][0]*r[0][0], vx[1][0], vy[1][0], r[1][0]*r[1][0], vx[2][0], vy[2][0], r[2][0]*r[2][0], vx[3][0], vy[3][0], f[3][0]*r[3][0]);
	return;
}

#endif

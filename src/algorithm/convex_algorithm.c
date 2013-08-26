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

#include "util/util_circle.c"
#include "util/util_vector.c"
/*
 * This function finds a point, if existing, in the intersection of
 * the no_anchors given open disks with center points
 * (vx,vy) and radii r.
 * Runtime is in O(no_anchors^3). Used space is in O(1).
 */	
static inline size_t __attribute__((__always_inline__,__gnu_inline__,__nonnull__,__artificial__))
convex_run(const VECTOR* vx, const VECTOR* vy, const VECTOR *restrict r,
         size_t no_anchors,
         int width __attribute__((__unused__)),
         int height __attribute__((__unused__)),
         VECTOR *restrict resx,
         VECTOR *restrict resy)
{
	size_t no_res = 0;
	size_t found = 0;
	float resultx[2], resulty[2];
	float resx1[2], resy1[2];
	resx1[0] = 0.0f;
	resy1[0] = 0.0f;
	resx1[1] = 0.0f;
	resy1[1] = 0.0f;
	float dist = 0.0f;

	// For each circle find the intersection points with the other circles.
 for (size_t i = 0; i < no_anchors && found < 2; i++){
   for(size_t j = i+1; j < no_anchors && found < 2; j++){
		       no_res = circle_get_intersection(vx[i][0],vy[i][0],r[i][0],vx[j][0],vy[j][0],r[j][0],resx1,resy1);
		       // if no intersection is found
		       // check why
		       if(no_res == 0){
			       dist = distance_s(vx[i][0],vy[i][0],vx[j][0],vy[j][0]);
			       // case 1:
			       // if circles don't intersect
			       // there is no solution
			       if(dist > r[i][0] + r[j][0]){
				       return 0;
			       }
			       // case 2:
			       // if one circle is inside of the other
			       // then the the outer circle can be removed
			       // because it has not impact to solution
			       // TODO begin
			       // can be checked, if the circle is inside
			       // of all other circles, then its center is a solution
			       // TODO end
			       // for now it will just be continued
			       // to the next circle
			       if(dist < fabs(r[i][0] - r[j][0])){
				       continue;
			       }
			       // case 3:
			       // if circles are coincident
			       // TODO begin
			       // this should not happen
			       // this case implies either case 2
			       // or the circles are identical and one can be omitted
			       // TODO end
			       // will be treated like case 2 for now
			       if(dist == 0.0f){
				       continue;
			       }
		       }
			   // if two circles only touch
			   // there is no intersection
			   // of the open disks
			   // hence no solution exists
			   if(no_res == 1)
					return 0;
		       // check if found intersection point
		       // is inside of all other closed disks
			   size_t k = 0;
		       for(k = 0; k < no_anchors; k++){
			       if(k == i || k == j)
				       continue;
			       dist = distance_s(vx[k][0],vy[k][0],resx1[0],resy1[0]);
			       if(dist > r[k][0]){
				       break;
				   }
		       }
			   
			   // if the loop ran to the end
			   // then there was a point inside all
			   // other disks
			   if(k == no_anchors){
					resultx[found] = resx1[0];
					resulty[found] = resy1[0];
					found++;
					if(found >= 2){
						goto end;
					}
			   }
			   
			   // run again for second intersection point
		       for(k = 0; k < no_anchors; k++){
			       if(k == i || k == j)
				       continue;
			       dist = distance_s(vx[k][0],vy[k][0],resx1[1],resy1[1]);
			       if(dist > r[k][0]){
				       break;
				   }
		       }
			   
			   if(k == no_anchors){
					resultx[found] = resx1[1];
					resulty[found] = resy1[1];
					found++;
					if(found >= 2){
						goto end;
					}
		       }
	       }
       }
	   end:
       (*resx)[0] = (float)(resultx[0] + 0.5f*(resultx[1] - resultx[0]));
       (*resy)[0] = (float)(resulty[0] + 0.5f*(resulty[1] - resulty[0]));
       //TODO: sanity check
	   //possibly unnecessary
	   
       //TODO: write results!
       return 1;
}
#endif

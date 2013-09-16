#ifndef TEST_ALGORITHM_C_INCLUDED
#define TEST_ALGORITHM_C_INCLUDED 1

#include "algorithm/convex_algorithm.c"

extern FILE * convex_f;

static inline void
__attribute__((__always_inline__,__gnu_inline__,__artificial__,__nonnull__))
test_run(const VECTOR* vx, const VECTOR* vy, const VECTOR *restrict r,
size_t no_anchors,
int width,
int height,
VECTOR *restrict resx, VECTOR *restrict resy)
{
	VECTOR sx, sy;
	convex_run(vx, vy, r, no_anchors, width, height, &sx, &sy);

        for(size_t i = 0; i<no_anchors; i++){
                fprintf(convex_f,"%f\t%f\t%f\t", vx[i][0],vy[i][0],r[i][0]);
        }
        fprintf(convex_f,"%f\t%f\r\n", sx[0],sy[0]);
        resx[0][0] = 500.0f;
        resy[0][0] = 500.0f;

	
	double mydist;
	if(!isnan(sx[0]) && !isnan(sy[0])){
		for(size_t i = 0; i < no_anchors; i++){
			mydist = distance_squared_sf(vx[i][0], vy[i][0], sx[0], sy[0]);
			if(mydist >= r[i][0]*r[i][0]){
				fprintf(convex_f, "found a bad point\r\n");
			}
		}
	}
	else{
		convex_pixel(vx,vy,r,no_anchors,width,height,&sx,&sy);
		if(!isnan(sx[0]))
			fprintf(convex_f, "faulty -nan returned\r\n");
	}
}
#endif

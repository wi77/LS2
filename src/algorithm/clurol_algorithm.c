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

/* @algorithm_name: CLUROL */

/*******************************************************************
 ***
 ***   CLUROL;
 *** 
 *******************************************************************/

#ifndef CLUROL_ALGORITHM_C_INCLUDED
#define CLUROL_ALGORITHM_C_INCLUDED 1

#if HAVE_CONFIG_H
#  include "ls2/ls2-config.h"
#endif

#if HAVE_POPT_H
#  include <popt.h>
#endif
#include <math.h>
#include <string.h>
#include "algorithm/nllsq_algorithm.c"

#define min(x,y) (x<=y?x:y)

typedef struct Point2dC{
    double x;
    double y;
    int cluster;
} Point2dC;

static inline void
__attribute__((__always_inline__,__gnu_inline__,__nonnull__,__artificial__))
center_of_mass_eqwp(const Point2dC *pts,const int count,
              Point2dC *restrict ret)
{
    double x = 0, y = 0;
    for (int i=0; i < count; i++) {
        x += pts[i].x;
        y += pts[i].y;
    }
    ret->x = x/(double)count;
    ret->y = y/(double)count;
}

    /**
     * Pairwise distance tuple class.
     */
    typedef struct PairwiseDistanceTuple {
        /** The first point */
        Point2dC *x;
        /** The second point */
        Point2dC *y;
        /** The distance between both points */
        double d;
    } PairwiseDistanceTuple;

        /**
         * Creates a new instance of {@code PairwiseDistanceTuple}.
         * <p>
         * The distance between {@code x} and {@code y} is automatically set.
         *
         * @param x The first point.
         * @param y The second point.
         */
    static inline void 
    __attribute__((__always_inline__,__gnu_inline__,__nonnull__,__artificial__))
    PairwiseDistanceTuple_new(Point2dC *x, Point2dC *y, PairwiseDistanceTuple *this) {
         this->x = x;
         this->y = y;
         this->d = distance_sf(x->x,x->y,y->x,y->y);
    }

    static inline int 
    __attribute__((__always_inline__,__gnu_inline__,__nonnull__,__artificial__))
    PairwiseDistanceTuple_compareTo(const void *o_, const void *this_) {
        PairwiseDistanceTuple *o = (PairwiseDistanceTuple*) o_;
        PairwiseDistanceTuple *this = (PairwiseDistanceTuple*) this_;
        return o->d < this->d ? -1 : o->d > this->d;            
    }

    static inline int 
    __attribute__((__always_inline__,__gnu_inline__,__nonnull__,__artificial__))
    clusterSize(Point2dC* points, int plength, int cluster){
        int size = 0;
        for (int i=0; i < plength; i++) {
            if (points[i].cluster == cluster)
              size++;
        }
        return size;
    }

    static inline void 
    __attribute__((__always_inline__,__gnu_inline__,__nonnull__,__artificial__))
    cluster2array(Point2dC* points, int plength, int cluster, Point2dC *ret){
        int count = 0;
        for (int i=0; i < plength; i++) {
            if (points[i].cluster == cluster) {
                ret[count]=points[i];
                count++;
            }
        }
    }

    static inline int 
    __attribute__((__always_inline__,__gnu_inline__,__nonnull__,__artificial__))
    findMaxCluster(PairwiseDistanceTuple *D, int dlength, double dth, Point2dC* points, int plength) {
        // Build initial cluster
        int maxCluster=1;
        D[0].x->cluster = maxCluster;
        D[0].y->cluster = maxCluster;

        for (int i = 1; i < dlength; i++) {
            Point2dC *x = D[i].x;
            Point2dC *y = D[i].y;
            double dist = D[i].d;

            // test first condition, paper lines 5..7
            if (x->cluster == 0 && y->cluster == 0) {
                // x and y not in any cluster, add to a new cluster
                maxCluster++;
                x->cluster= maxCluster;
                y->cluster= maxCluster;
                continue;
            }

            // test second condition, paper lines 8..12
            if (x->cluster != 0 && y->cluster == 0) {
                // x in cluster and y does not belong to any cluster
                if (dist <= dth) {
                    y->cluster = x->cluster;
                }
                // continue with for-loop
                continue;
            }

            // test third condition, paper lines 13..17
            if (x->cluster == 0 && y->cluster != 0) {
                // y in cluster and x does not belong to any cluster
                if (dist <= dth) {
                    x->cluster = y->cluster;
                }
                // continue with for-loop
                continue;
            }

            // test fourth condition, paper lines 18..26
            if (x->cluster != 0 && y->cluster != 0 && x->cluster != y->cluster) {
                // need to check if Cx and Cy can be merged
                int XCount = clusterSize(points,plength,x->cluster);
                int YCount = clusterSize(points,plength,y->cluster);
                Point2dC XArray[XCount];
                Point2dC YArray[YCount];
                cluster2array(points,plength,x->cluster,XArray);
                cluster2array(points,plength,y->cluster,YArray);
                Point2dC centroidCx, centroidCy;
                
                center_of_mass_eqwp(XArray,XCount,&centroidCx);
                center_of_mass_eqwp(YArray,YCount,&centroidCy);
                if (distance_sf(centroidCx.x,centroidCx.y,centroidCy.x,centroidCy.y) <= dth) {
                    // Merge both sets, remove one set from list
                    int yc = y->cluster;
                    for (int ic=0; ic < plength; ic++){
                        if (points[ic].cluster == yc) {
                            points[ic].cluster = x->cluster;
                        }
                    }
                }
            }
        }

        // return cluster with maximum cardinality as result
        int maxC = 1;
        int maxClustersize = clusterSize(points,plength,1);
        int x;
        for (int i = 2; i < maxCluster; i++) {
            if ((x = clusterSize(points,plength,i)) > maxClustersize) {
                maxC = i;
                maxClustersize = x;
            }
        }
        return maxC;
    }



static inline void
__attribute__((__always_inline__,__gnu_inline__,__artificial__,__nonnull__))
clurol_run(const VECTOR* vx, const VECTOR* vy, const VECTOR *restrict r,
              size_t no_anchors, int width, int height,
              VECTOR *restrict resx, VECTOR *restrict resy)
{
if (width==height){};

for (int ii = 0; ii < VECTOR_OPS; ii++) {
        // step 1: initialize variables

        int n = (int)no_anchors;
        int ancStatusTaken = 0;
        int ancStatus[n];
        memset(ancStatus, 0 , no_anchors * sizeof(int));
        
        // step 2: calculate circle intersections
        int bino = binom(n, 2);
        double intersections_x[bino*2];
        double intersections_y[bino*2];
        int num = 0;
        int is = 0;
        for (int i = 0; i < n-1; i++) {
            for (int j = i+1; j < n; j++) {
                // Berechne Schnittpunkte mit aktueller Permutation
                is += (int)circle_get_intersectionf(vx[i][ii],vy[i][ii],vx[j][ii],vy[j][ii],r[i][ii],r[j][ii],&intersections_x[is],&intersections_y[is]);
            }
        }

        // step 3: copy intersections into one array
        Point2dC points[is];
        for (int i = 0; i < is; i++) {
            points[i] = (Point2dC) {intersections_x[i], intersections_y[i], 0};
        }

        if (is < 2) {
            // no sense to do CluRoL, return NLLS here
            nllsq_run(vx, vy, r, no_anchors, width, height, resx, resy);
            continue;
        }

        // step 4: build all pairwise distance tuples
        num = 0;
        bino = binom(is, 2);
        PairwiseDistanceTuple D[bino];
        for (int i = 0; i < is-1; i++) {
            for (int j = i+1; j < is; j++) {
                PairwiseDistanceTuple_new(&points[i], &points[j], &D[num]);
                num++;
            }
        }

        // step 5: sort D in ascending order of the pairwise distances
        qsort(D,(size_t)bino,sizeof(PairwiseDistanceTuple),PairwiseDistanceTuple_compareTo);

        // step 6: calculate distance threshold as n-th percentile tuple’s
        //         pairwise distance value
        int alpha = binom((int)(ceil((double)n/2.0) + 2), 2);
        int beta = 2 * binom(n, 2);
        double nth = (binom(alpha, 2)/(double)binom(beta, 2));
        int nthPercentile = (int)lround((nth * (double)bino) + 0.5);
        nthPercentile = min(nthPercentile, bino);
        double dth = D[nthPercentile-1].d;

        // step 7: call findMaxCluster subroutine
        int cMax = findMaxCluster(D, bino, dth , points, is);

        // step 8: determine anchors for Minimum Squared Error (MSE) method
        double dMax = 0.05;
        double boundConst = (1 + dMax) * (1 + dMax);
        for (int cm=0;cm<is;cm++) {
            if (points[cm].cluster==cMax){
                for (int i = 0; i < n; i++) {
                    double ub = r[i][ii] * boundConst;
                    double lb = r[i][ii] / boundConst;
                    double dist = distance_sf(points[cm].x,points[cm].y,vx[i][ii],vy[i][ii]);
                    if (dist <= ub && dist >= lb) {
                        if (!ancStatus[i]) {
                            ancStatusTaken++;
                        }
                        ancStatus[i] = 1;
                    }
                }
            }
        }

        // step 9: copy anchors and ranges into new array and localize
        num = 0;
        VECTOR anchorsTaken_x[ancStatusTaken],anchorsTaken_y[ancStatusTaken];
        VECTOR rangesTaken[ancStatusTaken];
        VECTOR retx, rety;
        for (int i = 0; i < n; i++) {
            if (ancStatus[i]) {
                anchorsTaken_x[num][ii] = vx[i][ii];
                anchorsTaken_y[num][ii] = vy[i][ii];
                rangesTaken[num][ii] = r[i][ii];
                num++;
            }
        }
        if (num<3) {
            (*resx)[ii]=NAN;
            (*resy)[ii]=NAN;
        } else {
            nllsq_run(anchorsTaken_x,anchorsTaken_y,rangesTaken,(size_t)num,width,height,&retx,&rety);
            (*resx)[ii]=retx[ii];
            (*resy)[ii]=rety[ii];
        }
    }
}


#endif

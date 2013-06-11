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

/* @algorithm_name: ICLA */

/*******************************************************************
 ***
 ***   ICLA;
 *** 
 *******************************************************************/

#ifndef ICLA_ALGORITHM_C_INCLUDED
#define ICLA_ALGORITHM_C_INCLUDED 1

#if HAVE_CONFIG_H
#  include "ls2/ls2-config.h"
#endif

#if HAVE_POPT_H
#  include <popt.h>
#endif

//#include "util/util_math.c"

typedef struct Point2d{
    double x;
    double y;
} Point2d;

typedef struct IcmPoint {
    int id;
    int weight;
    int merged;
    int movingDirection;
    Point2d intersection;
    Point2d currentLocation;
    struct IcmPoint *nodesInRange[90];
    int nin_fill ;
    double attractingBoundary;
    struct IcmPoint *mergeList[90];
    int ml_fill;
} IcmPoint;

static int MAPPING_TABLE[3][3] = {{6, 7, 8},{5, 9, 1},{4, 3, 2}};
static Point2d MOVING_DIRECTIONS[10] = { {0,0}, {0, 1}, {1, 1}, {1, 0}, {1, -1}, {0, -1}, {-1, -1}, {-1, 0}, {-1, 1}, {0, 0} };

static inline void icmp_new(double x, double y, int id, IcmPoint *this) {
            this->id = id;
            this->intersection = (Point2d){x,y};
            this->ml_fill = 0;
            this->currentLocation = (Point2d){x, y};
            this->merged = 0;
            this->weight = 1;
            this->nin_fill = 0;
        }

static inline int icmp_getAttractingForceDirection(IcmPoint *p, IcmPoint *this) {
            double d = distance_sf(p->currentLocation.x,p->currentLocation.y,this->currentLocation.x,this->currentLocation.y);
            long int x = lround((double)(p->currentLocation.x - this->currentLocation.x) / d);
            long int y = lround((double)(p->currentLocation.y - this->currentLocation.y) / d);
            return MAPPING_TABLE[x+1][y+1];
        }

static inline void icmp_move(double step, IcmPoint *this) {
            Point2d dir = MOVING_DIRECTIONS[this->movingDirection];
            double dirLength = sqrt(dir.x * dir.x + dir.y * dir.y);
            if (dirLength == 0.0) return; // no movement
            this->currentLocation.x = this->currentLocation.x + (step / dirLength) * dir.x;
            this->currentLocation.y = this->currentLocation.y + (step / dirLength) * dir.y;
        }

static inline void icmp_merge(IcmPoint *p, IcmPoint *this) {
            p->merged = 1;
            this->mergeList[this->ml_fill] = p;
            this->ml_fill++;
            for (int i = 0; i < p->ml_fill; i++) {
                this->mergeList[this->ml_fill] = p->mergeList[i];
                this->ml_fill++;
            }
            this->weight += p->weight;
            this->attractingBoundary = fmax(this->attractingBoundary, p->attractingBoundary);
        }

// computes initial radius of each points' attracting boundary as longest
// distance from other points. Uses brute force method (may be optimized)!
static inline void computeInitialAttractingBoundary(IcmPoint *points, int count) {
        for (int i = 0; i < count; i++) {
            double lDist = 0;
            for (int j = 0; j < count; j++) {
                if (i != j) {
                    double d = distance_sf(points[i].intersection.x,points[i].intersection.y,points[j].intersection.x,points[j].intersection.y);
                    if (d > lDist) {
                        lDist = d;
                    }
                }
            }
            points[i].attractingBoundary = lDist;
        }
    }

// compute moving direction of each point
static inline void computeMovingDirection(IcmPoint *points, int count) {
        for (int i = 0; i < count; i++) {
            if (points[i].merged) continue;
            int forceVector[9]; // initial to 0's
            memset(forceVector,0,sizeof(int)*9);
            for (int j = 0; j < count; j++) {
                if (i != j && !points[j].merged && distance_sf(points[i].currentLocation.x,points[i].currentLocation.y,points[j].currentLocation.x,points[j].currentLocation.y) <= points[i].attractingBoundary) {
                    int idx = icmp_getAttractingForceDirection(&points[j],&points[i]);
                    forceVector[idx-1] += points[j].weight;
                }
            }
            int maxWeightIdx = 0;
            for (int j = 1; j < 9; j++) {
                if (forceVector[j] > forceVector[maxWeightIdx]) {
                    maxWeightIdx = j;
                }
            }
            points[i].movingDirection = maxWeightIdx + 1;
        }
    }

    static inline int getNodesInRange(IcmPoint *points, int self, int pcount, IcmPoint **result) {
        int count = 0;
        for (int i = 0; i < pcount; i++) {
            if (points[i].merged) continue;
            if (i != self && distance_sf(points[i].currentLocation.x,points[i].currentLocation.y,points[self].currentLocation.x,points[self].currentLocation.y) <= points[self].attractingBoundary) {
                count++;
            }
        }
        if (count > 0) {
            count = 0;
            for (int i = 0; i < pcount; i++) {
                if (points[i].merged) continue;
                if (i != self && distance_sf(points[i].currentLocation.x,points[i].currentLocation.y, points[self].currentLocation.x,points[self].currentLocation.y) <= points[self].attractingBoundary) {
                    result[count] = &points[i];
                    count++;
                }
            }
        }
        return count;
    }

static inline void
__attribute__((__always_inline__,__gnu_inline__,__artificial__,__nonnull__))
icla_run(const VECTOR* vx, const VECTOR* vy, const VECTOR *restrict r,
              size_t no_anchors, int width, int height,
              VECTOR *restrict resx, VECTOR *restrict resy)
{
    static const double alpha = 1.5;
    static const double moveStep = 25;

    if (width==height){};
    for (int ii = 0; ii < VECTOR_OPS; ii++) {
        // step 1: calculate circle intersections
        int n = (int)no_anchors;
        int k = 2;
        int p[k];
        int bino = binom(n, k);
        double intersections_x[bino*2];
        double intersections_y[bino*2];
        int icount = 0;

        // initialisation for calculating k-permutations
        for (int i = 0; i < k; i++) {
            p[i] = i;
        }

        // build all k-permutations
        for (int i = 0; i < bino; i++) {
            int is;
            is = (int)circle_get_intersectionf(vx[p[0]][ii],vy[p[0]][ii],vx[p[1]][ii],vy[p[1]][ii],r[p[0]][ii],r[p[1]][ii],&intersections_x[icount],&intersections_y[icount]);
            icount += is;
            // build next permutation
            if (i == bino - 1) break;
            int j = k - 1;
            while (j >= 0) {
                if (!incCounter(p, j, n, k)) break;
                j--;
            }
            for (int l = j+1; l < k; l++) {
                p[l] = p[l-1] + 1;
            }
        }

        if (icount == 0) continue;
        IcmPoint points[icount];
        for (int i = 0; i < icount; i++) {
                icmp_new(intersections_x[i], intersections_y[i], i, &points[i]);
        }

        // step 2: adapt iterative clustering model (ICM)
        int iterate = 1;

        // ICM step 1: define initiation range
        computeInitialAttractingBoundary(points,icount);

        while (iterate) {
            // ICM step 2: determine moving direction
            computeMovingDirection(points,icount);
            // ICM step 3: move all points according to current direction on step forward

            for (int i = 0; i < icount; i++) {
                if (points[i].merged) continue;
                points[i].nin_fill = getNodesInRange(points, i,icount, points[i].nodesInRange);
            }
            for (int i = 0; i < icount; i++) {
                if (points[i].merged || points[i].attractingBoundary == 0) continue;
                icmp_move(moveStep,&points[i]);
            }
 
            // ICM step 4: if merging condition is true, merge points
            for (int i = 0; i < icount; i++) {
                if (points[i].merged) continue;
                if (points[i].nin_fill == 0) continue;
                // if only one node in range and not in this nodes' attracting
                // boundary => kick this node from list, set attracting boundary
                // to zero to stop from moving
                if (points[i].nin_fill == 1) {
                    double d = distance_sf(points[i].currentLocation.x,points[i].currentLocation.y,points[i].nodesInRange[0]->currentLocation.x,points[i].nodesInRange[0]->currentLocation.y);
                    if (d > points[i].nodesInRange[0]->attractingBoundary) {
                        // kick point
                        points[i].nin_fill = 0;
                        points[i].attractingBoundary = 0;
                    } else {
                        // merge points
                        if (!points[i].nodesInRange[0]->merged) {
                            icmp_merge(points[i].nodesInRange[0], &points[i]);
                        } else {
                            // merged with other point before we could merge
                            points[i].nin_fill = 0;
                            points[i].attractingBoundary = 0;
                        }
                    }
                } else {
                    // test if points can be merged
                    double dShort = DBL_MAX;
                    IcmPoint *pShort = NULL;
                    for (int j = 0; j < points[i].nin_fill; j++) {
                        if (points[i].nodesInRange[j]->merged) continue;
                        double d = distance_sf(points[i].currentLocation.x,points[i].currentLocation.y,points[i].nodesInRange[j]->currentLocation.x,points[i].nodesInRange[j]->currentLocation.y);
                        if (d <= moveStep*sqrtf(2.0)) {
                             icmp_merge(points[i].nodesInRange[j],&points[i]);
                        }
                        if (d < dShort) {
                            dShort = d;
                            pShort = points[i].nodesInRange[j];
                        }
                    }
                    // merge points when distance between them is the
                    // shortest in both points attracting boundary
                    if (pShort != NULL && !pShort->merged && pShort->nin_fill != 0) {
                        double dShort2 = DBL_MAX;
                        IcmPoint *pShort2 = NULL;
                        for (int l = 0; l < pShort->nin_fill; l++) {
                            if (pShort->nodesInRange[l]->merged) continue;
                            double d = distance_sf(pShort->currentLocation.x,pShort->currentLocation.y,pShort->nodesInRange[l]->currentLocation.x,pShort->nodesInRange[l]->currentLocation.y);
                            if (d < dShort2) {
                                dShort2 = d;
                                pShort2 = pShort->nodesInRange[l];
                            }
                        }
                        if (pShort2 != NULL && pShort2 == &points[i]) {
                            //points[i].merge(pShort);
                        }
                    }
                }
            }

            // ICM step 6: Update ranges of points whose range radii are not zero
            for (int i = 0; i < icount; i++) {
                double dMax = 0;
                double dMin = DBL_MAX;
                if (points[i].merged) continue;
                if (points[i].attractingBoundary == 0) continue;
                if (points[i].nin_fill != 0) {
                    for (int j = 0; j < points[i].nin_fill; j++) {
                        if (points[i].nodesInRange[j]->merged) continue;
                        double d = distance_sf(points[i].currentLocation.x,points[i].currentLocation.y,points[i].nodesInRange[j]->currentLocation.x,points[i].nodesInRange[j]->currentLocation.y);
                        dMin = (d < dMin) ? d : dMin;
                        dMax = (d > dMax) ? d : dMax;
                    }
                }
                points[i].attractingBoundary = dMin == DBL_MAX ? 0.0 : dMin + ((dMax - dMin) / alpha);
            }

            // ICM step 5: check if we can terminate
            iterate = 0;
            for (int i = 0; i < icount; i++) {
                if (points[i].merged) continue;
                if (points[i].attractingBoundary != 0) {
                    iterate = 1;
                    break;
                }
            }
        }

        // step 3: return centroid of selected intersection points
        int maxCluster = -1;
        for (int i = 0; i < icount; i++) {
            if (!points[i].merged) {
                if (maxCluster == -1) {
                    maxCluster = i;
                } else {
                    if (points[i].ml_fill > points[maxCluster].ml_fill) {
                        maxCluster = i;
                    }
                }
            }
        }
        int lcount = points[maxCluster].ml_fill+1;
        float pCenterOfMass_x[lcount];
        float pCenterOfMass_y[lcount];
        float mass[lcount];
        pCenterOfMass_x[0] = (float)points[maxCluster].intersection.x;
        pCenterOfMass_y[0] = (float)points[maxCluster].intersection.y;
        mass[0] = 1.0f;
        for (int i = 1; i < lcount; i++) {
            pCenterOfMass_x[i] = (float)points[maxCluster].mergeList[i-1]->intersection.x;
            pCenterOfMass_y[i] = (float)points[maxCluster].mergeList[i-1]->intersection.y;
            mass[i]=1.0f;
        }
        center_of_mass(lcount, pCenterOfMass_x, pCenterOfMass_y, mass, &((*resx)[ii]), &((*resy)[ii]));
    }
}


#endif

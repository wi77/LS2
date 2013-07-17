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

/* @algorithm_name: Maximum Likelihood Estimator (Gauss) */

/*******************************************************************
 ***
 ***   Maximum Likelihood Estimator (Gauss);
 ***
 *******************************************************************/

#ifndef MLE_GAUSS_ALGORITHM_C_INCLUDED
#define MLE_GAUSS_ALGORITHM_C_INCLUDED 1

#if HAVE_CONFIG_H
#  include "ls2/ls2-config.h"
#endif

#if HAVE_POPT_H
#  include <popt.h>
#endif

#include <assert.h>
#include <gsl/gsl_multimin.h>

#include "algorithm/llsq_algorithm.c"

#ifndef MLE_GAUSS_DEFAULT_MEAN
#define MLE_GAUSS_DEFAULT_MEAN 50.0
#endif

#ifndef MLE_GAUSS_DEFAULT_DEVIATION
#define MLE_GAUSS_DEFAULT_DEVIATION 20.0
#endif

#ifndef MLE_GAUSS_DEFAULT_EPSILON
#define MLE_GAUSS_DEFAULT_EPSILON 1e-3
#endif

#ifndef MLE_GAUSS_DEFAULT_ITERATIONS
#define MLE_GAUSS_DEFAULT_ITERATIONS 100
#endif

static double mle_gauss_mean       = MLE_GAUSS_DEFAULT_MEAN;
static double mle_gauss_deviation  = MLE_GAUSS_DEFAULT_DEVIATION;
static double mle_gauss_epsilon    = MLE_GAUSS_DEFAULT_EPSILON;
static    int mle_gauss_iterations = MLE_GAUSS_DEFAULT_ITERATIONS;


struct poptOption mle_gauss_arguments[] = {
        { "mle-gauss-deviation", 0, POPT_ARG_DOUBLE | POPT_ARGFLAG_SHOW_DEFAULT,
          &mle_gauss_deviation, 0,
          "deviation of the gauss distribution", NULL },
        { "mle-gauss-mean", 0, POPT_ARG_DOUBLE | POPT_ARGFLAG_SHOW_DEFAULT,
          &mle_gauss_mean, 0,
          "mean of the gauss distribution", NULL },
        { "mle-gauss-epsilon", 0, POPT_ARG_DOUBLE | POPT_ARGFLAG_SHOW_DEFAULT,
          &mle_gauss_epsilon, 0,
          "maximum size of the simplex for termination", NULL },
        { "mle-gauss-iterations", 0, POPT_ARG_INT | POPT_ARGFLAG_SHOW_DEFAULT,
          &mle_gauss_iterations, 0,
          "maximum number of iterations before termination", NULL },
        POPT_TABLEEND
};


struct mle_gauss_point2d {
    double x;
    double y;
};

/* This structure holds the parameters to the likelihood function. */
struct mle_gauss_params {
    struct mle_gauss_point2d *anchors;
    double *ranges;
    size_t no_anchors;
};


static double
__attribute__((__nonnull__))
mle_gauss_likelihood_function(const gsl_vector *X, void *restrict params)
{
    double result = 0.0;
    const struct mle_gauss_params *const p = params;
    const double mean2 = mle_gauss_mean * mle_gauss_mean;
    const double variance = mle_gauss_deviation * mle_gauss_deviation;
    const double thetaX = gsl_vector_get(X, 0);
    const double thetaY = gsl_vector_get(X, 1);

    for (size_t j = 0; j < p->no_anchors; j++) {
        const double d2 = (thetaX - p->anchors[j].x) * (thetaX - p->anchors[j].x) +
                          (thetaY - p->anchors[j].y) * (thetaY - p->anchors[j].y);
        const double d = sqrt(d2);
        const double z = p->ranges[j];
        double logProbability = z * z ;
        logProbability += (-2 * mle_gauss_mean - 2 * d) * z;
        logProbability += variance * log(M_PI);
        logProbability += mean2;
        logProbability += 2 * d * mle_gauss_mean;
        logProbability += 2 * variance * log(mle_gauss_deviation); 
        logProbability += log(2) * mle_gauss_deviation;
        logProbability += d2; 
        logProbability /= 2 * variance;
        result += logProbability;
    }
    return result;
}





static void
__attribute__((__nonnull__))
mle_gauss_likelihood_gradient(const gsl_vector *X, void *restrict params,
                              gsl_vector *restrict g)
{
    const struct mle_gauss_params *const p = params;
    const double variance = mle_gauss_deviation * mle_gauss_deviation;
    const double thetaX = gsl_vector_get(X, 0);
    const double thetaY = gsl_vector_get(X, 1);
    double gradX = 0.0, gradY = 0.0;

    for (size_t j = 0; j < p->no_anchors; j++) {
        const double d2 = (thetaX - p->anchors[j].x) * (thetaX - p->anchors[j].x) +
                          (thetaY - p->anchors[j].y) * (thetaY - p->anchors[j].y);
        const double d = sqrt(d2);
        const double z = p->ranges[j];
        double x = ((2 * mle_gauss_mean * (thetaX - p->anchors[j].x) -
                     2 * z * (thetaX - p->anchors[j].x)) / d +
                    2 * (thetaX - p->anchors[j].x)) / (2 * variance);
        double y = ((2 * mle_gauss_mean * (thetaY - p->anchors[j].y) -
                     2 * z * (thetaY - p->anchors[j].y)) / d +
                    2 * (thetaY - p->anchors[j].y)) / (2 * variance);
        gradX += x;
        gradY += y;
    }
    gsl_vector_set(g, 0, gradX);
    gsl_vector_set(g, 1, gradY);
}




static void
__attribute__((__nonnull__,__flatten__,__hot__))
mle_gauss_likelihood_fdf(const gsl_vector *X, void *restrict params,
                         double *restrict y, gsl_vector *restrict g)
{
    *y = mle_gauss_likelihood_function(X, params);
    mle_gauss_likelihood_gradient(X, params, g);
}




static inline void
__attribute__((__always_inline__,__gnu_inline__,__artificial__,__nonnull__))
mle_gauss_run(const VECTOR* vx, const VECTOR* vy, const VECTOR *restrict r,
              size_t no_anchors, int width, int height,
              VECTOR *restrict resx, VECTOR *restrict resy)
{
    /* Step 0: Set up the likelihood function. */
    struct mle_gauss_params p;
    struct mle_gauss_point2d anchors[no_anchors];
    double ranges[no_anchors];
    p.no_anchors = no_anchors;
    p.anchors = anchors;
    for (size_t i = 0; i < no_anchors; i++) {
        p.anchors[i].x = vx[i][0];
        p.anchors[i].y = vy[i][0];
    }
    p.ranges = ranges;

    gsl_multimin_function_fdf fdf;
    fdf.n = 2u;
    fdf.f = &mle_gauss_likelihood_function;
    fdf.df = &mle_gauss_likelihood_gradient;
    fdf.fdf = &mle_gauss_likelihood_fdf;
    fdf.params = (void*) &p;

    /* Step 1: Calculate an initial estimate. */
    VECTOR sx, sy;
    llsq_run(vx, vy, r, no_anchors, width, height, &sx, &sy);

    /* Step 2: Call the optimiser. */
    const gsl_multimin_fdfminimizer_type *T = gsl_multimin_fdfminimizer_vector_bfgs2;
    gsl_multimin_fdfminimizer *s = gsl_multimin_fdfminimizer_alloc(T, 2);
    gsl_vector *x;
    x = gsl_vector_alloc(2);
    
    for (int i = 0; i < VECTOR_OPS; i++) {
        /* Step 2a: Initialize the parameters. */
        for (size_t j = 0; j < no_anchors; j++) {
            p.ranges[j] = r[j][i];
        }
        gsl_vector_set(x, 0, sx[i]);
        gsl_vector_set(x, 1, sy[i]);

        int iter = 0;
        int status;

        gsl_multimin_fdfminimizer_set(s, &fdf, x, 1e-2, 1e-4);

        /* Step 2b: Iterate the minimization algorithm. */
        do {
            iter++;
            status = gsl_multimin_fdfminimizer_iterate(s);
            if (__builtin_expect(status, 0))
                break;
            status =
                gsl_multimin_test_gradient (s->gradient, mle_gauss_epsilon);
        } while (status == GSL_CONTINUE && iter < mle_gauss_iterations);

        /* Step 2c: Store the result. */
        (*resx)[i] = (float) gsl_vector_get(s->x, 0);
        (*resy)[i] = (float) gsl_vector_get(s->x, 1);
    }

    /* Step 3: Clean up. */
    gsl_multimin_fdfminimizer_free(s);
    gsl_vector_free(x);
}

#endif

/*
  This file is part of LS² - the Localization Simulation Engine of FU Berlin.

  Copyright 2011-2013   Heiko Will, Marcel Kyas, Thomas Hillebrandt,
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

#ifndef _GNU_SOURCE
#  define _GNU_SOURCE
#endif

#if HAVE_CONFIG_H
# include "ls2/ls2-config.h"
#endif

#include <immintrin.h>

#include <assert.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <string.h>
#include <inttypes.h>
#include <float.h>
#include <math.h>
#include <sys/time.h>
#include <sys/resource.h>

#ifdef HAVE_POPT_H
# include <popt.h>
#endif

#include "ls2/ls2.h"
#include "ls2/library.h"
#include "vector_shooter.h"


/*******************************************************************
 *******************************************************************
 ***
 ***   Helpers
 ***
 *******************************************************************
 *******************************************************************/

#include "util/util_median.c"
#include "util/util_random.c"
#include "util/util_vector.c"
#include "util/util_matrix.c"
#include "util/util_circle.c"
#include "util/util_vcircle.c"
#include "util/util_triangle.c"
#include "util/util_points.c"
#include "util/util_misc.c"
#include "util/util_points_opt.c"

#if defined(STAND_ALONE)
#  include INCLUDE_ALG(ALGORITHM)
#  include INCLUDE_EM(ERRORMODEL)
#else
#  include "library.c"
#endif


/*******************************************************************
 *******************************************************************
 ***
 ***   Thread management
 ***
 *******************************************************************
 *******************************************************************/

/*! Whether the threads have been cancelled. */
static volatile bool cancelled;

/*! The number of still running threads. */
static volatile int running;

/*! An array of thread identifiers. */
static pthread_t *ls2_thread;

/*! The total number of threads that are running. */
static volatile size_t ls2_num_threads;



int
cancel_running(void)
{
    if (running <= 0)
        return 0;
    for(size_t t = 0; t < ls2_num_threads; t++)
    	pthread_cancel(ls2_thread[t]);
    cancelled = true;
    return 1;
}


/*******************************************************************
 *******************************************************************
 ***
 ***   Simulation engine: Location based
 ***
 *******************************************************************
 *******************************************************************/

/*! Array of the view names. */
static const char * ls2_result_view_name[] = {
#undef  LS2OUT_VARIANT
#define LS2OUT_VARIANT(tag, name, h5code) name,
#include "ls2/output-variants.h"
    NULL
};



const char * const *
get_result_view_names(void)
{
    return ls2_result_view_name;
}



int __attribute__((__const__))
get_number_of_result_views(void)
{
    return NUM_VARIANTS;
}



ls2_output_variant
ls2_get_view_by_name(const char *name)
{
    ls2_output_variant var;
    for (var = AVERAGE_ERROR; var < NUM_VARIANTS; var++) {
        if (strcasecmp(name, ls2_result_view_name[var]) == 0)
            return var;
    }
    return NUM_VARIANTS;
}



/*! Whether to collect statistics about this thread */
int ls2_verbose = 0;


/*! Parameters to the location-based simulator. */

typedef struct locbased_runparams_t {
    size_t id;
    unsigned int seed;
    vector2 const *anchors;
    size_t no_anchors;
    float * restrict * restrict results;
    uint16_t width;
    uint16_t height;
    size_t from;
    size_t count;
    int_fast64_t runs;
    algorithm_t algorithm;
    error_model_t error_model;
} locbased_runparams_t;


/* The following two arrays are used to store the results.
 * The beginning of the array starts on a cache line, if the cache line
 * size is 64 bytes large.
 *
 * It would be neat if SIZE is divisible by cache line, then each image
 * line would start on a cache line. Alas, default SIZE is 8 * 125.
 * SIZE = 1024 might be much better.
 */
static void*
__attribute__((__nonnull__,__hot__,__flatten__))
ls2_shooter_run(void *rr)
{
    assert(rr != NULL);
    locbased_runparams_t *params = (locbased_runparams_t *) rr;

    VECTOR vx[MAX_ANCHORS];
    VECTOR vy[MAX_ANCHORS];
    VECTOR r[MAX_ANCHORS]; 
    VECTOR distances[MAX_ANCHORS];

    __m128i seed;
    do {
        int seed0 = rand_r(&(params->seed));
        int seed1 = rand_r(&(params->seed));
        int seed2 = rand_r(&(params->seed));
        int seed3 = rand_r(&(params->seed));
        seed = _mm_set_epi32(seed0, seed1, seed2, seed3);
    } while (0);


    // Precalculate Values
    for (size_t i = 0; i < params->no_anchors; i++) {
        vx[i] = VECTOR_BROADCASTF(params->anchors[i].x);
        vy[i] = VECTOR_BROADCASTF(params->anchors[i].y);
    }

    // Precalculate whether we are in the common case of running for the
    // common case.
    const long shortcut =
        (params->results[AVERAGE_ERROR] ||
          params->results[STANDARD_DEVIATION]) &&
          (params->results[ROOT_MEAN_SQUARED_ERROR] ||
          params->results[AVERAGE_X_ERROR] ||
          params->results[STANDARD_DEVIATION_X_ERROR] ||
          params->results[AVERAGE_Y_ERROR] ||
          params->results[STANDARD_DEVIATION_Y_ERROR]);

    // Calculation for every pixel
    for (size_t j = params->from; j < params->from + params->count; j++) {
	const uint16_t x = (uint16_t) (j % params->width);
	const uint16_t y = (uint16_t) (j / params->width);
        const VECTOR tagx = VECTOR_BROADCASTF((float) x);
        const VECTOR tagy = VECTOR_BROADCASTF((float) y);

    	// precalculate real distances
        for (size_t k = 0; k < params->no_anchors; k++) {
            distances[k] = distance(vx[k], vy[k], tagx, tagy);
        }

        float M = 0.0F, M_old, S = 0.0F, cnt = 0.0F, mse = 0.0F;
        float M_X = 0.0F, M_X_old, S_X = 0.0F, C_X = 0.0F;
        float M_Y = 0.0F, M_Y_old, S_Y = 0.0F, C_Y = 0.0F;
        VECTOR min_error = VECTOR_BROADCASTF(FLT_MAX),
               max_error = VECTOR_BROADCASTF(0.0F);

    	// Calculate every pixel runs times 
    	for (int_fast64_t i = 0; i < params->runs; i += VECTOR_OPS) {
            // The results of the algorithm
            VECTOR resx, resy;
            
#if defined(STAND_ALONE)
            EMFUNCTION(error)(&seed, params->no_anchors, distances,
			      vx, vy, tagx, tagy, r);
	    ALGORITHM_RUN(params->no_anchors, vx, vy, r, &resx, &resy);
#else
	    pthread_testcancel();   // Check whether this thread is cancelled.

	    error_model(params->error_model, &seed, distances, vx, vy,
                        params->no_anchors, tagx, tagy, r);
	    algorithm(params->algorithm, vx, vy, r, params->no_anchors,
                      params->width, params->height, &resx, &resy);
#endif

            // Get Errors
	    VECTOR errors = distance(resx, resy, tagx, tagy);

            max_error = VECTOR_MAX(errors, max_error);
	    min_error = VECTOR_MIN(errors, min_error);

            if (params->results[AVERAGE_ERROR] != NULL ||
                params->results[STANDARD_DEVIATION] != NULL) {
                for (int k = 0; k < VECTOR_OPS; k++) {
                    if (!isnan(errors[k])) {
                        cnt += 1.0F;
                        M_old = M;
                        M += (errors[k] - M) / cnt;
                        if (params->results[STANDARD_DEVIATION] != NULL)
                            S += (errors[k] - M) * (errors[k] - M_old);
                     }
                }
            }

            // The common case is to compute the average error, so we
            // optimise for this case by not testing all cases below.
            if (shortcut)
                continue;

            if (params->results[ROOT_MEAN_SQUARED_ERROR] != NULL) {
                VECTOR tmp = errors;
                tmp *= tmp;

	        tmp = VECTOR_HADD(tmp, tmp);
	        tmp = VECTOR_HADD(tmp, tmp);
#  ifdef __AVX__
	        tmp = VECTOR_HADD(tmp, tmp);
#  endif
	        mse += tmp[0];
            }

            if (params->results[AVERAGE_X_ERROR] != NULL ||
                params->results[STANDARD_DEVIATION_X_ERROR] != NULL) {
                for (int k = 0; k < VECTOR_OPS; k++) {
                    if (!isnan(resx[k])) {
                        C_X += 1.0F;
                        M_X_old = M_X;
			const float dx = resx[k] - x;
                        M_X += (dx - M_X_old) / C_X;
                        if (params->results[STANDARD_DEVIATION_X_ERROR] != NULL)
                            S_X += (dx - M_X) * (dx - M_X_old);
                     }
                }
            }
            if (params->results[AVERAGE_Y_ERROR] != NULL ||
                params->results[STANDARD_DEVIATION_Y_ERROR] != NULL) {
                for (int k = 0; k < VECTOR_OPS; k++) {
                    if (!isnan(resy[k])) {
                        C_Y += 1.0F;
                        M_Y_old = M_Y;
			const float dy = resy[k] - y;
                        M_Y += (dy - M_Y_old) / C_Y;
                        if (params->results[STANDARD_DEVIATION_Y_ERROR] != NULL)
                            S_Y += (dy - M_X) * (dy - M_Y_old);
                     }
                }
            }
	}

	const size_t pos = (size_t) (x +  y * params->width);
        if (params->results[AVERAGE_ERROR] != NULL) {
	    params->results[AVERAGE_ERROR][pos] = M;
        }
        if (params->results[STANDARD_DEVIATION] != NULL) {
            params->results[STANDARD_DEVIATION][pos] = sqrtf(S / (cnt - 1.0F));
        }
        if (params->results[MAXIMUM_ERROR] != NULL) {
	    params->results[MAXIMUM_ERROR][pos] =
                vector_max_ps(max_error, 0.0F);
        }
        if (params->results[MINIMUM_ERROR] != NULL) {
	    params->results[MINIMUM_ERROR][pos] =
                vector_min_ps(min_error, FLT_MAX);
        }
        if (params->results[ROOT_MEAN_SQUARED_ERROR] != NULL) {
	    params->results[ROOT_MEAN_SQUARED_ERROR][pos] =
                sqrtf(mse / (float) params->runs);
        }
        if (params->results[AVERAGE_X_ERROR] != NULL) {
	    params->results[AVERAGE_X_ERROR][pos] = M_X;
        }
        if (params->results[STANDARD_DEVIATION_X_ERROR] != NULL) {
	    params->results[STANDARD_DEVIATION_X_ERROR][pos] = sqrtf(S_X / (C_X - 1.0F));
        }
        if (params->results[AVERAGE_Y_ERROR] != NULL) {
	    params->results[AVERAGE_Y_ERROR][pos] = M_Y;
        }
        if (params->results[STANDARD_DEVIATION_Y_ERROR] != NULL) {
	    params->results[STANDARD_DEVIATION_Y_ERROR][pos] = sqrtf(S_Y / (C_Y - 1.0F)); 
        }
    }
    running--;

    if (ls2_verbose >= 2) {
        struct rusage resources;
        getrusage(RUSAGE_THREAD, &resources);
        fprintf(stderr, "Thread %zu: %d.%06d sec.\n", params->id,
                (int)resources.ru_utime.tv_sec, (int)resources.ru_utime.tv_usec);
        fflush(stderr);
    }
    return NULL;
}



/************************************************************************
 *****
 ***** Start threads and distribute work to them
 *****
 ************************************************************************/

/*!
 * \brief Estimates the position for each place on the playing field.
 *
 * \param[in] alg        A number that indicates the position estimation
 *                       algorithm.
 * \param[in] em         A number that indicates the error model.
 * \param[in] no_anchors The number of anchor nodes to use.
 * \param[in] anchors    Array of anchors nodes of length [no_anchors].
 * \param[in] width      Width of the playing field.
 * \param[in] height     Height of the playing field.
 */
void __attribute__((__nonnull__))
ls2_distribute_work_shooter(const int alg, const int em,
                            const int num_threads, const int64_t runs,
                            const long seed,
                            const vector2* anchors, const size_t no_anchors,
			    float *results[NUM_VARIANTS],
                            const int width, const int height)
{
    ls2_num_threads = (size_t) num_threads;
    const size_t slice = (size_t) (width * height) / ls2_num_threads;
    locbased_runparams_t *params;

    running = 0;

#if defined(STAND_ALONE)
    EMFUNCTION(setup)(anchors, no_anchors);
#else
    error_model_setup(em, anchors, no_anchors);
#endif

    params = (locbased_runparams_t *) calloc(ls2_num_threads, sizeof(locbased_runparams_t));
    if (params == NULL) {
        perror("calloc()");
        exit(EXIT_FAILURE);
    }

    ls2_thread = (pthread_t *) calloc(ls2_num_threads, sizeof(pthread_t));
    if (ls2_thread == NULL) {
        perror("calloc()");
        exit(EXIT_FAILURE);
    }

    srand((unsigned int) seed);

    // distribute work to threads
    for (size_t t = 0; t < ls2_num_threads; t++) {
        params[t].id = t;
        params[t].seed = (unsigned int) rand() + (unsigned int) t;
        params[t].no_anchors = (size_t)no_anchors;
        params[t].anchors = anchors;
        params[t].width = (uint16_t) width;
        params[t].height = (uint16_t) height;
        params[t].results = results;
        params[t].from = (uint32_t) (t * slice);
        params[t].count = (uint32_t) slice;
        params[t].runs = runs;
        params[t].algorithm = alg;
        params[t].error_model = em;

        if (pthread_create(&ls2_thread[t], NULL, ls2_shooter_run, &params[t])) {
            perror("pthread_create()");
            free(ls2_thread);
            exit(EXIT_FAILURE);
        }
        running++;
    }


    // Sync Threads if work has been done
    for(size_t t = 0; t < ls2_num_threads; t++)
        pthread_join(ls2_thread[t], NULL);

    free(ls2_thread);
    free(params);
}





extern int
compute_locbased(const int alg, const int em,
                 const int num_threads, const int64_t runs,
                 const float *anchor_x, const float *anchor_y,
                 const int no_anchors, float* results[NUM_VARIANTS],
                 const int width, const int height)
{
    vector2 *anchors;

    cancelled = false;

    // parse and normalize arguments
    anchors = calloc((size_t) no_anchors, sizeof(vector2));
    if (anchors == NULL) {
        perror("compute_locbased(): calloc()");
        return -1;
    }
 
    for (int i = 0; i < no_anchors; i++) {
	    anchors[i].x = anchor_x[i];
	    anchors[i].y = anchor_y[i];
    }

    ls2_distribute_work_shooter(alg, em, num_threads, runs, time(NULL),
                                anchors, (size_t) no_anchors,
                                results, width, height);

    free(anchors);

    if (cancelled)
	return -1;
    else
        return 0;
}





/*******************************************************************
 *******************************************************************
 ***
 ***   Simulation engine: Inverted
 ***
 *******************************************************************
 *******************************************************************/

typedef struct inverted_runparams_t {
    int id;
    unsigned int seed;
    float tag_x, tag_y;
    vector2 const *anchors;
    size_t no_anchors;
    int_fast64_t runs;
    uint_fast64_t *result;
    float cx, cy, cn;
    float *results;
    int width, height;
    algorithm_t algorithm;
    error_model_t error_model;
} inverted_runparams_t;


static void* ls2_inverse_run(void *rr)
{
    inverted_runparams_t *params = (inverted_runparams_t *) rr;
    __m128i seed;
    const int_fast64_t runs =
        params->runs / (int_fast64_t) (ls2_num_threads * VECTOR_OPS);

    do { // Limit life time of temporaries.
        int seed0 = rand_r(&(params->seed));
        int seed1 = rand_r(&(params->seed));
        int seed2 = rand_r(&(params->seed));
        int seed3 = rand_r(&(params->seed));
        seed = _mm_set_epi32(seed0, seed1, seed2, seed3);
    } while (0);

    VECTOR vx[MAX_ANCHORS];
    VECTOR vy[MAX_ANCHORS];
    VECTOR r[MAX_ANCHORS]; 
    VECTOR distances[MAX_ANCHORS];    

    // Precalculate Values
    const VECTOR tagx = VECTOR_BROADCASTF(params->tag_x);
    const VECTOR tagy = VECTOR_BROADCASTF(params->tag_y);

    // precalculate real distances
    for (size_t i = 0; i < params->no_anchors; i++) {
        vx[i] = VECTOR_BROADCASTF(params->anchors[i].x);
        vy[i] = VECTOR_BROADCASTF(params->anchors[i].y);
        distances[i] = distance(vx[i], vy[i], tagx, tagy);
    }

    // Calculation for every pixel
    for (int_fast64_t j = 0; j < runs; j++) {
        VECTOR resx, resy;

#if defined(STAND_ALONE)
	EMFUNCTION(error)(&seed, params->no_anchors, distances,
			  vx, vy, tagx, tagy, r);
	ALGORITHM_RUN(params->no_anchors, vx, vy, r, &resx, &resy);
#else
        pthread_testcancel();
	error_model(params->error_model, &seed, distances, vx, vy,
                    params->no_anchors, tagx, tagy, r);
	algorithm(params->algorithm, vx, vy, r, params->no_anchors,
                  params->width, params->height, &resx, &resy);
#endif

        // errors[j] = distance(resx[j], resy[j], tagx, tagy);
        for (int k = 0; k < VECTOR_OPS; ++k) {
            if (!isnan(resx[k]) && !isnan(resy[k])) {
                const int x = (int) roundf(resx[k]);
                const int y = (int) roundf(resy[k]);		
		if (0 <= x && x < params->width && 0 <= y && y < params->height) {
		    params->result[x + params->width * y] += 1;
		}
                params->cx += resx[k];
                params->cy += resy[k];
                params->cn += 1.0F;
            }
        }
    }
    running--;

    if (ls2_verbose >= 2) {
        struct rusage resources;
        getrusage(RUSAGE_THREAD, &resources);
        fprintf(stderr, "Thread %d: %d.%06d sec.\n", params->id,
                (int)resources.ru_utime.tv_sec, (int)resources.ru_utime.tv_usec);
        fflush(stderr);
    }
    return NULL;
}



/************************************************************************
 *****
 ***** Start threads and distribute work to them
 *****
 ************************************************************************/

void
ls2_distribute_work_inverted(const int alg, const int em,
			     const int num_threads, const int64_t runs,
                             const long seed,
                             const float tag_x, const float tag_y,
			     const vector2* anchors, const size_t no_anchors,
			     float *results, const int width, const int height,
			     float *center_x, float *center_y)
{
    running = 0;

#if defined(STAND_ALONE)
    EMFUNCTION(setup)(anchors, no_anchors);
#else
    error_model_setup(em, anchors, no_anchors);
#endif

    // distribute work to threads
    inverted_runparams_t *params;
    const size_t sz =
	((size_t) width) * ((size_t) height) * sizeof(uint_fast64_t);

    ls2_num_threads = (size_t) num_threads;

    params = (inverted_runparams_t *) calloc(ls2_num_threads, sizeof(inverted_runparams_t));
    if (params == NULL) {
        perror("calloc()");
        exit(EXIT_FAILURE);
    }

    srand((unsigned int) seed);

    for (int t = 0; t < num_threads; t++) {
        params[t].id = t;
        params[t].seed = (unsigned int) rand() + (unsigned int) t;
	params[t].tag_x = tag_x;
	params[t].tag_y = tag_y;
	params[t].anchors = anchors;
	params[t].no_anchors = no_anchors;
	params[t].results = results;
	params[t].width = width;
	params[t].height = height;
	params[t].runs = runs;
	params[t].algorithm = alg;
	params[t].error_model = em;
        if (posix_memalign((void **) &(params[t].result), ALIGNMENT, sz) != 0) {
            fprintf(stderr, "allocate result\n");
            exit(EXIT_FAILURE);
        }
        memset(params[t].result, 0, sz);
    }

    ls2_thread = (pthread_t *) calloc((size_t) num_threads, sizeof(pthread_t));
    if (ls2_thread == NULL) {
        perror("calloc()");
        exit(EXIT_FAILURE);
    }

    for (int t = 0; t < num_threads; t++) {
        if (pthread_create(&ls2_thread[t], NULL, ls2_inverse_run, &params[t])) {
            perror("pthread_create()");
            free(ls2_thread);
            exit(EXIT_FAILURE);
        }
        running++;
    }

    // Sync Threads if work has been done
    for (int t = 0; t < num_threads; t++)
        pthread_join(ls2_thread[t], NULL);

    free(ls2_thread);

    /*
     * Evaluate the results.
     */

    /* Accumulate all results and store them in the first thread's image. */
    uint_fast64_t max_val = 0;
    for (int t = 1; t < num_threads; t++) {
	params[0].cx += params[t].cx;
	params[0].cy += params[t].cy;
	params[0].cn += params[t].cn;
        for (int i = 0; i < width * height; i++) {
            params[0].result[i] += params[t].result[i];
            if (params[0].result[i] > max_val) max_val = params[0].result[i];
        }
    }

    fprintf(stderr, "Maximum value: %" PRIuFAST64 "\n", max_val);

    /* After having counted the hits per positions, convert them to an
     * array of doubles from [0; 1.0], indicating the relative number of
     * hits in a location.
     */
    const double f_max_val = (double) max_val;
    for (int y = 0; y < params->height; y++) {
        for (int x = 0; x < params->width; x++) {
	    const int pos = x + params->width * y;
            const double samples = (double) params[0].result[pos];
            if (max_val > 0) {
	        params->results[pos] = (float) (samples / f_max_val);
	    } else {
		params->results[pos] = 0.0;
	    }
	}
    }
    *center_x = params[0].cx / params[0].cn;
    *center_y = params[0].cy / params[0].cn;

    for (int t = 0; t < num_threads; t++) {
	free(params[t].result);
    }
    free(params);
}





extern int
compute_inverse(const int alg, const int em,
		const int num_threads, const int64_t runs,
		const float *anchor_x, const float *anchor_y,
		const int no_anchors, const float tag_x, const float tag_y,
		float* results, const int width, const int height,
		float *center_x, float *center_y)
{
    vector2 *anchors;

    cancelled = false;

    // parse and normalize arguments
    anchors = calloc((size_t) no_anchors, sizeof(vector2));
    if (anchors == NULL) {
        perror("compute_locbased(): calloc()");
        return -1;
    }
 
    for (int i = 0; i < no_anchors; i++) {
	    anchors[i].x = anchor_x[i];
	    anchors[i].y = anchor_y[i];
    }

    ls2_distribute_work_inverted(alg, em, num_threads, runs, time(NULL),
				 tag_x, tag_y,
				 anchors, (size_t) no_anchors, results,
                                 width, height, center_x, center_y);

    free(anchors);

    if (cancelled)
	return -1;
    else
        return 0;
}





/**********************************************************************
 **********************************************************************
 ***
 *** Estimator
 ***
 ***
 **********************************************************************
 **********************************************************************/


typedef struct estimator_runparams_t {
    size_t id;
    vector2 const *anchors;
    size_t no_anchors;
    float **results;
    uint16_t width;
    uint16_t height;
    size_t from;
    size_t count;
    estimator_t estimator;
} estimator_runparams_t;


/* The following two arrays are used to store the results.
 * The beginning of the array starts on a cache line, if the cache line
 * size is 64 bytes large.
 */
static void*
ls2_estimator_run(void *rr)
{
    estimator_runparams_t *params = (estimator_runparams_t *) rr;
    assert(params != NULL);

    // Calculation for every pixel
    for (size_t pos = params->from; pos < params->from + params->count; pos++) {
        vector2 location;
        float result;

        location.x = (float) (pos % params->width);
        location.y = (float) (pos / params->width);

	pthread_testcancel();   // Check whether this thread is cancelled.

	result = estimate(params->estimator, params->anchors,
                          params->no_anchors, &location);

        if (params->results[ROOT_MEAN_SQUARED_ERROR] != NULL) {
	    params->results[ROOT_MEAN_SQUARED_ERROR][pos] = sqrtf(result);
        }
    }
    running--;

    if (ls2_verbose >= 2) {
        struct rusage resources;
        getrusage(RUSAGE_THREAD, &resources);
        fprintf(stderr, "Thread %zu: %d.%06d sec.\n", params->id,
                (int)resources.ru_utime.tv_sec, (int)resources.ru_utime.tv_usec);
        fflush(stderr);
    }
    return NULL;
}




/************************************************************************
 *****
 ***** Start threads and distribute work to them
 *****
 ************************************************************************/

/*!
 * \brief Estimates the position for each place on the playing field.
 *
 * \param[in] est        A number that indicates the variance estimation
 *                       algorithm.
 * \param[in] no_anchors The number of anchor nodes to use.
 * \param[in] anchors    Array of anchors nodes of length [no_anchors].
 * \param[out] results   Array of playing fields, only MEAN_SQUARED_ERROR
 *                       must be not NULL.
 * \param[in] width      Width of the playing field.
 * \param[in] height     Height of the playing field.
 */
void
ls2_distribute_work_estimator(const int est, const int num_threads,
			      const vector2* anchors, const size_t no_anchors,
			      float *results[NUM_VARIANTS],
			      const int width, const int height)
{
    ls2_num_threads = (size_t) num_threads;
    const size_t slice = (size_t) (width * height) / ls2_num_threads;
    estimator_runparams_t *params;

    running = 0;

    params = (estimator_runparams_t *) calloc(ls2_num_threads, sizeof(estimator_runparams_t));
    if (params == NULL) {
        perror("calloc()");
        exit(EXIT_FAILURE);
    }

    ls2_thread = (pthread_t *) calloc(ls2_num_threads, sizeof(pthread_t));
    if (ls2_thread == NULL) {
        perror("calloc()");
        exit(EXIT_FAILURE);
    }

    // distribute work to threads
    for (size_t t = 0; t < ls2_num_threads; t++) {
        params[t].id = t;
        params[t].anchors = anchors;
        params[t].no_anchors = (size_t)no_anchors;
        params[t].results = results;
        params[t].width = (uint16_t) width;
        params[t].height = (uint16_t) height;
        params[t].from = (uint32_t) (t * slice);
        params[t].count = (uint32_t) slice;
        params[t].estimator = est;

        if (pthread_create(&ls2_thread[t], NULL, ls2_estimator_run, &params[t])) {
            perror("pthread_create()");
            free(ls2_thread);
            exit(EXIT_FAILURE);
        }
        running++;
    }


    // Sync Threads if work has been done
    for(size_t t = 0; t < ls2_num_threads; t++)
        pthread_join(ls2_thread[t], NULL);

    free(ls2_thread);
    free(params);
}





int
compute_estimates(const int est, const int num_threads,
		  const float *anchor_x, const float *anchor_y,
		  const int no_anchors,
		  float* results[NUM_VARIANTS], const int width,
		  const int height)
{
    vector2 *anchors;

    cancelled = false;

    // parse and normalize arguments
    anchors = calloc((size_t) no_anchors, sizeof(vector2));
    if (anchors == NULL) {
        perror("compute_locbased(): calloc()");
        return -1;
    }
 
    for (int i = 0; i < no_anchors; i++) {
	    anchors[i].x = anchor_x[i];
	    anchors[i].y = anchor_y[i];
    }

    ls2_distribute_work_estimator(est, num_threads,
				  anchors, (size_t) no_anchors,
				  results, width, height);

    free(anchors);

    if (cancelled)
	return -1;
    else
        return 0;
}

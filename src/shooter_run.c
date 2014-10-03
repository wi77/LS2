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

#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <string.h>
#include <inttypes.h>
#include <float.h>
#include <math.h>
#include <signal.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <time.h>
#include <unistd.h>

#include <glib.h>

#include "ls2/library.h"
#include "ls2/output.h"
#include "ls2/ls2.h"
#include "ls2/progress.h"
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

#include "library.c"


#ifndef DEFAULT_RUNS
#define DEFAULT_RUNS  0x8000U
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
volatile size_t ls2_running;

/*! An array of thread identifiers. */
static pthread_t *ls2_thread;

/*! The total number of threads that are running. */
volatile size_t ls2_num_threads;



int
cancel_running(void)
{
    if (ls2_running <= 0)
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


/* From progress.c */
extern volatile size_t ls2_progress_total;



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
    uint_fast64_t runs;
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
    g_assert(rr != NULL);
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

        float M = 0.0F, M_old, S = 0.0F, cnt = 0.0F;
        float MSE = 0.0F, MSE_old, C_MSE = 0.0F;
        float M_X = 0.0F, M_X_old, S_X = 0.0F, C_X = 0.0F;
        float M_Y = 1.0F, M_Y_old, S_Y = 0.0F, C_Y = 0.0F;
        uint_fast64_t failures = 0U; // How often did it fail (nan)?

        VECTOR min_error = VECTOR_BROADCASTF(FLT_MAX),
               max_error = VECTOR_BROADCASTF(0.0F);

    	// Calculate every pixel runs times 
    	for (uint_fast64_t i = 0; i < params->runs; i += VECTOR_OPS) {
            // The results of the algorithm
            VECTOR resx, resy;
            
	    pthread_testcancel();   // Check whether this thread is cancelled.

            if (__builtin_expect(ls2_progress_total > 0, 0)) {
                const uint_fast64_t step =
                    (j - params->from) * params->runs + i;
                if (__builtin_expect((step & (DEFAULT_RUNS - 1U)) == 0, 0)) {
                    ls2_update_progress_bar(DEFAULT_RUNS);
                }
            }

	    error_model(params->error_model, &seed, distances, vx, vy,
                        params->no_anchors, tagx, tagy, r);
	    algorithm(params->algorithm, vx, vy, r, params->no_anchors,
                      params->width, params->height, &resx, &resy);

            // Calculate errors and update statistics.
	    const VECTOR errors = distance(resx, resy, tagx, tagy);
            const VECTOR sqerror = errors * errors;

            max_error = VECTOR_MAX(errors, max_error);
	    min_error = VECTOR_MIN(errors, min_error);

            for (int k = 0; k < VECTOR_OPS; k++) {
                if (isnan(errors[k]) == 0) {
                    cnt += 1.0F;
                    M_old = M;
                    M += (errors[k] - M) / cnt;
                    if (params->results[STANDARD_DEVIATION] != NULL)
                        S += (errors[k] - M) * (errors[k] - M_old);
                } else {
                    failures += 1;
                }

                if (isnan(sqerror[k]) == 0) {
                    C_MSE += 1.0F;
                    MSE_old = MSE;
                    MSE += (sqerror[k] - MSE_old) / C_MSE;
                }

                if (isnan(resx[k]) == 0 && isnan(resy[k]) == 0) {
                    C_X += 1.0F;
                    M_X_old = M_X;
		    const float dx = resx[k] - x;
                    M_X += (dx - M_X_old) / C_X;
                    S_X += (dx - M_X) * (dx - M_X_old);

                    C_Y += 1.0F;
                    M_Y_old = M_Y;
		    const float dy = resy[k] - y;
                    M_Y += (dy - M_Y_old) / C_Y;
                    S_Y += (dy - M_Y) * (dy - M_Y_old);
                }
            }
	}

        /* Enter the results into the result arrays. */
	const size_t pos = (size_t) (x +  y * params->width);

        /* Set to NAN explicitely, if we are never successful.  */
	params->results[AVERAGE_ERROR][pos] = (cnt > 0.0F) ? M : NAN;

        /* Set to NAN explicitely, if we are never successful.  */
        params->results[STANDARD_DEVIATION][pos] =
             (cnt > 1.0F) ? sqrtf(S / (cnt - 1.0F)) : NAN;

	params->results[MAXIMUM_ERROR][pos] = vector_max_ps(max_error, 0.0F);
	params->results[MINIMUM_ERROR][pos] = vector_min_ps(min_error, FLT_MAX);

	params->results[FAILURES][pos] =
            ((float) failures) / ((float) params->runs);
        if (__builtin_expect(ls2_verbose > 0, 0)) {
            if (__builtin_expect(params->results[FAILURES][pos] > 0.0, 0)) {
                g_warning("%" PRIuFAST64 " of %" PRIuFAST64
                          " runs failed at (%d, %d)\n",
                          failures, params->runs, x, y);
            }
        }

	params->results[ROOT_MEAN_SQUARED_ERROR][pos] = sqrtf(MSE);
	params->results[AVERAGE_X_ERROR][pos] = (C_X > 0.0F) ? M_X : NAN;
	params->results[STANDARD_DEVIATION_X_ERROR][pos] =
            (C_X > 1.0F) ? sqrtf(S_X / (C_X - 1.0F)) : NAN;
	params->results[AVERAGE_Y_ERROR][pos] = (C_Y > 0.0F) ? M_Y : NAN;
	params->results[STANDARD_DEVIATION_Y_ERROR][pos] =
            (C_Y > 1.0F) ? sqrtf(S_Y / (C_Y - 1.0F)) : NAN;
    }
    ls2_running--;

    if (__builtin_expect(ls2_verbose >= 2, 0)) {
        struct rusage resources;
        getrusage(RUSAGE_THREAD, &resources);
        g_message("Thread %zu: %d.%06d sec.\n", params->id,
                  (int)resources.ru_utime.tv_sec, (int)resources.ru_utime.tv_usec);
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
void __attribute__((__nonnull__(8,10)))
ls2_distribute_work_shooter(void *alg_params __attribute__((__unused__)), const algorithm_t alg,
                            void *em_params __attribute__((__unused__)), const error_model_t em,
                            const int num_threads, const int64_t runs,
                            const long seed,
                            const vector2* anchors, const size_t no_anchors,
			    float *results[NUM_VARIANTS],
                            const size_t width, const size_t height)
{
    ls2_num_threads = (size_t) num_threads;
    const size_t slice = (size_t) (width * height) / ls2_num_threads;
    locbased_runparams_t *params;

    ls2_running = 0;

    error_model_setup(em, anchors, no_anchors);

    params = g_new(locbased_runparams_t, ls2_num_threads);
    ls2_thread = g_new(pthread_t, ls2_num_threads);

    srand((unsigned int) seed);

    // Set up the parameters.
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
        params[t].runs = (uint_fast64_t) runs;
        params[t].algorithm = alg;
        params[t].error_model = em;
    }

    /* Create the threads. */
    if (ls2_num_threads > 1) {
        for (size_t t = 0; t < ls2_num_threads; t++) {
            if (pthread_create(&ls2_thread[t], NULL, ls2_shooter_run,
                               &params[t])) {
                perror("pthread_create()");
                exit(EXIT_FAILURE);
            }
            ls2_running++;
        }


        // Sync Threads if work has been done
        for(size_t t = 0; t < ls2_num_threads; t++) {
            pthread_join(ls2_thread[t], NULL);
        }
    } else {
        /* Since there is only one thread, we call the run code directly.
         * This should make debugging simpler.
         */
        ls2_running = 1;
        ls2_shooter_run(&(params[0]));
    }
    g_free(ls2_thread);
    g_free(params);
}





/*
 * This function should only be called by tha Java api.
 */
extern int __attribute__((__nonnull__(7,8,10)))
ls2_compute_locbased(void *alg_params, const algorithm_t alg,
                     void *em_params, const error_model_t em,
                     const int num_threads, const int64_t runs,
                     const float *anchor_x, const float *anchor_y,
                     const int no_anchors, float* results[NUM_VARIANTS],
                     const int width, const int height)
{
    vector2 *anchors;

    cancelled = false;

    ls2_reset_progress_bar((size_t)(runs * width * height), NULL);

    // parse and normalize arguments
    anchors = g_new(vector2, (size_t) no_anchors);
 
    for (int i = 0; i < no_anchors; i++) {
	    anchors[i].x = anchor_x[i];
	    anchors[i].y = anchor_y[i];
    }

    ls2_distribute_work_shooter(alg_params, alg, em_params, em, num_threads,
                                runs, time(NULL), anchors, (size_t) no_anchors,
                                results, (size_t) width, (size_t) height);

    g_free(anchors);

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
    float cx, sx, cy, sy, cn;
    size_t width, height;
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

    float M_X = 0.0F, M_X_old, S_X = 0.0F, N = 0.0F,
          M_Y = 0.0F, M_Y_old, S_Y = 0.0F;

    // Calculation for every pixel
    for (int_fast64_t j = 0; j < runs; j++) {
        VECTOR resx, resy;

        pthread_testcancel();
        if (__builtin_expect(ls2_progress_total > 0, 0)) {
	    if (__builtin_expect(((j + 1u) & (DEFAULT_RUNS/VECTOR_OPS-1U)) == 0, 0)) {
                ls2_update_progress_bar(DEFAULT_RUNS);
            }
        }

	error_model(params->error_model, &seed, distances, vx, vy,
                    params->no_anchors, tagx, tagy, r);
	algorithm(params->algorithm, vx, vy, r, params->no_anchors,
                  params->width, params->height, &resx, &resy);

        // errors[j] = distance(resx[j], resy[j], tagx, tagy);
        for (int k = 0; k < VECTOR_OPS; ++k) {
            if (!isnan(resx[k]) && !isnan(resy[k])) {
                const int x = (int) roundf(resx[k]);
                const int y = (int) roundf(resy[k]);		
		if (0 <= x && x < (int) params->width && 0 <= y && y < (int) params->height) {
		    params->result[x + (int) params->width * y] += 1;
		}
                N   += 1.0F;
                M_X_old = M_X;
                M_X += (resx[k] - M_X_old) / N;
                S_X += (resx[k] - M_X) * (resx[k] - M_X_old);
                M_Y_old = M_Y;
                M_Y += (resy[k] - M_Y_old) / N;
                S_Y += (resy[k] - M_Y) * (resy[k] - M_Y_old);
            }
        }

        params->cn = N;
        params->cx = M_X;
        params->sx = S_X / N;
        params->cy = M_Y;
        params->sy = S_Y / N;
    }
    ls2_running--;

    if (ls2_verbose >= 2) {
        struct rusage resources;
        getrusage(RUSAGE_THREAD, &resources);
        g_message("Thread %d: %d.%06d sec.\n", params->id,
                  (int)resources.ru_utime.tv_sec,
                  (int)resources.ru_utime.tv_usec);
    }
    return NULL;
}



/************************************************************************
 *****
 ***** Start threads and distribute work to them
 *****
 ************************************************************************/

void
ls2_distribute_work_inverted(void *alg_params, const algorithm_t alg,
                             void *em_params, const error_model_t em,
			     const int num_threads, const int64_t runs,
                             const long seed,
                             const float tag_x, const float tag_y,
			     const vector2 *restrict anchors,
                             const size_t no_anchors,
			     uint64_t *restrict results,
                             const size_t width, const size_t height,
			     float *restrict center_x, float *restrict sdev_x,
                             float *center_y, float *restrict sdev_y)
{
    ls2_running = 0;

    if (alg_params) { /* silence warning */ }
    if (em_params) { /* silence warning */ }

    error_model_setup(em, anchors, no_anchors);

    // distribute work to threads
    inverted_runparams_t *params;
    const size_t sz =
	((size_t) width) * ((size_t) height) * sizeof(uint_fast64_t);

    ls2_num_threads = (size_t) num_threads;

    params = g_new(inverted_runparams_t, ls2_num_threads);

    srand((unsigned int) seed);

    for (int t = 0; t < num_threads; t++) {
        params[t].id = t;
        params[t].seed = (unsigned int) rand() + (unsigned int) t;
	params[t].tag_x = tag_x;
	params[t].tag_y = tag_y;
	params[t].anchors = anchors;
	params[t].no_anchors = no_anchors;
	params[t].width = width;
	params[t].height = height;
	params[t].runs = runs;
	params[t].algorithm = alg;
	params[t].error_model = em;
        if (posix_memalign((void **) &(params[t].result), ALIGNMENT, sz) != 0) {
            g_error("allocate result\n");
            exit(EXIT_FAILURE);
        }
        memset(params[t].result, 0, sz);
    }

    ls2_thread = g_new(pthread_t, (size_t) num_threads);

    for (int t = 0; t < num_threads; t++) {
        if (pthread_create(&ls2_thread[t], NULL, ls2_inverse_run, &params[t])) {
            perror("pthread_create()");
            g_free(ls2_thread);
            exit(EXIT_FAILURE);
        }
        ls2_running++;
    }

    // Sync Threads if work has been done
    for (int t = 0; t < num_threads; t++)
        pthread_join(ls2_thread[t], NULL);

    g_free(ls2_thread);

    /*
     * Evaluate the results.
     */

   *center_x = params[0].cx;
   *sdev_x   = params[0].sx;
   *center_y = params[0].cy;
   *sdev_y   = params[0].sy;

    /* Accumulate all results and store them in the first thread's image. */
    for (int t = 1; t < num_threads; t++) {
        for (size_t i = 0; i < width * height; i++) {
            params[0].result[i] += params[t].result[i];
        }
	*center_x += params[t].cx;
        *sdev_x   += params[t].sx;
	*center_y += params[t].cy;
        *sdev_y   += params[t].sy;
    }
    *center_x /= ((float) num_threads);
    *sdev_x = sqrtf(*sdev_x / (float) num_threads);
    *center_y /= ((float) num_threads);
    *sdev_y = sqrtf(*sdev_y / (float) num_threads);

    for (size_t y = 0; y < params->height; y++) {
        for (size_t x = 0; x < params->width; x++) {
	    const size_t pos = x + params->width * y;
            const uint64_t samples = (uint64_t) params[0].result[pos];
	    results[pos] = samples;
	}
    }
    for (int t = 0; t < num_threads; t++) {
	free(params[t].result);
    }
    g_free(params);
}





extern int
ls2_compute_inverse(void *alg_params, const algorithm_t alg,
                    void *em_params, const error_model_t em,
		    const int num_threads,
                    const int64_t runs, const float *restrict anchor_x,
                    const float *restrict anchor_y, const int no_anchors,
                    const float tag_x, const float tag_y,
		    uint64_t *restrict result,
                    const int width, const int height,
		    float *restrict center_x, float *restrict sdev_x,
                    float *restrict center_y, float *restrict sdev_y)
{
    vector2 *anchors;

    cancelled = false;

    // parse and normalize arguments
    anchors = g_new(vector2, (size_t) no_anchors);
 
    for (int i = 0; i < no_anchors; i++) {
	    anchors[i].x = anchor_x[i];
	    anchors[i].y = anchor_y[i];
    }

    ls2_distribute_work_inverted(alg_params, alg, em_params, em, num_threads,
                                 runs, time(NULL), tag_x, tag_y, anchors,
                                 (size_t) no_anchors,
                                 result, (size_t) width, (size_t) height,
                                 center_x, sdev_x, center_y, sdev_y);

    g_free(anchors);

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
    g_assert(params != NULL);

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
    ls2_running--;

    if (ls2_verbose >= 2) {
        struct rusage resources;
        getrusage(RUSAGE_THREAD, &resources);
        g_message("Thread %zu: %d.%06d sec.\n", params->id,
                  (int)resources.ru_utime.tv_sec,
                  (int)resources.ru_utime.tv_usec);
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
ls2_distribute_work_estimator(void *est_params, const estimator_t est,
                              const int num_threads,
			      const vector2* anchors, const size_t no_anchors,
			      float *results[NUM_VARIANTS],
			      const int width, const int height)
{
    ls2_num_threads = (size_t) num_threads;
    const size_t slice = (size_t) (width * height) / ls2_num_threads;
    estimator_runparams_t *params;

    ls2_running = 0;

    if (est_params) { /* silence warning */ }
    params = g_new(estimator_runparams_t, ls2_num_threads);
    ls2_thread = g_new(pthread_t, ls2_num_threads);

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
            g_free(ls2_thread);
            exit(EXIT_FAILURE);
        }
        ls2_running++;
    }


    // Sync Threads if work has been done
    for(size_t t = 0; t < ls2_num_threads; t++)
        pthread_join(ls2_thread[t], NULL);

    g_free(ls2_thread);
    g_free(params);
}





int
ls2_compute_estimates(void *est_params, const estimator_t est,
                      const int num_threads,
		      const float *anchor_x, const float *anchor_y,
		      const size_t no_anchors,
		      float* results[NUM_VARIANTS], const int width,
		      const int height)
{
    vector2 *anchors;

    cancelled = false;

    // parse and normalize arguments
    anchors = g_new(vector2, no_anchors);
 
    for (size_t i = 0; i < no_anchors; i++) {
	    anchors[i].x = anchor_x[i];
	    anchors[i].y = anchor_y[i];
    }

    ls2_distribute_work_estimator(est_params, est, num_threads,
				  anchors, no_anchors,
				  results, width, height);

    g_free(anchors);

    if (cancelled)
	return -1;
    else
        return 0;
}

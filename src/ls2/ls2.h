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


 * SVN revision information:
 * @version $Revision$:
 * @author  $Author$:
 * @date    $Date$:
 * @URL     $URL$:
 */

#ifndef INCLUDED_LS2_H
#define INCLUDED_LS2_H

// The number of anchors which are supported, usually defined via makefile.
// more anchors mean more memory consumption
#ifndef MAX_ANCHORS
#  define MAX_ANCHORS 16
#endif

/*! Enumerates the different variants of information gathered by
 * shooter_run().
 */
typedef enum ls2_output_variant {
#undef  LS2OUT_VARIANT
#define LS2OUT_VARIANT(tag, name, h5code) tag,
#include "ls2/output-variants.h"
    NUM_VARIANTS
} ls2_output_variant;

typedef struct vector2 {
    float x;
    float y;
} vector2;

/*! Whether to collect statistics about this thread */
extern int ls2_verbose;


extern int
get_algorithm_by_name(const char *)  __attribute__((__const__));

extern const char * const *
get_algorithms(void) __attribute__((__const__));

extern int
get_estimator_by_name(const char *)  __attribute__((__const__));

extern const char * const *
get_estimators(void) __attribute__((__const__));

extern int
get_error_model_by_name(const char *) __attribute__((__const__));

extern const char * const *
get_error_models(void) __attribute__((__const__));

extern const char * const *
get_result_view_names(void) __attribute__((__const__));

extern int
get_number_of_result_views(void) __attribute__((__const__));

extern ls2_output_variant
ls2_get_view_by_name(const char *name);




/*!
 * \brief Estimates the position for each place on the playing field.
 *
 * \param[in] alg        A number that indicates the position estimation
 *                       algorithm.
 * \param[in] em         A number that indicates the error model.
 * \param[in] num_threads  Number of threads to use.
 * \param[in] runs      Number of runs per location on the discrete grid.
 * \param[in] anchors    Array of anchors nodes of length [no_anchors].
 * \param[in] no_anchors The number of anchor nodes to use.
 * \param[out] results   Array of result arrays, indexed by ls2_output_variant.
 *                       If an entry is \c NULL, no information will be
 *                       collected for this ls2_output_variant. Otherwise it must
 *                       be an array that can hold width * height float
 *                       values.
 * \param[in] width      Width of the playing field.
 * \param[in] height     Height of the playing field.
 */
extern void __attribute__((__nonnull__))
ls2_distribute_work_shooter(int alg, int em,
                            const int num_threads, const int64_t runs,
                            const long seed,
                            const vector2* anchors, const size_t no_anchors,
			    float *results[NUM_VARIANTS],
                            const int width, const int height);

/*!
 * Perform a simulation based on locations.
 *
 * \param[in] alg  The algorithm to use. Use any value of algorithm_t.
 * \param[in] em   The error model to use. Use any value of error_model_t.
 * \param[in] num_threads  Number of threads to use.
 * \param[in] runs      Number of runs per location on the discrete grid.
 * \param[in] anchor_x  Array of X coordinates of the anchors.
 * \param[in] anchor_y  Array of Y coordinates of the anchors.
 * \param[in] no_anchors  Number of anchors
 * \param[out] results   An array of arrays holding the results for ls2_output_variant.
 *                       If an entry is \c NULL, no information will be
 *                       collected for this ls2_output_variant. Otherwise it must
 *                       be an array that can hold width * height float
 *                       values.
 * \param[in] width  Width of the playing field.
 * \param[in] height Height of the playing field.
 */
extern int __attribute__((__nonnull__))
compute_locbased(const int alg, const int em, const int num_threads,
                 const int64_t runs, const float *anchor_x,
                 const float *anchor_y, const int no_anchors,
                 float* results[NUM_VARIANTS], const int width,
                 const int height);

extern void __attribute__((__nonnull__))
ls2_distribute_work_inverted(const int alg, const int em,
			     const int num_threads, const int64_t runs,
                             const long seed, const float tag_x,
                             const float tag_y,
			     const vector2 *restrict anchors, const size_t no_anchors,
			     float *restrict results, const int width, const int height,
			     float *restrict center_x, float *restrict sdev_x,
                             float *restrict center_y, float *restrict sdev_y);

/*!
 * Perform a simulation with a fixed location.
 *
 * \arg[in] alg  The algorithm to use. Use any value of algorithm_t.
 * \arg[in] em   The error model to use. Use any value of error_model_t.
 * \arg[in] num_threads  Number of threads to use.
 * \arg[in] runs      Number of runs per location on the discrete grid.
 * \arg[in] anchor_x  Array of X coordinates of the anchors.
 * \arg[in] anchor_y  Array of Y coordinates of the anchors.
 * \arg[in] no_anchors  Number of anchors
 * \arg{in] tag_x  X coordinate of the tag node.
 * \arg{in] tag_y  Y coordinate of the tag node.
 * \arg[out] results    An array holding the average and maximal errors per
 *                      position.
 * \arg[in] width  Width of the playing field.
 * \arg[in] height Height of the playing field.
 * \arg[out] center_x  X coordinate of the center of mass of all hits.
 * \arg[out] center_y  Y coordinate of the center of mass of all hits.
 */
extern int __attribute__((__nonnull__))
compute_inverse(const int alg, const int em, const int num_threads,
		const int64_t runs, const float *restrict anchor_x,
		const float *restrict anchor_y, const int no_anchors,
		const float tag_x, const float tag_y,
                float *restrict result, const int width, const int height,
		float *restrict center_x, float *restrict sdev_x,
                float *restrict center_y, float *restrict sdev_y);

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
extern void __attribute__((__nonnull__))
ls2_distribute_work_estimator(const int est, const int num_threads,
			      const vector2* anchors, const size_t no_anchors,
			      float *results[NUM_VARIANTS],
			      const int width, const int height);

extern int __attribute__((__nonnull__))
compute_estimates(const int est, const int num_threads,
		  const float *anchor_x, const float *anchor_y,
		  const int no_anchors,
		  float* results[NUM_VARIANTS], const int width,
		  const int height);

/*! Cancel a computation. */
extern int
cancel_running(void);



#endif

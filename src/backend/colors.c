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

extern double ls2_backend_steps;

static inline double ls2_gradation(double value)
{
	if (ls2_backend_steps > 0.0) {
	    double t = round(value * ls2_backend_steps) / ls2_backend_steps;
	    return t;
        } else {
	    return value;
        }
}


/*!
 *
 */
static inline void
__attribute__((__nonnull__,__gnu_inline__,__always_inline__,__const__))
ls2_pick_color_locbased(const float sample, double *r, double *g,
			double *b, double __attribute__((__unused__)) *a)
{
    const float good_color = 50.0F;
    const float bad_color = 250.0F;
    if (isnan(sample)) {
        // Mark not-a-number in magenta.
        *r = 1.0; *g = 0.0; *b = 1.0;
    } else if (sample < good_color) {
        // Use a very good color.
        const double t = ls2_gradation(sample / 50.0);
        *r = t; *g = 1.0; *b = t;
    } else if (sample < bad_color) {
        // Use a good color.
        const double t = ls2_gradation(1.0 - (sample - good_color) / (bad_color - good_color));
        *r = t; *g = t; *b = t;
    } else {
        // Error too large
        *r = 0.0; *g = 0.0; *b = 1.0;
    }
}




/*!
 *
 */
static inline void
__attribute__((__nonnull__,__gnu_inline__,__always_inline__,__const__))
ls2_pick_color_density(const float sample, double *restrict hue,
                       double *restrict saturation, double *restrict lightness)
{
    if (isnan(sample)) {
        // Mark not-a-number in magenta.
        *hue = 300.0;
        *saturation = 1.0;
        *lightness = 0.5;
    } else {
	g_assert(0.0 <= sample && sample <= 1.0);
        const double t = 0.8 * cbrt(1.0 - sample);
        if (sample <= 0.25) {
            *hue = 150.0;
        } else if (0.25 < sample && sample <= 0.5) {
            *hue = 120.0;
        } else if (0.5 < sample && sample <= 0.75) {
            *hue = 90.0;
        } else if (0.75 < sample && sample <= 0.95) {
            *hue = 60.0;
        } else if (0.95 < sample && sample <= 0.98) {
            *hue = 30.0;
        } else if (0.98 < sample && sample <= 0.99) {
            *hue = 0.0;
        } else {
            *hue = 300.0;
        }
        *saturation = 1.0;
        *lightness = t;
    }
}





/*!
 *
 */
static inline void
__attribute__((__nonnull__,__gnu_inline__,__always_inline__,__const__))
ls2_pick_color_diff(const double sample, const double similar,
		    const double dynamic, double *restrict hue,
		    double *restrict saturation, double *restrict lightness)
{
    static const double min_saturation = 0.125;
    static const double max_saturation = 0.875;

    if (sample < -similar && isinf(sample) == 0) {
        /* First was better                    */
	/* The color is green and gets darker. */
	g_assert(sample < 0);
        if (sample < -dynamic / 2.0) {
	    *hue = 150.0;
        } else {
	    *hue = 120.0;
        }
	*saturation = MAX(min_saturation, 0.5 + sample / dynamic);
	g_assert(min_saturation <= *saturation && *saturation <= 0.5);
	*lightness = MAX(0.0, 0.5 + sample / dynamic);
	g_assert(0 <= *lightness && *lightness <= 0.5);
    } else if (-similar <= sample && sample <= similar) {    // Similar
	*hue = 60.0;  // This is the color yellow
	*lightness = 0.5;
	*saturation = 0.5;
    } else if (similar < sample && isinf(sample) == 0) {
        /* Second was better                   */
	/* We start with red and get brighter. */
	g_assert (0 < sample);
	if (sample < dynamic / 2.0 ) {
            *hue = 15.0;
        } else {
            *hue = 0.0;
        }
	*saturation = MIN(0.5 + sample / dynamic, max_saturation);
	g_assert(0.5 <= *saturation && *saturation <= max_saturation);
	*lightness = MIN(0.5 + sample / dynamic, 1.0);
	g_assert(0.5 <= *lightness && *lightness <= 1.0);
    } else if (isnan(sample) != 0) {
        // Mark not-a-number in magenta.
        *hue = 300.0;
        *saturation = 1.0;
        *lightness = 0.5;
    } else {
	g_assert(false);
    }
}





#define INVERTED_THRESHOLD 0.04
#define INVERTED_CUBED_THRESHHOLD (INVERTED_THRESHOLD * INVERTED_THRESHOLD * INVERTED_THRESHOLD)

/*!
 *
 */
static inline void
__attribute__((__nonnull__,__gnu_inline__,__always_inline__,__const__))
ls2_pick_color_inverted(const double sample, double *restrict h,
			double *restrict s, double *restrict l)
{
    *h = 0.0;  /* Red. */
    if (sample > 0.0) {
        // constant used for compression
        static const double alpha = INVERTED_CUBED_THRESHHOLD;
        double t = sample;
        t = (1.0 - alpha) * t + alpha;
        *l = 1.0 - cbrt(t);
    } else {
        *l = 1.0;
    }
    if (sample < 0.97) {
        *s = 0.0;
    } else {
        *s = 0.8;
    }
}

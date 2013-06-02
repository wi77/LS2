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
        const double t = sample / 50.0;
        *r = t; *g = 1.0; *b = t;
    } else if (sample < bad_color) {
        // Use a good color.
        const double t = 1.0 - (sample - good_color) / (bad_color - good_color);
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
ls2_pick_color_diff(const double sample, const double similar,
		    const double dynamic, double *restrict hue,
		    double *restrict saturation, double *restrict lightness)
{
    if (sample < -similar) {                     // First was better
	/* The color is green and gets darker. */
	assert(sample < 0);
	*hue = 120.0;
	*saturation = MAX(0.125, 0.5 + sample / dynamic);
	assert(0.125 <= *saturation && *saturation <= 0.5);
	*lightness = MAX(0.0, 0.5 + sample / dynamic);
	assert(0 <= *lightness && *lightness <= 0.5);
    } else if (-similar <= sample && sample <= similar) {    // Similar
	*hue = 60.0;  // This is the color yellow
	*lightness = 0.5;
	*saturation = 0.5;
    } else if (similar < sample) {              // Second was better
	/* We start with red and get brighter. */
	assert (0 < sample);
	*hue = 0.0;
	*saturation = MIN(0.5 + sample / dynamic, 0.875);
	assert(0.5 <= *saturation && *saturation <= 0.875);
	*lightness = MIN(0.5 + sample / dynamic, 0.875);
	assert(0.5 <= *lightness && *lightness <= 0.875);
    } else {
	assert(false);
    }
}





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
        static const double ___a = 0.04;
        static const double alpha = ___a * ___a * ___a; 
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

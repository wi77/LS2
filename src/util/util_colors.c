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

/********************************************************************
 **
 **  This file is made only for including in the lib_lat project
 **  and not desired for stand alone usage!
 **
 ********************************************************************/

/*******************************************************************
 ***
 *** Selecting colors.
 ***
 *******************************************************************/

static inline void  __attribute__((__always_inline__))
hsl_to_rgb(const double h, const double s, const double l,
	   double *r, double *g, double *b)
{
    const double c = (1.0 - fabs(2.0 * l - 1.0)) * s;
    const double hi = h / 60.0;
    const double x = c * (1.0 - (double)fabs(fmod(floor(hi), 2.0) - 1.0));
    const double m = l - 0.5 * c;

    if (0.0 <= hi && hi < 1.0) {
	*r = c;
	*g = x;
	*b = 0.0;
    } else if (1.0 <= hi && hi < 2.0) {
	*r = x;
	*g = c;
	*b = 0.0;
    } else if (2.0 <= hi && hi < 3.0) {
	*r = 0.0;
	*g = c;
	*b = x;
    } else if (3.0 <= hi && hi < 4.0) {
	*r = 0.0;
	*g = x;
	*b = c;
    } else if (4.0 <= hi && hi < 5.0) {
	*r = x;
	*g = 0.0;
	*b = c;	  
    } else if (5.0 <= hi && hi < 6.0) {
	*r = c;
	*g = 0.0;
	*b = x;	  
    } else {
	*r = 0.0;
	*g = 0.0;
	*b = 0.0;
    }
    *r += m;
    *g += m;
    *b += m;
}

static inline void  __attribute__((__always_inline__))
hsl_to_rgbf(const float h, const float s, const float l,
	    float *r, float *g, float *b)
{
    const float c = (1.0f - (float) fabs(2.0f * l - 1.0f)) * s;
    const float hi = h / 60.0f;
    const float x = c * (1.0f - (float)fabs(fmod(floor(hi), 2.0) - 1.0));
    const float m = l - 0.5f * c;

    if (0.0f <= hi && hi < 1.0f) {
	*r = c;
	*g = x;
	*b = 0.0f;
    } else if (1.0f <= hi && hi < 2.0f) {
	*r = x;
	*g = c;
	*b = 0.0f;
    } else if (2.0f <= hi && hi < 3.0f) {
	*r = 0.0f;
	*g = c;
	*b = x;
    } else if (3.0f <= hi && hi < 4.0f) {
	*r = 0.0f;
	*g = x;
	*b = c;
    } else if (4.0f <= hi && hi < 5.0f) {
	*r = x;
	*g = 0.0f;
	*b = c;	  
    } else if (5.0f <= hi && hi < 6.0f) {
	*r = c;
	*g = 0.0f;
	*b = x;
    } else {
	*r = 0.0f;
	*g = 0.0f;
	*b = 0.0f;
    }
    *r += m;
    *g += m;
    *b += m;
}

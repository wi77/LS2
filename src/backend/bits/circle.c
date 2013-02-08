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

#define DECLARE_DRAW_CIRCLE(name, img_t, ...)				\
     void								\
     name(img_t image, const uint16_t width, const uint16_t height,	\
	  const uint16_t x0, const uint16_t y0, const uint16_t radius,	\
          __VA_ARGS__)							\
     {									\
     /* Use Bresenham's algorithm to plot the circles. */		\
     int x = 0;								\
     int y = radius;							\
     int f = 1 - radius;						\
     int dx = 0;							\
     int dy = -2 * radius;						\
									\
     DRAW_PIXEL(x0, y0 + y);						\
     DRAW_PIXEL(x0, y0 - y);						\
     DRAW_PIXEL(x0 + y, y0);						\
     DRAW_PIXEL(x0 - y, y0);						\
									\
     while (x < y) {							\
	  if (f >= 0) {							\
	       y--;							\
	       dy += 2;							\
	       f += dy;							\
	  }								\
	  x++;								\
	  dx += 2;							\
	  f += dx + 1;							\
									\
	  DRAW_PIXEL(x0 + x, y0 + y);					\
	  DRAW_PIXEL(x0 - x, y0 + y);					\
	  DRAW_PIXEL(x0 - x, y0 - y);					\
	  DRAW_PIXEL(x0 + x, y0 - y);					\
	  DRAW_PIXEL(x0 + y, y0 + x);					\
	  DRAW_PIXEL(x0 - y, y0 + x);					\
	  DRAW_PIXEL(x0 - y, y0 - x);					\
	  DRAW_PIXEL(x0 + y, y0 - x);					\
     }									\
     }

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

#define DECLARE_DRAW_RECTANGLE(name, img_t, ...)			\
  void									\
  name(img_t image, uint16_t width, uint16_t height, const uint16_t x0, \
       const uint16_t y0, const uint16_t xs, const uint16_t ys,		\
       __VA_ARGS__)							\
  {									\
       for (uint16_t y = y0; y < y0 + ys; y++) {			\
            for (uint16_t x = x0; x < x0 + xs; x++) {			\
	         DRAW_PIXEL(x, y);					\
	    }								\
       }								\
  }

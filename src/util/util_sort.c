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

#ifndef UTIL_SORT_C_INCLUDED
#define UTIL_SORT_C_INCLUDED 1

/* An implementation of Comb sort by Wlodzimierz Dobosiewicz (1980)
 *
 * It is a variant opf shell sort and improves bubblesort.
 */
#define SORT_TEMPLATE(T, NAME) \
static inline void \
__attribute__((__always_inline__,__gnu_inline__,__nonnull__,__artificial__)) \
NAME(T *values, const size_t length) \
{ \
    size_t gap = length; \
    bool swapped = false; \
    \
    while ((gap > 1) || swapped) { \
        if (gap > 1) \
            gap = (gap * 13u) / 10u; \
        swapped = false; \
        for (size_t i = 0; i + gap < length; i += gap) { \
            if (values[i] > values[i + gap]) { \
                T t = values[i]; \
                values[i] = values[i + gap]; \
                values[i + gap] = t; \
                swapped = true; \
            } \
        } \
    } \
}

SORT_TEMPLATE(int, isort)

SORT_TEMPLATE(float, fsort)

SORT_TEMPLATE(double, dsort)

#endif

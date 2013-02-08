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

#include <math.h>
#include <float.h>
#include <stdlib.h>

#include "ls2/util.h"

void __attribute__((__nonnull__))
ls2_statistics(const float * restrict values, const size_t size,
	       float * restrict mu, float * restrict sigma,
               float * restrict min, float * restrict max)
{
     float cnt = 1.0F, M_old, M = 0.0F, S = 0.0F;
     *min = values[0];
     *max = values[0];

     for (size_t i = 0; i < size; i++) {
	  const float tmp = values[i];
	  if(isnan(tmp)) continue;
	  M_old = M;
	  M += (tmp - M_old) / cnt;
	  S += (tmp - M_old) * (tmp - M);
	  cnt += 1.0F;
	  *min = (*min < tmp) ? *min : tmp;
	  *max = (*max > tmp) ? *max : tmp;
     }
     *mu = M;
     *sigma = sqrtf(S/(cnt - 1.0f));
}

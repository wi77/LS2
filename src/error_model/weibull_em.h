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
  along with LS². If not, see <http://www.gnu.org/licenses/>.

 */

#ifndef INCLUDED_WEIBULL_H
#define INCLUDED_WEIBULL_H

extern void weibull_setup(const vector2 *vv, size_t num);

#if HAVE_POPT_H
extern struct poptOption weibull_arguments[];
#endif

#if defined(STAND_ALONE)
#  define ERROR_MODEL_NAME "Weibull noise"
#  define ERROR_MODEL_ARGUMENTS weibull_arguments
#endif

#endif

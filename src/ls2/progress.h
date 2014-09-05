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

#ifndef INCLUDED_LS2_PROGRESS_H
#define INCLUDED_LS2_PROGRESS_H

#ifdef __cplusplus
extern "C"
{
#endif




extern void
ls2_initialize_progress_bar(size_t __total, const char *__algorithm);


extern void
ls2_stop_progress_bar(void);



extern void
ls2_reset_progress_bar(size_t total, const char *name);



/*!
 * Returns the progress information.
 *
 * \param[out] threads  The number of currently active threads.
 * \return     Fraction of progress between 0 and 1.
 */
extern double
ls2_get_progress(int *threads);


extern void
ls2_update_progress_bar(size_t value);


#ifdef __cplusplus
}
#endif

#endif

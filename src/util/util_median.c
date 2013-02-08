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

/* Find the k-smallest element in an array.
 *
 * Adapted from Wirth, Niklaus.  * "Algorithms + data structures = programs",
 * Englewood Cliffs: Prentice-Hall, 1976.
 */
static inline float
__attribute__((__always_inline__,__gnu_inline__,__pure__,__nonnull__,__artificial__))
select_s(const size_t length, float const * const values, const size_t k)
{
    if (length <= 0) return 0.0F;
    float a[length];
    memcpy(a, values, length * sizeof(float));

    size_t i, j, l = 0, m = length - 1;
    float x;

    while (l < m) {
        x = a[k];
        i = l;
        j = m;
        do {
            while(a[i] < x) i++;
            while(x < a[j]) j--;
            if (i <= j) {
                float t = a[i];
                a[i] = a[j];
                a[j] = t;
                i++;
                j--;
            }
        } while (i <= j);
        if (j < k) l = i;
        if (k < i) m = j;
    }
    return a[k];
}

/* Compute the media. */
static inline float
__attribute__((__always_inline__,__gnu_inline__,__pure__,__nonnull__,__artificial__))
median_s(const size_t length, float const * const values)
{
    const size_t k = length / 2 - ((length & 1) ? 0 : 1);
    return select_s(length, values, k);
}

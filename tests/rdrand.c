/*
   This file is part of LS² - the Localization Simulation Engine of FU Berlin.
 
   Copyright 2011-2013   Heiko Will, Marcel Kyas, Thomas Hillebrandt,
   Stefan Adler, Malte Rohde, Jonathan Gunthermann
 
   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.
 
   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.
 
   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.
 
*/

#include <immintrin.h>
#include <stdio.h>
#include <pthread.h>


void *run(void * param)
{
    unsigned long long user = 0;
    for (unsigned long long k = 0; k < 1000000000; k++) {
        unsigned long long result;
        int carry;
        do {
            carry = _rdrand64_step(&result);
        } while (carry == 0);
        user += result;
    }
    return (void*) user;
}

int
main(int argc, char **argv)
{
    const int no_threads = 8;
    pthread_t threads[no_threads];

    for (int t = 0; t < no_threads; t++) {
        pthread_create(&(threads[t]), NULL, run, NULL);
    }

    for (int t = 0; t < no_threads; t++) {
        pthread_join(threads[t], NULL);
    }
}

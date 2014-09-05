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


#ifndef _GNU_SOURCE
#  define _GNU_SOURCE
#endif

#if HAVE_CONFIG_H
# include "ls2/ls2-config.h"
#endif

#include <immintrin.h>

#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <string.h>
#include <inttypes.h>
#include <float.h>
#include <math.h>
#include <signal.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <time.h>
#include <unistd.h>

#include <glib.h>

#include "ls2/library.h"
#include "ls2/ls2.h"
#include "ls2/progress.h"
#include "vector_shooter.h"


/*******************************************************************
 *******************************************************************
 ***
 ***   Progress bar.
 ***
 *******************************************************************
 *******************************************************************/

static char const* ls2_display_name;
volatile size_t ls2_progress_total;
static volatile unsigned int ls2_spinner;
static volatile size_t ls2_progress_current;
static volatile size_t ls2_progress_last;
extern volatile size_t ls2_running;
extern volatile size_t ls2_num_threads;


static timer_t timer_id;

#define DEFAULT_WIDTH 80

#define DEFAULT_NAME  24

#define DEFAULT_STEPS 32U

static pthread_mutex_t progress_bar_mutex = PTHREAD_MUTEX_INITIALIZER;

extern void
ls2_update_progress_bar(size_t value)
{
        pthread_mutex_lock(&progress_bar_mutex);
        ls2_progress_current =
                MIN(ls2_progress_current + value, ls2_progress_total);
        pthread_mutex_unlock(&progress_bar_mutex);
}



/*
 *
 */
extern double
ls2_get_progress(int *threads)
{
    double result;
    if (threads != NULL)
        *threads = (int) ls2_running;
    if (ls2_progress_total > 0)
        result = (double) ls2_progress_current / (double) ls2_progress_total;
    else
        result = 0.0;
    return result;
}



/*
 * Handle a progress event and draw the bar to the console.
 */
static void
ls2_handle_progress_bar(int signal __attribute__((__unused__)),
                        siginfo_t *si __attribute__((__unused__)),
                        void *uc __attribute__((__unused__)))
{
    static const char spinner_char[4] = { '|', '/', '-', '\\' };
    char buffer [DEFAULT_WIDTH + 1];
    int pos = 0;

    if (isatty(STDERR_FILENO)) {
        if (ls2_display_name != NULL) {
            strncpy(buffer, ls2_display_name, (size_t)(DEFAULT_NAME - 1));
            pos = (int) strlen(buffer);
        } else {
            pos = 0;
        }
        while (pos < DEFAULT_NAME + 2)
            buffer[pos++] = ' ';
        buffer[pos++] = '|';

        const int ratio =
          (int) ((DEFAULT_STEPS * ls2_progress_current) / ls2_progress_total);
        while (pos < ratio + DEFAULT_NAME + 2) {
            buffer[pos++] = '=';
        }
        if (ratio < (int) DEFAULT_STEPS)
            buffer[pos++] = '>';
        while (pos < DEFAULT_NAME + (int) DEFAULT_STEPS + 2) {
            buffer[pos++] = ' ';
        }

        // Turn the spinner if not finished
        if (ls2_progress_current < ls2_progress_total) {
            buffer[pos++] = spinner_char[ls2_spinner];
            if (ls2_progress_current != ls2_progress_last)
                ls2_spinner = (ls2_spinner + 1U) & 0x3U;
        } else {
            buffer[pos++] = '|';
        }

        float progress =
            ((float) ls2_progress_current) * 100.0f / ((float) ls2_progress_total);
        if (ls2_num_threads < 100) {
            pos += snprintf(buffer + pos, (size_t) (DEFAULT_WIDTH - pos),
                            " %5.1f%% %2zu/%2zu thr.", progress,
                            ls2_running, ls2_num_threads);
        } else {
            pos += snprintf(buffer + pos, (size_t) (DEFAULT_WIDTH - pos),
                            " %5.1f%% %4zu thr.", progress, ls2_running);
        }
        pos = MIN(DEFAULT_WIDTH - 3, pos);
        buffer[pos++] = '\r';
        buffer[pos] = '\0';
        if (write(STDERR_FILENO, buffer, (size_t) pos)) {}
    } else { // Not a tty, just write the percent percentage.
        float progress =
            ((float) ls2_progress_current) * 100.0f / ((float) ls2_progress_total);
        int s = snprintf(buffer, sizeof(buffer), " %5.1f%% %4zu\n",
                         progress, ls2_running);
        if (s > 0) {
            if (write(STDERR_FILENO, buffer, (size_t)s)) {
                // Do nothing.
            }
        }
    }
    fdatasync(STDERR_FILENO);
    ls2_progress_last = ls2_progress_current;
}




/*
 * Register a handler for updating progress.
 */
static void
ls2_setup_progress_handler(void (*handler)(int, siginfo_t *, void*))
{
    struct sigevent sev;
    struct itimerspec its;
    sigset_t mask;
    struct sigaction sa;

    sa.sa_flags = SA_SIGINFO;
    sa.sa_sigaction = handler;
    sigemptyset(&sa.sa_mask);
    if (sigaction(SIGRTMIN, &sa, NULL) == -1) {
        perror("sigaction");
        exit(EXIT_FAILURE);
    }

    sigemptyset(&mask);
    sigaddset(&mask, SIGRTMIN);
    if (sigprocmask(SIG_SETMASK, &mask, NULL) == -1) {
        perror("sigprocmask");
        exit(EXIT_FAILURE);
    }

    sev.sigev_notify = SIGEV_SIGNAL;
    sev.sigev_signo = SIGRTMIN;
    sev.sigev_value.sival_ptr = &timer_id;
    if (timer_create(CLOCK_REALTIME, &sev, &timer_id) == -1) {
        perror("timer_create");
        exit(EXIT_FAILURE);
    }

    /* Start the timer */

    its.it_value.tv_sec = 0;
    its.it_value.tv_nsec = 500000000;
    its.it_interval.tv_sec = its.it_value.tv_sec;
    its.it_interval.tv_nsec = its.it_value.tv_nsec;

    if (timer_settime(timer_id, 0, &its, NULL) == -1) {
        perror("timer_settime");
        exit(EXIT_FAILURE);
    }

    if (sigprocmask(SIG_UNBLOCK, &mask, NULL) == -1) {
        perror("sigprocmask");
        exit(EXIT_FAILURE);
    }
}




void
ls2_reset_progress_bar(size_t total, const char *name)
{
        ls2_spinner          = 0U;
        ls2_progress_current = 0U;
        ls2_progress_last    = 0U;
        ls2_progress_total   = total;
        ls2_display_name     = name;   
}



void
ls2_initialize_progress_bar(size_t total, const char *name)
{
        ls2_reset_progress_bar(total, name);
        ls2_setup_progress_handler(ls2_handle_progress_bar);
        ls2_handle_progress_bar(SIGRTMIN, NULL, NULL);
}


static void
ls2_teardown_progress_handler(void)
{
    timer_delete(timer_id);
    signal(SIGRTMIN, SIG_IGN);
}

void
ls2_stop_progress_bar(void)
{
    ls2_teardown_progress_handler();
    ls2_handle_progress_bar(SIGRTMIN, NULL, NULL);
    if (write(STDERR_FILENO, "\n", 1u)) {}
    fdatasync(STDERR_FILENO);
}

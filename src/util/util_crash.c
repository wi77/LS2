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

/*******************************************************************
 ***
 *** Handle a crash inside the program.
 ***
 *******************************************************************/

#define _GNU_SOURCE 1

#if HAVE_CONFIG_H
#  include "ls2/ls2-config.h"
#endif

#include <errno.h>
#include <pthread.h>
#include <signal.h>
#include <sys/time.h>
#include <sys/types.h>
#include <unistd.h>
#include <execinfo.h>
#include <ucontext.h>
#include <stdio.h>
#include <stdlib.h>

#include "ls2/util.h"


static void
sigsegv_handler(int signal __attribute__((__unused__)),
                siginfo_t * siginfo __attribute__((__unused__)),
                void *context __attribute__((__unused__)))
{
    int ret;
    static const char message1[] =
        "\n" PACKAGE " crashed with signal SEGV\n"
        "Please report this incident to " PACKAGE_BUGREPORT "\n\n"
        "================================= Backtrace: =================================\n\n";
    static const char message2[] =
        "\n=============================== Backtrace End. ===============================\n\n";

    if (write(STDERR_FILENO, message1, sizeof(message1))) {}

    /* Get the backtrace. */
    void *array[12];
    ret = backtrace(array, 12);

#if defined(__i386__)
    array[3] = (void *) (((ucontext_t*) context)->uc_mcontext.gregs[REG_EIP]);
#endif

    backtrace_symbols_fd(array, ret, STDERR_FILENO);
    if (write(STDERR_FILENO, message2, sizeof(message2))) {}
    fdatasync(STDERR_FILENO);

    abort();
}



void
register_sigsegv_handler(void)
{
    struct sigaction sa;
    sa.sa_flags = SA_SIGINFO;
    sa.sa_sigaction = sigsegv_handler;
    sigemptyset(&sa.sa_mask);
    if (sigaction(SIGSEGV, &sa, NULL) == -1) {
        perror("sigaction");
        exit(EXIT_FAILURE);
    }
}

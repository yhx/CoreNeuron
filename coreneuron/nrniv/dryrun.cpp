/*
   Copyright (c) 2018, Blue Brain Project
   All rights reserved.

   Redistribution and use in source and binary forms, with or without modification,
   are permitted provided that the following conditions are met:
   1. Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
   2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
   3. Neither the name of the copyright holder nor the names of its contributors
   may be used to endorse or promote products derived from this software
   without specific prior written permission.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
   AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
   IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
   ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
   LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
   CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
   SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
   INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
   CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
   ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
   THE POSSIBILITY OF SUCH DAMAGE.
   */
#include "coreneuron/nrniv/dryrun.h"
#include "coreneuron/nrnoc/membfunc.h"
#include "nrnoc/multicore.h"
#include <sys/time.h>
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <vector>
static inline double elapsed_in_micro_sec (timeval& start, timeval& stop) {
    return (stop.tv_sec - start.tv_sec)*10.0e6 + (stop.tv_usec - start.tv_usec);
}

static std::vector<double> channel_measures;
static std::vector<int>    channel_counters;

static inline void profile_mechanism_function (mod_f_t f, NrnThread* cur_nt, NrnThreadMembList* tml) {
    timeval start, stop;
    gettimeofday(&start,NULL);
    (*f)(cur_nt, tml->ml, tml->index);
    gettimeofday(&stop,NULL);
    channel_measures[tml->index] += elapsed_in_micro_sec (start, stop);
    channel_counters[tml->index] +=1;

}
extern "C" int execute_dryrun () {
    NrnThreadMembList* tml;
    double elapsed;
    std::fstream mcomplex_file;
    channel_measures.clear();
    channel_counters.clear();
    channel_measures.resize(get_nb_mechs_avail()); 
    channel_counters.resize(get_nb_mechs_avail());
    for (int t_idx = 0 ; t_idx < nrn_nthread; t_idx++) {
        NrnThread* cur_nt = &nrn_threads[t_idx];
        for (tml = cur_nt->tml; tml; tml = tml->next) {
            if (memb_func[tml->index].state) {
                profile_mechanism_function(memb_func[tml->index].state, cur_nt, tml);
            }
        }
        for (tml = cur_nt->tml; tml; tml = tml->next) {
            if (memb_func[tml->index].current) {
                profile_mechanism_function(memb_func[tml->index].current, cur_nt, tml);
            }
        }
        for (tml = cur_nt->tml; tml; tml = tml->next) {
            if (memb_func[tml->index].jacob) {
                profile_mechanism_function(memb_func[tml->index].jacob, cur_nt, tml);
            }
        }
    }
    for (int i = 0 ; i < channel_measures.size(); i++) {
        if (channel_counters[i]) {
            channel_measures[i] /= (double) channel_counters[i];
        }
    }
    mcomplex_file.open("mcomplex.dat", std::ios_base::out);
    for (int i = 2 ; i < channel_measures.size(); i++) {
        if (memb_func[i].sym) {
            mcomplex_file << channel_measures[i] << " " << memb_func[i].sym << std::endl;
        }
    }
    mcomplex_file.close();
    return 0;
}

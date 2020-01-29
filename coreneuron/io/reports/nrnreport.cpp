/*
Copyright (c) 2016, Blue Brain Project
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
#include <iostream>
#include <vector>
#include <algorithm>
#include <map>
#include <set>
#include <cmath>

#include "coreneuron/network/netcon.hpp"
#include "coreneuron/utils/nrn_assert.h"
#include "coreneuron/network/netcvode.hpp"
#include "coreneuron/sim/multicore.hpp"
#include "coreneuron/io/reports/nrnreport.hpp"
#include "coreneuron/io/nrnsection_mapping.hpp"
#include "coreneuron/mechanism/mech_mapping.hpp"
#include "coreneuron/mechanism/membfunc.hpp"
#ifdef ENABLE_REPORTING
#include "bbp/sonata/reports.h"
#include "reportinglib/Records.h"
#endif

namespace coreneuron {

// Size in MB of the report buffer
static int size_report_buffer = 4;

void nrn_flush_reports(double t) {
#ifdef ENABLE_REPORTING
    // flush before buffer is full
    sonata_end_iteration(t);
    records_end_iteration(t);
#endif
}

/** in the current implementation, we call flush during every spike exchange
 *  interval. Hence there should be sufficient buffer to hold all reports
 *  for the duration of mindelay interval. In the below call we specify the
 *  number of timesteps that we have to buffer.
 *  TODO: revisit this because spike exchange can happen few steps before/after
 *  mindelay interval and hence adding two extra timesteps to buffer.
 */
void setup_report_engine(double dt_report, double mindelay) {
#ifdef ENABLE_REPORTING
    /** reportinglib setup */
    int min_steps_to_record = static_cast<int>(std::round(mindelay / dt_report));
    sonata_set_min_steps_to_record(min_steps_to_record);
    sonata_setup_communicators();
    sonata_prepare_datasets();

    records_set_min_steps_to_record(min_steps_to_record);
    records_setup_communicator();
    records_finish_and_share();
#endif  // ENABLE_REPORTING
}

// Size in MB of the report buffers
void set_report_buffer_size(int n) {
    size_report_buffer = n;
#ifdef ENABLE_REPORTING
    sonata_set_max_buffer_size_hint(size_report_buffer);
    records_set_max_buffer_size_hint(size_report_buffer);
#endif
}

void finalize_report() {
#ifdef ENABLE_REPORTING
    sonata_flush(nrn_threads[0]._t);
    records_flush(nrn_threads[0]._t);
#endif
}
}

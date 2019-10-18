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

#include <stdlib.h>
#include <vector>

#include "coreneuron/nrnconf.h"
#include "coreneuron/nrnoc/nrnpthread.h"
#include "coreneuron/nrnoc/multicore.h"
#include "coreneuron/nrnmpi/nrnmpi.h"
#include "coreneuron/nrnoc/nrnoc_decl.h"
#include "coreneuron/nrniv/memory.h"
/*
Now that threads have taken over the actual_v, v_node, etc, it might
be a good time to regularize the method of freeing, allocating, and
updating those arrays. To recapitulate the history, Node used to
be the structure that held direct values for v, area, d, rhs, etc.
That continued to hold for the cray vectorization project which
introduced v_node, v_parent, memb_list. Cache efficiency introduced
actual_v, actual_area, actual_d, etc and the Node started pointing
into those arrays. Additional nodes after allocation required updating
pointers to v and area since those arrays were freed and reallocated.
Now, the threads hold all these arrays and we want to update them
properly under the circumstances of changing topology, changing
number of threads, and changing distribution of cells on threads.
Note there are no longer global versions of any of these arrays.
We do not want to update merely due to a change in area. Recently
we have dealt with diam, area, ri on a section basis. We generally
desire an update just before a simulation when the efficient
structures are necessary. This is reasonably well handled by the
v_structure_change flag which historically freed and reallocated
v_node and v_parent and, just before this comment,
ended up setting the NrnThread tml. This makes most of the old
memb_list vestigial and we now got rid of it except for
the artificial cells (and it is possibly not really necessary there).
Switching between sparse and tree matrix just cause freeing and
reallocation of actual_rhs.

If we can get the freeing, reallocation, and pointer update correct
for _actual_v, I am guessing everything else can be dragged along with
it. We have two major cases, call to pc.nthread and change in
model structure. We want to use Node* as much as possible and defer
the handling of v_structure_change as long as possible.
*/

namespace coreneuron {

int nrn_nthread = 0;
NrnThread* nrn_threads = nullptr;
void (*nrn_mk_transfer_thread_data_)();

static int table_check_cnt_;
static ThreadDatum* table_check_;

int nrn_inthread_;

void nrn_threads_create(int n) {
    int i, j;
    if (nrn_nthread != n) {
        /*printf("sizeof(NrnThread)=%d   sizeof(Memb_list)=%d\n", sizeof(NrnThread),
         * sizeof(Memb_list));*/

        nrn_threads = nullptr;
        nrn_nthread = n;
        if (n > 0) {
            nrn_threads = (NrnThread*)emalloc_align(n * sizeof(NrnThread));

            for (i = 0; i < nrn_nthread; ++i) {
                NrnThread& nt = nrn_threads[i];
                nt._t = 0.;
                nt._dt = -1e9;
                nt.id = i;
                nt._stop_stepping = 0;
                nt.n_vecplay = 0;
                nt._vecplay = nullptr;
                nt.tml = nullptr;
                nt._ml_list = nullptr;
                nt.pntprocs = nullptr;
                nt.presyns = nullptr;
                nt.presyns_helper = nullptr;
                nt.pnt2presyn_ix = nullptr;
                nt.netcons = nullptr;
                nt.weights = nullptr;
                nt.n_pntproc = 0;
                nt.n_presyn = 0;
                nt.n_netcon = 0;
                nt.n_weight = 0;
                nt._ndata = 0;
                nt._nidata = 0;
                nt._nvdata = 0;
                nt._data = nullptr;
                nt._idata = nullptr;
                nt._vdata = nullptr;
                nt.ncell = 0;
                nt.end = 0;
                for (j = 0; j < BEFORE_AFTER_SIZE; ++j) {
                    nt.tbl[j] = nullptr;
                }
                nt._actual_rhs = nullptr;
                nt._actual_d = nullptr;
                nt._actual_a = nullptr;
                nt._actual_b = nullptr;
                nt._actual_v = nullptr;
                nt._actual_area = nullptr;
                nt._actual_diam = nullptr;
                nt._v_parent_index = nullptr;
                nt._permute = nullptr;
                nt._shadow_rhs = nullptr;
                nt._shadow_d = nullptr;
                nt._ecell_memb_list = nullptr;
                nt._sp13mat = nullptr;
                nt._ctime = 0.0;

                nt._net_send_buffer_size = 0;
                nt._net_send_buffer = nullptr;
                nt._net_send_buffer_cnt = 0;
                nt._watch_types = nullptr;
                nt.mapping = nullptr;
                nt.trajec_requests = nullptr;

                nt.nrn_fast_imem = nullptr;
            }
        }
        v_structure_change = 1;
        diam_changed = 1;
    }
    /*printf("nrn_threads_create %d %d\n", nrn_nthread, nrn_thread_parallel_);*/
}

void nrn_threads_free() {
    if (nrn_nthread) {
        free_memory((void*)nrn_threads);
        nrn_threads = nullptr;
        nrn_nthread = 0;
    }
}

void nrn_mk_table_check() {
    if (table_check_) {
        free((void*)table_check_);
        table_check_ = nullptr;
    }

    // Allocate int array of size of mechanism types
    auto ix = std::vector<int>(memb_func.size(), -1);
    table_check_cnt_ = 0;
    for (int id = 0; id < nrn_nthread; ++id) {
        auto& nt = nrn_threads[id];
        for (NrnThreadMembList* tml = nt.tml; tml; tml = tml->next) {
            int index = tml->index;
            if (memb_func[index].thread_table_check_ && ix[index] == -1) {
                ix[index] = id;
                table_check_cnt_ += 2;
            }
        }
    }
    if (table_check_cnt_) {
        table_check_ = (ThreadDatum*)emalloc(table_check_cnt_ * sizeof(ThreadDatum));
    }
    int i = 0;
    for (int id = 0; id < nrn_nthread; ++id) {
        auto& nt = nrn_threads[id];
        NrnThreadMembList* tml;
        for (tml = nt.tml; tml; tml = tml->next) {
            int index = tml->index;
            if (memb_func[index].thread_table_check_ && ix[index] == id) {
                table_check_[i++].i = id;
                table_check_[i++]._pvoid = (void*)tml;
            }
        }
    }
}

void nrn_thread_table_check() {
    for (int i = 0; i < table_check_cnt_; i += 2) {
        auto& nt = nrn_threads[table_check_[i].i];
        auto tml = static_cast<NrnThreadMembList*>(table_check_[i + 1]._pvoid);
        Memb_list* ml = tml->ml;
        (*memb_func[tml->index].thread_table_check_)(0, ml->_nodecount_padded, ml->data, ml->pdata,
                                                     ml->_thread, &nt, tml->index);
    }
}

void nrn_multithread_job(void* (*job)(NrnThread*)) {
    int i;
#if defined(_OPENMP)
// clang-format off
    #pragma omp parallel for private(i) shared(nrn_threads, job, nrn_nthread, \
                                           nrnmpi_myid) schedule(static, 1)
    for (i = 0; i < nrn_nthread; ++i) {
        (*job)(nrn_threads + i);
    }
// clang-format on
#else
    for (i = 1; i < nrn_nthread; ++i) {
        (*job)(nrn_threads + i);
    }
    (*job)(nrn_threads);
#endif
}
}  // namespace coreneuron

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
#ifndef _H_NRNCHECKPOINT_
#define _H_NRNCHECKPOINT_

/*
 * functional protoype of writting checkpoint file in same format of input data
 * */

class NrnThread;
class FileHandler;

void write_checkpoint(NrnThread* nt,
                      int nb_threads,
                      const char* dir,
                      bool swap_bytes_order = false);

void checkpoint_restore_tqueue(NrnThread&, FileHandler&);

int* inverse_permute(int* p, int n);
void nrn_inverse_i_layout(int i, int& icnt, int cnt, int& isz, int sz, int layout);

/* return true if special checkpoint initialization carried out and
   one should not do finitialize
*/
bool checkpoint_initialize();

/** return time to start simulation : if restore_dir provided
 *  then tries to read time.dat file otherwise returns 0
 */
double restore_time(const char* restore_path);

extern int patstimtype;

#ifndef CHKPNTDEBUG
#define CHKPNTDEBUG 1
#endif
#if CHKPNTDEBUG
// Factored out from checkpoint changes to nrnoc/multicore.h and nrnoc/nrnoc_ml.h
// Put here to avoid potential issues with gpu transfer and to allow
// debugging comparison with respect to checkpoint writing to verify that
// data is same as on reading when inverse transforming SoA and permutations.
// Following is a mixture of substantive information which is lost during
// nrn_setup.cpp and debugging only information which is retrievable from
// NrnThread and Memb_list. Ideally, this should all go away

typedef struct Memb_list_ckpnt {
    double* data_not_permuted;  // FIXME temporary store to understand data layout
    Datum* pdata_not_permuted;  // FIXME temporary store to undertsand data layout
    size_t data_offset;

    int* nodeindices;
} Memb_list_chkpnt;

typedef struct NrnThreadChkpnt {
    int* src_gids;     // FIXME temporary struct to store netcon_srcgid from file phase1 (nrn_setup.cpp:278)
    int* output_gids;  // We keep it as current version of coreNeuron dont keep Artificial Gids as output when they appears in Phase1 file
    int nmech;         // Size of linked list tml
    int n_outputgids;  // FIXME temp..
    int ndata_unpadded;        // FIXME temp..
    int* output_vindex;        // FIXME temp..
    double* output_threshold;  // FIXME temp..
    int* pnttype;              // FIXME temp..
    int* pntindex;             // FIXME temp..
    double* delay;             // FIXME temp..
    int npnt;                  // FIXME temp..
    int* icnt;                 // FIXME temp..
    int* dcnt;                 // FIXME temp..
    int* mtype;                // FIXME temp..
    int* vtype;                // FIXME temp..
    int* type;                 // FIXME temp..
    int* vecplay_ix;           // FIXME temp..
    int* vecplay_sz;           // FIXME temp..
    double** vecplay_yvec;     // FIXME temp..
    double** vecplay_tvec;     // FIXME temp..
    Memb_list_chkpnt** mlmap;  // parallel to NrnThread._ml_list
    int file_id;       /* File Id of this NrnThread */

    double* area;
    int* parent;
} NrnThreadChkpnt;

extern NrnThreadChkpnt* nrnthread_chkpnt; // nrn_nthread of these

#endif //CHKPNTDEBUG

#endif

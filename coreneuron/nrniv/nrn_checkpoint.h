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

#include "coreneuron/nrniv/nrn_filehandler.h"

class NrnThread;
void write_checkpoint(NrnThread* nt,
                      int nb_threads,
                      const char* dir,
                      bool swap_bytes_order = false);

void checkpoint_restore_tqueue(NrnThread&, FileHandler&);

template <typename T>
T* chkpnt_soa2aos(T* data, int cnt, int sz, int layout, int* permute) {
  // inverse of F -> data. Just a copy if layout=1. If SoA, original file order depends on
  // padding and permutation.
  // Good for a, b, area, v, diam, Memb_list.data, or anywhere values do not change.
  T* d = new T[cnt * sz];
  if (layout == 1) { /* AoS */
    for (int i=0; i < cnt*sz; ++i) {
      d[i] = data[i];
    }
  }else if (layout == 0) { /* SoA */
    int align_cnt = nrn_soa_padded_size(cnt, layout);
    for (int i=0; i < cnt; ++i) {
      int ip = i;
      if (permute) { ip = permute[i]; }
      for (int j = 0; j < sz; ++j) {
        d[i*sz + j] = data[ip + j*align_cnt];
      }
    }
  }
  return d;
}
     
template <typename T>
void chkpnt_data_write(FileHandler& F, T* data, int cnt, int sz, int layout, int* permute) {
  T* d = chkpnt_soa2aos(data, cnt, sz, layout, permute);
  F.write_array<T>(d, cnt * sz);
  delete [] d;
}

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

#endif

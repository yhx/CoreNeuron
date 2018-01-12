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
#include "coreneuron/nrnoc/multicore.h"
#include "coreneuron/nrniv/nrniv_decl.h"
#include "coreneuron/nrnoc/nrnoc_decl.h"
#include "coreneuron/nrniv/nrn_filehandler.h"
#include "coreneuron/nrniv/nrn_checkpoint.h"
#include "coreneuron/nrniv/netcvode.h"
#include "coreneuron/nrniv/vrecitem.h"
#include "coreneuron/mech/mod2c_core_thread.h"
#include "coreneuron/utils/file_utils.h"
#include "coreneuron/nrniv/node_permute.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <cassert>
#include <stdio.h> // only needed for debugging printf

#define UseFileHandlerWrap 0

#if UseFileHandlerWrap
#include <iomanip>
class FileHandlerWrap {
public:
  FileHandler F;
  std::fstream G;
  FileHandlerWrap(){};

  void open(const char* filename, bool reorder, std::ios::openmode mode = std::ios::in) {
    F.open(filename, reorder, mode);
    std::ostringstream fname;
    fname << filename << ".txt";
    G.open(fname.str().c_str(), mode);
  }

  void close() {
    F.close();
    G.close();
  }

  void checkpoint(int c) {
    F.checkpoint(c);
  }

  template <typename T>
  void write_array(T* p, size_t nb_elements) {
    // G first before chkpnt is incremented
    G << "chkpnt " << F.checkpoint() << std::endl;
    for (size_t i = 0; i < nb_elements; ++i) {
      G << std::setprecision(8) << p[i] << std::endl;
    }
    F.write_array(p, nb_elements); // chkpnt incremented
  }

  template <typename T>
  FileHandlerWrap& operator<<(const T& scalar) {
    F << scalar;
    G << scalar;
    return *this;
  }
};

#else

#define FileHandlerWrap FileHandler

#endif // UseFileHandlerWrap

template <typename T>
void chkpnt_data_write(FileHandlerWrap& F, T* data, int cnt, int sz, int layout, int* permute) {
  T* d = chkpnt_soa2aos(data, cnt, sz, layout, permute);
  F.write_array<T>(d, cnt * sz);
  delete [] d;
}

NrnThreadChkpnt* nrnthread_chkpnt;

int patstimtype;

#ifndef LAYOUT
#define LAYOUT 1
#endif
/*
 * LAYOUT = 0 => SoA, LAYOUT = 1 => AoS
 *
 * */
static const char* output_dir;  // output directory to write simple checkpoint
static bool swap_bytes;

// todo : only keep phase2 as rest (phase1, gan and 3) are constant
static void write_phase2(NrnThread& nt, FileHandlerWrap& file_handle);
static void write_phase3(NrnThread& nt, FileHandlerWrap& file_handle);
static void write_tqueue(NrnThread& nt, FileHandlerWrap& file_handle);
static void write_time(const char *dir);

void write_checkpoint(NrnThread* nt, int nb_threads, const char* dir, bool swap_bytes_order) {
    if (!strlen(dir))
        return;  // empty directory means the option is not enabled
    output_dir = dir;

    // todo : mpi barrier required
    if (nrnmpi_myid == 0) {
        mkdir_p(output_dir);
    }
#if NRNMPI
    nrnmpi_barrier();
#endif

    int i;
    swap_bytes = swap_bytes_order;
    /*
    #if defined(_OPENMP)
      #pragma omp parallel for private(i) shared(nt, nb_threads) schedule(runtime)
    #endif
    */
    FileHandlerWrap f;
    for (i = 0; i < nb_threads; i++) {
      if (nt[i].ncell) {
        write_phase2(nt[i], f);
        write_phase3(nt[i], f);
      }
    }
    if(nrnmpi_myid == 0) {
        write_time(output_dir);
    }
}

static void write_phase2(NrnThread& nt, FileHandlerWrap& fh) {
#if CHKPNTDEBUG
    NrnThreadChkpnt& ntc = nrnthread_chkpnt[nt.id];
#endif
    std::ostringstream filename;
#if CHKPNTDEBUG
    filename << output_dir << "/" << ntc.file_id << "_2.dat";
#endif
    fh.open(filename.str().c_str(), swap_bytes, std::ios::out);
    fh.checkpoint(2);
#if CHKPNTDEBUG
    fh << ntc.n_outputgids << " ngid\n";
#endif
    fh << nt.ncell << " n_real_gid\n";
    fh << nt.end << " nnode\n";
    fh << ((nt._actual_diam == NULL) ? 0 : nt.end) << " ndiam\n";
#if CHKPNTDEBUG
    fh << ntc.nmech << " nmech\n";
#endif

    for (NrnThreadMembList* current_tml = nt.tml; current_tml; current_tml = current_tml->next) {
        if (current_tml->index == patstimtype) { continue; }
        fh << current_tml->index << "\n";
        fh << current_tml->ml->nodecount << "\n";
    }

#if CHKPNTDEBUG
    fh << ntc.ndata_unpadded << " ndata\n";
#endif
    fh << nt._nidata << " nidata\n";
    fh << nt._nvdata << " nvdata\n";
    fh << nt.n_weight << " nweight\n";

    // see comment about parent in node_permute.cpp
    int* pinv_nt = NULL;
    if (nt._permute) {
      int* d = new int[nt.end];
      pinv_nt = inverse_permute(nt._permute, nt.end);
      for (int i=0; i < nt.end; ++i) {
        int x = nt._v_parent_index[nt._permute[i]];
        if (x >= 0) {
          d[i] = pinv_nt[x];
        }else{
          d[i] = 0; // really should be -1;
        }
      }
#if CHKPNTDEBUG
      for (int i=0; i < nt.end; ++i) {
        assert(d[i] == ntc.parent[i]);
      }
#endif
      fh.write_array<int>(d, nt.end);
      delete [] d;
    }else{
#if CHKPNTDEBUG
      for (int i=0; i < nt.end; ++i) {
        assert(nt._v_parent_index[i] == ntc.parent[i]);
      }
#endif
      fh.write_array<int>(nt._v_parent_index, nt.end);
      pinv_nt = new int[nt.end];
      for (int i=0; i < nt.end; ++i) { pinv_nt[i] = i; }
    }

    chkpnt_data_write(fh, nt._actual_a, nt.end, 1, 0, nt._permute);
    chkpnt_data_write(fh, nt._actual_b, nt.end, 1, 0, nt._permute);
#if CHKPNTDEBUG
      for (int i=0; i < nt.end; ++i) {
        assert(nt._actual_area[i] == ntc.area[pinv_nt[i]]);
      }
#endif
    chkpnt_data_write(fh, nt._actual_area, nt.end, 1, 0, nt._permute);
    chkpnt_data_write(fh, nt._actual_v, nt.end, 1, 0, nt._permute);

    if (nt._actual_diam)
        chkpnt_data_write(fh, nt._actual_diam, nt.end, 1, 0, nt._permute);

    // Will need the ml_pinv inverse permutation of ml._permute for ions
    int** ml_pinv = (int**)ecalloc(n_memb_func, sizeof(int*));

    for (NrnThreadMembList* current_tml = nt.tml; current_tml; current_tml = current_tml->next) {
	Memb_list* ml = current_tml->ml;
        int type = current_tml->index;
        if (type == patstimtype) { continue; }
        int cnt = ml->nodecount;
        int sz = nrn_prop_param_size_[type];
        int layout = nrn_mech_data_layout_[type];
        int* semantics = memb_func[type].dparam_semantics;

        if (!nrn_is_artificial_[type]) {
            // ml->nodeindices values are permuted according to nt._permute
            // and locations according to ml._permute
            // i.e. according to comment in node_permute.cpp
            // nodelist[p_m[i]] = p[nodelist_original[i]
            // so pinv[nodelist[p_m[i]] = nodelist_original[i]
            int* nd_ix = new int[cnt];
            for (int i=0; i < cnt; ++i) {
                int ip = ml->_permute ? ml->_permute[i] : i;
                int ipval = ml->nodeindices[ip];
                nd_ix[i] = pinv_nt[ipval];
            }
            fh.write_array<int>(nd_ix, cnt);
            delete[] nd_ix;
        }

        chkpnt_data_write(fh, ml->data, cnt, sz, layout, ml->_permute);

        sz = nrn_prop_dparam_size_[type];
        if (sz) {
            int* d = chkpnt_soa2aos(ml->pdata, cnt, sz, layout, ml->_permute);
            // need to update some values according to Datum semantics.
            if (!nrn_is_artificial_[type]) for (int i_instance = 0; i_instance < cnt; ++i_instance) {
                for (int i = 0; i < sz; ++i) {
                    int ix = i_instance*sz + i;
                    int s = semantics[i];
                    if (s == -1) { // area
                        int p = pinv_nt[d[ix] - (nt._actual_area - nt._data)];
                        d[ix] = p + (5*nt.end); // unpadded area offset
                    }else if (s == -9) { // diam
                        int p = pinv_nt[d[ix] - (nt._actual_diam -nt._data)];

                        d[ix] = p + (6*nt.end);
                    }else if (s == -5) { // Assume pointer to membrane voltage
                        int p = pinv_nt[d[ix] - (nt._actual_v - nt._data)];
                        d[ix] = p + (4*nt.end);
                    }else if (s >= 0 && s < 1000) { // ion
                        // determine ei_instance and ei
                        int etype = s;
                        Memb_list* eml = nt._ml_list[etype];
                        int ecnt = eml->nodecount;
                        int esz = nrn_prop_param_size_[etype];
                        int elayout = nrn_mech_data_layout_[etype];
                        // current index into eml->data is a  function
                        // of elayout, eml._permute, ei_instance, ei, and
                        // eml padding.
                        int p = d[ix] - (eml->data - nt._data);
                        int ei_instance, ei;
                        nrn_inverse_i_layout(p, ei_instance, ecnt, ei, esz, elayout);
                        if (elayout == 0) {
                            if (eml->_permute) {
                                if (!ml_pinv[etype]) {
                                    ml_pinv[etype] = inverse_permute(eml->_permute, eml->nodecount);
                                }
                                ei_instance = ml_pinv[etype][ei_instance];
                            }
                        }
                        int offset = ntc.mlmap[etype]->data_offset;
                        d[ix] = ei_instance*esz + ei + offset;
                    }
#if CHKPNTDEBUG
                    if (s != -8) { // WATCH values change
                      assert(d[ix] == ntc.mlmap[type]->pdata_not_permuted[i_instance*sz + i]);
                    }
#endif
                }
            }
            fh.write_array<int>(d, cnt*sz);
            delete [] d;
        }
    }

    for (int i=0; i < n_memb_func; ++i) {
        if (ml_pinv[i]) { delete [] ml_pinv[i]; }
    }
    free(ml_pinv);
    delete [] pinv_nt;

    int nnetcon = nt.n_netcon - nrn_setup_extracon;
#if CHKPNTDEBUG
    fh.write_array<int>(ntc.output_vindex, nt.n_presyn);
    fh.write_array<double>(ntc.output_threshold, nt.ncell);
    fh.write_array<int>(ntc.pnttype, nnetcon);
    fh.write_array<int>(ntc.pntindex, nnetcon);
#endif
    fh.write_array<double>(nt.weights, nt.n_weight);
#if CHKPNTDEBUG
    fh.write_array<double>(ntc.delay, nnetcon);
#endif
#if CHKPNTDEBUG
    int npnt = ntc.npnt;
#endif
    fh << npnt << " bbcorepointer\n";

    int* iArray = NULL;
    double* dArray = NULL;

    for (int i = 0; i < npnt; i++) {
        fh << ntc.type[i] << "\n";
        fh << ntc.icnt[i] << "\n";
        fh << ntc.dcnt[i] << "\n";

        iArray = new int[ntc.icnt[i]];
        dArray = new double[ntc.dcnt[i]];
        Memb_list* ml = nt._ml_list[ntc.type[i]];
        int dsz = nrn_prop_param_size_[ntc.type[i]];
        int pdsz = nrn_prop_dparam_size_[ntc.type[i]];

        if (nrn_bbcore_write_[ntc.type[i]] && (ntc.icnt[i] || ntc.dcnt[i])) {
            int d_offset = 0;
            int i_offset = 0;
            double* d = ml->data;
            Datum* pd = ml->pdata;
            int layout = nrn_mech_data_layout_[ntc.type[i]];

            for (int j = 0; j < ml->nodecount; j++) {
                int jp = j;

                if (ml->_permute) {
                    jp = ml->_permute[j];
                }

                d = ml->data + nrn_i_layout(jp, ml->nodecount, 0, dsz, layout);
                pd = ml->pdata + nrn_i_layout(jp, ml->nodecount, 0, pdsz, layout);
                int aln_cntml = nrn_soa_padded_size(ml->nodecount, layout);

                // extra parameters after i_offset dont seems to be used
                (*nrn_bbcore_write_[ntc.type[i]])(dArray, iArray, &d_offset, &i_offset, 0, aln_cntml,
                                                 d, pd, ml->_thread, &nt, 0.0);
            }
        }

        /// if there is data from bbcore pointer but bbcore_write is null means
        /// mod2c generated c file hasn't registered a callback
        if((ntc.icnt[i] || ntc.dcnt[i]) && nrn_bbcore_write_[ntc.type[i]] == NULL) {
            std::cerr << " WARNING: bbcore_write not registered for type : " << nrn_get_mechname(ntc.type[i]) <<  std::endl;
        }

        if (ntc.icnt[i]) {
            fh.write_array<int>(iArray, ntc.icnt[i]);
        }

        if (ntc.dcnt[i]) {
            fh.write_array<double>(dArray, ntc.dcnt[i]);
        }
        delete[] iArray;
        delete[] dArray;
    }

    fh << nt.n_vecplay << " VecPlay instances\n";
    for (int i = 0; i < nt.n_vecplay; i++) {
        fh << ntc.vtype[i] << "\n";
        fh << ntc.mtype[i] << "\n";
        fh << ntc.vecplay_ix[i] << "\n";
        fh << ntc.vecplay_sz[i] << "\n";
        fh.write_array<double>(ntc.vecplay_yvec[i], ntc.vecplay_sz[i]);
        fh.write_array<double>(ntc.vecplay_tvec[i], ntc.vecplay_sz[i]);
    }
    write_tqueue(nt, fh);
    fh.close();
}

static void write_phase3(NrnThread&, FileHandlerWrap&) {
}

static void write_time(const char *output_dir) {
    std::ostringstream filename;
    FileHandler f;
    filename << output_dir << "/" << "time.dat";
    f.open(filename.str().c_str(), swap_bytes, std::ios::out);
    f.write_array(&t, 1);
    f.close();
}

/// todo : need to broadcast this rather than all reading a double
double restore_time(const char *restore_dir) {
    double rtime = 0;
    if (strlen(restore_dir)) {
        std::ostringstream filename;
        FileHandler f;

        filename << restore_dir << "/" << "time.dat";
        f.open(filename.str().c_str(), swap_bytes, std::ios::in);
        f.read_array(&rtime, 1);
        f.close();
    }
    return rtime;
}

static void write_tqueue(TQItem* q, NrnThread& nt, FileHandlerWrap& fh) {
    DiscreteEvent* d = (DiscreteEvent*)q->data_;
    //printf("  p %.20g %d\n", q->t_, d->type());
    //d->pr("", q->t_, net_cvode_instance);

    fh << d->type() << "\n";
    fh.write_array(&q->t_, 1);

    switch (d->type()) {
      case NetConType: {
        NetCon* nc = (NetCon*)d;
        assert (nc >= nt.netcons && (nc < (nt.netcons + nt.n_netcon)));
        fh << (nc - nt.netcons) << "\n";
        break;
      }
      case SelfEventType: {
        SelfEvent* se = (SelfEvent*)d;
        fh << int(se->target_->_type) << "\n";
        fh << se->target_ - nt.pntprocs << "\n"; // index of nrnthread.pntprocs
        fh << se->target_->_i_instance << "\n"; // not needed except for assert check
        fh.write_array(&se->flag_, 1);
        fh << (se->movable_ - nt._vdata) << "\n"; // DANGEROUS?
        fh << se->weight_index_ << "\n";
        //printf("    %d %ld %d %g %ld %d\n", se->target_->_type, se->target_ - nt.pntprocs, se->target_->_i_instance, se->flag_, se->movable_ - nt._vdata, se->weight_index_);
        break;
      }
      case PreSynType: {
        PreSyn* ps = (PreSyn*)d;
        assert(ps >= nt.presyns && (ps < (nt.presyns + nt.n_presyn)));
        fh << (ps - nt.presyns) << "\n";
        break;
      }
      case NetParEventType: {
        // nothing extra to write
        break;
      }
      case PlayRecordEventType: {
        PlayRecord* pr = ((PlayRecordEvent*)d)->plr_;
        fh << pr->type() << "\n";
	if (pr->type() == VecPlayContinuousType) {
          VecPlayContinuous* vpc = (VecPlayContinuous*)pr;
          int ix = -1;
          for (int i = 0; i < nt.n_vecplay; ++i) {
            // if too many for fast search, put ix in the instance
            if (nt._vecplay[i] == (void*)vpc) {
              ix = i;
              break;
            }
          }
          assert(ix >= 0);
          fh << ix << "\n";
        }else{
          assert(0);
        }
        break;
      }
      default: {
        // In particular, InputPreSyn does not appear in tqueue as it
        // immediately fans out to NetCon.
        assert(0);
        break;
      }
    }    
}

static int patstim_index;
static double patstim_te;

static void checkpoint_restore_tqitem(int type, NrnThread& nt, FileHandler& fh) {
    double te;
    fh.read_array(&te, 1);
    //printf("restore tqitem type=%d te=%.20g\n", type, te);

    switch (type) {
      case NetConType: {
        int ncindex = fh.read_int();
        //printf("  NetCon %d\n", ncindex);
        NetCon* nc = nt.netcons + ncindex;
        nc->send(te, net_cvode_instance, &nt);
        break;
      }
      case SelfEventType: {
        int target_type = fh.read_int(); // not really needed (except for assert below)
        int pinstance = fh.read_int();
        int target_instance = fh.read_int();
        double flag;
        fh.read_array(&flag, 1);
        int movable = fh.read_int();
        int weight_index = fh.read_int();
        if (target_type == patstimtype) {
          if (nt.id == 0) {
            patstim_te = te;
          }
          break;
        }
        Point_process* pnt = nt.pntprocs + pinstance;
        //printf("  SelfEvent %d %d %d %g %d %d\n", target_type, pinstance, target_instance, flag, movable, weight_index);
        assert(target_instance == pnt->_i_instance);
        assert(target_type == pnt->_type);
        net_send(nt._vdata + movable, weight_index, pnt, te, flag);
        break;
      }
      case PreSynType: {
        int psindex = fh.read_int();
        //printf("  PreSyn %d\n", psindex);
        PreSyn* ps = nt.presyns + psindex;
        int gid = ps->output_index_;
        ps->output_index_ = -1;
        ps->send(te, net_cvode_instance, &nt);
        ps->output_index_ = gid;
        break;
      }
      case NetParEventType: {
        // nothing extra to read
        // printf("  NetParEvent\n");
        break;
      }
      case PlayRecordEventType: {
        int prtype = fh.read_int();
        if (prtype == VecPlayContinuousType) {
          VecPlayContinuous* vpc = (VecPlayContinuous*)(nt._vecplay[fh.read_int()]);
          vpc->e_->send(te, net_cvode_instance, &nt);
        }else{
          assert(0);
        }
        break;
      }
      default: {
        assert(0);
        break;
      }
    }    
}

extern "C" {
extern int checkpoint_save_patternstim(_threadargsproto_);
extern void checkpoint_restore_patternstim(int, double, _threadargsproto_);
}

static void write_tqueue(NrnThread& nt, FileHandlerWrap& fh) {
    // VecPlayContinuous
    fh << nt.n_vecplay << " VecPlayContinuous state\n";
    for (int i=0; i < nt.n_vecplay; ++i) {
      VecPlayContinuous* vpc = (VecPlayContinuous*)nt._vecplay[i];
      fh << vpc->last_index_ << "\n";
      fh << vpc->discon_index_ << "\n";
      fh << vpc->ubound_index_ << "\n";
    }

    // PatternStim
    int patstim_index = -1;
    for (NrnThreadMembList* tml = nrn_threads[0].tml; tml; tml = tml->next) {
      if (tml->index == patstimtype) {
        Memb_list* ml = tml->ml;
        patstim_index = checkpoint_save_patternstim(
          /* below correct only for AoS */
          0, ml->nodecount, ml->data, ml->pdata, ml->_thread, nrn_threads, 0.0);
        break;
      }
    }
    fh << patstim_index << " PatternStim\n";

    // Avoid extra spikes due to some presyn voltages above threshold
    fh << -1 << " Presyn ConditionEvent flags\n";
    for (int i=0; i < nt.n_presyn; ++i) {
      // PreSyn.flag_ not used. HPC memory utilizes PreSynHelper.flag_ array
      fh << nt.presyns_helper[i].flag_ << "\n";
    }

    NetCvodeThreadData& ntd = net_cvode_instance->p[nt.id];
    // printf("write_tqueue %d %p\n", nt.id, ndt.tqe_);
    TQueue<QTYPE>* tqe = ntd.tqe_;
    TQItem* q;

    fh << -1 << " TQItems from atomic_dq\n";
    while((q = tqe->atomic_dq(1e20)) != NULL) {
      write_tqueue(q, nt, fh);
    }
    fh << 0 << "\n";

    fh << -1 << " TQItemsfrom binq_\n";
    for (q = tqe->binq_->first(); q; q = tqe->binq_->next(q)) {
      write_tqueue(q, nt, fh);
    }
    fh << 0 << "\n";
}

static bool checkpoint_restored_ = false;

void checkpoint_restore_tqueue(NrnThread& nt, FileHandler& fh) {
    int type;
    checkpoint_restored_ = true;

    // VecPlayContinuous
    assert(fh.read_int() == nt.n_vecplay); // VecPlayContinuous state
    for (int i=0; i < nt.n_vecplay; ++i) {
      VecPlayContinuous* vpc = (VecPlayContinuous*)nt._vecplay[i];
      vpc->last_index_ = fh.read_int();
      vpc->discon_index_ = fh.read_int();
      vpc->ubound_index_ = fh.read_int();
    }

    // PatternStim
    patstim_index = fh.read_int(); // PatternStim
    if (nt.id == 0) {
      patstim_te = -1.0; // changed if relevant SelfEvent;
    }

    assert(fh.read_int() == -1); // -1 PreSyn ConditionEvent flags
    for (int i=0; i < nt.n_presyn; ++i) {
      nt.presyns_helper[i].flag_ = fh.read_int();
    }

    assert(fh.read_int() == -1); // -1 TQItems from atomic_dq
    while ((type = fh.read_int()) != 0) {
      checkpoint_restore_tqitem(type, nt, fh);
    }

    assert(fh.read_int() == -1); // -1 TQItems from binq_
    while ((type - fh.read_int()) != 0) {
      checkpoint_restore_tqitem(type, nt, fh);
    }
}

// A call to finitialize must be avoided after restoring the checkpoint
// as that would change all states to a voltage clamp initialization.
// Nevertheless t and some spike exchange and other computer state needs to
// be initialized.
// Consult finitialize.c to help decide what should be here
bool checkpoint_initialize() {
    dt2thread(-1.);
    nrn_thread_table_check();
    nrn_spike_exchange_init();

    // if PatternStim exists, needs initialization
    for (NrnThreadMembList* tml = nrn_threads[0].tml; tml; tml = tml->next) {
      if (tml->index == patstimtype && patstim_index >= 0 && patstim_te > 0.0) {
        Memb_list* ml = tml->ml;
        checkpoint_restore_patternstim(patstim_index, patstim_te,
          /* below correct only for AoS */
          0, ml->nodecount, ml->data, ml->pdata, ml->_thread, nrn_threads, 0.0);
        break;
      }
    }
    return checkpoint_restored_;
}

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
#include "coreneuron/nrniv/nrn_checkpoint.h"
#include "coreneuron/nrnoc/multicore.h"
#include "coreneuron/nrniv/nrniv_decl.h"
#include "coreneuron/nrnoc/nrnoc_decl.h"
#include "coreneuron/nrniv/nrn_filehandler.h"
#include "coreneuron/nrniv/netcvode.h"
#include "coreneuron/nrniv/vrecitem.h"
#include "coreneuron/mech/mod2c_core_thread.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <cassert>
#include <stdio.h> // only needed for debugging printf

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
static void write_phase2(NrnThread& nt, FileHandler& file_handle);
static void write_phase3(NrnThread& nt, FileHandler& file_handle);
static void write_tqueue(NrnThread& nt, FileHandler& file_handle);
static void write_time(const char *dir);

void write_checkpoint(NrnThread* nt, int nb_threads, const char* dir, bool swap_bytes_order) {
    if (!strlen(dir))
        return;  // empty directory means the option is not enabled
    output_dir = dir;
    int i;
    swap_bytes = swap_bytes_order;
    /*
    #if defined(_OPENMP)
      #pragma omp parallel for private(i) shared(nt, nb_threads) schedule(runtime)
    #endif
    */
    FileHandler f;
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

static void write_phase2(NrnThread& nt, FileHandler& file_handle) {
    std::ostringstream filename;
    filename << output_dir << "/" << nt.file_id << "_2.dat";
    file_handle.open(filename.str().c_str(), swap_bytes, std::ios::out);
    file_handle.checkpoint(2);
    file_handle << nt.n_outputgids << " ngid\n";
    file_handle << nt.ncell << " n_real_gid\n";
    file_handle << nt.end << " nnode\n";
    file_handle << ((nt._actual_diam == NULL) ? 0 : nt.end) << " ndiam\n";
    file_handle << nt.nmech << " nmech\n";

    for (NrnThreadMembList* current_tml = nt.tml; current_tml; current_tml = current_tml->next) {
        if (current_tml->index == patstimtype) { continue; }
        file_handle << current_tml->index << "\n";
        file_handle << current_tml->ml->nodecount << "\n";
    }

    file_handle << nt.ndata_unpadded << " ndata\n";
    file_handle << nt._nidata << " nidata\n";
    file_handle << nt._nvdata << " nvdata\n";
    file_handle << nt.n_weight << " nweight\n";
    file_handle.write_array<int>(nt._v_parent_index, nt.end);
    file_handle.write_array<double>(nt._actual_a, nt.end);
    file_handle.write_array<double>(nt._actual_b, nt.end);
    file_handle.write_array<double>(nt._actual_area, nt.end);
    file_handle.write_array<double>(nt._actual_v, nt.end);

    if (nt._actual_diam)
        file_handle.write_array<double>(nt._actual_diam, nt.end);

    for (NrnThreadMembList* current_tml = nt.tml; current_tml; current_tml = current_tml->next) {
        int type = current_tml->index;
        if (type == patstimtype) { continue; }
        int nb_nodes = current_tml->ml->nodecount;
        int size_of_line_data = nrn_soa_padded_size(nb_nodes, nrn_mech_data_layout_[type]);

        if (!nrn_is_artificial_[type]) {
            file_handle.write_array<int>(current_tml->ml->nodeindices, nb_nodes);
        }

        // if LAYOUT is SoA: we need to transpose to the structure to write in file format order
        file_handle.write_array<double>(current_tml->ml->data, nb_nodes, size_of_line_data,
                                        nrn_prop_param_size_[type], !LAYOUT);

        if (nrn_prop_dparam_size_[type]) {
            // if LAYOUT is SoA: we need to transpose to the structure to write in file format order
            file_handle.write_array<int>(current_tml->ml->pdata, nb_nodes,
                                         size_of_line_data, nrn_prop_dparam_size_[type], !LAYOUT);
        }
    }

    int nnetcon = nt.n_netcon - nrn_setup_extracon;
    file_handle.write_array<int>(nt.output_vindex, nt.n_presyn);
    file_handle.write_array<double>(nt.output_threshold, nt.ncell);
    file_handle.write_array<int>(nt.pnttype, nnetcon);
    file_handle.write_array<int>(nt.pntindex, nnetcon);
    file_handle.write_array<double>(nt.weights, nt.n_weight);
    file_handle.write_array<double>(nt.delay, nnetcon);
    file_handle << nt.npnt << " bbcorepointer\n";

    int* iArray = NULL;
    double* dArray = NULL;

    for (int i = 0; i < nt.npnt; i++) {
        file_handle << nt.type[i] << "\n";
        file_handle << nt.icnt[i] << "\n";
        file_handle << nt.dcnt[i] << "\n";

        iArray = new int[nt.icnt[i]];
        dArray = new double[nt.dcnt[i]];
        Memb_list* ml = nt.mlmap[nt.type[i]];
        int dsz = nrn_prop_param_size_[nt.type[i]];
        int pdsz = nrn_prop_dparam_size_[nt.type[i]];

        if (nrn_bbcore_write_[nt.type[i]] && (nt.icnt[i] || nt.dcnt[i])) {
            int d_offset = 0;
            int i_offset = 0;
            double* d = ml->data;
            Datum* pd = ml->pdata;
            int layout = nrn_mech_data_layout_[nt.type[i]];

            for (int j = 0; j < ml->nodecount; j++) {
                int jp = j;

                if (ml->_permute) {
                    jp = ml->_permute[j];
                }

                d = ml->data + nrn_i_layout(jp, ml->nodecount, 0, dsz, layout);
                pd = ml->pdata + nrn_i_layout(jp, ml->nodecount, 0, pdsz, layout);
                int aln_cntml = nrn_soa_padded_size(ml->nodecount, layout);

                // extra parameters after i_offset dont seems to be used
                (*nrn_bbcore_write_[nt.type[i]])(dArray, iArray, &d_offset, &i_offset, 0, aln_cntml,
                                                 d, pd, ml->_thread, &nt, 0.0);
            }
        }

        /// if there is data from bbcore pointer but bbcore_write is null means
        /// mod2c generated c file hasn't registered a callback
        if((nt.icnt[i] || nt.dcnt[i]) && nrn_bbcore_write_[nt.type[i]] == NULL) {
            std::cerr << " WARNING: bbcore_write not registered for type : " << nrn_get_mechname(nt.type[i]) <<  std::endl;
        }

        if (nt.icnt[i]) {
            file_handle.write_array<int>(iArray, nt.icnt[i]);
        }

        if (nt.dcnt[i]) {
            file_handle.write_array<double>(dArray, nt.dcnt[i]);
        }
        delete[] iArray;
        delete[] dArray;
    }

    file_handle << nt.n_vecplay << " VecPlay instances\n";
    for (int i = 0; i < nt.n_vecplay; i++) {
        file_handle << nt.vtype[i] << "\n";
        file_handle << nt.mtype[i] << "\n";
        file_handle << nt.vecplay_ix[i] << "\n";
        file_handle << nt.vecplay_sz[i] << "\n";
        file_handle.write_array<double>(nt.vecplay_yvec[i], nt.vecplay_sz[i]);
        file_handle.write_array<double>(nt.vecplay_tvec[i], nt.vecplay_sz[i]);
    }
    write_tqueue(nt, file_handle);
    file_handle.close();
}

static void write_phase3(NrnThread&, FileHandler&) {
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

static void write_tqueue(TQItem* q, NrnThread& nt, FileHandler& fh) {
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

static void write_tqueue(NrnThread& nt, FileHandler& fh) {
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

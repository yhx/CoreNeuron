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

#include "coreneuron/nrniv/nrn_file_reader.h"
#include <sstream>
namespace coreneuron {
    FileReader::FileReader (const char* path_w, const char* restore_path_w, int* gidgroups, bool reoder):
        path_w (path_w), restore_path_w(restore_path_w), gidgroups_w(gidgroups), reorder(reorder) {

    }
     FileReader::~FileReader(){

    }

     void 
        FileReader::mkmech_info(std::ostream& ){

        }

     void*
        FileReader::get_global_dbl_item(void*, const char*& name, int& size, double*& val){

        }

     int 
        FileReader::get_global_int_item(const char* name){

        }

     void 
        FileReader::get_partrans_setup_info(int tid, int& ntar, int& nsrc,
                int& type, int& ix_vpre, int*& sid_target, int*& sid_src, int*& v_indices){

        }
     int 
        FileReader::get_dat1_(int tid, int& n_presyn, int& n_netcon,
                int*& output_gid, int*& netcon_srcgid, int nrn_setup_extracon){
            std::ostringstream filename;
            filename << path_w << "/" << gidgroups_w[tid] <<"_1.dat";
            F.open (filename.str().c_str(), reorder);
            n_presyn  = F.read_int();
            n_netcon = F.read_int();
            output_gid    = new int [n_presyn];
            netcon_srcgid = new int [n_netcon + nrn_setup_extracon];
            F.read_array<int>(output_gid, n_presyn);
            F.read_array<int>(netcon_srcgid, n_netcon);
            F.close();
        }

     int 
        FileReader::get_dat2_1(int tid, int& ngid, int& n_real_gid, int& nnode, int& ndiam,
                int& nmech, int*& tml_index, int*& ml_nodecount, int& nidata, int& nvdata, int& nweight){
            std::ostringstream filename;
            filename << restore_path_w << "/" << gidgroups_w[tid] <<"_2.dat";
            F.open (filename.str().c_str(), reorder);
            int n_outputgid = F.read_int();  // is it used ? shall we check it is == tid ?
            ngid        = F.read_int();
            n_real_gid  = F.read_int();
            ndiam       = F.read_int();  // 0 if not needed, else nt->end
            nmech       = F.read_int();
            tml_index   = new int[nmech];
            ml_nodecount = new int[nmech];
            for (int i=0; i < nmech; ++i) {
                tml_index[i] = F.read_int();
                ml_nodecount[i] = F.read_int();
            }
            nidata   = F.read_int();
            nvdata   = F.read_int();
            nweight = F.read_int();
        }

     int 
        FileReader::get_dat2_2(int tid, int*& v_parent_index, double*& a, double*& b, //FIXME needs to pass ndiam and nt_end
                double*& area, double*& v, double*& diamvec, int ndiam, int nt_end){
            F.read_array<int>(v_parent_index, nt_end);
            F.read_array<double>(a, nt_end);
            F.read_array<double>(b, nt_end);
            F.read_array<double>(area, nt_end);
            F.read_array<double>(v, nt_end);
            if (ndiam) {
                F.read_array<double>(diamvec, nt_end);
            }
        }

     int 
        FileReader::get_dat2_mech(int tid, size_t i, int dsz_inst, int*& nodeindices, //FIXME needs to pass ml_nodecount & is_art
                double*& data, int*& pdata, int ml_nodecount, int is_art){
            if (is_art) {
                F.read_array<int>(nodeindices, ml_nodecount);
            }
        }

     int 
        FileReader::get_dat2_3(int tid, int nweight, int*& output_vindex, double*& output_threshold,
                int*& netcon_pnttype, int*& netcon_pntindex, double*& weights, double*& delays, int nt_n_presyn){ //FIXME needs to pass nt_n_presyn
                output_vindex = F.read_array<int>(nt_n_presyn);
        }

     int 
        FileReader::get_dat2_corepointer(int tid, int& n){
            n = F.read_int();
        }

     int 
        FileReader::get_dat2_corepointer_mech(int tid, int type,
                int& icnt, int& dcnt, int*& iarray, double*& darray){
            type = F.read_int();
            icnt = F.read_int();
            dcnt = F.read_int();
            if (icnt) {
                iarray = F.read_array<int>(iarray, icnt);
            }
            if (dcnt) {
                darray = F.read_array<double>(darray, dcnt);
            }
        }

     int 
        FileReader::get_dat2_vecplay(int tid, int& n){
            n = F.read_int();
        }

     int 
        FileReader::get_dat2_vecplay_inst(int tid, int i, int& vptype, int& mtype,
                int& ix, int& sz, double*& yvec, double*& tvec){
            vptype = F.read_int();
            mtype = F.read_int();
            ix    = F.read_int();
            sz    = F.read_int();
            F.read_array<double>(yvec, sz);
            F.read_array<double>(tvec, sz);
        }
} //namespace coreneuron


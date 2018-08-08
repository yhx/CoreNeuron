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
#ifndef nrn_neuron_reader_h
#define nrn_neuron_reader_h

#include <iostream>
#include <fstream>
#include <vector>
#include "coreneuron/utils/endianness.h"
#include "coreneuron/utils/swap_endian.h"
#include "coreneuron/nrniv/nrn_assert.h"
namespace coreneuron {


    /**
     *  Base class to interface inputs reading for nrn_setup
     */
    class NeuronReader {
        public:
            NeuronReader() {};
            virtual ~NeuronReader() = 0;

            virtual void 
                mkmech_info(std::ostream& ) = 0;

            virtual void*
                get_global_dbl_item(void*, const char*& name, int& size, double*& val) = 0;

            virtual int 
                get_global_int_item(const char* name) = 0;

            virtual void 
                get_partrans_setup_info(int tid, int& ntar, int& nsrc,
                        int& type, int& ix_vpre, int*& sid_target, int*& sid_src, int*& v_indices) = 0;
            virtual int 
                get_dat1_(int tid, int& n_presyn, int& n_netcon,
                        int*& output_gid, int*& netcon_srcgid) = 0;

            virtual int get_dat2_1(int tid, int& ngid, int& n_real_gid, int& nnode, int& ndiam,
                    int& nmech, int*& tml_index, int*& ml_nodecount, int& nidata, int& nvdata, int& nweight) = 0;


            virtual int 
                get_dat2_2(int tid, int*& v_parent_index, double*& a, double*& b,
                        double*& area, double*& v, double*& diamvec)=0;

            virtual int 
                get_dat2_mech(int tid, size_t i, int dsz_inst, int*& nodeindices,
                        double*& data, int*& pdata)=0;

            virtual int 
                get_dat2_3(int tid, int nweight, int*& output_vindex, double*& output_threshold,
                        int*& netcon_pnttype, int*& netcon_pntindex, double*& weights, double*& delays)=0;

            virtual int 
                get_dat2_corepointer(int tid, int& n)=0;

            virtual int 
                get_dat2_corepointer_mech(int tid, int type,
                        int& icnt, int& dcnt, int*& iarray, double*& darray)=0;

            virtual int 
                get_dat2_vecplay(int tid, int& n)=0;

            virtual int 
                get_dat2_vecplay_inst(int tid, int i, int& vptype, int& mtype,
                        int& ix, int& sz, double*& yvec, double*& tvec)=0;
    };

} // namespace coreneuron
#endif // nrn_neuron_reader_h

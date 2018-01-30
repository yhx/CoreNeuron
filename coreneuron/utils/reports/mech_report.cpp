/*
Copyright (c) 2017, Blue Brain Project
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

#include <string.h>
#include <sstream>
#include <iostream>
#include <stdexcept>

#include "coreneuron/nrnconf.h"
#include "coreneuron/nrniv/nrniv_decl.h"
#include "coreneuron/nrnmpi/nrnmpi.h"
#include "coreneuron/nrnmpi/nrnmpi_impl.h"
#include "coreneuron/nrnmpi/nrnmpidec.h"
#include "coreneuron/nrnoc/multicore.h"

#if NRNMPI

void print_mech_count(std::vector<unsigned long long>& mech_count) {
    std::vector<unsigned long long> global_mech_count(n_memb_func);
    MPI_Allreduce(&mech_count[0], &global_mech_count[0], mech_count.size(), MPI_UNSIGNED_LONG_LONG, MPI_SUM, MPI_COMM_WORLD);
    if(nrnmpi_myid == 0) {
        printf("\n Mechanism counts by type : \n");
        for(int i = 0; i < mech_count.size(); i++) {
            printf("Id, name, count :: %4d. %20s %10lld\n", i, nrn_get_mechname(i), global_mech_count[i]);
        }
    }
}

/** Write mechanism count for every gid using MPI I/O */
void write_mech_report(std::string path, int ngroups, int* gidgroups) {
	std::string fname(path);
    fname += "/mech.stats";

    /// remove if file already exist
    if(nrnmpi_myid == 0) {
        remove(fname.c_str());
    }
    nrnmpi_barrier();

    /// count no of gids/threads to allocate buffer
    int non_empty_thread_count = 0;
    for(int i=0; i<nrn_nthread; i++) {
        NrnThread& nt = nrn_threads[i];
        if(nt.ncell) {
            non_empty_thread_count++;
        }
    }


    // each column is 64 chars chars (gid followed by all mechanism name)
    const int RECORD_LEN = 64;
    unsigned num_records = 1 + n_memb_func;
    int header = 0;

    /// header is written by first rank only
    if(nrnmpi_myid == 0) {
        header = 1;
    }

    /// allocate memoty
    unsigned num_bytes = (sizeof(char) * num_records * RECORD_LEN * (non_empty_thread_count+header));
    char *record_data = (char*) malloc(num_bytes);
    if(record_data == NULL) {
        throw std::runtime_error("Memory allocation error while writing mech stats");
    }

    strcpy(record_data, "");
    char record_entry[RECORD_LEN];

    /// prepare first row as header
    if(nrnmpi_myid == 0) {
        strcat(record_data, "gid,");
        for(int i=0; i<n_memb_func; i++) {
            const char*name = nrn_get_mechname(i);
            if(name == NULL)
                name = "";
            snprintf(record_entry, RECORD_LEN, "%s", name);
            strcat(record_data, record_entry);
            if((i+1) < n_memb_func) {
                strcat(record_data, ",");
            }
        }
        strcat(record_data, "\n");
    }

    std::vector<unsigned long long> total_mech_count(n_memb_func, 0);

    /// each gid record goes on separate row, only check non-empty threads
	for(int i=0; i<non_empty_thread_count; i++) {
		NrnThread& nt = nrn_threads[i];
		std::vector<int> mech_count(n_memb_func, 0);
		NrnThreadMembList* tml;
		for (tml = nt.tml; tml; tml = tml->next) {
			int type = tml->index;
			Memb_list* ml = tml->ml;
			mech_count[type] = ml->nodecount;
            total_mech_count[type] +=  ml->nodecount;
		}

        /// copy to buffer
		snprintf(record_entry, RECORD_LEN, "%d,", gidgroups[i]);
		strcat(record_data, record_entry);
		for(int i=0; i<n_memb_func; i++) {
			snprintf(record_entry, RECORD_LEN, "%d", mech_count[i]);
			strcat(record_data, record_entry);
            if((i+1) < n_memb_func) {
                strcat(record_data, ",");
            }
		}
		strcat(record_data, "\n");
	}

    /// print global stats
    print_mech_count(total_mech_count);

    // calculate offset into global file. note that we don't write
    // all num_bytes but only "populated" buffer
    unsigned long num_chars = strlen(record_data);
    unsigned long offset = 0;

    // global offset into file
    MPI_Exscan(&num_chars, &offset, 1, MPI_UNSIGNED_LONG, MPI_SUM, MPI_COMM_WORLD);

    // write to file using parallel mpi i/o
    MPI_File fh;
    MPI_Status status;

    // ibm mpi (bg-q) expects char* instead of const char* (even though it's standard)
    int op_status = MPI_File_open(MPI_COMM_WORLD, (char*) fname.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);
    if(op_status != MPI_SUCCESS && nrnmpi_myid == 0) {
        throw std::runtime_error("Error while opening mech stat file " + fname);
    }

    op_status = MPI_File_write_at(fh, offset, record_data, num_chars, MPI_BYTE, &status);
    if(op_status != MPI_SUCCESS && nrnmpi_myid == 0) {
        throw std::runtime_error("Error while writing mech stats");
    }
    MPI_File_close(&fh);
}

#else

void output_spikes_parallel(const char* outpath) {
    std::cout << "Build and run with MPI to get mechanism stats \n";
}

#endif


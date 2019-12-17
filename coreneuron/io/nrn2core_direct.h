#ifndef nrn2core_direct_h
#define nrn2core_direct_h

#include <iostream>

extern "C" {
// The callbacks into nrn/src/nrniv/nrnbbcore_write.cpp to get
// data directly instead of via files.

extern int corenrn_embedded;
extern int corenrn_embedded_nthread;

extern void (*nrn2core_group_ids_)(int*);

extern void (*nrn2core_mkmech_info_)(std::ostream&);

extern void* (*nrn2core_get_global_dbl_item_)(void*, const char*& name, int& size, double*& val);
extern int (*nrn2core_get_global_int_item_)(const char* name);

extern void (*nrn2core_get_partrans_setup_info_)(int tid,
                                                 int& ntar,
                                                 int& nsrc,
                                                 int& type,
                                                 int& ix_vpre,
                                                 int*& sid_target,
                                                 int*& sid_src,
                                                 int*& v_indices);

extern int (*nrn2core_get_dat1_)(int tid,
                                 int& n_presyn,
                                 int& n_netcon,
                                 int*& output_gid,
                                 int*& netcon_srcgid);

extern int (*nrn2core_get_dat2_1_)(int tid,
                                   int& ngid,
                                   int& n_real_gid,
                                   int& nnode,
                                   int& ndiam,
                                   int& nmech,
                                   int*& tml_index,
                                   int*& ml_nodecount,
                                   int& nidata,
                                   int& nvdata,
                                   int& nweight);

extern int (*nrn2core_get_dat2_2_)(int tid,
                                   int*& v_parent_index,
                                   double*& a,
                                   double*& b,
                                   double*& area,
                                   double*& v,
                                   double*& diamvec);

extern int (*nrn2core_get_dat2_mech_)(int tid,
                                      size_t i,
                                      int dsz_inst,
                                      int*& nodeindices,
                                      double*& data,
                                      int*& pdata);

extern int (*nrn2core_get_dat2_3_)(int tid,
                                   int nweight,
                                   int*& output_vindex,
                                   double*& output_threshold,
                                   int*& netcon_pnttype,
                                   int*& netcon_pntindex,
                                   double*& weights,
                                   double*& delays);

extern int (*nrn2core_get_dat2_corepointer_)(int tid, int& n);

extern int (*nrn2core_get_dat2_corepointer_mech_)(int tid,
                                                  int type,
                                                  int& icnt,
                                                  int& dcnt,
                                                  int*& iarray,
                                                  double*& darray);

extern int (*nrn2core_get_dat2_vecplay_)(int tid, int& n);

extern int (*nrn2core_get_dat2_vecplay_inst_)(int tid,
                                              int i,
                                              int& vptype,
                                              int& mtype,
                                              int& ix,
                                              int& sz,
                                              double*& yvec,
                                              double*& tvec);

extern void (*nrn2core_part2_clean_)();

/* what variables to send back to NEURON on each time step */
extern void (*nrn2core_get_trajectory_requests_)(int tid,
                                                 int& bsize,
                                                 int& n_pr,
                                                 void**& vpr,
                                                 int& n_trajec,
                                                 int*& types,
                                                 int*& indices,
                                                 double**& pvars,
                                                 double**& varrays);

/* send values to NEURON on each time step */
extern void (*nrn2core_trajectory_values_)(int tid, int n_pr, void** vpr, double t);

/* Filled the Vector data arrays and send back the sizes at end of run */
extern void (*nrn2core_trajectory_return_)(int tid, int n_pr, int vecsz, void** vpr, double t);
}

#endif /* nrn2core_direct_h */
#ifndef nrn2core_direct_h
#define nrn2core_direct_h

#include <iostream>

extern int corenrn_embedded;

extern void (*nrn2core_mkmech_info_)(std::ostream&);

extern void* (*nrn2core_get_global_dbl_item_)(void*, const char*& name, int& size, double*& val);
extern int (*nrn2core_get_global_int_item_)(const char* name);

extern void (*nrn2core_get_partrans_setup_info_)(int tid, int& ntar, int& nsrc,
  int& type, int& ix_vpre, int*& sid_target, int*& sid_src, int*& v_indices);

extern void (*nrn2core_get_dat1_)(int tid, int& n_presyn, int& n_netcon,
  int*& output_gid, int*& netcon_srcgid);

#endif /* nrn2core_direct_h */

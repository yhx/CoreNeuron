#ifndef nrnbbcore_direct_h
#define nrnbbcore_direct_h

#include <iostream>

extern void (*nrnbbcore_write_mkmech_info_)(std::ostream&);
extern void* (*nrnbbcore_get_global_dbl_item_)(void*, const char*& name, int& size, double*& val);
extern int (*nrnbbcoe_get_global_int_item_)(const char* name);

#endif /* nrnbbcore_direct_h */

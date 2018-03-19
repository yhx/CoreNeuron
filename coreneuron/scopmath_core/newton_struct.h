#ifndef newton_struct_h
#define newton_struct_h



#include "coreneuron/mech/mod2c_core_thread.h"
#include "coreneuron/public/newton.hpp"
namespace coreneuron {
    /* avoid incessant alloc/free memory */
    typedef struct NewtonSpace {
        int n;
        int n_instance;
        double* delta_x;
        double** jacobian;
        int* perm;
        double* high_value;
        double* low_value;
        double* rowmax;
    } NewtonSpace;
}

#endif //newton_struct_h

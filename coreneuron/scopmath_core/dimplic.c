/*
  hopefully a temporary expedient to work around the inability to
  pass function pointers as arguments
*/

#include "coreneuron/mech/cfile/scoplib.h"
#include "coreneuron/mech/mod2c_core_thread.h"
#include "_kinderiv.h"

int derivimplicit_thread(int n, int* slist, int* dlist, DIFUN fun, _threadargsproto_) {
    difun(fun);
    return 0;
}

int nrn_derivimplicit_steer(int fun, _threadargsproto_) {
    switch (fun) { _NRN_DERIVIMPLICIT_CASES }
    return 0;
}

int nrn_newton_steer(int fun, _threadargsproto_) {
    switch (fun) { _NRN_DERIVIMPLICIT_NEWTON_CASES }
    return 0;
}

int nrn_kinetic_steer(int fun, SparseObj* so, double* rhs, _threadargsproto_) {
    switch (fun) { _NRN_KINETIC_CASES }
    return 0;
}

// Derived from nrn/src/scopmath/euler.c

#define der_(arg) _p[der[arg]*_STRIDE]
#define var_(arg) _p[var[arg]*_STRIDE]

int euler_thread(int neqn,
                 int* var,
                 int* der,
                 DIFUN fun,
                 _threadargsproto_) {
    /* calculate the derivatives */
     difun(fun);

    double dt = _nt->_dt;
    int i;

    /* update dependent variables */
    for (i = 0; i < neqn; i++)
        var_(i) += dt * (der_(i));

    return 0;
}

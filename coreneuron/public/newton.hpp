#ifndef CORENEURON_EXTENSION_NEWTON_HPP
#define CORENEURON_EXTENSION_NEWTON_HPP

namespace coreneuron {
struct NewtonSpace;
typedef int NEWTFUN;

#pragma acc routine seq
int nrn_crout_thread(NewtonSpace* ns, int n, 
        double** a, int* perm, 
        int _iml, int _cntml_padded, 
        double *_p, Datum *_ppvar, 
        ThreadDatum *_thread, 
        NrnThread *_nt, double _v);

#pragma acc routine seq
void nrn_scopmath_solve_thread(int n,
        double** a, double* value,
        int* perm,
        double* delta_x, int* s,
        int _iml, int _cntml_padded, 
        double *_p, Datum *_ppvar, 
        ThreadDatum *_thread, 
        NrnThread *_nt, double _v);

#pragma acc routine seq
int nrn_newton_thread(NewtonSpace* ns,
        int n,
        int* s,
        NEWTFUN pfunc,
        double* value,
        int _iml, int _cntml_padded, 
        double *_p, Datum *_ppvar, 
        ThreadDatum *_thread, NrnThread *_nt,
        double _v);

#pragma acc routine seq
void nrn_buildjacobian_thread(NewtonSpace* ns,
        int n,
        int* s,
        NEWTFUN pfunc,
        double* value,
        double** jacobian,
        int _iml, int _cntml_padded, 
        double *_p, Datum *_ppvar, 
        ThreadDatum *_thread, NrnThread *_nt,
        double _v);

NewtonSpace* nrn_cons_newtonspace(int n, int n_instance);

void nrn_destroy_newtonspace(NewtonSpace* ns);

void nrn_newtonspace_copyto_device(NewtonSpace* ns);
}
#endif // CORENEURON_EXTENSION_NEWTON_HPP


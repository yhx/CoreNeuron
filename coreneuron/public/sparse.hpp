#ifndef CORENEURON_SPARSE_HPP
#define CORENEURON_SPARSE_HPP
    
struct SparseObj;

namespace coreneuron {
    typedef int SPFUN;
    #pragma acc routine seq
    double* _nrn_thread_getelm(SparseObj* so, int row, int col, int _iml);
    void* nrn_cons_sparseobj(SPFUN, int, Memb_list*, 
              int _iml, int _cntml_padded, 
              double *_p, Datum *_ppvar, 
              ThreadDatum *_thread, NrnThread *_nt,
              double _v);

    void _nrn_destroy_sparseobj_thread(SparseObj* so);

    #pragma acc routine seq
    int sparse_thread(SparseObj*, int, int*, int*, double*, 
            double, SPFUN, int, int , int , double *, 
            Datum*, ThreadDatum*, NrnThread*,
            double _v);

} //end coreneuron namespace
#endif

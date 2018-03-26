#ifndef register_mech_h
#define register_mech_h

namespace coreneuron {
    void hoc_reg_bbcore_read(int type, bbcore_read_t f);
    void hoc_reg_bbcore_write(int type, bbcore_write_t f);

} // namespace coreneuron

#endif


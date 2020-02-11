#pragma once
#include <algorithm>
#include <map>
#include <vector>
#include "coreneuron/network/netcon.hpp"
#include "coreneuron/network/netcvode.hpp"

namespace coreneuron {

#if defined(ENABLE_REPORTINGLIB) || defined(ENABLE_SONATA_REPORTS)
struct VarWithMapping {
    int id;
    double* var_value;
    VarWithMapping(int id_, double* v_) : id(id_), var_value(v_) {
    }
};

// mapping the set of variables pointers to report to its gid
typedef std::map<int, std::vector<VarWithMapping> > VarsToReport;

class ReportEvent : public DiscreteEvent {
  public:
    ReportEvent(double dt, double tstart, const VarsToReport& filtered_gids, const char* name);

    /** on deliver, call ReportingLib and setup next event */
    virtual void deliver(double t, NetCvode* nc, NrnThread* nt);
    virtual bool require_checkpoint();

  private:
    double dt;
    double step;
    std::string report_path;
    std::vector<int> gids_to_report;
    double tstart;
};
#endif  // defined(ENABLE_REPORTINGLIB) || defined(ENABLE_SONATA_REPORTS)

}  // Namespace coreneuron

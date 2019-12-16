#pragma once
#include <map>
#include <algorithm>
#include <vector>
#include "coreneuron/nrniv/netcon.h"
#include "coreneuron/nrniv/netcvode.h"

#define MAX_REPORT_NAME_LEN 256
#define MAX_FILEPATH_LEN 4096

namespace coreneuron {

#ifdef ENABLE_REPORTING
struct VarWithMapping {
    int id;
    double *var_value;
    VarWithMapping(int id_, double *v_) : id(id_), var_value(v_) {}
};

// mapping the set of variables pointers to report to its gid
typedef std::map<int, std::vector<VarWithMapping> > VarsToReport;

class ReportEvent : public DiscreteEvent {
  public:
    ReportEvent(double dt, double tstart, const VarsToReport &filtered_gids, const char *name);

    /** on deliver, call ReportingLib and setup next event */
    virtual void deliver(double t, NetCvode *nc, NrnThread *nt);
    virtual bool require_checkpoint();

  private:
    double dt;
    double step;
    char report_path[MAX_FILEPATH_LEN];
    std::vector<int> gids_to_report;
    double tstart;
};
#endif  // ENABLE_REPORTING

} // Namespace coreneuron

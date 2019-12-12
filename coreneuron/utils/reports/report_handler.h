#pragma once
#include <vector>
#include <memory>
#include "coreneuron/utils/reports/nrnreport.h"
#include "coreneuron/utils/reports/report_event.h"
#include "coreneuron/nrnoc/multicore.h"

namespace coreneuron {

class ReportHandler {
  public:
    ReportHandler(ReportConfiguration& config): m_report_config(config) {};
    virtual ~ReportHandler() = default;

    virtual void create_report(double dt, double tstop, double delay);
#ifdef ENABLE_REPORTING
    virtual void register_soma_report(NrnThread& nt, ReportConfiguration& config, VarsToReport& vars_to_report) = 0;
    virtual void register_compartment_report(NrnThread& nt, ReportConfiguration& config, VarsToReport& vars_to_report) = 0;
    virtual void register_custom_report(NrnThread& nt, ReportConfiguration& config, VarsToReport& vars_to_report) = 0;

    VarsToReport get_soma_vars_to_report(NrnThread& nt, std::set<int>& target, double* report_variable);
    VarsToReport get_compartment_vars_to_report(NrnThread& nt, std::set<int>& target, double* report_variable);
    VarsToReport get_custom_vars_to_report(NrnThread& nt, ReportConfiguration& report, std::vector<int>& nodes_to_gids);
    std::vector<int> map_gids(NrnThread& nt);
#endif  // ENABLE_REPORTING
  protected:
    ReportConfiguration m_report_config;
#ifdef ENABLE_REPORTING
    std::unique_ptr<ReportEvent> m_report_event;
#endif  // ENABLE_REPORTING
};

} // namespace coreneuron

#pragma once
#include <vector>
#include <memory>
#include "coreneuron/utils/reports/nrnreport.h"
#include "coreneuron/utils/reports/report_event.h"
#include "coreneuron/nrnoc/multicore.h"

namespace coreneuron {

class ReportHandler {
  public:
    ReportHandler(ReportConfiguration& config) : m_report_config(config){};
    virtual ~ReportHandler() = default;

    virtual void create_report(double dt, double tstop, double delay);
#ifdef ENABLE_REPORTING
    virtual void register_soma_report(const NrnThread& nt,
                                      ReportConfiguration& config,
                                      const VarsToReport& vars_to_report) = 0;
    virtual void register_compartment_report(const NrnThread& nt,
                                             ReportConfiguration& config,
                                             const VarsToReport& vars_to_report) = 0;
    virtual void register_custom_report(const NrnThread& nt,
                                        ReportConfiguration& config,
                                        const VarsToReport& vars_to_report) = 0;

    VarsToReport get_soma_vars_to_report(const NrnThread& nt,
                                         const std::set<int>& target,
                                         double* report_variable) const;
    VarsToReport get_compartment_vars_to_report(const NrnThread& nt,
                                                const std::set<int>& target,
                                                double* report_variable) const;
    VarsToReport get_custom_vars_to_report(const NrnThread& nt,
                                           ReportConfiguration& report,
                                           const std::vector<int>& nodes_to_gids) const;
    std::vector<int> map_gids(const NrnThread& nt) const;
#endif  // ENABLE_REPORTING
  protected:
    ReportConfiguration m_report_config;
#ifdef ENABLE_REPORTING
    std::vector<std::unique_ptr<ReportEvent> > m_report_events;
#endif  // ENABLE_REPORTING
};

}  // namespace coreneuron

#pragma once
#include <vector>
#include <memory>
#include "coreneuron/utils/reports/report_handler.h"

namespace coreneuron {

class SonataReportHandler: public ReportHandler {
  public:
    SonataReportHandler(ReportConfiguration& config): ReportHandler(config) {}

    void create_report(double dt, double tstop, double delay) override;
#ifdef ENABLE_REPORTING
    void register_soma_report(NrnThread& nt, ReportConfiguration& config, const VarsToReport& vars_to_report) override;
    void register_compartment_report(NrnThread& nt, ReportConfiguration& config, const VarsToReport& vars_to_report) override;
    void register_custom_report(NrnThread& nt, ReportConfiguration& config, const VarsToReport& vars_to_report) override;

  private:
    void register_report(NrnThread& nt, ReportConfiguration& config, const VarsToReport& vars_to_report);
#endif  // ENABLE_REPORTING
};

} // namespace coreneuron

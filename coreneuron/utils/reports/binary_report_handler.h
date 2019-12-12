#pragma once
#include <vector>
#include <memory>
#include "coreneuron/utils/reports/report_handler.h"

namespace coreneuron {

class BinaryReportHandler: public ReportHandler {
  public:
    BinaryReportHandler(ReportConfiguration& config): ReportHandler(config) {}
    ~BinaryReportHandler() = default;

    void create_report(double dt, double tstop, double delay) override;
#ifdef ENABLE_REPORTING
    void register_soma_report(NrnThread& nt, ReportConfiguration& config, VarsToReport& vars_to_report) override;
    void register_compartment_report(NrnThread& nt, ReportConfiguration& config, VarsToReport& vars_to_report) override;
    void register_custom_report(NrnThread& nt, ReportConfiguration& config, VarsToReport& vars_to_report) override;
#endif  // ENABLE_REPORTING
};

} // namespace coreneuron

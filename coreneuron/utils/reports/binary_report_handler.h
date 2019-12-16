#pragma once
#include <vector>
#include <memory>
#include <functional>
#include "coreneuron/utils/reports/report_handler.h"
#include "coreneuron/utils/reports/nrnsection_mapping.h"

namespace coreneuron {

class BinaryReportHandler: public ReportHandler {
  public:
    BinaryReportHandler(ReportConfiguration& config): ReportHandler(config) {}

    void create_report(double dt, double tstop, double delay) override;
#ifdef ENABLE_REPORTING
    void register_soma_report(NrnThread& nt, ReportConfiguration& config, const VarsToReport& vars_to_report) override;
    void register_compartment_report(NrnThread& nt, ReportConfiguration& config, const VarsToReport& vars_to_report) override;
    void register_custom_report(NrnThread& nt, ReportConfiguration& config, const VarsToReport& vars_to_report) override;

  private:
    void register_report(NrnThread& nt, ReportConfiguration& config, const VarsToReport& vars_to_report,
            std::function<void(CellMapping* mapping, std::array<int, 5>& extra)>& create_extra);
#endif  // ENABLE_REPORTING
};

} // namespace coreneuron

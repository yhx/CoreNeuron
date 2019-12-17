#pragma once
#include <functional>
#include <memory>
#include <vector>
#include "coreneuron/utils/reports/report_handler.h"
#include "coreneuron/utils/reports/nrnsection_mapping.h"

namespace coreneuron {

class BinaryReportHandler : public ReportHandler {
  public:
    BinaryReportHandler(ReportConfiguration& config) : ReportHandler(config) {
    }

    void create_report(double dt, double tstop, double delay) override;
#ifdef ENABLE_REPORTING
    void register_soma_report(const NrnThread& nt,
                              ReportConfiguration& config,
                              const VarsToReport& vars_to_report) override;
    void register_compartment_report(const NrnThread& nt,
                                     ReportConfiguration& config,
                                     const VarsToReport& vars_to_report) override;
    void register_custom_report(const NrnThread& nt,
                                ReportConfiguration& config,
                                const VarsToReport& vars_to_report) override;

  private:
    using create_extra_func = std::function<void(const CellMapping&, std::array<int, 5>&)>;
    void register_report(NrnThread& nt,
                         ReportConfiguration& config,
                         const VarsToReport& vars_to_report,
                         create_extra_func& create_extra);
#endif  // ENABLE_REPORTING
};

}  // namespace coreneuron

#pragma once
#include <vector>
#include <memory>
#include "coreneuron/utils/reports/report_handler.h"

namespace coreneuron {

    class BinaryReportHandler: public ReportHandler {
    public:
        BinaryReportHandler();
        ~BinaryReportHandler() = default;

        void register_report(double dt, double tstop, double delay, ReportConfiguration &config) override;
        void setup_report_engine(double dt_report, double mindelay) override;
        void set_report_buffer_size(int n) override;
        void update_reports(double t) override;
        void finalize_reports() override;
    };

} // namespace coreneuron
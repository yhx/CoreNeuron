#pragma once
#include <vector>
#include <memory>
#include "coreneuron/utils/reports/nrnreport.h"

namespace coreneuron {

    class ReportHandler {
    public:
        ReportHandler(std::vector<ReportConfiguration> report_configs): m_report_configs(report_configs) {}
        virtual ~ReportHandler() = default;

        virtual void register_report(double dt, double tstop, double delay, ReportConfiguration &config) = 0;
        virtual void setup_report_engine(double dt_report, double mindelay) = 0;
        virtual void set_report_buffer_size(int n) = 0;
        virtual void update_reports(double t) = 0;
        virtual void finalize_reports() = 0;

    protected:
        std::vector<ReportConfiguration> m_report_configs;
    };

} // namespace coreneuron

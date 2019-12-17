#include "sonata_report_handler.h"
#include "coreneuron/nrniv/netcvode.h"
#include "coreneuron/nrniv/netcon.h"
#include "coreneuron/utils/reports/nrnsection_mapping.h"
#include "coreneuron/nrnoc/mech_mapping.hpp"
#ifdef ENABLE_REPORTING
#include "bbp/sonata/reports.h"
#endif  // ENABLE_REPORTING

namespace coreneuron {

void SonataReportHandler::create_report(double dt, double tstop, double delay) {
#ifdef ENABLE_REPORTING
    sonata_set_atomic_step(dt);
#endif  // ENABLE_REPORTING
    ReportHandler::create_report(dt, tstop, delay);
}
#ifdef ENABLE_REPORTING
void SonataReportHandler::register_soma_report(const NrnThread& nt,
                                               ReportConfiguration& config,
                                               const VarsToReport& vars_to_report) {
    register_report(nt, config, vars_to_report);
}

void SonataReportHandler::register_compartment_report(const NrnThread& nt,
                                                      ReportConfiguration& config,
                                                      const VarsToReport& vars_to_report) {
    register_report(nt, config, vars_to_report);
}

void SonataReportHandler::register_custom_report(const NrnThread& nt,
                                                 ReportConfiguration& config,
                                                 const VarsToReport& vars_to_report) {
    register_report(nt, config, vars_to_report);
}

void SonataReportHandler::register_report(NrnThread& nt,
                                          ReportConfiguration& config,
                                          const VarsToReport& vars_to_report) {
    sonata_create_report(config.output_path, config.start, config.stop, config.report_dt,
                         config.type_str);
    sonata_set_report_max_buffer_size_hint(config.output_path, config.buffer_size);

    for (const auto& kv : vars_to_report) {
        int gid = kv.first;
        const std::vector<VarWithMapping>& vars = kv.second;
        if (!vars.size())
            continue;

        sonata_add_node(config.output_path, gid);
        sonata_set_report_max_buffer_size_hint(config.output_path, config.buffer_size);
        for (const auto& variable : vars) {
            sonata_add_element(config.output_path, gid, variable.id, variable.var_value);
        }
    }
}
#endif  // ENABLE_REPORTING
}  // Namespace coreneuron

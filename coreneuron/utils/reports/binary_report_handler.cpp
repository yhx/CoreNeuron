#include "binary_report_handler.h"
#ifdef ENABLE_REPORTING
#include "reportinglib/Records.h"
#endif  // ENABLE_REPORTING

namespace coreneuron {

void BinaryReportHandler::create_report(double dt, double tstop, double delay) {
#ifdef ENABLE_REPORTING
    records_set_atomic_step(dt);
#endif  // ENABLE_REPORTING
    ReportHandler::create_report(dt, tstop, delay);
}

#ifdef ENABLE_REPORTING

void create_soma_extra(CellMapping* mapping, std::array<int, 5>& extra) {
    extra = {1, 0, 0, 0, 0};
    /* report extra "mask" all infos not written in report: here only soma count is reported */
    extra[1] = mapping->get_seclist_segment_count("soma");
}

void create_compartment_extra(CellMapping* mapping, std::array<int, 5>& extra) {
    extra = {1, 0, 0, 0, 1};
    extra[1] = mapping->get_seclist_section_count("soma");
    extra[2] = mapping->get_seclist_section_count("axon");
    extra[3] = mapping->get_seclist_section_count("dend");
    extra[4] = mapping->get_seclist_section_count("apic");
    extra[0] = std::accumulate(extra.begin() + 1, extra.end(), 0);
}

void create_custom_extra(CellMapping* mapping, std::array<int, 5>& extra) {
    extra = {1, 0, 0, 0, 1};
    extra[1] = mapping->get_seclist_section_count("soma");
    // todo: axon seems masked on neurodamus side, need to check
    // extra[2] and extra[3]
    extra[4] = mapping->get_seclist_section_count("apic");
    extra[0] = std::accumulate(extra.begin() + 1, extra.end(), 0);
}

void BinaryReportHandler::register_soma_report(NrnThread& nt,
                                               ReportConfiguration& config,
                                               const VarsToReport& vars_to_report) {
    std::function<void(CellMapping * mapping, std::array<int, 5> & extra)> create_extra =
        create_soma_extra;
    register_report(nt, config, vars_to_report, create_extra);
}

void BinaryReportHandler::register_compartment_report(NrnThread& nt,
                                                      ReportConfiguration& config,
                                                      const VarsToReport& vars_to_report) {
    std::function<void(CellMapping * mapping, std::array<int, 5> & extra)> create_extra =
        create_compartment_extra;
    register_report(nt, config, vars_to_report, create_extra);
}

void BinaryReportHandler::register_custom_report(NrnThread& nt,
                                                 ReportConfiguration& config,
                                                 const VarsToReport& vars_to_report) {
    std::function<void(CellMapping * mapping, std::array<int, 5> & extra)> create_extra =
        create_custom_extra;
    register_report(nt, config, vars_to_report, create_extra);
}

void BinaryReportHandler::register_report(
    NrnThread& nt,
    ReportConfiguration& config,
    const VarsToReport& vars_to_report,
    std::function<void(CellMapping* mapping, std::array<int, 5>& extra)>& create_extra) {
    int sizemapping = 1;
    int extramapping = 5;
    std::array<int, 1> mapping = {0};
    std::array<int, 5> extra;
    for (const auto& var : vars_to_report) {
        int gid = var.first;
        auto& vars = var.second;
        if (vars.empty()) {
            continue;
        }
        const auto* mapinfo = static_cast<NrnThreadMappingInfo*>(nt.mapping);
        CellMapping* m = mapinfo->get_cell_mapping(gid);
        extra[0] = vars.size();
        create_extra(m, extra);
        records_add_report(config.output_path, gid, gid, gid, config.start, config.stop,
                           config.report_dt, sizemapping, config.type_str, extramapping,
                           config.unit);

        records_set_report_max_buffer_size_hint(config.output_path, config.buffer_size);
        /** add extra mapping : @todo api changes in reportinglib*/
        records_extra_mapping(config.output_path, gid, 5, extra.data());
        for (const auto& var : vars) {
            mapping[0] = var.id;
            records_add_var_with_mapping(config.output_path, gid, var.var_value, sizemapping,
                                         mapping.data());
        }
    }
}

#endif  // ENABLE_REPORTING

}  // Namespace coreneuron

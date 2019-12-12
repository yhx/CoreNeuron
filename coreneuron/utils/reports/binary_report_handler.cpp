#include "binary_report_handler.h"
#include "coreneuron/utils/reports/nrnsection_mapping.h"
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
void BinaryReportHandler::register_soma_report(NrnThread& nt, ReportConfiguration& config,
                                                VarsToReport& vars_to_report) {
    int sizemapping = 1;
    int extramapping = 5;
    // first column i.e. section numbers
    int mapping[] = {0};
    // first row, from 2nd value (skip gid)
    int extra[] = {1, 0, 0, 0, 0};

    VarsToReport::iterator it;
    for (it = vars_to_report.begin(); it != vars_to_report.end(); ++it) {
        int gid = it->first;
        std::vector<VarWithMapping>& vars = it->second;
        if (!vars.size())
            continue;
        NrnThreadMappingInfo* mapinfo = (NrnThreadMappingInfo*)nt.mapping;
        CellMapping* m = mapinfo->get_cell_mapping(gid);

        /* report extra "mask" all infos not written in report: here only soma count is reported */
        extra[0] = vars.size();
        extra[1] = m->get_seclist_segment_count("soma");

        /** for this gid, get mapping information */
        records_add_report((char*)config.output_path, gid, gid, gid, config.start, config.stop,
                           config.report_dt, sizemapping, (char*)config.type_str, extramapping,
                           (char*)config.unit);

        records_set_report_max_buffer_size_hint((char*)config.output_path, config.buffer_size);
        /** add extra mapping */
        records_extra_mapping(config.output_path, gid, 5, extra);
        for (int var_idx = 0; var_idx < vars.size(); ++var_idx) {
            /** 1st key is section-id and 1st value is segment of soma */
            mapping[0] = vars[var_idx].id;
            records_add_var_with_mapping(config.output_path, gid, vars[var_idx].var_value, sizemapping, mapping);
        }
    }
}

void BinaryReportHandler::register_compartment_report(NrnThread& nt, ReportConfiguration& config,
                                                        VarsToReport& vars_to_report) {
    int sizemapping = 1;
    int extramapping = 5;
    int mapping[] = {0};
    int extra[] = {1, 0, 0, 0, 1};

    VarsToReport::iterator it;
    for (it = vars_to_report.begin(); it != vars_to_report.end(); ++it) {
        int gid = it->first;
        std::vector<VarWithMapping>& vars = it->second;
        if (!vars.size())
            continue;
        NrnThreadMappingInfo* mapinfo = (NrnThreadMappingInfo*)nt.mapping;
        CellMapping* m = mapinfo->get_cell_mapping(gid);
        extra[1] = m->get_seclist_section_count("soma");
        extra[2] = m->get_seclist_section_count("axon");
        extra[3] = m->get_seclist_section_count("dend");
        extra[4] = m->get_seclist_section_count("apic");
        extra[0] = extra[1] + extra[2] + extra[3] + extra[4];

        records_add_report((char*)config.output_path, gid, gid, gid, config.start, config.stop,
                           config.report_dt, sizemapping, (char*)config.type_str, extramapping,
                           (char*)config.unit);

        records_set_report_max_buffer_size_hint((char*)config.output_path, config.buffer_size);
        /** add extra mapping */
        records_extra_mapping(config.output_path, gid, 5, extra);
        for (int var_idx = 0; var_idx < vars.size(); ++var_idx) {
            mapping[0] = vars[var_idx].id;
            records_add_var_with_mapping(config.output_path, gid, vars[var_idx].var_value, sizemapping, mapping);
        }
    }
}

void BinaryReportHandler::register_custom_report(NrnThread& nt, ReportConfiguration& config,
                                                    VarsToReport& vars_to_report) {
    int sizemapping = 1;
    int extramapping = 5;
    int mapping[] = {0};
    int extra[] = {1, 0, 0, 0, 1};
    int segment_count = 0;
    int section_count = 0;

    VarsToReport::iterator it;
    for (it = vars_to_report.begin(); it != vars_to_report.end(); ++it) {
        int gid = it->first;
        std::vector<VarWithMapping>& vars = it->second;
        if (!vars.size())
            continue;
        NrnThreadMappingInfo* mapinfo = (NrnThreadMappingInfo*)nt.mapping;
        CellMapping* m = mapinfo->get_cell_mapping(gid);
        extra[1] = m->get_seclist_section_count("soma");
        // todo: axon seems masked on neurodamus side, need to check
        // extra[2] and extra[3]
        extra[4] = m->get_seclist_section_count("apic");
        extra[0] = extra[1] + extra[2] + extra[3] + extra[4];

        records_add_report((char*)config.output_path, gid, gid, gid, config.start, config.stop,
                           config.report_dt, sizemapping, (char*)config.type_str, extramapping,
                           (char*)config.unit);

        records_set_report_max_buffer_size_hint((char*)config.output_path, config.buffer_size);
        /** add extra mapping : @todo api changes in reportinglib*/
        records_extra_mapping((char*)config.output_path, gid, 5, extra);
        for (int var_idx = 0; var_idx < vars.size(); ++var_idx) {
            mapping[0] = vars[var_idx].id;
            records_add_var_with_mapping(config.output_path, gid, vars[var_idx].var_value, sizemapping, mapping);
        }
    }
}
#endif  // ENABLE_REPORTING

} // Namespace coreneuron

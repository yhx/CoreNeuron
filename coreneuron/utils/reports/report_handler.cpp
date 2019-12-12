#include "report_handler.h"
#include "coreneuron/utils/reports/nrnsection_mapping.h"
#include "coreneuron/nrnoc/mech_mapping.hpp"

namespace coreneuron {

void ReportHandler::create_report(double dt, double tstop, double delay) {
#ifdef ENABLE_REPORTING
    if (m_report_config.start < t) {
        m_report_config.start = t;
    }
    m_report_config.stop = std::min(m_report_config.stop, tstop);

    m_report_config.mech_id = nrn_get_mechtype(m_report_config.mech_name);
    if (m_report_config.type == SynapseReport && m_report_config.mech_id == -1) {
        std::cerr << "[ERROR] mechanism to report: " << m_report_config.mech_name
                  << " is not mapped in this simulation, cannot report on it" << std::endl;
        nrn_abort(1);
    }
    for (int ith = 0; ith < nrn_nthread; ++ith) {
        NrnThread &nt = nrn_threads[ith];
        if (!nt.ncell)
            continue;
        std::vector<int> nodes_to_gid = map_gids(nt);
        VarsToReport vars_to_report;
        switch (m_report_config.type) {
            case SomaReport:
                vars_to_report = get_soma_vars_to_report(nt, m_report_config.target, nt._actual_v);
                register_soma_report(nt, m_report_config, vars_to_report);
                break;
            case CompartmentReport:
                vars_to_report = get_compartment_vars_to_report(nt, m_report_config.target, nt._actual_v);
                register_compartment_report(nt, m_report_config, vars_to_report);
                break;
            case IMembraneReport:
                vars_to_report = get_compartment_vars_to_report(nt, m_report_config.target, nt.nrn_fast_imem->nrn_sav_rhs);
                register_compartment_report(nt, m_report_config, vars_to_report);
                break;
            default:
                vars_to_report = get_custom_vars_to_report(nt, m_report_config, nodes_to_gid);
                register_custom_report(nt, m_report_config, vars_to_report);
        }
        if (vars_to_report.size()) {
            m_report_event = std::make_unique<ReportEvent>(dt, t, vars_to_report, m_report_config.output_path);
            m_report_event->send(t, net_cvode_instance, &nt);
        }
    }
#else
    if (nrnmpi_myid == 0) {
        printf("\n WARNING! : Can't enable reports, recompile with ReportingLib! \n");
    }
#endif  // ENABLE_REPORTING
}

#ifdef ENABLE_REPORTING
VarsToReport ReportHandler::get_soma_vars_to_report(NrnThread &nt, std::set<int> &target, double *report_variable) {
    VarsToReport vars_to_report;
    NrnThreadMappingInfo *mapinfo = (NrnThreadMappingInfo *) nt.mapping;

    if (!mapinfo) {
        std::cout << "[SOMA] Error : mapping information is missing for a Cell group " << std::endl;
        nrn_abort(1);
    }

    for (int i = 0; i < nt.ncell; i++) {
        int gid = nt.presyns[i].gid_;
        // only one element for each gid in this case
        std::vector<VarWithMapping> to_report;
        if (target.find(gid) != target.end()) {
            CellMapping *m = mapinfo->get_cell_mapping(gid);
            if (m == nullptr) {
                std::cout << "[SOMA] Error : mapping information is missing for Soma Report! \n";
                nrn_abort(1);
            }
            /** get  section list mapping for soma */
            SecMapping *s = m->get_seclist_mapping("soma");
            /** 1st key is section-id and 1st value is segment of soma */
            int section_id = s->secmap.begin()->first;
            int idx = s->secmap.begin()->second.front();
            double *variable = report_variable + idx;
            to_report.push_back(VarWithMapping(section_id, variable));
            vars_to_report[gid] = to_report;
        }
    }
    return vars_to_report;
}

VarsToReport
ReportHandler::get_compartment_vars_to_report(NrnThread &nt, std::set<int> &target, double *report_variable) {
    VarsToReport vars_to_report;
    NrnThreadMappingInfo *mapinfo = (NrnThreadMappingInfo *) nt.mapping;
    if (!mapinfo) {
        std::cout << "[COMPARTMENTS] Error : mapping information is missing for a Cell group "
                  << nt.ncell << std::endl;
        nrn_abort(1);
    }

    for (int i = 0; i < nt.ncell; i++) {
        int gid = nt.presyns[i].gid_;
        if (target.find(gid) != target.end()) {
            CellMapping *m = mapinfo->get_cell_mapping(gid);
            if (m == nullptr) {
                std::cout << "[COMPARTMENTS] Error : Compartment mapping information is missing! \n";
                nrn_abort(1);
            }
            std::vector<VarWithMapping> to_report;
            to_report.reserve(m->size());
            for (int j = 0; j < m->size(); j++) {
                SecMapping *s = m->secmapvec[j];
                for (secseg_it_type iterator = s->secmap.begin(); iterator != s->secmap.end();
                     iterator++) {
                    int compartment_id = iterator->first;
                    segvec_type &vec = iterator->second;
                    for (size_t k = 0; k < vec.size(); k++) {
                        int idx = vec[k];
                        /** corresponding voltage in coreneuron voltage array */
                        double *variable = report_variable + idx;
                        to_report.push_back(VarWithMapping(compartment_id, variable));
                    }
                }
            }
            vars_to_report[gid] = to_report;
        }
    }
    return vars_to_report;
}

VarsToReport ReportHandler::get_custom_vars_to_report(NrnThread &nt,
                                                      ReportConfiguration &report,
                                                      std::vector<int> &nodes_to_gids) {
    VarsToReport vars_to_report;
    for (int i = 0; i < nt.ncell; i++) {
        int gid = nt.presyns[i].gid_;
        if (report.target.find(gid) == report.target.end())
            continue;
        Memb_list *ml = nt._ml_list[report.mech_id];
        if (!ml)
            continue;
        std::vector<VarWithMapping> to_report;
        to_report.reserve(ml->nodecount);

        for (int j = 0; j < ml->nodecount; j++) {
            double *is_selected =
                    get_var_location_from_var_name(report.mech_id, SELECTED_VAR_MOD_NAME, ml, j);
            bool report_variable = false;

            /// if there is no variable in mod file then report on every compartment
            /// otherwise check the flag set in mod file
            if (is_selected == nullptr)
                report_variable = true;
            else
                report_variable = *is_selected;

            if ((nodes_to_gids[ml->nodeindices[j]] == gid) && report_variable) {
                double *var_value =
                        get_var_location_from_var_name(report.mech_id, report.var_name, ml, j);
                double *synapse_id =
                        get_var_location_from_var_name(report.mech_id, SYNAPSE_ID_MOD_NAME, ml, j);
                assert(synapse_id && var_value);
                to_report.push_back(VarWithMapping((int) *synapse_id, var_value));
            }
        }
        if (to_report.size()) {
            vars_to_report[gid] = to_report;
        }
    }
    return vars_to_report;
}

// map GIDs of every compartment, it consist in a backward sweep then forward sweep algorithm
std::vector<int> ReportHandler::map_gids(NrnThread &nt) {
    std::vector<int> nodes_gid(nt.end, -1);
    // backward sweep: from presyn compartment propagate back GID to parent
    for (int i = 0; i < nt.n_presyn; i++) {
        int gid = nt.presyns[i].gid_;
        int thvar_index = nt.presyns[i].thvar_index_;
        // only for non artificial cells
        if (thvar_index >= 0) {
            // setting all roots gids of the presyns nodes,
            // index 0 have parent set to 0, so we must stop at j > 0
            // also 0 is the parent of all, so it is an error to attribute a GID to it.
            nodes_gid[thvar_index] = gid;
            for (int j = thvar_index; j > 0; j = nt._v_parent_index[j]) {
                nodes_gid[nt._v_parent_index[j]] = gid;
            }
        }
    }
    // forward sweep: setting all compartements nodes to the GID of its root
    //  already sets on above loop. This is working only because compartments are stored in order
    //  parents follow by childrens
    for (int i = nt.ncell + 1; i < nt.end; i++) {
        nodes_gid[i] = nodes_gid[nt._v_parent_index[i]];
    }
    return nodes_gid;
}
#endif  // ENABLE_REPORTING

} // Namespace coreneuron

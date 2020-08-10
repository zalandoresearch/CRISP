#include "pybind11/pybind11.h"
#include "pybind11/stl.h"
#include "pybind11/numpy.h"
namespace py = pybind11;

#include "crisp.hpp"
#include "gibbs.hpp"
#include "lbp.hpp"

PYBIND11_MODULE(crisp, m) {
    py::class_<Outcome>(m, "Outcome")
        .def(py::init<const tuple<int, int, int>& > ());

    py::class_<Contact>(m, "Contact")
        .def(py::init<const tuple<int, int, int, int>& > ());

    py::class_<Distribution>(m, "Distribution")
        .def(py::init<vector<double>& /*pdf*/>())
        .def("get_log_p", &Distribution::getLogP)
        .def("get_log_p_tail", &Distribution::getLogPTail)
        .def("get_max_outcome_value", &Distribution::getMaxOutcomeValue)
        .def("get_min_outcome_value", &Distribution::getMinOutcomeValue);

    py::class_<InfectionTrace>(m, "InfectionTrace")
        .def(py::init<const int /*t0=0*/,
                      const int /*dE=0*/,
                      const int /*dI=0*/>())
        .def("get_t0", &InfectionTrace::getT0)
        .def("get_dE", &InfectionTrace::getDI)
        .def("get_dI", &InfectionTrace::getDE)
        .def("array", [] (const InfectionTrace &i) { 
                            return (py::array)py::cast(vector<int>({i.getT0(), i.getDE(), i.getDI()}));
                       }, 
            py::return_value_policy::copy);


    py::class_<GibbsPopulationInfectionStatus>(m, "GibbsPIS")
        .def(py::init<int /*S*/,
                      int /*T*/,
                      const vector<ContactTuple>& /*contacts*/,
                      const vector<OutcomeTuple>& /*outcomes*/,
                      Distribution& /*qE*/,
                      Distribution& /*qI*/,
                      double /*alpha*/,
                      double /*beta*/,
                      double /*p0*/,
                      double /*p1*/,
                      bool /*patientZero*/>())

        .def("get_marginals", 
            [] (GibbsPopulationInfectionStatus &g, int N, int burnin, int skip) {
                return (py::array)py::cast(g.getMarginals(N, burnin, skip)); 
            },
            py::arg("N"), 
            py::arg("burnin")=0, 
            py::arg("skip")=0, 
            py::return_value_policy::move)
        
        .def("sample", 
            [] (GibbsPopulationInfectionStatus &g, int N, int burnin, int skip) {
                return (py::array)py::cast(g.sample(N, burnin, skip)); 
            },
            py::arg("N"), 
            py::arg("burnin")=0, 
            py::arg("skip")=0, 
            py::return_value_policy::move)
        
        .def("gibbs_sample", 
            [] (GibbsPopulationInfectionStatus &g, int N, int burnin, int skip) {
                return (py::array)py::cast(g.gibbsSample(N, burnin, skip)); 
            },
            py::arg("N"), 
            py::arg("burnin")=0, 
            py::arg("skip")=0, 
            py::return_value_policy::move)
        
        .def("gibbs_sample_u", 
            [] (GibbsPopulationInfectionStatus &g, int u, int N) {
                return (py::array)py::cast(g.gibbsSampleU(u, N)); 
            },
            py::arg("u"), 
            py::arg("N") = 1, 
            py::return_value_policy::move)
        
        .def("advance", 
            &GibbsPopulationInfectionStatus::advance,
            py::arg("contacts"), 
            py::arg("outcomes"))
        
        .def("get_individual_traces", 
            [] (const GibbsPopulationInfectionStatus &g) {
                return (py::array)py::cast(g.getIndividualTrace()); 
            },
            py::return_value_policy::move)
        
        .def("get_infection_status", 
            [] (GibbsPopulationInfectionStatus &g, int N, int burnin, int skip) {
                return (py::array)py::cast(g.getInfectionStatus(N, burnin, skip)); 
            },
            py::arg("N")=0, 
            py::arg("burnin")=0, 
            py::arg("skip")=0, 
            py::return_value_policy::move);

    py::class_<LBPPopulationInfectionStatus>(m, "LBPPIS")
        .def(py::init<int /*S*/,
                      int /*T*/,
                      const vector<ContactTuple>& /*contacts*/,
                      const vector<OutcomeTuple>& /*outcomes*/,
                      Distribution& /*qE*/,
                      Distribution& /*qI*/,
                      double /*alpha*/,
                      double /*beta*/,
                      double /*p0*/,
                      double /*p1*/,
                      bool /*patientZero*/>())
        .def("propagate",
             [] (LBPPopulationInfectionStatus &g, int N, std::string pt) {
                static const std::map<std::string, LBPPopulationInfectionStatus::PropType> option_strings {
                        { "full",       LBPPopulationInfectionStatus::full },
                        { "forward",    LBPPopulationInfectionStatus::forward },
                        { "baum_welch", LBPPopulationInfectionStatus::baum_welch },
                };                 
                return g.propagate(N, option_strings.at(pt));
            },
            py::arg("N"),
            py::arg("prop_type")="full")

        .def("reset",
            &LBPPopulationInfectionStatus::reset)

        .def("get_marginals", 
            [] (LBPPopulationInfectionStatus &g, int N, int burnin, int skip) {
                return (py::array)py::cast(g.getMarginals(N, burnin, skip)); 
            },
            py::arg("N")=1, 
            py::arg("burnin")=0, 
            py::arg("skip")=0, 
            py::return_value_policy::move)
        
        .def("sample", 
            [] (LBPPopulationInfectionStatus &g, int N, int burnin, int skip) {
                return (py::array)py::cast(g.sample(N, burnin, skip)); 
            },
            py::arg("N")=1, 
            py::arg("burnin")=0, 
            py::arg("skip")=0, 
            py::return_value_policy::move)
                
        .def("advance", 
            &LBPPopulationInfectionStatus::advance,
            py::arg("contacts"), 
            py::arg("outcomes") )
            
        .def("get_infection_status", 
            [] (LBPPopulationInfectionStatus &g, int N, int burnin, int skip) {
                return (py::array)py::cast(g.getInfectionStatus(N, burnin, skip)); 
            },
            py::arg("N")=0, 
            py::arg("burnin")=0, 
            py::arg("skip")=0, 
            py::return_value_policy::move);
        
}

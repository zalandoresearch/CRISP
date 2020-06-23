#include "pybind11/pybind11.h"
#include "pybind11/stl.h"
#include "pybind11/numpy.h"
namespace py = pybind11;

#include "crisp.hpp"
#include "gibbs.hpp"




// // Returns all three numbers in an array
// py::array_t<int> array(void) {
//     int *ret = new int[3];
//     ret[0] = getT0();
//     ret[1] = getDE();
//     ret[2] = getDI();
//     return py::array_t<int>(3, ret);
// }



template <class T>
py::array_t<T> array_from_vector( size_t size[], const vector<T> &vec) {
    return py::array_t(size, vec.data());
}



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


    py::class_<GibbsPopulationInfectionStatus>(m, "PopulationInfectionStatus")
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
        .def(py::init<const GibbsPopulationInfectionStatus &>())

        .def("gibbs_sample", 
            [] (GibbsPopulationInfectionStatus &g, int N, bool burnin, int skip) {
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
            py::arg("outcomes"), 
            py::arg("ignore_tests"))
        
        .def("get_individual_traces", 
            [] (const GibbsPopulationInfectionStatus &g) {
                return (py::array)py::cast(g.getIndividualTrace()); 
            },
            py::return_value_policy::move)
        
        .def("get_infection_status", 
            [] (GibbsPopulationInfectionStatus &g, int N, bool burnin, int skip) {
                return (py::array)py::cast(g.getInfectionStatus(N, burnin, skip)); 
            },
            py::arg("N")=0, 
            py::arg("burnin")=0, 
            py::arg("skip")=0, 
            py::return_value_policy::move);
}

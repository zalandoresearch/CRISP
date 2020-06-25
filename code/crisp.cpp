#include "crisp.hpp"


std::ostream &operator<<(std::ostream &os, Contact const &c) {
    return os << "(" << c.getTargetIndividual() << "," << c.getSourceIndividual() << "," << c.getTime() << "," << c.getCount() <<")";
}


PopulationInfectionStatus::PopulationInfectionStatus(int S, int T,
                const vector<ContactTuple>& contacts, const vector<OutcomeTuple>& outcomes,
                Distribution& qE, Distribution& qI,
                double alpha, double beta, double p0, double p1,
                bool patientZero):

        _noIndividuals(S),
        _noTimeSteps(T),
        _gen(_rd()),
        _contacts(_noIndividuals),
        _outcomes(_noIndividuals),
        _qS(p0),
        _qE(qE),
        _qI(qI),
        _alpha(alpha),
        _beta(beta),
        _p0(p0),
        _p1(p1),
        _log1MinusP1(log(1.0-p1)),
        _minExposure(qE.getMinOutcomeValue()),
        _minInfectious(qI.getMinOutcomeValue()),
        _maxExposure(qE.getMaxOutcomeValue()),
        _maxInfectious(qI.getMaxOutcomeValue()) 
{}

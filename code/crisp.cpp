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
        _noTimeSteps(1),
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

PopulationInfectionStatus::PopulationInfectionStatus( const PopulationInfectionStatus& other) :
    _noIndividuals(other._noIndividuals),
    _noTimeSteps(other._noTimeSteps),
    _gen(other._gen),
    _contacts(other._contacts),
    _outcomes(other._outcomes),
    _qS(other._qS),
    _qE(other._qE),
    _qI(other._qI),
    _alpha(other._alpha),
    _beta(other._beta),
    _p0(other._p0),
    _p1(other._p1),
    _log1MinusP1(other._log1MinusP1),
    _minExposure(other._minExposure),
    _minInfectious(other._minInfectious),
    _maxExposure(other._maxExposure),
    _maxInfectious(other._maxInfectious)
{}
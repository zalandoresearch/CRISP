#ifndef _LBP_HPP
#define _LBP_HPP

#include "factor.hpp"
#include "seir.hpp"
#include "crisp.hpp"

#include <vector>
#include <memory>

using namespace std;


class LBPPopulationInfectionStatus: public PopulationInfectionStatus {

    const SEIRStateSpace &_states;
    vector<vector< unique_ptr<SEIRNode>>> _nodes;
    vector<unique_ptr<Factor>> _factors;

    void propagate(int N);

    bool _forward; // model is a forward message passing model

protected:
    // advance the whole model by one time step, adding new contacts and tests
    virtual void _advance(const vector<ContactTuple>& contacts, const vector<OutcomeTuple>& outcomes, bool updatePrior);

public:
    LBPPopulationInfectionStatus(int S, int T,
                const vector<ContactTuple>& contacts, const vector<OutcomeTuple>& outcomes,
                Distribution& qE, Distribution& qI,
                double alpha, double beta, double p0, double p1,
                bool forward,
                bool patientZero=false);
    ~LBPPopulationInfectionStatus() { delete &_states;}

    // get the posterior marginal distributions P(z_{u,T}|D_{contact}, D_{test})
    virtual vector<vector<double>> getInfectionStatus(int N=0, int burnIn=0, int skip=0);

    // get the posterior marginal distributions P(z_{u,t}|D_{contact}, D_{test})
    virtual array3<double> getMarginals(int N=0, int burnIn=0, int skip=0);

    // sample posterior mariginals $P_{u,t}(z_{u,t})$
    virtual array3<int> sample( int N, int burnIn=0, int skip=0);       


};



#endif
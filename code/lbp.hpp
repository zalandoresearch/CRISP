#ifndef _LBP_HPP
#define _LBP_HPP

#include <vector>
#include <memory>
#include <map>
#include "factor.hpp"
#include "seir.hpp"
#include "crisp.hpp"

using namespace std;

namespace std {
  template <>
  struct hash<tuple<int,int>>
  {
    size_t operator()(const tuple<int,int>& k) const
    {
      using std::size_t;
      using std::hash;

      // Compute individual hash values for first,
      // second and third and combine them using XOR
      // and bit shifting:

      // return ( (hash<int>()(get<0>(k)) ^ (hash<int>()(get<1>(k)) << 1)) );
      return hash<int>()(get<0>(k)*10000 + get<1>(k)) ;
    }
  };
}

// A class for loopy belief propagation over a whole population of individuals
class LBPPopulationInfectionStatus: public PopulationInfectionStatus {
    const SEIRStateSpace &_states;                      // Complete set of detailed SEIR states
    vector<vector<unique_ptr<SEIRNode>>> _nodes;        // 2D array of all variable node pointers over individuals and time
    vector<unique_ptr<Factor>> _factors;                // Pointer over all factors

    map<tuple<int,int>, vector<int>> _contact_map;      // A map for all individuals u at time t that maps to all contacts she had
protected:
    // Advance the whole model by one time step, adding new contacts and tests
    virtual void _advance(const vector<ContactTuple>& contacts, const vector<OutcomeTuple>& outcomes, bool updatePrior);
    // Adds all contacts to the contact map
    void _contact_helper(const vector<ContactTuple>& contacts);
public:
    // Standard constructor for the loopy belief propagation
    LBPPopulationInfectionStatus(int S, int T,
                const vector<ContactTuple>& contacts, const vector<OutcomeTuple>& outcomes,
                Distribution& qE, Distribution& qI,
                double alpha, double beta, double p0, double p1,
                bool patientZero=false);
    virtual ~LBPPopulationInfectionStatus() { delete &_states; }

    // Different types of message update schedules
    enum PropType { forward, baum_welch, full};
    // Runs the message passing schedule for N loops
    void propagate(int N, PropType prop_type = full);
    // Resets the entire set of (cached) message and marginals
    void reset();

    // Get the posterior marginal distributions P(z_{u,T}|D_{contact}, D_{test})
    virtual array2<double> getInfectionStatus(int N=0, int burnIn=0, int skip=0);

    // Get the posterior marginal distributions P(z_{u,t}|D_{contact}, D_{test})
    virtual array3<double> getMarginals(int N=0, int burnIn=0, int skip=0);

    // Sample posterior mariginals $P_{u,t}(z_{u,t})$
    virtual array3<int> sample(int N, int burnIn=0, int skip=0);       
};

#endif
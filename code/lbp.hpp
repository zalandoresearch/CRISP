#ifndef _LBP_HPP
#define _LBP_HPP

#include "factor.hpp"
#include "seir.hpp"
#include "crisp.hpp"

#include <vector>
#include <memory>
#include <map>

using namespace std;

namespace std {

  template <>
  struct hash<tuple<int,int>>
  {
    std::size_t operator()(const tuple<int,int>& k) const
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


class LBPPopulationInfectionStatus: public PopulationInfectionStatus {

    const SEIRStateSpace &_states;
    vector<vector< unique_ptr<SEIRNode>>> _nodes;
    vector<unique_ptr<Factor>> _factors;

    map< tuple<int,int>, vector<int>> _contact_map;

protected:
    // advance the whole model by one time step, adding new contacts and tests
    virtual void _advance(const vector<ContactTuple>& contacts, const vector<OutcomeTuple>& outcomes, bool updatePrior);

    void _contact_helper(const vector<ContactTuple>& contacts);

public:
    LBPPopulationInfectionStatus(int S, int T,
                const vector<ContactTuple>& contacts, const vector<OutcomeTuple>& outcomes,
                Distribution& qE, Distribution& qI,
                double alpha, double beta, double p0, double p1,
                bool patientZero=false);
    virtual ~LBPPopulationInfectionStatus() { delete &_states;}

    enum PropType { forward, baum_welch, full};

    void propagate(int N, PropType prop_type = full);

    void reset();

    // get the posterior marginal distributions P(z_{u,T}|D_{contact}, D_{test})
    virtual vector<vector<double>> getInfectionStatus(int N=0, int burnIn=0, int skip=0);

    // get the posterior marginal distributions P(z_{u,t}|D_{contact}, D_{test})
    virtual array3<double> getMarginals(int N=0, int burnIn=0, int skip=0);

    // sample posterior mariginals $P_{u,t}(z_{u,t})$
    virtual array3<int> sample( int N, int burnIn=0, int skip=0);       


};



#endif
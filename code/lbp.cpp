#include "lbp.hpp"




LBPPopulationInfectionStatus::LBPPopulationInfectionStatus(int S, int T,
                const vector<ContactTuple>& contacts, const vector<OutcomeTuple>& outcomes,
                Distribution& qE, Distribution& qI,
                double alpha, double beta, double p0, double p1,
                bool forward,
                bool patientZero) :
    PopulationInfectionStatus( S, T, contacts, outcomes, qE, qI, alpha, beta, p0, p1, patientZero),
    _states(* new SEIRStateSpace( qE.getMaxOutcomeValue(), qI.getMaxOutcomeValue()) ),
    _nodes( S),
    _forward(forward)
{
    for(int u=0; u<_noIndividuals; u++) 
        for( int t=0; t<_noTimeSteps; t++)
            _nodes[u].emplace_back( new SEIRNode(_states));

    for( size_t u=0; u<_noIndividuals; u++) 
        _factors.emplace_back(new SEIRInitFactor(*_nodes[u][0],patientZero && u==0));
    
    for( size_t t=1; t<_noTimeSteps; t++) {
        for( size_t u=0; u<_noIndividuals; u++) {

            vector<SEIRNode *> contact_nodes;
            for(auto c = contacts.begin(); c != contacts.end();++c) {
                Contact contact(*c);
                const int u_ = contact.getTargetIndividual();
                const int v_ = contact.getSourceIndividual();
                const int t_ = contact.getTime();
                if(u_==u && t_==t-1) {
                    contact_nodes.push_back(_nodes[v_][t_].get());
                }
            }
            if(_forward)
                _factors.emplace_back(new SEIRFactor(_qE, _qI, _p0, _p1, *_nodes[u][t-1], *_nodes[u][t], contact_nodes, {_nodes[u][t].get()}) );
            else
                _factors.emplace_back(new SEIRFactor(_qE, _qI, _p0, _p1, *_nodes[u][t-1], *_nodes[u][t], contact_nodes) );

            for(auto o = outcomes.begin(); o != outcomes.end();++o) {
                Outcome outcome(*o);
                if(outcome.getIndividual()==u && outcome.getTime()==t) {
                    _factors.emplace_back( new SEIRTestFactor( *_nodes[u][t], outcome.getOutcome(), _alpha, _beta) );
                }
            }


        }
    }  
    propagate(2);
}


void LBPPopulationInfectionStatus::_advance(const vector<ContactTuple>& contacts, const vector<OutcomeTuple>& outcomes, bool /* updatePrior */) {
     
    // 0. increase _noTimeSteps
    _noTimeSteps++;

    int t = _noTimeSteps-1;
    for( size_t u=0; u<_noIndividuals; u++) {
        // 1. append to all nodes
        _nodes[u].emplace_back(new SEIRNode(_states));

        // 2. add forward and contact factors
        vector<SEIRNode *> contact_nodes;
        for(auto c = contacts.begin(); c != contacts.end();++c) {
            Contact contact(*c);
            const int u_ = contact.getTargetIndividual();
            const int v_ = contact.getSourceIndividual();
            const int t_ = contact.getTime();
            if(u_==u && t_==t-1) {
                contact_nodes.push_back(_nodes[v_][t_].get());
            }
        }
            if(_forward)
                _factors.emplace_back(new SEIRFactor(_qE, _qI, _p0, _p1, *_nodes[u][t-1], *_nodes[u][t], contact_nodes, {_nodes[u][t].get()}) );
            else
                _factors.emplace_back(new SEIRFactor(_qE, _qI, _p0, _p1, *_nodes[u][t-1], *_nodes[u][t], contact_nodes) );
  
        // 3. add test outcome factors
        for(auto o = outcomes.begin(); o != outcomes.end();++o) {
            Outcome outcome(*o);
            if(outcome.getIndividual()==u && outcome.getTime()==t) {
                _factors.emplace_back( new SEIRTestFactor( *_nodes[u][t], outcome.getOutcome(), _alpha, _beta) );
            }
        }

        _nodes[u][t]->update();
    }
}


void LBPPopulationInfectionStatus::propagate(int N) {
    
    for( int n=0; n<N; n++) {
        for( int t=0; t<_noTimeSteps; t++)
            for( int u=0; u<_noIndividuals; u++) 
                _nodes[u][t]->update();
            
        for( int t=_noTimeSteps-1; t>=0; t--)
            for( int u=_noIndividuals-1; u>=0; u--) 
                _nodes[u][t]->update();
    }
}


vector<vector<double>> LBPPopulationInfectionStatus::getInfectionStatus(int N, int burnIn, int skip) {

    vector<vector<double>> res(_noIndividuals, vector<double> (4));
    for( int u=0; u<_noIndividuals; u++)
        res[u] = normalize(basic_states(*_nodes[u][_noTimeSteps-1]->message_to(), _states));

    return res;
}

// get the posterior marginal distributions P(z_{u,t}|D_{contact}, D_{test})
array3<double> LBPPopulationInfectionStatus::getMarginals(int N, int burnIn, int skip) {

    array3<double> res(_noIndividuals, array2<double>( _noTimeSteps, array1<double>( 4, 0.0)));
    for( int u=0; u<_noIndividuals; u++)
        for(int t=0; t<_noTimeSteps; t++)   
            res[u][t] = normalize(basic_states(*_nodes[u][t]->message_to(), _states));

    return res;

}

// sample posterior mariginals $P_{u,t}(z_{u,t})$
array3<int> LBPPopulationInfectionStatus::sample( int N, int burnIn, int skip) {

    array3<int> Z(N, array2<int>( _noIndividuals, array1<int>( _noTimeSteps, 0)));
    for( int u=0; u<_noIndividuals; u++)
        for(int t=0; t<_noTimeSteps; t++)
            for( int n=0; n<N; n++) {   
                auto p = normalize(basic_states(*_nodes[u][t]->message_to(), _states));
                double r = (double)rand() / RAND_MAX;
                for( auto pi: p ) {
                    r -=pi;
                    if( r<=0) break;
                    Z[n][u][t]++;
                }
            }
    return Z;
}       


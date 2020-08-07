#include <chrono>
#include "lbp.hpp"

// Standard constructor for the loopy belief propagation
LBPPopulationInfectionStatus::LBPPopulationInfectionStatus(int S, int T,
                const vector<ContactTuple>& contacts, const vector<OutcomeTuple>& outcomes,
                Distribution& qE, Distribution& qI,
                double alpha, double beta, double p0, double p1,
                bool patientZero) :
    PopulationInfectionStatus(S, T, contacts, outcomes, qE, qI, alpha, beta, p0, p1, patientZero),
    _states(*(new SEIRStateSpace(qE.getMaxOutcomeValue(), qI.getMaxOutcomeValue()))),
    _nodes(S)
{
    // Allocate all the variables nodes
    for(int u=0; u<_noIndividuals; ++u) {
        for(int t=0; t<_noTimeSteps; ++t) {
            _nodes[u].emplace_back(new SEIRNode(_states));
        }
    }

    // Allocate enough space for all factors
    _factors.reserve(_noIndividuals*_noTimeSteps + outcomes.size());

    // Setup the prior nodes (if patient_zero == true, then u=0 gets a different prior)
    for(int u=0; u<_noIndividuals; ++u) {
        _factors.emplace_back(new SEIRInitFactor(*(_nodes[u][0]),patientZero && u==0, _p0));
    }

    // Prepare the contact map
    _contact_helper(contacts);

    // Setup all the infecton spread and test outcome nodes
    for(int t=1; t<_noTimeSteps; ++t) {
        for(int u=0; u<_noIndividuals; ++u) {
            vector<SEIRNode*> contact_nodes;
            auto cmi = _contact_map.find(make_tuple(u,t-1));
            if (cmi != _contact_map.end()) {
                for(auto v:  cmi->second) {
                    contact_nodes.push_back(_nodes[v][t-1].get());
                }
            }
            _factors.emplace_back(new SEIRFactor(_qE, _qI, _p0, _p1, *(_nodes[u][t-1]), *(_nodes[u][t]), contact_nodes));

            // Loop through all test outcomes to setup the test outcome(s) for individual u at time t
            // TODO: This is *very* inefficient and assumes that we have a very small number of tests only!
            for(auto o = outcomes.begin(); o != outcomes.end(); ++o) {
                Outcome outcome(*o);
                if(outcome.getIndividual()==u && outcome.getTime()==t) {
                    _factors.emplace_back( new SEIRTestFactor( *_nodes[u][t], outcome.getOutcome(), _alpha, _beta));
                }
            }
        }
    }
}

// Advance the whole model by one time step, adding new contacts and tests
void LBPPopulationInfectionStatus::_advance(const vector<ContactTuple>& contacts, 
                                            const vector<OutcomeTuple>& outcomes, 
                                            bool /* updatePrior */) {
    // (0) Increase _noTimeSteps
    _noTimeSteps++;

    int t = _noTimeSteps-1;
    for(int u=0; u<_noIndividuals; ++u) {
        // (1) Append to all new variable nodes
        _nodes[u].emplace_back(new SEIRNode(_states));

        // (2) Add forward and contact factors
        _contact_helper(contacts);
        vector<SEIRNode*> contact_nodes;
        auto cmi = _contact_map.find(make_tuple(u,t-1));
        if (cmi != _contact_map.end()) {
            for( auto v:  cmi->second) {
                contact_nodes.push_back(_nodes[v][t-1].get());
            }
        }
        _factors.emplace_back(new SEIRFactor(_qE, _qI, _p0, _p1, *_nodes[u][t-1], *_nodes[u][t], contact_nodes));

        // (3) Add test outcome factors
        // TODO: This is *very* inefficient and assumes that we have a very small number of tests only!
        for(auto o = outcomes.begin(); o != outcomes.end();++o) {
            Outcome outcome(*o);
            if(outcome.getIndividual()==u && outcome.getTime()==t) {
                _factors.emplace_back( new SEIRTestFactor( *_nodes[u][t], outcome.getOutcome(), _alpha, _beta));
            }
        }
    }

    // TODO: We should get away from fixed time steps in the past
    for(; t>=_noTimeSteps-5 && t>=0; --t) {
        for(int u=0; u<_noIndividuals; u++)
            _nodes[u][t]->update(SEIRNode::full);
    }
    t+=2;
    for(; t<_noTimeSteps; ++t) {
        for(int u=0; u<_noIndividuals; u++)
            _nodes[u][t]->update(SEIRNode::full);
    }

    return;
}

// Adds all contacts to the contact map
void LBPPopulationInfectionStatus::_contact_helper(const vector<ContactTuple>& contacts) {
    for(auto c = contacts.begin(); c != contacts.end(); ++c) {
        Contact contact(*c);
        const int u = contact.getTargetIndividual();
        const int v = contact.getSourceIndividual();
        const int t = contact.getTime();

        _contact_map[make_tuple(u,t)].push_back(v);
    }
    return;
}

// Runs the message passing schedule for N loops
void LBPPopulationInfectionStatus::propagate(int N, PropType prop_type) {
    Node::message_update_counter = 0;
    cerr << _factors.size() << " factors, N=" << N << endl;
    auto tic = chrono::steady_clock::now();

    for(int n=0; n<N; ++n) {
        // First, update all variables in either forward or full mode forward in time
        for(int t=0; t<_noTimeSteps; ++t) {
            for(int u=0; u<_noIndividuals; ++u) {
                switch(prop_type) {
                    case forward:
                    case baum_welch:
                        _nodes[u][t]->update(SEIRNode::forward);
                        break;
                    case full:
                        _nodes[u][t]->update(SEIRNode::full);
                        break;
                }
            }

            cerr << n << ". " << t;
            auto toc = chrono::steady_clock::now();
            cerr << " " << Node::message_update_counter << " message updates in";
            double d = chrono::duration_cast<chrono::nanoseconds>(toc-tic).count();
            d /= 1e9;
            cerr << " (" << (Node::message_update_counter/d) << " message updates/s)" << endl;
        }

        // In case of Baum-Welch or Full, now update all variables in backward in time
        if(prop_type==baum_welch || prop_type==full) {
            for(int t=_noTimeSteps-1; t>=0; --t) {
                for(int u=_noIndividuals-1; u>=0; --u) {
                    switch(prop_type) {
                        case forward:
                            // We will never get here because of the if clause
                            break;
                        case baum_welch:
                            _nodes[u][t]->update(SEIRNode::backward);
                            break;
                        case full:
                            _nodes[u][t]->update(SEIRNode::full);
                            break;
                    }
                }
                cerr << n << ". " << t;
                auto toc = chrono::steady_clock::now();
                cerr << " " << Node::message_update_counter << " message updates in";
                double d = chrono::duration_cast<chrono::nanoseconds>(toc-tic).count();
                d /= 1e9;
                cerr << " (" << (Node::message_update_counter/d) << " message updates/s)" << endl;
            }
        }
    }

    return;
}

// Resets the entire set of (cached) message and marginals
void LBPPopulationInfectionStatus::reset() {
    // Loop over all u
    for(auto &nu: _nodes) {
        // Loop over all t
        for(auto &nut: nu) {
            nut->reset();
        }
    }
    return;
}

// Gets the posterior probability by summing over all detailed state outcomes
vector<vector<double>> LBPPopulationInfectionStatus::getInfectionStatus(int /*N*/, int /*burnIn*/, int /*skip*/) {
    vector<vector<double>> res(_noIndividuals, vector<double> (4));
    for(int u=0; u<_noIndividuals; ++u) {
        res[u] = basic_states(*_nodes[u][_noTimeSteps-1]->normalized_message_to_factor(), _states);
    }

    return res;
}

// Get the posterior marginal distributions P(z_{u,t}|D_{contact}, D_{test})
array3<double> LBPPopulationInfectionStatus::getMarginals(int /*N*/, int /*burnIn*/, int /*skip*/) {

    array3<double> res(_noIndividuals, array2<double>( _noTimeSteps, array1<double>( 4, 0.0)));
    for(int u=0; u<_noIndividuals; ++u) {
        for(int t=0; t<_noTimeSteps; ++t) {
            res[u][t] = basic_states(*_nodes[u][t]->normalized_message_to_factor(), _states);
        }
    }

    return res;
}

// Sample posterior mariginals $P_{u,t}(z_{u,t})$
array3<int> LBPPopulationInfectionStatus::sample(int N, int /*burnIn*/, int /*skip*/) {

    array3<int> Z(N, array2<int>(_noIndividuals, array1<int>(_noTimeSteps, 0)));
    for(int u=0; u<_noIndividuals; ++u) {
        for(int t=0; t<_noTimeSteps; ++t) {
            for( int n=0; n<N; n++) {
                auto p = basic_states(*_nodes[u][t]->normalized_message_to_factor(), _states);
                double r = (double)rand() / RAND_MAX;
                for(auto pi: p) {
                    r -=pi;
                    if(r<=0) break;
                    Z[n][u][t]++;
                }
            }
        }
    }
    return Z;
}
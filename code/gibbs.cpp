#include <iostream>
#include <map>
#include <random>
#include <cmath>
#include <cassert>
#include "gibbs.hpp"

using namespace std;

InfectionTrace::InfectionTrace(const int t0, const int dE, const int dI) :
    _t0(t0),
    _t0dE(t0 + dE),
    _t0dEdI(t0+dE+dI) { }

// Overloading [] operator to access elements in array style
SEIRState InfectionTrace::operator[](int index) const {
    if (index < _t0) {
        return SEIRState::S;
    } else if (index < _t0dE) {
        return SEIRState(SEIRState::E,index-_t0dE);
    } else if (index < _t0dEdI) {
        return SEIRState(SEIRState::I, index-_t0dEdI);
    }
    return SEIRState::R;
}


std::ostream &operator<<(std::ostream &os, InfectionTrace const &it) {
    return os << "(" << it.getT0() << "," << it.getDE() << "," << it.getDI() <<")";
}


// determine whether a contact is valid or not
bool GibbsPopulationInfectionStatus::valid(int u, int v, int t) {
    return((u>=0) && (u<_noIndividuals) && (v>=0) && (v<_noIndividuals) && (t>=0) && (t<_noTimeSteps));
}

// compute the list of all possible infection traces (t0, dE, dI) together with their a priori probabilities
void GibbsPopulationInfectionStatus::initPrior() {

    _logPrior.clear();
    double sumP = 0;
    for(int t0 = 0; t0 <= _noTimeSteps; ++t0) {
        if(t0 == _noTimeSteps) {
            InfectionTrace z(t0,0,0);
            double logP = _qS.getLogPTail(t0);
            // store infection trace (t0 = _noTimeSteps, dE = 0, dI = 0) and logP = log l0(t0)
            if (isfinite(logP))
                _logPrior.push_back(tuple<InfectionTrace,double> (z, logP));
            sumP += exp(logP);
        } else {
            for(int dE = 1; dE <= _maxExposure; ++dE) {
                if(t0+dE == _noTimeSteps) {
                    InfectionTrace z(t0,_noTimeSteps-t0,0);
                    double logP = _qS.getLogP(t0) + _qE.getLogPTail(dE);
                    // store infection trace (t0, dE, dI = 0) and logP = log l0(t0) + log lE(dE)
                    if (isfinite(logP))
                        _logPrior.push_back(tuple<InfectionTrace,double> (z, logP));
                    sumP += exp(logP);
                    break;
                } else {
                    for (int dI = 1; dI <= _maxInfectious; ++dI) {
                        InfectionTrace z(t0,dE,dI);
                        if (t0+dE+dI == _noTimeSteps) {
                            double logP = _qS.getLogP(t0) + _qE.getLogP(dE) + _qI.getLogPTail(dI);
                            // store infection trace (t0, dE, dI) and logP = log l0(t0) + log lE(dE) + log lI(dI)
                            if (isfinite(logP))
                                _logPrior.push_back(tuple<InfectionTrace,double> (z, logP));
                            sumP += exp(logP);
                            break;
                        } else {
                            double logP = _qS.getLogP(t0) + _qE.getLogP(dE) + _qI.getLogP(dI);
                            // store infection trace (t0, dE, dI) and logP = log l0(t0) + log lE(dE) + log lI(dI)
                            if (isfinite(logP))
                                _logPrior.push_back(tuple<InfectionTrace,double> (z, logP));
                            sumP += exp(logP);
                        }
                    }
                }
            }
        }
    }
/*
    // check if prior probabilities sum up to 1, keep this for debugging purposes
    for(auto z=_logPrior.begin(); z!=_logPrior.end(); ++z) {
        cerr << get<0>(*z) <<","<< get<1>(*z) << endl;
    }
    cerr << _logPrior.size() << " traces, sum(p_i)=" << sumP << endl;
*/

}


// updates _infectiousContactCounts with contribution from zu removed (before) or added (after)
void GibbsPopulationInfectionStatus::updateInfectiousContactCounts(int u, const InfectionTrace &zu, bool before) {
    for(int t=0; t<_noTimeSteps-1; ++t) {
        if(zu[t] == SEIRState::I) {
            for(auto fc: _futureContact(u,t)) {
                int v = get<0>(fc);
                int x = get<1>(fc);
                _infectiousContactCounts[v][t+1] += before ? -x : x;
            }
        }
    }

    return;
}

// Prints an infection trace on the screen
void GibbsPopulationInfectionStatus::print(InfectionTrace z) {
    for(int t=0; t < _noTimeSteps; ++t) {
        switch(z[t].phase()) {
            case SEIRState::S: cout << "S"; break;
            case SEIRState::E: cout << "E"; break;
            case SEIRState::I: cout << "I"; break;
            case SEIRState::R: cout << "R"; break;
            default: cout << "*"; break;
        }
    }
}

// advance the whole model by one time step, adding new contacts and tests
void GibbsPopulationInfectionStatus::_advance(const vector<ContactTuple>& contacts, const vector<OutcomeTuple>& outcomes, bool updatePrior) {
    // 0. increase _noTimeSteps
    _noTimeSteps++;
    //int T=_noTimeSteps-1;

    // 1. add contacts and tests
        for(int i = 0; i < _noIndividuals; ++i) {
        _contacts[i].resize(_noTimeSteps);
    }
    for(auto c = contacts.begin(); c != contacts.end();++c) {
        Contact contact(*c);
        const int u = contact.getTargetIndividual();
        const int v = contact.getSourceIndividual();
        const int t = contact.getTime();
        if(valid(u,v,t) && t==_noTimeSteps-1) {
            _contacts[u][t].push_back(make_tuple(v,contact.getCount()));
        }
    }
    for(auto o = outcomes.begin(); o != outcomes.end();++o) {
        Outcome outcome(*o);
        if(outcome.getTime() == _noTimeSteps-1) {
            _outcomes[outcome.getIndividual()].push_back(outcome);
        }
    }

    // 2. reinitialize prior
    if (updatePrior) initPrior();

    // 3. update _infectiousContactCounts - stores the number of infectious contacts that individual u has at each time step t
    for(int u=0; u<_noIndividuals; ++u) {
        _infectiousContactCounts[u].resize(_noTimeSteps);
        _infectiousContactCounts[u][_noTimeSteps-1] = 0;
        for(auto pc: _pastContact(u,_noTimeSteps-1)) {
            int v = get<0>(pc);
            int x = get<1>(pc);
            assert(u!=v);
            if (_individualTrace[v][_noTimeSteps-2].phase()==SEIRState::I) {
                _infectiousContactCounts[u][_noTimeSteps-1] += x;
            }
        }
    }

    // 4. append to _individualTraces
    uniform_real_distribution<> rng(0.0, 1.0);
    for(int u=0; u<_noIndividuals; ++u) {
        auto zu = _individualTrace[u];
        double rnd =  rng(_gen);
        switch(zu[_noTimeSteps-2].phase()) {
            case SEIRState::S: {
                // Compute f(u, t, Zt) using Equation (5) in paper
                double p = (1.0-_p0)*pow(1-_p1, _infectiousContactCounts[u][_noTimeSteps-1]);
                // Transition from S to S with probability f(u, t, Zt), and from S to E with probability 1 - f(u, t, Zt) - see Equation (4) in paper
                if (rnd<p) {
                    zu.setT0(zu.getT0()+1);
                } else {
                    zu.setDE(zu.getDE()+1);
                } }
                break;

            case SEIRState::E: {
                int dE = zu.getDE();
                // Compute g(u, t, Zut) using Equation (7) in paper
                double pEE = 1.0-exp(_qE.getLogP(dE)-_qE.getLogPTail(dE));
                // Transition from E to E with probability 1 - g(u, t, Zut), and from E to I with probability g(u, t, Zut) - see Equation (4) in paper
                if (rnd<pEE){
                    zu.setDE(zu.getDE()+1);
                } else {
                    zu.setDI(zu.getDI()+1);
                } }
                break;

            case SEIRState::I: {
                int dI = zu.getDI();
                // Compute h(u, t, Zut) using Equation (8) in paper
                double pII = 1.0-exp(_qI.getLogP(dI)-_qI.getLogPTail(dI));
                // Transition from I to I with probability 1 - h(u, t, Zut), and from I to R with probability h(u, t, Zut) - see Equation (4) in paper
                if (rnd<pII){
                    zu.setDI(zu.getDI()+1);
                } }
                break;

            case SEIRState::R:
                break;
        }
        _individualTrace[u] = zu;
    }
}

// Initializes a population of S individuals over T time steps with contact and outcomes
GibbsPopulationInfectionStatus::GibbsPopulationInfectionStatus(int S, int T,
                const vector<ContactTuple>& contacts, const vector<OutcomeTuple>& outcomes,
                Distribution& qE, Distribution& qI,
                double alpha, double beta, double p0, double p1,
                bool patientZero) :
    PopulationInfectionStatus( S, 1, contacts, outcomes, qE, qI, alpha, beta, p0, p1, patientZero),
    _individualTrace(S, InfectionTrace(1,0,0)),
    _infectiousContactCounts(_noIndividuals) {

    uniform_real_distribution<> rng(0.0, 1.0);
    for(int u=0; u<_noIndividuals; ++u) {
        if(u==0 && patientZero) {
            _individualTrace[u] = InfectionTrace(0,1,0);
        }
        else {
            _individualTrace[u] = rng(_gen)<_p0 ? InfectionTrace(0,1,0) : InfectionTrace(1,0,0);
        }
    }

    for(int t=1; t<T; ++t) {
        _advance(contacts, outcomes,  /*updatePrior =*/ false);
    }
    initPrior();
}

// implements Gibbs sampling over the whole population
vector<vector<vector<int>>> GibbsPopulationInfectionStatus::gibbsSample(int N, int burnIn, int skip) {

    _precisionWarning1Issued = false;

    //vector<int> sample(N*_noIndividuals*3);
    vector<vector<vector<int>>> sample(N, vector<vector<int>>(_noIndividuals, vector<int>(3, 0)));

    for(int n=0; n<burnIn; ++n) {
        for(int u=0; u < _noIndividuals; ++u) {
            gibbsSampleU(u);
        }
    }

    for(int n=0; n<N; ++n) {

        for(int s=0; s<skip; ++s) {
            for(int u=0; u< _noIndividuals; ++u) {
                gibbsSampleU(u);
            }
        }

        for(int u=0; u< _noIndividuals; ++u) {
            if(n>0) gibbsSampleU(u);
            sample[n][u][0] = _individualTrace[u].getT0();
            sample[n][u][1] = _individualTrace[u].getDE();
            sample[n][u][2] = _individualTrace[u].getDI();
        }
    }

    return sample;
}

// Implements a single Gibbs sampling step according to Algorithm 1 in the paper
vector<vector<int>> GibbsPopulationInfectionStatus::gibbsSampleU(int u, int N) {

    vector<double> active(_noTimeSteps, 0.0);
    vector<double> passiveSS(_noTimeSteps, 0.0);
    vector<double> passiveSE(_noTimeSteps, 0.0);
    vector<double> non_infectious(_noTimeSteps, 0.0);

    updateInfectiousContactCounts(u, _individualTrace[u], /*before =*/ true);


    for(int t=0; t<_noTimeSteps; ++t) {
        double logP = 0;
        for(const auto &fc: _futureContact(u,t)) {
            int v = get<0>(fc);
            int x = get<1>(fc);
            assert(v!=u);

            // Compute log f(v, t, Zt^I) - log f(v, t, Zt^{neg I}) for v in CS(u, t) in the paper
            if(_individualTrace[v][t] == SEIRState::S && _individualTrace[v][t+1] == SEIRState::S) {
                logP += x * _log1MinusP1;
            }

                // Compute log f(v, t, Zt^I) - log f(v, t, Zt^{neg I}) for v in CE(u, t) in the paper
            if(_individualTrace[v][t] == SEIRState::S && _individualTrace[v][t+1] == SEIRState::E) {
                int parent_contacts = _infectiousContactCounts[v][t+1];

                double tmp = 0;
                tmp += log(1.0-(1.0-_p0)*pow(1.0-_p1, parent_contacts+x));
                tmp -= log(1.0-(1.0-_p0)*pow(1.0-_p1, parent_contacts));
                if (tmp<_precisionThreshold1 && !_precisionWarning1Issued) {
                    cerr << "WARNING: too large R0, sampling precision may be low" << endl;
                    _precisionWarning1Issued = true;
                }
                logP += tmp;
            }

        }
        // Pre-compute and store in active[t] the quantity (log B(zu, t, Zt^I) - log B(zu, t, Zt^{neg I})) using Equation (25) in paper
        active[t] = logP;

        int parent_contacts = 0;
        for(auto pc = _pastContact(u,t).begin(); pc != _pastContact(u,t).end(); ++pc) {
            int v = get<0>(*pc);
            int x = get<1>(*pc);
            if (_individualTrace[v][t-1]==SEIRState::I) {
                parent_contacts += x;
            }
        }
        // Pre-compute and store in passiveSS[t] the quantity log l_infected(t) using Equation (22) in paper
        passiveSS[t] = parent_contacts * _log1MinusP1;
        // Pre-compute and store in passiveSE[t] the quantity log l'_infected(t) using Equation (23) in paper
        passiveSE[t] = log(1.0-(1.0-_p0)*pow(1.0-_p1, parent_contacts)) - log(_p0);
    }

    // Computes the delta log factor due to the test outcomes
    for(auto test = _outcomes[u].begin(); test != _outcomes[u].end(); ++test) {
        const int t = test->getTime();
        const TestOutcome o = test->getOutcome();

        if (o == Positive) {
                active[t] += log(1.0 - _alpha);
                non_infectious[t] += log(_beta);
            } else {
                active[t] += log(_alpha);
                non_infectious[t] += log(1.0-_beta);
            }
    }

    vector<double> w(_logPrior.size());
    double max_w = -INFINITY;
    // Compute (un-normalized log) conditional probabilities in w[i] for each infection trace i = (t0, dE, dI)
    for(unsigned long i=0; i<_logPrior.size(); ++i) {
            const InfectionTrace &z = get<0>(_logPrior[i]);
            // Initialize w[i] to log l_0(t_0) + log l_E(d_E) + log l_I(d_I)
            w[i] = get<1>(_logPrior[i]);

            for(int t=0; t<_noTimeSteps; ++t) {
                switch (z[t].phase()) {
                    case SEIRState::I:
                        // Add terms in log B(zu) - see Equation (26) in paper
                        w[i] += active[t];
                        break;
                    case SEIRState::S:
                        // Add terms log l_infected(t) - see Equation (24) in paper
                        w[i] += passiveSS[t];
                        w[i] += non_infectious[t];
                        break;
                    default:
                        w[i] += non_infectious[t];
                }
            }
            // Add term log l'_infected(t0) - see Equation (24) in paper
            if (z.getT0()<_noTimeSteps)
                w[i] += passiveSE[z.getT0()];
            max_w = max(max_w, w[i]);
    }

    double sum_w = 0;
    for(unsigned long i=0; i<_logPrior.size(); ++i) {
        w[i] = exp(w[i]-max_w);
        sum_w += w[i];
        }
    // Normalize (un-normalized) conditional probabilities
    for(unsigned long i=0; i<_logPrior.size(); ++i) {
        w[i] /= sum_w;
    }

    discrete_distribution<> dist(w.begin(),w.end());

    vector<vector<int>> sample(N, vector<int>(3,0));
    for(int n=0; n<N; ++n) {
        // Sample i = (t0, dE, dI) with probability proportional to w[i]
        const int sampleIndex = dist(_gen);
        InfectionTrace sample_trace = get<0>(_logPrior[sampleIndex]);

        if (n==0) {
            // Update log-likelihood and sample trace with the first sample
            _individualTrace[u] = sample_trace;
            updateInfectiousContactCounts(u, sample_trace, /*before =*/ false);
            }

            sample[n][0]   = sample_trace.getT0();
            sample[n][1] = sample_trace.getDE();
            sample[n][2] = sample_trace.getDI();
    }
    return sample;
}

vector<vector<int>> GibbsPopulationInfectionStatus::getIndividualTrace() const {
    vector<vector<int>> res(_noIndividuals,vector<int>(3) );

    for(int i=0; i<_noIndividuals; ++i) {
        res[i][0] = _individualTrace[i].getT0();
        res[i][1] = _individualTrace[i].getDE();
        res[i][2] = _individualTrace[i].getDI();
    }
    return res;
}

vector<vector<double>> GibbsPopulationInfectionStatus::getInfectionStatus(int N, int burnIn, int skip) {
    _precisionWarning1Issued = false;

    vector<vector<double>> P_z_uT(_noIndividuals,vector<double>(4,0.0));

    if(N>0) {
        for(int n=0; n<burnIn; ++n) { for(int u=0; u < _noIndividuals; ++u) { gibbsSampleU(u); } }

        for(int n=0; n<N; ++n) {
            for(int s=0; s<skip; ++s) { for(int u=0; u< _noIndividuals; ++u) { gibbsSampleU(u); } }
            for(int u=0; u< _noIndividuals; ++u) {
                gibbsSampleU(u);
                P_z_uT[u][(int)(_individualTrace[u][_noTimeSteps-1].phase())] += 1;
            }
        }
        for(int u=0; u< _noIndividuals; ++u) {
            for(int i=0; i<4; ++i) {
                P_z_uT[u][i] /= N;
            }
        }
    }
    else {
            for(int u=0; u< _noIndividuals; ++u) {
                P_z_uT[u][(int)(_individualTrace[u][_noTimeSteps-1].phase())] = 1;
            }
    }

    return P_z_uT;
}


// get the posterior marginal distributions P(z_{u,t}|D_{contact}, D_{test})
array3<double> GibbsPopulationInfectionStatus::getMarginals(int N, int burnIn, int skip) {

    auto Z = gibbsSample(N, burnIn, skip);
    array3<double> p(_noIndividuals, array2<double>( _noTimeSteps, array1<double>( 4, 0.0)));

    for( int n=0; n<N; n++)
        for( int u=0; u<_noIndividuals; u++) {
            int t=0;
            for( int t_=0; t_<Z[n][u][0]; t_++) p[u][t++][0] += 1.0/N;
            for( int t_=0; t_<Z[n][u][1]; t_++) p[u][t++][1] += 1.0/N;
            for( int t_=0; t_<Z[n][u][2]; t_++) p[u][t++][2] += 1.0/N;
            while(t<_noTimeSteps) p[u][t++][3] += 1.0/N;
        }

    return p;
}


// sample posterior mariginals $P_{u,t}(z_{u,t})$
array3<int> GibbsPopulationInfectionStatus::sample( int N, int burnIn, int skip) {

    auto Z = gibbsSample(N, burnIn, skip);
    array3<int> Z_(N, array2<int>( _noIndividuals, array1<int>( _noTimeSteps, 0)));

    for( int n=0; n<N; n++)
        for( int u=0; u<_noIndividuals; u++) {
            int t=0;
            for( int t_=0; t_<Z[n][u][0]; t_++) Z_[n][u][t++] = 0;
            for( int t_=0; t_<Z[n][u][1]; t_++) Z_[n][u][t++] = 1;
            for( int t_=0; t_<Z[n][u][2]; t_++) Z_[n][u][t++] = 2;
            while(t<_noTimeSteps) Z_[n][u][t++] = 3;
        }

    return Z_;
}

#include <iostream>
#include <map>
#include <random>
#include <cmath>

#include "pybind11/pybind11.h"
#include "pybind11/stl.h"
#include "pybind11/numpy.h"

using namespace std;
namespace py = pybind11;

// A class for a distribution over durations
class Distribution {
    private:
        double* _logPdf;           // The logarithm of the PDF of the distribution
        double* _logPdfTail;       // The logarithm of 1 minus the CDF of the distribution (equivalent to the log-pdf of the truncated distribution)
        int     _maxOutcomeValue;  // The maximum outcome value
        int     _minOutcomeValue;  // The minimum outcome value
    public:
        // Constructs the distribution with zero entries
        Distribution() {
            _maxOutcomeValue = 0;
            _minOutcomeValue = 0;
            _logPdf = new double[_maxOutcomeValue+1];
            _logPdfTail = new double[_maxOutcomeValue+1];
            _logPdf[0] = 0.0;
            _logPdfTail[0] = 0.0;
        }

        // Copy constructor
        Distribution(const Distribution& d) {
            _maxOutcomeValue = d._maxOutcomeValue;
            _minOutcomeValue = d._minOutcomeValue;
            _logPdf = new double[_maxOutcomeValue+1];
            _logPdfTail = new double[_maxOutcomeValue+1];
            for(int i = 0; i <= _maxOutcomeValue; ++i) {
                _logPdf[i] = d._logPdf[i];
                _logPdfTail[i] = d._logPdfTail[i];
            }
        }

        // Constructs the distribution
        Distribution(vector<double>& pdf) {
            _maxOutcomeValue = pdf.size()-1;
            _minOutcomeValue = -1;
            _logPdf = new double[_maxOutcomeValue+1];
            _logPdfTail = new double[_maxOutcomeValue+1];
            for(int i = 0; i <= _maxOutcomeValue; ++i) {
                if(_minOutcomeValue == -1 && pdf[i] > 0.0) {
                    _minOutcomeValue = i;
                }
                _logPdf[i] = pdf[i];
                _logPdfTail[i] = (i==0)?1.0:_logPdfTail[i-1] - _logPdf[i-1];
            }
            // now take the log of all values
            for(int i = 0; i <= _maxOutcomeValue; ++i) {
                _logPdf[i] = log(_logPdf[i]);
                _logPdfTail[i] = log(_logPdfTail[i]);
            }
        }

        // Returns the log-probability mass function at value k
        double getLogP(int k) {
            return((k <= _maxOutcomeValue)?_logPdf[k]:-INFINITY);
        }

        // Returns the log probability mass function of the truncated distribution at value k
        double getLogPTail(int k) {
            return((k <= _maxOutcomeValue)?_logPdfTail[k]:-INFINITY);
        }

        // Get the maximum outcome value
        int getMaxOutcomeValue() { return _maxOutcomeValue; }

        // Get the minimum outcome value
        int getMinOutcomeValue() { return _minOutcomeValue; }

        // Destructs the distribution
        ~Distribution() {
            delete[] _logPdf;
            delete[] _logPdfTail;
        }

};


class Geometric {
    private:
        double _p0;

    public:
        Geometric(double p0): _p0(p0) {}

        // Returns the log-probability mass function at value k
        double getLogP(int k) {
            return((k<0)?(-INFINITY):((k)*log(1.0-_p0)+log(_p0)));
        }

        // Returns the log probability mass function of the truncated distribution at value k
        double getLogPTail(int k) {
            return((k<0)?0:((k)*log(1.0-_p0)));
        }
};


// This enum captures the four states of infection that an individual can be in
enum SEIRState 
{   Susceptible = 0, 
    Exposed     = 1, 
    Infectious  = 2,
    Recovered   = 3
};

// This enum captures the two test outcomes
enum TestOutcome 
{   Negative = 0, 
    Positive = 1
};

// This class stores all the outcome information
typedef std::tuple<int,int,int> OutcomeTuple;
class Outcome {
    private:
        int         _individual;
        int         _time;
        TestOutcome _outcome;
    public:
        Outcome(int u, int t, TestOutcome o) : 
            _individual(u),
            _time(t),
            _outcome(o) {}
        Outcome(const OutcomeTuple &o): Outcome(get<0>(o), get<1>(o), get<2>(o)==0?Negative:Positive) {}
        int getIndividual() { return(_individual); }
        int getTime() { return(_time); }
        TestOutcome getOutcome() { return(_outcome); }
};

// This class stores all the contact information
typedef std::tuple<int,int,int,int> ContactTuple;
class Contact {
    private:
        int  _fromIndividual;
        int  _toIndividual;
        int  _time;
        int  _count;
    public:
        Contact(int u, int v, int t, int count) : 
            _fromIndividual(u),
            _toIndividual(v),
            _time(t), 
            _count(count) {}; 
        Contact(const ContactTuple &c): Contact(get<0>(c), get<1>(c), get<2>(c), get<3>(c) ) {}
        int getTargetIndividual() const { return(_toIndividual); }
        int getSourceIndividual() const { return(_fromIndividual); }
        int getTime() const { return(_time); }
        int getCount() const { return(_count); }
};

std::ostream &operator<<(std::ostream &os, Contact const &c) {
    return os << "(" << c.getTargetIndividual() << "," << c.getSourceIndividual() << "," << c.getTime() << "," << c.getCount() <<")";
}


// This class captures an infection trace for a single individual (called \mathbf{z}_u in the paper) 
class InfectionTrace {
    private:
    // Every infection trace is unique characterized by three numbers:
    //     (1) The time step until when the patient is susceptible (_t0)
    //     (2) The time step until when the patient is exposed (_t0dE)
    //     (3) The time step until when the patient is infectious (_t0dEdI)
    // After this time step, the patient is recovered.
        int     _t0;
        int     _t0dE;
        int     _t0dEdI;

    public:
        InfectionTrace(const int t0=0, const int dE=0, const int dI=0) : 
            _t0(t0), 
            _t0dE(t0 + dE),
            _t0dEdI(t0+dE+dI) { }

        // Overloading [] operator to access elements in array style 
        SEIRState operator[](int index) const {
            if (index < _t0) {
                return Susceptible;
            } else if (index < _t0dE) {
                return Exposed;
            } else if (index < _t0dEdI) {
                return Infectious;
            }
            return Recovered;
        }

        // Returns the number of susceptible days
        int getT0(void) const { return _t0; }

        // Returns the number of exposed days
        int getDE(void) const { return _t0dE-_t0; }

        // Returns the number of infectious days
        int getDI(void) const { return _t0dEdI-_t0dE; }

        // Sets the number of susceptible days
        inline void setT0(int t0) {
            _t0dEdI += - getT0() + t0;
            _t0dE   += - getT0() + t0;
            _t0 = t0;
        }
        // Sets the number of exposed days
        inline void setDE(int dE) {
            _t0dEdI +=  -getDE() + dE;
            _t0dE   +=  -getDE() + dE;
        }
        // Sets the number of infectious days
        inline void setDI(int dI) {
            _t0dEdI +=  -getDI() + dI;
        }

        // Returns all three numbers in an array
        py::array_t<int> array(void) {
            int *ret = new int[3];
            ret[0] = getT0();
            ret[1] = getDE();
            ret[2] = getDI();
            return py::array_t<int>(3, ret);
        }
};

std::ostream &operator<<(std::ostream &os, InfectionTrace const &it) {
    return os << "(" << it.getT0() << "," << it.getDE() << "," << it.getDI() <<")";
}


// This class captures the infection status of an entire population (variable Z in the paper)
class PopulationInfectionStatus {
    private:
        // Infection trace for every single individual in the simulation
        vector<InfectionTrace> _individualTrace;

        // Total number of individuals
        int _noIndividuals;
        // Total number of time steps
        int _noTimeSteps;

        // Random number generator
        random_device _rd;
        mt19937 _gen;

        // Contacts data
        vector<vector<vector<tuple<int,int>>>> _contacts;
        vector<tuple<int,int>> _empty;
        inline vector<tuple<int,int>>& _futureContact(int u, int t) { return t>=0 && t<_noTimeSteps ? _contacts[u][t] : _empty ; }
        inline vector<tuple<int,int>>& _pastContact(int u, int t) { return t>=1 && t<=_noTimeSteps ? _contacts[u][t-1] : _empty ; }

        inline const vector<tuple<int,int>>& _futureContact(int u, int t) const { return t>=0 && t<_noTimeSteps ? _contacts[u][t] : _empty ; }
        inline const vector<tuple<int,int>>& _pastContact(int u, int t) const { return t>=1 && t<=_noTimeSteps ? _contacts[u][t-1] : _empty ; }

        // Test outcomes for all people
        vector<vector<Outcome> >  _outcomes;

        // Distribution of the length of the susceptible phase
        Geometric _qS;

        // Distribution of the duration of exposure
        Distribution _qE;

        // Distribution of the duration of infectiouness
        Distribution _qI;

        // False-Negative rate of the test outcome
        double _alpha;

        // False-Positive rate of the test outcome
        double _beta;

        // cached value of p0
        double _p0;

        // cached value of p1
        double _p1;

        // cached value of log(1-_p1)
        double _log1MinusP1;

        // Maximum & minimum duration of exposure and infectiousness (depends on the discrete distribution qE and qI)
        int _minExposure;
        int _minInfectious;
        int _maxExposure;
        int _maxInfectious;

        // warning metrics for numerical that can occur at (I -> S-E)-edges
        bool _precisionWarning1Issued;
        double _precisionThreshold1 = 1e-6;

        // determine whether a contact is valid or not
        bool valid(int u, int v, int t) {
            return((u>=0) && (u<_noIndividuals) && (v>=0) && (v<_noIndividuals) && (t>=0) && (t<_noTimeSteps));
        }

        //  list of tuples (infection trace, probability) of all possible infection traces
        vector< tuple< InfectionTrace, double> > _logPrior;

        // compute the list of all possible infection traces (t0, dE, dI) together with their a priori probabilities
        void initPrior() {

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

        // an array that stores the number of infectious past contacts for all u and t
        vector<vector<int>> _infectiousContactCounts;

        // updates _infectiousContactCounts with contribution from zu removed (before) or added (after)
        void updateInfectiousContactCounts(int u, const InfectionTrace &zu, bool before) {
            for(int t=0; t<_noTimeSteps-1; ++t) {
                if(zu[t] == Infectious) {
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
        void print(InfectionTrace z) {
            for(int t=0; t < _noTimeSteps; ++t) {
                switch(z[t]) {
                    case Susceptible: cout << "S"; break;
                    case Exposed: cout << "E"; break;
                    case Infectious: cout << "I"; break;
                    case Recovered: cout << "R"; break;
                    default: cout << "*"; break;
                }
            }
        }

        // advance the whole model by one time step, adding new contacts and tests
        void _advance(const vector<ContactTuple>& contacts, const vector<OutcomeTuple>& outcomes, bool ignore_tests, bool updatePrior) {
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
                    if (_individualTrace[v][_noTimeSteps-2]==Infectious) {
                        _infectiousContactCounts[u][_noTimeSteps-1] += x;
                    }
                }
            }

            // 4. append to _individualTraces
            uniform_real_distribution<> rng(0.0, 1.0);
            for(int u=0; u<_noIndividuals; ++u) {
                auto zu = _individualTrace[u];
                bool tested = false;
                TestOutcome o;
                if(!ignore_tests) {
                    for(auto test: _outcomes[u]) {
                        if (test.getTime()==_noTimeSteps-1) {
                            tested = true;
                            o = test.getOutcome();
                            break;
                        }
                    }
                }
                double rnd =  rng(_gen);
                switch(zu[_noTimeSteps-2]) {
                    case Susceptible: {
			            // Compute f(u, t, Zt) using Equation (5) in paper
                        double p = (1.0-_p0)*pow(1-_p1, _infectiousContactCounts[u][_noTimeSteps-1]);
			            // Transition from S to S with probability f(u, t, Zt), and from S to E with probability 1 - f(u, t, Zt) - see Equation (4) in paper
                        if (rnd<p) {
                            zu.setT0(zu.getT0()+1);
                        } else {
                            zu.setDE(zu.getDE()+1);
                        } }
                        break;

                    case Exposed: {
                        int dE = zu.getDE();
			            // Compute g(u, t, Zut) using Equation (7) in paper
                        double pEE = 1.0-exp(_qE.getLogP(dE)-_qE.getLogPTail(dE));
			            // Adjust with probability of test outcome - see Equation (9) in paper
                        if(tested) {
                            double pEI = 1.0 - pEE;
                            pEE *= o==Positive ? _beta : (1.0-_beta);
                            pEI *= o==Positive ? (1.0-_alpha) : _alpha;
                            pEE /= (pEE+pEI);
                        }
			            // Transition from E to E with probability 1 - g(u, t, Zut), and from E to I with probability g(u, t, Zut) - see Equation (4) in paper
                        if (rnd<pEE){
                            zu.setDE(zu.getDE()+1);
                        } else {
                            zu.setDI(zu.getDI()+1);
                        } }
                        break;

                    case Infectious: {
                        int dI = zu.getDI();
			            // Compute h(u, t, Zut) using Equation (8) in paper
                        double pII = 1.0-exp(_qI.getLogP(dI)-_qI.getLogPTail(dI));
			            // Adjust with probability of test outcome - see Equation (9) in paper
                        if(tested) {
                            double pIR = 1.0 - pII;
                            pII *= o==Positive ? (1.0-_alpha) : _alpha;
                            pIR *= o==Positive ? _beta : (1.0-_beta);
                        }
			            // Transition from I to I with probability 1 - h(u, t, Zut), and from I to R with probability h(u, t, Zut) - see Equation (4) in paper
                        if (rnd<pII){
                            zu.setDI(zu.getDI()+1);
                        } }
                        break;

                    case Recovered:
                        break;
                }
                _individualTrace[u] = zu;
            }
        }


    public:
        // Initializes a population of S individuals over T time steps with contact and outcomes
        PopulationInfectionStatus(int S, int T,
                        const vector<ContactTuple>& contacts, const vector<OutcomeTuple>& outcomes,
                        Distribution& qE, Distribution& qI,
                        double alpha, double beta, double p0, double p1,
                        bool patientZero=false) :
            _individualTrace(S, InfectionTrace(1,0,0)),
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
            _maxInfectious(qI.getMaxOutcomeValue()),
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
                _advance(contacts, outcomes, /*ignore_tests =*/ true, /*updatePrior =*/ false);
            }
            initPrior();
        }

        PopulationInfectionStatus( const PopulationInfectionStatus &other) :
            _individualTrace(other._individualTrace),
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
            _maxInfectious(other._maxInfectious),
            _infectiousContactCounts(other._infectiousContactCounts) {

                    initPrior();
        }

        // implements Gibbs sampling over the whole population
        py::array_t<int> gibbsSample(int N=1, int burnIn=0, int skip=0) {

            _precisionWarning1Issued = false;

            int* sample = new int[N*_noIndividuals*3];

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
                    gibbsSampleU(u);
                    sample[n*_noIndividuals*3 +u*3 + 0] = _individualTrace[u].getT0();
                    sample[n*_noIndividuals*3 +u*3 + 1] = _individualTrace[u].getDE();
                    sample[n*_noIndividuals*3 +u*3 + 2] = _individualTrace[u].getDI();
                }
            }

            return(py::array_t<int>({N, _noIndividuals, 3}, sample));
        }

        // Implements a single Gibbs sampling step according to Algorithm 1 in the paper
        py::array_t<int> gibbsSampleU(int u, int N=1) {

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
                    if(_individualTrace[v][t] == Susceptible && _individualTrace[v][t+1] == Susceptible) {
                        logP += x * _log1MinusP1;
                    }

		             // Compute log f(v, t, Zt^I) - log f(v, t, Zt^{neg I}) for v in CE(u, t) in the paper
                    if(_individualTrace[v][t] == Susceptible && _individualTrace[v][t+1] == Exposed) {
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
                    if (_individualTrace[v][t-1]==Infectious) {
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
                        switch (z[t]) {
                            case Infectious:
				                // Add terms in log B(zu) - see Equation (26) in paper
                                w[i] += active[t];
                                break;
                            case Susceptible:
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

            int* sample = new int [3*N];
            for(int n=0; n<N; ++n) {
		        // Sample i = (t0, dE, dI) with probability proportional to w[i]
                const int sampleIndex = dist(_gen);
                InfectionTrace sample_trace = get<0>(_logPrior[sampleIndex]);

                if (n==0) {
                    // Update log-likelihood and sample trace with the first sample
                    _individualTrace[u] = sample_trace;
                    updateInfectiousContactCounts(u, sample_trace, /*before =*/ false);
                 }

                 sample[3*n]   = sample_trace.getT0();
                 sample[3*n+1] = sample_trace.getDE();
                 sample[3*n+2] = sample_trace.getDI();
            }
            return(py::array_t<int>({N,3}, sample));
        }

        // advance the whole model by one time step, adding new contacts and tests
        void advance(const vector<ContactTuple>& contacts, const vector<OutcomeTuple>& outcomes, bool ignore_tests) {
             return _advance(contacts, outcomes, ignore_tests, true);
        }



        py::array_t<int> getIndividualTrace() const {
            int *res = new int[_noIndividuals*3];
            for(int i=0; i<_noIndividuals; ++i) {
                res[3*i    ] = _individualTrace[i].getT0();
                res[3*i + 1] = _individualTrace[i].getDE();
                res[3*i + 2] = _individualTrace[i].getDI();
            }
            return py::array_t<int>({_noIndividuals,3}, res);
        }

        py::array_t<double> getInfectionStatus(int N=0, int burnIn=0, int skip=0) {
            _precisionWarning1Issued = false;

            double* P_z_uT = new double[_noIndividuals*4]();

            if(N>0) {
                for(int n=0; n<burnIn; ++n) { for(int u=0; u < _noIndividuals; ++u) { gibbsSampleU(u); } }

                for(int n=0; n<N; ++n) {
                    for(int s=0; s<skip; ++s) { for(int u=0; u< _noIndividuals; ++u) { gibbsSampleU(u); } }
                    for(int u=0; u< _noIndividuals; ++u) {
                        gibbsSampleU(u);
                        P_z_uT[u*4 + (int)(_individualTrace[u][_noTimeSteps-1])] += 1;
                    }
                }
                for(int i=0; i<_noIndividuals*4; ++i) {
                    P_z_uT[i] /= N;
                }
            }
            else {
                 for(int u=0; u< _noIndividuals; ++u) {
                        P_z_uT[u*4 + (int)(_individualTrace[u][_noTimeSteps-1])] = 1;
                 }
            }

            return(py::array_t<double>({_noIndividuals,4}, P_z_uT));
        }

        ~PopulationInfectionStatus() {
        }
};

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
        .def("array",  &InfectionTrace::array, py::return_value_policy::move);


    py::class_<PopulationInfectionStatus>(m, "PopulationInfectionStatus")
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
        .def(py::init<const PopulationInfectionStatus &>())
        .def("gibbs_sample", &PopulationInfectionStatus::gibbsSample,
            py::arg("N"), py::arg("burnin")=0, py::arg("skip")=0, py::return_value_policy::move)
        .def("gibbs_sample_u", &PopulationInfectionStatus::gibbsSampleU,
            py::arg("u"), py::arg("N") = 1, py::return_value_policy::move)
        .def("advance", &PopulationInfectionStatus::advance,
                        py::arg("contacts"), py::arg("outcomes"), py::arg("ignore_tests"))
        .def("get_individual_traces", &PopulationInfectionStatus::getIndividualTrace,
                        py::return_value_policy::move)
        .def("get_infection_status", &PopulationInfectionStatus::getInfectionStatus,
            py::arg("N")=0, py::arg("burnin")=0, py::arg("skip")=0, py::return_value_policy::move);
}

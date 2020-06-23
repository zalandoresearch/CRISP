#ifndef _POPULATION_HPP
#define _POPULATION_HPP

#include <vector>
#include <iostream>
# include <random>

using namespace std;

#include "distribution.hpp"




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

std::ostream &operator<<(std::ostream &os, Contact const &c);




class PopulationInfectionStatus {

protected:
        // Total number of individuals S
        int _noIndividuals;

        // Total number of time steps T
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

    // advance the whole model by one time step, adding new contacts and tests
    virtual void _advance(const vector<ContactTuple>& contacts, const vector<OutcomeTuple>& outcomes, bool ignore_tests, bool updatePrior) = 0;


public:
    PopulationInfectionStatus(int S, int T,
                    const vector<ContactTuple>& contacts, const vector<OutcomeTuple>& outcomes,
                    Distribution& qE, Distribution& qI,
                    double alpha, double beta, double p0, double p1,
                    bool patientZero=false);

    PopulationInfectionStatus( const PopulationInfectionStatus& other);

    // advance the whole model by one time step, adding new contacts and tests
    void advance(const vector<ContactTuple>& contacts, const vector<OutcomeTuple>& outcomes, bool ignore_tests) {
            return _advance(contacts, outcomes, ignore_tests, true);
    }

    virtual vector<vector<int>> getIndividualTrace() const = 0;

    virtual vector<vector<double>> getInfectionStatus(int N=0, int burnIn=0, int skip=0) = 0;

};


#endif
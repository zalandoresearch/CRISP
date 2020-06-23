#ifndef _GIBBS_HPP
#define _GIBBS_HPP

#include "seir.hpp"
#include "crisp.hpp" 

// ToDo: make this a proper class header file 

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
        InfectionTrace(const int t0=0, const int dE=0, const int dI=0);  
        // Overloading [] operator to access elements in array style 
        SEIRState operator[](int index) const;

        // Returns the number of susceptible days
        inline int getT0(void) const { return _t0; }

        // Returns the number of exposed days
        inline int getDE(void) const { return _t0dE-_t0; }

        // Returns the number of infectious days
        inline int getDI(void) const { return _t0dEdI-_t0dE; }

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

};



// This class captures the infection status of an entire population (variable Z in the paper)
class GibbsPopulationInfectionStatus : public PopulationInfectionStatus {
    private:
        // Infection trace for every single individual in the simulation
        vector<InfectionTrace> _individualTrace;

        // warning metrics for numerical that can occur at (I -> S-E)-edges
        bool _precisionWarning1Issued;
        double _precisionThreshold1 = 1e-6;

        // determine whether a contact is valid or not
        bool valid(int u, int v, int t);

        //  list of tuples (infection trace, probability) of all possible infection traces
        vector< tuple< InfectionTrace, double> > _logPrior;

        // compute the list of all possible infection traces (t0, dE, dI) together with their a priori probabilities
        void initPrior();

        // an array that stores the number of infectious past contacts for all u and t
        vector<vector<int>> _infectiousContactCounts;

        // updates _infectiousContactCounts with contribution from zu removed (before) or added (after)
        void updateInfectiousContactCounts(int u, const InfectionTrace &zu, bool before);

        // Prints an infection trace on the screen
        void print(InfectionTrace z);

        // advance the whole model by one time step, adding new contacts and tests
        void _advance(const vector<ContactTuple>& contacts, const vector<OutcomeTuple>& outcomes, bool ignore_tests, bool updatePrior);

    public:
        // Initializes a population of S individuals over T time steps with contact and outcomes
        GibbsPopulationInfectionStatus(int S, int T,
                        const vector<ContactTuple>& contacts, const vector<OutcomeTuple>& outcomes,
                        Distribution& qE, Distribution& qI,
                        double alpha, double beta, double p0, double p1,
                        bool patientZero=false);

        GibbsPopulationInfectionStatus( const GibbsPopulationInfectionStatus &other);
        

        // implements Gibbs sampling over the whole population
        vector<vector<vector<int>>> gibbsSample(int N=1, int burnIn=0, int skip=0);

        // Implements a single Gibbs sampling step according to Algorithm 1 in the paper
        vector<vector<int>>  gibbsSampleU(int u, int N=1);

        vector<vector<int>>  getIndividualTrace() const;
        
        vector<vector<double>> getInfectionStatus(int N=0, int burnIn=0, int skip=0);
        
};


#endif
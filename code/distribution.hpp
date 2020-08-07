#ifndef _DISTRIBUTION_HPP
#define _DISTRIBUTION_HPP

#include <vector>
#include <cmath>

using namespace std;

// A class for a distribution over durations
class Distribution {
        double* _logPdf;           // The logarithm of the PDF of the distribution
        double* _logPdfTail;       // The logarithm of 1 minus the CDF of the distribution (equivalent to the log-pdf of the truncated distribution)
        int     _maxOutcomeValue;  // The maximum outcome value
        int     _minOutcomeValue;  // The minimum outcome value
    public:
        // Constructs the distribution with zero entries
        Distribution(); 
        // Copy constructor
        Distribution(const Distribution& d);
        // Constructs the distribution
        Distribution(const vector<double>& pdf);
        // Destructor
        ~Distribution();

        // Returns the log-probability mass function at value k
        double getLogP(int k) const {
            return((k <= _maxOutcomeValue)?_logPdf[k]:-INFINITY);
        }

        // Returns the log probability mass function of the truncated distribution at value k
        double getLogPTail(int k) const {
            return((k <= _maxOutcomeValue)?_logPdfTail[k]:-INFINITY);
        }

        // Get the maximum outcome value
        int getMaxOutcomeValue() const { return _maxOutcomeValue; }

        // Get the minimum outcome value
        int getMinOutcomeValue() const { return _minOutcomeValue; }
};

// A class for the geometric distribution (i.e., waiting time distribution for the first success of a coin toss)
class Geometric {
        double _p0;             // Success probability in each trial
    public:
        Geometric(double p0): _p0(p0) {}

        // Returns the log-probability mass function at value k
        double getLogP(int k) const {
            return((k<0)?(-INFINITY):((k)*log(1.0-_p0)+log(_p0)));
        }

        // Returns the log probability mass function of the truncated distribution at value k
        double getLogPTail(int k) const {
            return((k<0)?0:((k)*log(1.0-_p0)));
        }
};

#endif
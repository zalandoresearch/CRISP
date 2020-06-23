#ifndef _DISTRIBUTION_HPP
#define _DISTRIBUTION_HPP

#include <vector>
#include <cmath>
using namespace std;


// A class for a distribution over durations
class Distribution {
    private:
        double* _logPdf;           // The logarithm of the PDF of the distribution
        double* _logPdfTail;       // The logarithm of 1 minus the CDF of the distribution (equivalent to the log-pdf of the truncated distribution)
        int     _maxOutcomeValue;  // The maximum outcome value
        int     _minOutcomeValue;  // The minimum outcome value
    public:
        // Constructs the distribution with zero entries
        Distribution(); 

        // Copy constructor
        Distribution(const Distribution& d);

        // Destructs the distribution
        ~Distribution();

        // Constructs the distribution
        Distribution(vector<double>& pdf);

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



#endif
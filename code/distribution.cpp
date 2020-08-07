#include "distribution.hpp"

// Constructs the distribution with zero entries
Distribution::Distribution() {
    _maxOutcomeValue = 0;
    _minOutcomeValue = 0;
    _logPdf = new double[_maxOutcomeValue+1];
    _logPdfTail = new double[_maxOutcomeValue+1];
    _logPdf[0] = 0.0;
    _logPdfTail[0] = 0.0;
}

// Destructs the distribution
Distribution::~Distribution() {
    delete[] _logPdf;
    delete[] _logPdfTail;
}

// Copy constructor
Distribution::Distribution(const Distribution& d) {
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
Distribution::Distribution(const vector<double>& pdf) {
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